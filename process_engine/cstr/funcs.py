"""通用全混流反应器（CSTR）独立工具函数。

设计原则：
- 每个 public calc_* 只接收原始题目参数
- public calc_* 之间禁止互相调用
- 私有 _ 辅助函数和常量表可共享
- 每个函数只返回自己的最终结果

全程国际单位制。

入参约定：
- 所有接收 inlet_molar_flows + stoichiometric_coefficients 的函数，inlet_molar_flows
  必须包含 stoichiometric_coefficients 中出现的所有组分（反应物 + 产物），
  即使某产物在进料中流量为 0 也必须显式写入（如 {"产物": 0.0}），
  否则该产物不会出现在出口流量字典中。

修复记录（化工角度）：
  [FIX-1] _compute_outlet_molar_flows 增加出口流量非负校验：若转化率超过
           化学计量上限导致任一组分出口流量为负，抛出含定量信息的 ValueError。
  [FIX-2] _compute_reaction_extent 校验关键组分计量系数必须为负（反应物），
           误传产物为关键组分时给出明确错误。
  [FIX-3] _compute_archimedes_number 校验颗粒密度须大于流体密度（密度差 > 0），
           避免 sqrt 负数崩溃；_compute_minimum_fluidization_velocity 校验 umf > 0。
  [FIX-4] _compute_volumetric_flow_at_TP 液相后备路径逐组分 Vml 非 None 保护，
           任一组分 Vml 缺失时改为 ValueError 而非 TypeError。
  [FIX-5] _compute_fluid_density 气相路径改用 Mixture.MW（单次调用）替代
           逐组分循环建 Mixture 求混合摩尔质量，同时增加 MW 合理性校验。
  [FIX-6] calc_residence_time method='ghsv' 且无反应信息时禁止用进口流量
           近似出口体积流量：强制要求提供反应参数或明确报错，避免摩尔数
           变化体系（气相分解等）产生数倍误差。
  [FIX-7] calc_cstr_volume / calc_residence_time method='ghsv' 补充 STP 与
           NTP 的说明性注释，并允许调用方通过 ghsv_reference_T_K /
           ghsv_reference_P_Pa 指定参考状态（默认保持原 STP）。
"""

import os

_VERBOSE = os.getenv("FUNCS_VERBOSE", "0") == "1"


def _vprint(*args, **kwargs):
    if _VERBOSE:
        print(*args, **kwargs)


import math

from thermo import Mixture, Chemical
from thermo import ChemicalConstantsPackage
from thermo.flash import FlashVL
from thermo.phases import CEOSGas, CEOSLiquid
from thermo.eos_mix import PRMIX


# ─────────────────────────────────────────────
# 模块级常量
# ─────────────────────────────────────────────
IDEAL_GAS_CONSTANT_J_per_mol_K = 8.314
STP_TEMPERATURE_K = 273.15      # 标准状态温度（IUPAC STP）
STP_PRESSURE_Pa   = 101_325     # 标准状态压力（IUPAC STP）
# 注：工业 GHSV 有时用 NTP（20°C, 101325 Pa），与 STP 差约 3.7%；
# 可通过 calc_cstr_volume / calc_residence_time 的 ghsv_reference_T_K 参数指定。


# ─────────────────────────────────────────────
# PR EOS 辅助函数（私有，可共享）
# ─────────────────────────────────────────────
def _make_flasher(components):
    constants, correlations = ChemicalConstantsPackage.from_IDs(components)
    eos_kwargs = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas)
    gas = CEOSGas(PRMIX, HeatCapacityGases=correlations.HeatCapacityGases, eos_kwargs=eos_kwargs)
    liquid = CEOSLiquid(PRMIX, HeatCapacityGases=correlations.HeatCapacityGases, eos_kwargs=eos_kwargs)
    return FlashVL(constants, correlations, liquid=liquid, gas=gas)


def _calc_phase_pr(components, zs, T, P):
    flasher = _make_flasher(components)
    result = flasher.flash(T=T, P=P, zs=zs)
    return 'l' if result.VF < 0.5 else 'g'


# ─────────────────────────────────────────────
# 私有：核心计算（无 print，纯数学，可共享）
# ─────────────────────────────────────────────
def _compute_reaction_extent(inlet_molar_flows, stoichiometric_coefficients,
                              key_component, key_component_conversion):
    """计算反应进度 ξ (mol/s)。

    [FIX-2] 校验关键组分计量系数必须 < 0（即反应物），
    防止误传产物导致 ξ 方向错误而不报错。
    """
    if key_component not in inlet_molar_flows:
        raise ValueError(
            f"关键组分 '{key_component}' 不在 inlet_molar_flows 中: "
            f"{list(inlet_molar_flows.keys())}"
        )
    if key_component not in stoichiometric_coefficients:
        raise ValueError(
            f"关键组分 '{key_component}' 不在 stoichiometric_coefficients 中: "
            f"{list(stoichiometric_coefficients.keys())}"
        )
    coeff_key = stoichiometric_coefficients[key_component]

    # [FIX-2] 关键组分必须是反应物（计量系数 < 0）
    if coeff_key >= 0:
        raise ValueError(
            f"关键组分 '{key_component}' 的计量系数为 {coeff_key}（≥ 0），"
            f"关键组分必须是反应物（计量系数 < 0）。"
            f"请勿将产物指定为关键组分。"
        )
    if not (0.0 < key_component_conversion <= 1.0):
        raise ValueError(
            f"key_component_conversion 须在 (0, 1] 之间，得到 {key_component_conversion}"
        )

    inlet_key = inlet_molar_flows[key_component]
    if inlet_key <= 0:
        raise ValueError(
            f"关键组分 '{key_component}' 的进口摩尔流率须 > 0，"
            f"得到 {inlet_key} mol/s"
        )
    coeff_key_abs = abs(coeff_key)
    return (inlet_key * key_component_conversion) / coeff_key_abs


def _compute_outlet_molar_flows(inlet_molar_flows, stoichiometric_coefficients,
                                 key_component, key_component_conversion):
    """计算出口各组分摩尔流率 (mol/s)。

    [FIX-1] 检查所有出口流量 >= 0：若转化率超过化学计量上限，某反应物
    出口流量会变负，此时应报错而非传递负流量到下游物性计算。
    """
    extent = _compute_reaction_extent(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion
    )
    outlet = {
        name: flow + stoichiometric_coefficients.get(name, 0) * extent
        for name, flow in inlet_molar_flows.items()
    }

    # [FIX-1] 非负性校验（允许极小负值作数值误差容忍，截断为 0）
    _TOL_NEG = -1e-9  # mol/s，小于此绝对值视为数值误差
    for name, flow in outlet.items():
        if flow < _TOL_NEG:
            coeff = stoichiometric_coefficients.get(name, 0)
            max_conversion = (
                abs(inlet_molar_flows[name] / (coeff * (inlet_molar_flows[key_component]
                    / abs(stoichiometric_coefficients[key_component]))))
                if (coeff < 0 and inlet_molar_flows.get(name, 0) > 0) else float('inf')
            )
            raise ValueError(
                f"组分 '{name}' 出口摩尔流率为 {flow:.4e} mol/s（<0），"
                f"转化率 {key_component_conversion:.4f} 超过该组分的化学计量上限。"
                f"请降低转化率或检查计量系数。"
            )
        # 数值误差截断
        if _TOL_NEG <= flow < 0:
            outlet[name] = 0.0

    return outlet


def _compute_outlet_volumetric_flow(inlet_molar_flows, stoichiometric_coefficients,
                                     key_component, key_component_conversion,
                                     temperature_K, pressure_Pa):
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    return _compute_volumetric_flow_at_TP(outlet, temperature_K, pressure_Pa)


def _compute_volumetric_flow_at_TP(molar_flows, T, P):
    """通用：根据组成 + T/P 算体积流量。返回 (v, phase, Vml)。

    [FIX-4] 液相后备路径逐组分 Vml 增加非 None 保护。
    """
    components = list(molar_flows.keys())
    flows = list(molar_flows.values())
    active = [(c, f) for c, f in zip(components, flows) if f > 1e-12]
    if not active:
        raise ValueError("所有组分流量均为 0 或极小，无法计算体积流量")
    components = [c for c, _ in active]
    flows = [f for _, f in active]
    total = sum(flows)
    mole_fractions = [f / total for f in flows]

    if len(active) == 1:
        chem_single = Chemical(components[0], T=T, P=P)
        phase = chem_single.phase
    else:
        phase = _calc_phase_pr(components, mole_fractions, T, P)

    mixture = Mixture(IDs=components, zs=mole_fractions, T=T, P=P)

    if phase == 'l':
        Vml = mixture.Vml
        if Vml is None or Vml <= 0:
            # [FIX-4] 后备：逐组分加权，但必须校验每个 Vml 不为 None
            Vml_parts = []
            for c, z in zip(components, mole_fractions):
                single_mix = Mixture([c], zs=[1], T=T, P=P)
                vml_single = single_mix.Vml
                if vml_single is None or vml_single <= 0:
                    raise ValueError(
                        f"组分 '{c}' 在 T={T:.2f}K, P={P:.0f}Pa 下液相摩尔体积 "
                        f"Vml={vml_single}（无效），无法计算液相体积流量。"
                        f"请检查操作条件是否超出 thermo 数据范围，或该组分是否确为液相。"
                    )
                Vml_parts.append(z * vml_single)
            Vml = sum(Vml_parts)
        return total * Vml, 'l', Vml

    # 气相：理想气体
    v = (total * IDEAL_GAS_CONSTANT_J_per_mol_K * T) / P
    return v, 'g', None


def _compute_inlet_vol_flow_at_ref(inlet_molar_flows, ref_T_K, ref_P_Pa):
    """在参考状态 (ref_T_K, ref_P_Pa) 下计算进口体积流量。

    [FIX-7] 从原始 _compute_inlet_vol_flow_stp（固定 STP）改为接受参考状态参数，
    支持工业 NTP（293.15K）或用户自定义参考态。
    """
    v, _, _ = _compute_volumetric_flow_at_TP(inlet_molar_flows, ref_T_K, ref_P_Pa)
    return v


def _compute_fluid_density(inlet_molar_flows, temperature_K, pressure_Pa):
    """返回 (rho, phase)。

    [FIX-5] 气相混合摩尔质量改用 Mixture.MW（单次调用），
    不再逐组分循环建 Mixture，并增加 MW 合理性校验。
    """
    active = {c: f for c, f in inlet_molar_flows.items() if f > 0}
    if not active:
        raise ValueError("所有组分流量为 0，无法计算密度")
    components = list(active.keys())
    flows = list(active.values())
    total = sum(flows)
    mole_fractions = [f / total for f in flows]

    if len(active) == 1:
        chem_single = Chemical(components[0], T=temperature_K, P=pressure_Pa)
        phase = chem_single.phase
    else:
        phase = _calc_phase_pr(components, mole_fractions, temperature_K, pressure_Pa)

    mixture = Mixture(IDs=components, zs=mole_fractions, T=temperature_K, P=pressure_Pa)

    if phase == 'l':
        rho_l = mixture.rhol
        if rho_l is None or rho_l <= 0:
            raise ValueError(
                f"混合物在 T={temperature_K:.2f}K, P={pressure_Pa:.0f}Pa 下"
                f"液相密度 rhol={rho_l}（无效）。请检查操作条件或组分数据。"
            )
        return rho_l, 'l'

    # [FIX-5] 气相：用 Mixture.MW 一次性获取混合摩尔质量
    MW_mix_g_per_mol = mixture.MW  # g/mol
    if MW_mix_g_per_mol is None or MW_mix_g_per_mol <= 0:
        raise ValueError(
            f"无法获取混合物摩尔质量（MW={MW_mix_g_per_mol}），"
            f"请检查组分名称是否 thermo 可识别。"
        )
    MW_mix_kg_per_mol = MW_mix_g_per_mol / 1000.0
    rho = (pressure_Pa * MW_mix_kg_per_mol) / (IDEAL_GAS_CONSTANT_J_per_mol_K * temperature_K)
    return rho, 'g'


def _compute_fluid_viscosity(inlet_molar_flows, temperature_K, pressure_Pa):
    """返回 (mu, phase)。"""
    active = {c: f for c, f in inlet_molar_flows.items() if f > 0}
    if not active:
        raise ValueError("所有组分流量为 0，无法计算粘度")
    components = list(active.keys())
    flows = list(active.values())
    total = sum(flows)
    mole_fractions = [f / total for f in flows]

    if len(active) == 1:
        chem_single = Chemical(components[0], T=temperature_K, P=pressure_Pa)
        phase = chem_single.phase
    else:
        phase = _calc_phase_pr(components, mole_fractions, temperature_K, pressure_Pa)

    mixture = Mixture(IDs=components, zs=mole_fractions, T=temperature_K, P=pressure_Pa)
    if phase == 'l':
        mu = mixture.mul
        if mu is None or mu <= 0:
            raise ValueError(
                f"混合物在 T={temperature_K:.2f}K 下液相粘度 mul={mu}（无效）。"
            )
        return mu, 'l'
    mu = mixture.mug
    if mu is None or mu <= 0:
        raise ValueError(
            f"混合物在 T={temperature_K:.2f}K 下气相粘度 mug={mu}（无效）。"
        )
    return mu, 'g'


def _compute_archimedes_number(particle_diameter_m, gas_density_kg_per_m3,
                                particle_density_kg_per_m3, gas_viscosity_Pa_s,
                                gravity_acceleration_m_per_s2):
    """计算阿基米德数 Ar。

    [FIX-3] 校验颗粒密度必须大于流体密度（密度差 > 0），
    否则 Ar < 0，sqrt 会崩溃或 umf 无物理意义。
    """
    delta_rho = particle_density_kg_per_m3 - gas_density_kg_per_m3
    if delta_rho <= 0:
        raise ValueError(
            f"颗粒密度 ({particle_density_kg_per_m3:.2f} kg/m³) 必须大于"
            f"流体密度 ({gas_density_kg_per_m3:.4f} kg/m³)，"
            f"否则颗粒不会沉降，Ar < 0 无物理意义。"
            f"请检查颗粒密度参数或确认流体相态。"
        )
    if gas_viscosity_Pa_s <= 0:
        raise ValueError(
            f"流体粘度须 > 0，得到 {gas_viscosity_Pa_s} Pa·s"
        )
    if particle_diameter_m <= 0:
        raise ValueError(
            f"颗粒直径须 > 0，得到 {particle_diameter_m} m"
        )
    return (particle_diameter_m**3 * gas_density_kg_per_m3
            * delta_rho
            * gravity_acceleration_m_per_s2) / gas_viscosity_Pa_s**2


def _compute_minimum_fluidization_velocity(archimedes_number, gas_viscosity_Pa_s,
                                            gas_density_kg_per_m3, particle_diameter_m):
    """Wen-Yu 关联式计算 umf (m/s)。

    [FIX-3] 校验判别式 >= 0（理论上 Ar > 0 时恒成立，作防御性检查）
    并校验 umf > 0。
    """
    discriminant = 33.7**2 + 0.0408 * archimedes_number
    if discriminant < 0:
        raise ValueError(
            f"Wen-Yu 关联式判别式 < 0（discriminant={discriminant:.4f}，Ar={archimedes_number:.4e}），"
            f"无法计算 umf，请检查 Ar 值。"
        )
    umf = (gas_viscosity_Pa_s / (gas_density_kg_per_m3 * particle_diameter_m)) * \
          (math.sqrt(discriminant) - 33.7)
    if umf <= 0:
        raise ValueError(
            f"计算得 umf = {umf:.4e} m/s（≤ 0），物理上不合理。"
            f"Ar = {archimedes_number:.4e}，请检查颗粒与流体参数。"
        )
    return umf


# ═════════════════════════════════════════════
# 工具 1：反应进度
# ═════════════════════════════════════════════
def calc_reaction_extent(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
) -> dict:
    """计算 CSTR 全混流反应器的反应进度 ξ（单位 mol/s）。

    反应进度是描述反应推进程度的基准量，所有组分变化遵循统一计量关系：
    只需知道关键组分转化率，即可通过 ξ 推算所有组分变化。
    ξ = (关键组分进口流率 × 转化率) / |关键组分计量系数|。

    参数:
        inlet_molar_flows           — {组分名: 进口摩尔流率} (mol/s)
                                       必须包含 stoichiometric_coefficients 中的所有组分，
                                       包括进料流量为 0 的纯产物（如 {"产物": 0.0}）
        stoichiometric_coefficients — {组分名: 计量系数}，关键组分系数须 < 0（反应物）
        key_component               — 关键组分名称（必须是反应物，计量系数 < 0）
        key_component_conversion    — 关键组分转化率 (0, 1]

    返回:
        {"reaction_extent_mol_per_s": ξ (mol/s)}
    """
    extent = _compute_reaction_extent(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    _vprint(f"【calc_reaction_extent】ξ = {extent:.4f} mol/s")
    return {"reaction_extent_mol_per_s": extent}


# ═════════════════════════════════════════════
# 工具 2：出口各组分摩尔流率
# ═════════════════════════════════════════════
def calc_outlet_molar_flows(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
) -> dict:
    """计算 CSTR 各组分出口摩尔流率（单位 mol/s）。

    根据物料守恒：出口流率 = 进口流率 + 计量系数 × 反应进度。反应物（计量系数
    为负）出口减少，产物（为正）出口增加，惰性组分不变。本工具独立计算。

    参数:
        inlet_molar_flows           — {组分名: 进口摩尔流率} (mol/s)
                                       必须包含 stoichiometric_coefficients 中的所有组分，
                                       包括进料流量为 0 的纯产物（如 {"产物": 0.0}）
        stoichiometric_coefficients — {组分名: 计量系数}，关键组分系数须 < 0（反应物）
        key_component               — 关键组分名称
        key_component_conversion    — 关键组分转化率

    返回:
        {组分名: 出口摩尔流率 (mol/s)}
    """
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    _vprint("\n【calc_outlet_molar_flows】")
    for sp, f in outlet.items():
        _vprint(f"  F_{sp} = {f:.4f} mol/s")
    return outlet


# ═════════════════════════════════════════════
# 工具 3：出口总摩尔流率
# ═════════════════════════════════════════════
def calc_outlet_total_molar_flow(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
) -> dict:
    """计算 CSTR 出口总摩尔流率（单位 mol/s）。

    由于反应前后分子数可能变化，出口总流率可能与进口不同。总流率是计算体积
    流量和浓度的基础。本工具独立计算。

    返回:
        {"outlet_total_molar_flow_mol_per_s": F_total (mol/s)}
    """
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    total = sum(outlet.values())
    _vprint(f"\n【calc_outlet_total_molar_flow】F_total = {total:.4f} mol/s")
    return {"outlet_total_molar_flow_mol_per_s": total}


# ═════════════════════════════════════════════
# 工具 4：出口体积流量
# ═════════════════════════════════════════════
def calc_outlet_volumetric_flow(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
    temperature_K: float,
    pressure_Pa: float,
) -> dict:
    """计算 CSTR 出口体积流量（单位 m³/s）。

    根据出口组成、温度、压力自动判断相态：气相用理想气体方程 PV=nRT，
    液相用 thermo 库液体摩尔体积 Vml。本工具独立计算。

    返回:
        {"outlet_volumetric_flow_m3_per_s": v (m³/s)}
    """
    v, phase, Vml = _compute_outlet_volumetric_flow(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
        temperature_K, pressure_Pa)
    _vprint("\n【calc_outlet_volumetric_flow】")
    if phase == 'l':
        _vprint(f"  相态判定：液相，Vml = {Vml:.6e} m^3/mol")
    else:
        _vprint(f"  相态判定：气相，理想气体定律")
    _vprint(f"  v = {v:.6f} m^3/s")
    return {"outlet_volumetric_flow_m3_per_s": v}


# ═════════════════════════════════════════════
# 工具 5：出口浓度
# ═════════════════════════════════════════════
def calc_outlet_concentrations(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
    temperature_K: float,
    pressure_Pa: float,
) -> dict:
    """计算 CSTR 出口各组分浓度（单位 mol/m³）。

    浓度 = 摩尔流率 / 体积流量，反映反应后各组分在混合物中的密集程度。
    CSTR 全混流模型中，反应器内浓度等于出口浓度，本工具算出的即为反应器内浓度。
    本工具独立计算。

    返回:
        {组分名: 出口浓度 (mol/m³)}
    注意:
        inlet_molar_flows 必须包含 stoichiometric_coefficients 中的所有组分，
        包括进料流量为 0 的纯产物（如 {"产物": 0.0}）。
    """
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    v, _, _ = _compute_outlet_volumetric_flow(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
        temperature_K, pressure_Pa)
    concentrations = {name: flow / v for name, flow in outlet.items()}
    _vprint("\n【calc_outlet_concentrations】")
    for k, c in concentrations.items():
        _vprint(f"  {k} = {c:.4f} mol/m^3")
    return concentrations


# ═════════════════════════════════════════════
# 工具 6：流体密度
# ═════════════════════════════════════════════
def calc_fluid_density(
    inlet_molar_flows: dict,
    temperature_K: float,
    pressure_Pa: float,
) -> dict:
    """计算进料混合流体密度（单位 kg/m³）。

    自动判断相态：气相用理想气体定律 ρ=PM/(RT)，液相用 thermo 库液体密度。

    返回:
        {"fluid_density_kg_per_m3": ρ (kg/m³)}
    """
    rho, phase = _compute_fluid_density(inlet_molar_flows, temperature_K, pressure_Pa)
    _vprint(f"\n【calc_fluid_density】")
    _vprint(f"  相态判定：{'液相' if phase == 'l' else '气相'}")
    _vprint(f"  rho = {rho:.4f} kg/m^3")
    return {"fluid_density_kg_per_m3": rho}


# ═════════════════════════════════════════════
# 工具 7：流体粘度
# ═════════════════════════════════════════════
def calc_fluid_viscosity(
    inlet_molar_flows: dict,
    temperature_K: float,
    pressure_Pa: float,
) -> dict:
    """计算进料混合流体动力粘度（单位 Pa·s）。

    自动判断相态，用 thermo 库获取气相或液相粘度。

    返回:
        {"fluid_viscosity_Pa_s": μ (Pa·s)}
    """
    mu, phase = _compute_fluid_viscosity(inlet_molar_flows, temperature_K, pressure_Pa)
    _vprint(f"\n【calc_fluid_viscosity】")
    _vprint(f"  相态判定：{'液相' if phase == 'l' else '气相'}")
    _vprint(f"  mu = {mu:.6e} Pa·s")
    return {"fluid_viscosity_Pa_s": mu}


# ═════════════════════════════════════════════
# 工具 8：阿基米德数
# ═════════════════════════════════════════════
def calc_archimedes_number(
    particle_diameter_m: float,
    particle_density_kg_per_m3: float,
    inlet_molar_flows: dict,
    temperature_K: float,
    pressure_Pa: float,
    gravity_acceleration_m_per_s2: float = 9.81,
) -> dict:
    """计算阿基米德数 Ar（无量纲）。

    Ar 表征颗粒在流体中重力与粘性力的比值，是判断流化状态的关键准则数。
    Ar = dp³·ρg·(ρp−ρg)·g / μ²。

    参数:
        particle_diameter_m              — 颗粒直径 (m)，须 > 0
        particle_density_kg_per_m3       — 颗粒密度 (kg/m³)，须 > 流体密度
        inlet_molar_flows                — {组分名: 进口摩尔流率} (mol/s)
        temperature_K                    — 温度 (K)
        pressure_Pa                      — 压力 (Pa)
        gravity_acceleration_m_per_s2    — 重力加速度 (m/s²)，默认 9.81

    返回:
        {"archimedes_number": Ar}
    """
    rho, _ = _compute_fluid_density(inlet_molar_flows, temperature_K, pressure_Pa)
    mu, _ = _compute_fluid_viscosity(inlet_molar_flows, temperature_K, pressure_Pa)
    ar = _compute_archimedes_number(
        particle_diameter_m, rho, particle_density_kg_per_m3, mu, gravity_acceleration_m_per_s2)
    _vprint(f"\n【calc_archimedes_number】Ar = {ar:.4e}")
    return {"archimedes_number": ar}


# ═════════════════════════════════════════════
# 工具 9：最小流化速度
# ═════════════════════════════════════════════
def calc_minimum_fluidization_velocity(
    particle_diameter_m: float,
    particle_density_kg_per_m3: float,
    inlet_molar_flows: dict,
    temperature_K: float,
    pressure_Pa: float,
) -> dict:
    """用 Wen-Yu 关联式计算最小流化速度 umf（单位 m/s）。

    umf 是颗粒刚好被流体托起的临界速度，操作气速必须大于 umf 才能实现流化。
    本工具独立计算：内部完成密度、粘度、Ar 推算。

    参数:
        particle_diameter_m         — 颗粒直径 (m)，须 > 0
        particle_density_kg_per_m3  — 颗粒密度 (kg/m³)，须大于流体密度
        inlet_molar_flows           — {组分名: 进口摩尔流率} (mol/s)
        temperature_K               — 温度 (K)
        pressure_Pa                 — 压力 (Pa)

    返回:
        {"minimum_fluidization_velocity_m_per_s": umf (m/s)}
    """
    rho, _ = _compute_fluid_density(inlet_molar_flows, temperature_K, pressure_Pa)
    mu, _ = _compute_fluid_viscosity(inlet_molar_flows, temperature_K, pressure_Pa)
    ar = _compute_archimedes_number(particle_diameter_m, rho, particle_density_kg_per_m3, mu, 9.81)
    umf = _compute_minimum_fluidization_velocity(ar, mu, rho, particle_diameter_m)
    _vprint(f"\n【calc_minimum_fluidization_velocity】umf = {umf:.4f} m/s")
    return {"minimum_fluidization_velocity_m_per_s": umf}


# ═════════════════════════════════════════════
# 工具 10：反应速率
# ═════════════════════════════════════════════
def calc_reaction_rate(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
    temperature_K: float,
    pressure_Pa: float,
    arrhenius_pre_exponential_factor: float,
    activation_energy_J_mol: float,
    reaction_order: int = 1,
) -> dict:
    """计算关键组分反应速率 -rA（单位 mol/(m³·s)）。

    采用 Arrhenius 速率方程 k=A·exp(-Ea/RT)，结合出口浓度计算 n 级反应速率
    -rA = k·CA^n。

    CSTR 全混流模型中反应器内浓度 = 出口浓度，因此本函数用出口浓度代入速率
    方程是正确的。本工具独立计算：内部完成出口浓度推算。

    参数:
        inlet_molar_flows                — {组分名: 进口摩尔流率} (mol/s)
        stoichiometric_coefficients      — {组分名: 计量系数}，关键组分系数须 < 0
        key_component                    — 关键组分名称
        key_component_conversion         — 关键组分转化率 (0, 1]
        temperature_K                    — 反应温度 (K)
        pressure_Pa                      — 反应压力 (Pa)
        arrhenius_pre_exponential_factor — Arrhenius 指前因子（单位与 reaction_order 相关）
        activation_energy_J_mol          — 活化能 (J/mol)，须 >= 0
        reaction_order                   — 反应级数，默认 1，须 >= 0

    返回:
        {"reaction_rate_mol_per_m3_s": -rA (mol/(m³·s))}
    """
    if activation_energy_J_mol < 0:
        raise ValueError(
            f"活化能须 >= 0，得到 {activation_energy_J_mol} J/mol"
        )
    if reaction_order < 0:
        raise ValueError(
            f"反应级数须 >= 0，得到 {reaction_order}"
        )

    # 内部独立算出口浓度（= CSTR 反应器内浓度，全混流模型）
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    v, _, _ = _compute_outlet_volumetric_flow(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
        temperature_K, pressure_Pa)
    key_concentration = outlet[key_component] / v

    if key_concentration < 0:
        raise ValueError(
            f"关键组分 '{key_component}' 出口浓度为 {key_concentration:.4e} mol/m³（< 0），"
            f"请检查转化率和出口流量。"
        )

    rate_constant = arrhenius_pre_exponential_factor * math.exp(
        -activation_energy_J_mol / (IDEAL_GAS_CONSTANT_J_per_mol_K * temperature_K))
    rate = rate_constant * key_concentration**reaction_order

    _vprint(f"\n【calc_reaction_rate】")
    _vprint(f"  CA_out = {key_concentration:.4e} mol/m³（CSTR 全混流，等于反应器内浓度）")
    _vprint(f"  k = {rate_constant:.4e} (Arrhenius)")
    _vprint(f"  rA = {rate:.4e} mol/(m^3·s)")
    return {"reaction_rate_mol_per_m3_s": rate}


# ═════════════════════════════════════════════
# 工具 11：CSTR 体积
# ═════════════════════════════════════════════
def calc_cstr_volume(
    method: str = "kinetic",
    inlet_molar_flows: dict = None,
    stoichiometric_coefficients: dict = None,
    key_component: str = None,
    key_component_conversion: float = None,
    temperature_K: float = None,
    pressure_Pa: float = None,
    arrhenius_pre_exponential_factor: float = None,
    activation_energy_J_mol: float = None,
    reaction_order: int = 1,
    space_velocity_ghsv: float = None,
    ghsv_time_unit: str = "s",
    ghsv_reference_T_K: float = STP_TEMPERATURE_K,
    ghsv_reference_P_Pa: float = STP_PRESSURE_Pa,
) -> dict:
    """计算 CSTR 反应器体积（单位 m³）。

    支持两种方法：
    - method="kinetic"：CSTR 设计方程 V = FA0·X/(-rA)，内部独立计算反应速率。
      注：FA0 取 inlet_molar_flows[key_component]，转化率 X 即 key_component_conversion。
    - method="ghsv"：V = v_ref / GHSV，内部从 inlet_molar_flows 推算参考态体积流量。

    参数（method="ghsv" 新增）:
        ghsv_reference_T_K  — GHSV 参考态温度 (K)，默认 273.15K（IUPAC STP）。
                               工业习惯 NTP 请传 293.15K（20°C）。
        ghsv_reference_P_Pa — GHSV 参考态压力 (Pa)，默认 101325 Pa。

    返回:
        {"cstr_volume_m3": V (m³)}
    """
    if method == "kinetic":
        if (inlet_molar_flows is None or stoichiometric_coefficients is None
                or key_component is None or key_component_conversion is None
                or temperature_K is None or pressure_Pa is None
                or arrhenius_pre_exponential_factor is None
                or activation_energy_J_mol is None):
            raise ValueError("method='kinetic' 需要全部反应动力学参数")

        outlet = _compute_outlet_molar_flows(
            inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
        v_out, _, _ = _compute_outlet_volumetric_flow(
            inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
            temperature_K, pressure_Pa)
        key_conc = outlet[key_component] / v_out
        rate_constant = arrhenius_pre_exponential_factor * math.exp(
            -activation_energy_J_mol / (IDEAL_GAS_CONSTANT_J_per_mol_K * temperature_K))
        rA = rate_constant * key_conc**reaction_order

        FA0 = inlet_molar_flows[key_component]
        vol = (FA0 * key_component_conversion) / rA
        _vprint(f"\n【calc_cstr_volume】(kinetic) V = {vol:.4f} m^3")
        return {"cstr_volume_m3": vol}

    elif method == "ghsv":
        if inlet_molar_flows is None or space_velocity_ghsv is None:
            raise ValueError("method='ghsv' 需要 inlet_molar_flows 与 space_velocity_ghsv")
        if ghsv_time_unit == "h":
            ghsv_per_s = space_velocity_ghsv / 3600
        elif ghsv_time_unit == "s":
            ghsv_per_s = space_velocity_ghsv
        else:
            raise ValueError("ghsv_time_unit 须为 's' 或 'h'")

        # [FIX-7] 用参数化参考态替代固定 STP
        v_ref = _compute_inlet_vol_flow_at_ref(
            inlet_molar_flows, ghsv_reference_T_K, ghsv_reference_P_Pa
        )
        vol = v_ref / ghsv_per_s
        _vprint(
            f"\n【calc_cstr_volume】(ghsv) 参考态 T={ghsv_reference_T_K}K, "
            f"P={ghsv_reference_P_Pa}Pa, v_ref={v_ref:.4e} m³/s, V = {vol:.4f} m^3"
        )
        return {"cstr_volume_m3": vol}

    raise ValueError("method 须为 'kinetic' 或 'ghsv'")


# ═════════════════════════════════════════════
# 工具 12：停留时间
# ═════════════════════════════════════════════
def calc_residence_time(
    method: str = "kinetic",
    inlet_molar_flows: dict = None,
    stoichiometric_coefficients: dict = None,
    key_component: str = None,
    key_component_conversion: float = None,
    temperature_K: float = None,
    pressure_Pa: float = None,
    arrhenius_pre_exponential_factor: float = None,
    activation_energy_J_mol: float = None,
    reaction_order: int = 1,
    space_velocity_ghsv: float = None,
    ghsv_time_unit: str = "s",
    ghsv_reference_T_K: float = STP_TEMPERATURE_K,
    ghsv_reference_P_Pa: float = STP_PRESSURE_Pa,
) -> dict:
    """计算 CSTR 停留时间 τ（单位 s）。

    τ = V / v_out，表示物料在反应器内的平均停留时间。本工具独立计算：内部完成
    反应器体积与出口体积流量推算。

    [FIX-6] 当 method='ghsv' 时，τ = V/v_out 仍需出口体积流量，
    因此必须提供反应参数（stoichiometric_coefficients + key_component +
    key_component_conversion + temperature_K + pressure_Pa）。
    若上述参数缺失，不再用进口流量静默近似，而是抛出明确错误。

    参数同 calc_cstr_volume，新增:
        ghsv_reference_T_K  — GHSV 参考态温度 (K)，默认 STP（273.15K）
        ghsv_reference_P_Pa — GHSV 参考态压力 (Pa)，默认 101325 Pa

    返回:
        {"residence_time_s": τ (s)}
    """
    # 1) 内部独立算 V
    if method == "kinetic":
        if (inlet_molar_flows is None or stoichiometric_coefficients is None
                or key_component is None or key_component_conversion is None
                or temperature_K is None or pressure_Pa is None
                or arrhenius_pre_exponential_factor is None
                or activation_energy_J_mol is None):
            raise ValueError("method='kinetic' 需要全部反应动力学参数")
        outlet = _compute_outlet_molar_flows(
            inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
        v_out_for_rate, _, _ = _compute_outlet_volumetric_flow(
            inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
            temperature_K, pressure_Pa)
        key_conc = outlet[key_component] / v_out_for_rate
        rate_constant = arrhenius_pre_exponential_factor * math.exp(
            -activation_energy_J_mol / (IDEAL_GAS_CONSTANT_J_per_mol_K * temperature_K))
        rA = rate_constant * key_conc**reaction_order
        FA0 = inlet_molar_flows[key_component]
        V = (FA0 * key_component_conversion) / rA

    elif method == "ghsv":
        if inlet_molar_flows is None or space_velocity_ghsv is None:
            raise ValueError("method='ghsv' 需要 inlet_molar_flows 与 space_velocity_ghsv")
        if ghsv_time_unit == "h":
            ghsv_per_s = space_velocity_ghsv / 3600
        elif ghsv_time_unit == "s":
            ghsv_per_s = space_velocity_ghsv
        else:
            raise ValueError("ghsv_time_unit 须为 's' 或 'h'")
        v_ref = _compute_inlet_vol_flow_at_ref(
            inlet_molar_flows, ghsv_reference_T_K, ghsv_reference_P_Pa
        )
        V = v_ref / ghsv_per_s
    else:
        raise ValueError("method 须为 'kinetic' 或 'ghsv'")

    # 2) 内部独立算 v_out（出口操作体积流量）
    # [FIX-6] 要求反应参数完整，不允许用进口流量静默近似
    if (stoichiometric_coefficients is not None
            and key_component is not None
            and key_component_conversion is not None
            and temperature_K is not None
            and pressure_Pa is not None):
        v_out, _, _ = _compute_outlet_volumetric_flow(
            inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
            temperature_K, pressure_Pa)
    else:
        # [FIX-6] 不再静默用进口近似，明确报错
        raise ValueError(
            "calc_residence_time 需要提供完整反应参数（stoichiometric_coefficients、"
            "key_component、key_component_conversion、temperature_K、pressure_Pa）"
            "以计算出口体积流量。\n"
            "对于有显著摩尔数变化的反应（如气相分解 A→2B），进出口体积流量差异可达"
            "数倍，不可用进口流量近似。\n"
            "如确实不知出口组成，请先调用 calc_outlet_molar_flows 确定出口流率后"
            "再调用本函数。"
        )

    tau = V / v_out
    _vprint(f"\n【calc_residence_time】τ = {tau:.4f} s")
    return {"residence_time_s": tau}


if __name__ == "__main__":
    inlet_molar_flows = {
        'Cyclohexene': 100000 / 3600,
        'Acetic Acid': 130000 / 3600,
        'Cyclohexyl Acetate': 0.001,
    }
    stoichiometric_coefficients = {
        'Cyclohexene': -1, 'Acetic Acid': -1, 'Cyclohexyl Acetate': 1,
    }
    key_component = 'Cyclohexene'
    key_component_conversion = 0.6
    temperature_K = 473.15
    pressure_Pa = 300000.0
    arrhenius_pre_exponential_factor = 1.5e8
    activation_energy_J_mol = 75000.0

    calc_reaction_extent(inlet_molar_flows, stoichiometric_coefficients,
                          key_component, key_component_conversion)
    calc_outlet_molar_flows(inlet_molar_flows, stoichiometric_coefficients,
                             key_component, key_component_conversion)
    calc_outlet_total_molar_flow(inlet_molar_flows, stoichiometric_coefficients,
                                  key_component, key_component_conversion)
    calc_outlet_volumetric_flow(inlet_molar_flows, stoichiometric_coefficients,
                                 key_component, key_component_conversion,
                                 temperature_K, pressure_Pa)
    calc_outlet_concentrations(inlet_molar_flows, stoichiometric_coefficients,
                                key_component, key_component_conversion,
                                temperature_K, pressure_Pa)
    calc_fluid_density(inlet_molar_flows, temperature_K, pressure_Pa)
    calc_fluid_viscosity(inlet_molar_flows, temperature_K, pressure_Pa)
    calc_archimedes_number(0.0005, 2500.0, inlet_molar_flows, temperature_K, pressure_Pa)
    calc_minimum_fluidization_velocity(0.0005, 2500.0, inlet_molar_flows, temperature_K, pressure_Pa)
    calc_reaction_rate(inlet_molar_flows, stoichiometric_coefficients,
                        key_component, key_component_conversion,
                        temperature_K, pressure_Pa,
                        arrhenius_pre_exponential_factor, activation_energy_J_mol)
    calc_cstr_volume(method="kinetic",
                      inlet_molar_flows=inlet_molar_flows,
                      stoichiometric_coefficients=stoichiometric_coefficients,
                      key_component=key_component,
                      key_component_conversion=key_component_conversion,
                      temperature_K=temperature_K,
                      pressure_Pa=pressure_Pa,
                      arrhenius_pre_exponential_factor=arrhenius_pre_exponential_factor,
                      activation_energy_J_mol=activation_energy_J_mol)
    calc_residence_time(method="kinetic",
                         inlet_molar_flows=inlet_molar_flows,
                         stoichiometric_coefficients=stoichiometric_coefficients,
                         key_component=key_component,
                         key_component_conversion=key_component_conversion,
                         temperature_K=temperature_K,
                         pressure_Pa=pressure_Pa,
                         arrhenius_pre_exponential_factor=arrhenius_pre_exponential_factor,
                         activation_energy_J_mol=activation_energy_J_mol)