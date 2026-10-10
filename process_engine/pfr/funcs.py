"""平推流反应器（PFR）独立工具函数。

设计原则：
- 每个 public calc_* 只接收原始题目参数
- public calc_* 之间禁止互相调用
- 私有 _ 辅助函数和常量表可共享
- 每个函数只返回自己的最终结果

入参约定：
- 所有接收 inlet_molar_flows + stoichiometric_coefficients 的函数，inlet_molar_flows
  必须包含 stoichiometric_coefficients 中出现的所有组分（反应物 + 产物），
  即使某产物在进料中流量为 0 也必须显式写入（如 {"产物": 0.0}），
  否则该产物不会出现在出口流量字典中。

修复记录（化工角度）：
  [FIX-1] _compute_bed_volume_kinetic 积分被积函数：改用 _compute_volumetric_flow_at_TP
           自动判断气液相态，消除原先对液相 PFR 硬编码气相理想气体方程的根本性错误。
  [FIX-2] 积分被积函数中间步骤的出口流量负值保护：中间转化率 X 下若某非关键
           组分流量变负（超化学计量），截断为 0 而非传入物性计算崩溃；同时在积分
           上限处预校验不超过化学计量上限。
  [FIX-3] _compute_reaction_extent 校验关键组分计量系数必须 < 0（反应物），
           转化率须在 (0, 1] 之间，进料流量须 > 0。
  [FIX-4] 积分端点奇异性保护：quad 改用 limit 参数和 points 参数，并将积分下限
           从 X=0 替换为 X=ε（1e-9）以避免反应物初始浓度等于总进料时速率为 0
           导致的被积函数在 X→0 端点处的奇异。
  [FIX-5] _compute_outlet_molar_flows 增加出口流量非负校验，与 CSTR 版本保持一致。
  [FIX-6] _compute_volumetric_flow_at_TP 液相后备路径逐组分 Vml 非 None 保护。
  [FIX-7] calc_bed_length_single_tube 消除 A → D → A_section 的冗余往返转换，
           直接复用 A，消除浮点误差积累。
  [FIX-8] GHSV 参考态参数化：新增 ghsv_reference_T_K / ghsv_reference_P_Pa，
           默认保持原 STP（273.15K），支持工业 NTP（293.15K）或任意参考态。
  [FIX-9] 截面积设计基准改为进出口最大流量：calc_single_tube_diameter、
           calc_tube_number、calc_bed_length_multitube、calc_bed_length_single_tube
           均从进口和出口流量中取较大值作截面设计依据，确保出口段气速不超限。
"""

import os

_VERBOSE = os.getenv("FUNCS_VERBOSE", "0") == "1"


def _vprint(*args, **kwargs):
    if _VERBOSE:
        print(*args, **kwargs)


import math

from scipy.integrate import quad
from thermo import Chemical, Mixture
from thermo import ChemicalConstantsPackage
from thermo.flash import FlashVL
from thermo.phases import CEOSGas, CEOSLiquid
from thermo.eos_mix import PRMIX


# ─────────────────────────────────────────────
# 模块级常量
# ─────────────────────────────────────────────
IDEAL_GAS_CONSTANT_J_per_mol_K = 8.314
STP_TEMPERATURE_K = 273.15      # IUPAC STP
STP_PRESSURE_Pa   = 101_325
# 注：工业 GHSV 有时基于 NTP（293.15K, 101325Pa），
# 通过 ghsv_reference_T_K 参数指定，默认保持原 STP。

# 出口流量负值容差（mol/s）：小于此绝对值视为数值误差，截断为 0
_FLOW_NEG_TOL = -1e-9

# 积分下限偏移：避免 X=0 时某些反应被积函数奇异
_X_EPS = 1e-9


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
    # 单组分时 thermo PR flash 因 incipient_guesses N-1=0 触发 ZeroDivisionError
    if len(components) == 1:
        return Chemical(components[0], T=T, P=P).phase
    flasher = _make_flasher(components)
    result = flasher.flash(T=T, P=P, zs=zs)
    return 'l' if result.VF < 0.5 else 'g'


# ─────────────────────────────────────────────
# 私有：核心计算（无 print，纯数学，可共享）
# ─────────────────────────────────────────────
def _compute_reaction_extent(inlet_molar_flows, stoichiometric_coefficients,
                              key_component, key_component_conversion):
    """计算反应进度 ξ (mol/s)。

    [FIX-3] 校验：
      - key_component 在 inlet_molar_flows 和 stoichiometric_coefficients 中均存在
      - 关键组分计量系数 < 0（必须是反应物）
      - 转化率在 (0, 1] 内
      - 关键组分进料流率 > 0
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
    # [FIX-3] 关键组分必须是反应物
    if coeff_key >= 0:
        raise ValueError(
            f"关键组分 '{key_component}' 的计量系数为 {coeff_key}（≥ 0），"
            f"关键组分必须是反应物（计量系数 < 0），请勿将产物指定为关键组分。"
        )
    if not (0.0 < key_component_conversion <= 1.0):
        raise ValueError(
            f"key_component_conversion 须在 (0, 1] 之间，得到 {key_component_conversion}"
        )
    inlet_key = inlet_molar_flows[key_component]
    if inlet_key <= 0:
        raise ValueError(
            f"关键组分 '{key_component}' 的进口摩尔流率须 > 0，得到 {inlet_key} mol/s"
        )
    coeff_key_abs = abs(coeff_key)
    return (inlet_key * key_component_conversion) / coeff_key_abs


def _compute_outlet_molar_flows(inlet_molar_flows, stoichiometric_coefficients,
                                 key_component, key_component_conversion):
    """计算出口各组分摩尔流率 (mol/s)。

    [FIX-5] 出口流量非负校验：超化学计量上限时抛出 ValueError；
    极小负值（数值误差）截断为 0。
    """
    extent = _compute_reaction_extent(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion
    )
    outlet = {
        name: flow + stoichiometric_coefficients.get(name, 0) * extent
        for name, flow in inlet_molar_flows.items()
    }
    for name, flow in outlet.items():
        if flow < _FLOW_NEG_TOL:
            raise ValueError(
                f"组分 '{name}' 出口摩尔流率为 {flow:.4e} mol/s（< 0），"
                f"转化率 {key_component_conversion:.4f} 超过该组分的化学计量上限。"
                f"请降低转化率或检查计量系数。"
            )
        if _FLOW_NEG_TOL <= flow < 0:
            outlet[name] = 0.0  # 数值误差截断
    return outlet


def _clamp_intermediate_flows(intermediate_flows):
    """将中间转化率步骤中因数值误差产生的极小负流量截断为 0。

    [FIX-2] 积分被积函数内部调用，不抛出异常（积分过程中间点允许截断），
    仅截断 [_FLOW_NEG_TOL, 0) 范围内的值；超过 _FLOW_NEG_TOL 的负值说明
    该转化率在化学计量上不合法，返回 False 通知调用方。
    """
    for name, flow in intermediate_flows.items():
        if flow < _FLOW_NEG_TOL:
            return False  # 化学计量不合法
        if _FLOW_NEG_TOL <= flow < 0:
            intermediate_flows[name] = 0.0
    return True


def _compute_volumetric_flow_at_TP(molar_flows, T, P):
    """通用：根据组成 + T/P 算体积流量。返回 (v, phase, Vml)。

    [FIX-6] 液相后备路径逐组分 Vml 非 None / 非正 保护。
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
            # [FIX-6] 后备：逐组分加权，校验每个 Vml 有效
            Vml_parts = []
            for c, z in zip(components, mole_fractions):
                single_mix = Mixture([c], zs=[1], T=T, P=P)
                vml_single = single_mix.Vml
                if vml_single is None or vml_single <= 0:
                    raise ValueError(
                        f"组分 '{c}' 在 T={T:.2f}K, P={P:.0f}Pa 下液相摩尔体积 "
                        f"Vml={vml_single}（无效）。请检查操作条件是否超出数据范围，"
                        f"或该组分在此条件下是否确为液相。"
                    )
                Vml_parts.append(z * vml_single)
            Vml = sum(Vml_parts)
        return total * Vml, 'l', Vml

    # 气相：理想气体
    v = (total * IDEAL_GAS_CONSTANT_J_per_mol_K * T) / P
    return v, 'g', None


def _compute_inlet_vol_flow_at_ref(inlet_molar_flows, ref_T_K, ref_P_Pa):
    """在参考状态 (ref_T_K, ref_P_Pa) 下计算进口体积流量。

    [FIX-8] 参数化参考态，替代原固定 STP，支持 NTP 或任意参考态。
    """
    v, _, _ = _compute_volumetric_flow_at_TP(inlet_molar_flows, ref_T_K, ref_P_Pa)
    return v


def _compute_inlet_vol_flow_stp(inlet_molar_flows):
    """标准状态（STP）进口体积流量，保留为内部向后兼容入口。"""
    return _compute_inlet_vol_flow_at_ref(inlet_molar_flows, STP_TEMPERATURE_K, STP_PRESSURE_Pa)


def _compute_inlet_vol_flow_operating(inlet_molar_flows, temperature_K, pressure_Pa):
    v, _, _ = _compute_volumetric_flow_at_TP(inlet_molar_flows, temperature_K, pressure_Pa)
    return v


def _compute_max_vol_flow_operating(inlet_molar_flows, stoichiometric_coefficients,
                                    key_component, key_component_conversion,
                                    temperature_K, pressure_Pa):
    """[FIX-9] 计算进出口体积流量的较大值，用于截面积设计。

    对气相摩尔数增大的反应（如 A→2B），出口流量远大于进口，
    截面设计应取出口流量以确保整个床层气速不超限。
    """
    v_in, _, _ = _compute_volumetric_flow_at_TP(inlet_molar_flows, temperature_K, pressure_Pa)
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion
    )
    v_out, _, _ = _compute_volumetric_flow_at_TP(outlet, temperature_K, pressure_Pa)
    return max(v_in, v_out), v_in, v_out


def _compute_bed_volume_ghsv(inlet_molar_flows, space_velocity_ghsv, ghsv_time_unit="s",
                              ref_T_K=STP_TEMPERATURE_K, ref_P_Pa=STP_PRESSURE_Pa):
    """[FIX-8] 接受参考态参数。"""
    inlet_v_ref = _compute_inlet_vol_flow_at_ref(inlet_molar_flows, ref_T_K, ref_P_Pa)
    if ghsv_time_unit == "h":
        ghsv_per_s = space_velocity_ghsv / 3600
    elif ghsv_time_unit == "s":
        ghsv_per_s = space_velocity_ghsv
    else:
        raise ValueError("ghsv_time_unit 须为 's' 或 'h'")
    return inlet_v_ref / ghsv_per_s


def _compute_bed_volume_kinetic(inlet_molar_flows, stoichiometric_coefficients,
                                 key_component, key_component_conversion,
                                 temperature_K, pressure_Pa,
                                 arrhenius_pre_exponential_factor, activation_energy_J_mol,
                                 reaction_order):
    """PFR 设计方程数值积分 V = ∫FA0/(-rA) dX。

    [FIX-1] 积分被积函数中体积流量改用 _compute_volumetric_flow_at_TP，
            自动识别液相/气相，消除对液相 PFR 的硬编码气相假设。

    [FIX-2] 积分中间步骤出口流量负值保护：
            - 积分上限前预校验最终出口流量（调用 _compute_outlet_molar_flows）
            - 被积函数内截断极小负流量，化学计量不合法时返回极大值（等价于跳过）

    [FIX-4] 积分下限从 X=0 改为 X=_X_EPS，避免进料中产物为 0 时
            X=0 处某些速率方程（如 CA=0 且 n<1）产生奇异。
    """
    # [FIX-2] 预校验：最终转化率处所有组分出口流量非负
    # （_compute_outlet_molar_flows 内部已有 ValueError，此处自动触发）
    _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion
    )

    rate_constant = arrhenius_pre_exponential_factor * math.exp(
        -activation_energy_J_mol / (IDEAL_GAS_CONSTANT_J_per_mol_K * temperature_K))
    FA0 = inlet_molar_flows[key_component]
    coeff_key_abs = abs(stoichiometric_coefficients[key_component])

    def integrand(xA):
        extent = (FA0 * xA) / coeff_key_abs
        # 构造中间步骤出口流量
        outlet = {
            comp: inlet_molar_flows[comp] + stoichiometric_coefficients.get(comp, 0) * extent
            for comp in inlet_molar_flows
        }
        # [FIX-2] 截断中间步骤极小负流量
        if not _clamp_intermediate_flows(outlet):
            # 化学计量不合法（理论上已被预校验拦截，此为防御性保护）
            return float('inf')

        # [FIX-1] 自动判断相态，液相用 Vml，气相用理想气体
        try:
            v_flow, _, _ = _compute_volumetric_flow_at_TP(outlet, temperature_K, pressure_Pa)
        except Exception:
            return float('inf')

        if v_flow <= 0:
            return float('inf')

        c_key = outlet[key_component] / v_flow
        if c_key < 0:
            c_key = 0.0

        # 防止浓度为 0 时 n 级速率奇异（如 n<1 时 0^n 产生 NaN）
        if c_key == 0.0:
            return float('inf')

        rate = rate_constant * c_key**reaction_order
        if rate <= 0:
            return float('inf')

        return FA0 / rate

    # [FIX-4] 积分下限从 0 改为 _X_EPS，避免端点奇异
    x_lo = _X_EPS
    x_hi = key_component_conversion

    if x_hi <= x_lo:
        return 0.0

    bed_volume, _err = quad(
        integrand,
        x_lo, x_hi,
        limit=200,          # 增加自适应分段上限
        epsabs=1e-6,
        epsrel=1e-6,
    )
    return bed_volume


def _compute_bed_volume(method, inlet_molar_flows, stoichiometric_coefficients,
                        key_component, key_component_conversion,
                        temperature_K, pressure_Pa,
                        arrhenius_pre_exponential_factor, activation_energy_J_mol,
                        reaction_order, space_velocity_ghsv, ghsv_time_unit,
                        ghsv_reference_T_K=STP_TEMPERATURE_K,
                        ghsv_reference_P_Pa=STP_PRESSURE_Pa):
    if method == "ghsv":
        return _compute_bed_volume_ghsv(
            inlet_molar_flows, space_velocity_ghsv, ghsv_time_unit,
            ref_T_K=ghsv_reference_T_K, ref_P_Pa=ghsv_reference_P_Pa
        )
    if method == "kinetic":
        return _compute_bed_volume_kinetic(
            inlet_molar_flows, stoichiometric_coefficients,
            key_component, key_component_conversion,
            temperature_K, pressure_Pa,
            arrhenius_pre_exponential_factor, activation_energy_J_mol, reaction_order)
    raise ValueError("method 须为 'ghsv' 或 'kinetic'")


# ═════════════════════════════════════════════
# 工具 1：反应进度
# ═════════════════════════════════════════════
def calc_reaction_extent(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
) -> dict:
    """计算 PFR 平推流反应器的反应进度 ξ（单位 mol/s）。

    ξ = (关键组分进口流率 × 转化率) / |关键组分计量系数|。

    参数:
        inlet_molar_flows           — {组分名: 进口摩尔流率} (mol/s)
        stoichiometric_coefficients — {组分名: 计量系数}，关键组分系数须 < 0（反应物）
        key_component               — 关键组分名称（必须是反应物）
        key_component_conversion    — 关键组分转化率 (0, 1]

    返回:
        {"reaction_extent_mol_per_s": ξ (mol/s)}
    """
    extent = _compute_reaction_extent(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    _vprint(f"【calc_reaction_extent】ξ = {extent:.4f} mol/s")
    return {"reaction_extent_mol_per_s": extent}


# ═════════════════════════════════════════════
# 工具 2：进口总摩尔流率
# ═════════════════════════════════════════════
def calc_inlet_total_molar_flow(inlet_molar_flows: dict) -> dict:
    """计算 PFR 进口总摩尔流率（单位 mol/s）。

    返回:
        {"inlet_total_molar_flow_mol_per_s": F_total (mol/s)}
    """
    total = sum(inlet_molar_flows.values())
    _vprint(f"\n【calc_inlet_total_molar_flow】F_total = {total:.4f} mol/s")
    return {"inlet_total_molar_flow_mol_per_s": total}


# ═════════════════════════════════════════════
# 工具 3：进口体积流量（标准状态）
# ═════════════════════════════════════════════
def calc_inlet_volumetric_flow_stp(
    inlet_molar_flows: dict,
    reference_temperature_K: float = STP_TEMPERATURE_K,
    reference_pressure_Pa: float = STP_PRESSURE_Pa,
) -> dict:
    """计算 PFR 参考状态下进口体积流量（单位 m³/s）。

    [FIX-8] 新增 reference_temperature_K / reference_pressure_Pa 参数，
    默认 STP（273.15K, 101325Pa）；工业 NTP 传 293.15K。
    该流量用于 GHSV 法计算催化剂床层体积。

    返回:
        {"inlet_volumetric_flow_stp_m3_per_s": v (m^3/s)}
    """
    v = _compute_inlet_vol_flow_at_ref(inlet_molar_flows, reference_temperature_K, reference_pressure_Pa)
    _vprint(
        f"\n【calc_inlet_volumetric_flow_stp】"
        f"参考态 T={reference_temperature_K}K, P={reference_pressure_Pa}Pa, "
        f"v = {v:.6f} m^3/s"
    )
    return {"inlet_volumetric_flow_stp_m3_per_s": v}


# ═════════════════════════════════════════════
# 工具 4：进口体积流量（工况）
# ═════════════════════════════════════════════
def calc_inlet_volumetric_flow_operating(
    inlet_molar_flows: dict,
    temperature_K: float,
    pressure_Pa: float,
) -> dict:
    """计算 PFR 实际工况下进口体积流量（单位 m³/s）。

    返回:
        {"inlet_volumetric_flow_operating_m3_per_s": v (m^3/s)}
    """
    v = _compute_inlet_vol_flow_operating(inlet_molar_flows, temperature_K, pressure_Pa)
    _vprint(f"\n【calc_inlet_volumetric_flow_operating】v = {v:.6f} m^3/s")
    return {"inlet_volumetric_flow_operating_m3_per_s": v}


# ═════════════════════════════════════════════
# 工具 5：出口各组分摩尔流率
# ═════════════════════════════════════════════
def calc_outlet_molar_flows(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
) -> dict:
    """计算 PFR 各组分出口摩尔流率（单位 mol/s）。

    出口流率 = 进口流率 + 计量系数 × 反应进度。

    返回:
        {组分名: 出口摩尔流率 (mol/s)}
    """
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    _vprint("\n【calc_outlet_molar_flows】")
    for comp, flow in outlet.items():
        _vprint(f"  F_{comp} = {flow:.4f} mol/s")
    return outlet


# ═════════════════════════════════════════════
# 工具 6：出口总摩尔流率
# ═════════════════════════════════════════════
def calc_outlet_total_molar_flow(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
) -> dict:
    """计算 PFR 出口总摩尔流率（单位 mol/s）。

    返回:
        {"outlet_total_molar_flow_mol_per_s": F_total (mol/s)}
    """
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    total = sum(outlet.values())
    _vprint(f"\n【calc_outlet_total_molar_flow】F_total = {total:.4f} mol/s")
    return {"outlet_total_molar_flow_mol_per_s": total}


# ═════════════════════════════════════════════
# 工具 7：出口体积流量
# ═════════════════════════════════════════════
def calc_outlet_volumetric_flow(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
    temperature_K: float,
    pressure_Pa: float,
) -> dict:
    """计算 PFR 出口体积流量（单位 m³/s）。

    气相用理想气体方程 PV=nRT，液相用 thermo 库液体摩尔体积 Vml。

    返回:
        {"outlet_volumetric_flow_m3_per_s": v (m³/s)}
    """
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    v, phase, Vml = _compute_volumetric_flow_at_TP(outlet, temperature_K, pressure_Pa)
    _vprint("\n【calc_outlet_volumetric_flow】")
    if phase == 'l':
        _vprint(f"  相态判定：液相，Vml = {Vml:.6e} m^3/mol")
    else:
        _vprint(f"  相态判定：气相，理想气体定律")
    _vprint(f"  v = {v:.6f} m^3/s")
    return {"outlet_volumetric_flow_m3_per_s": v}


# ═════════════════════════════════════════════
# 工具 8：出口浓度
# ═════════════════════════════════════════════
def calc_outlet_concentrations(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
    temperature_K: float,
    pressure_Pa: float,
) -> dict:
    """计算 PFR 出口各组分浓度（单位 mol/m³）。

    返回:
        {组分名: 出口浓度 (mol/m³)}
    """
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    v, _, _ = _compute_volumetric_flow_at_TP(outlet, temperature_K, pressure_Pa)
    concentrations = {name: flow / v for name, flow in outlet.items()}
    _vprint("\n【calc_outlet_concentrations】")
    for k, c in concentrations.items():
        _vprint(f"  {k} = {c:.4f} mol/m^3")
    return concentrations


# ═════════════════════════════════════════════
# 工具 9：床层体积
# ═════════════════════════════════════════════
def calc_bed_volume(
    method: str = "ghsv",
    inlet_molar_flows: dict = None,
    space_velocity_ghsv: float = None,
    ghsv_time_unit: str = "s",
    ghsv_reference_T_K: float = STP_TEMPERATURE_K,
    ghsv_reference_P_Pa: float = STP_PRESSURE_Pa,
    stoichiometric_coefficients: dict = None,
    key_component: str = None,
    key_component_conversion: float = None,
    temperature_K: float = None,
    pressure_Pa: float = None,
    arrhenius_pre_exponential_factor: float = None,
    activation_energy_J_mol: float = None,
    reaction_order: int = 1,
) -> dict:
    """计算 PFR 催化剂床层体积（单位 m³）。

    支持两种方法：
    - method="ghsv"：V = v_ref/GHSV，v_ref 在参考态下计算。
      [FIX-8] 新增 ghsv_reference_T_K（默认 STP=273.15K）和
      ghsv_reference_P_Pa（默认 101325Pa）参数，支持工业 NTP（293.15K）。
    - method="kinetic"：积分 PFR 设计方程 V=∫(FA0/(-rA))dX。
      [FIX-1] 被积函数自动识别液气相态。
      [FIX-4] 积分下限从 X=0 改为 X=ε 避免端点奇异。

    返回:
        {"bed_volume_m3": V (m³)}
    """
    V = _compute_bed_volume(
        method, inlet_molar_flows, stoichiometric_coefficients,
        key_component, key_component_conversion, temperature_K, pressure_Pa,
        arrhenius_pre_exponential_factor, activation_energy_J_mol,
        reaction_order, space_velocity_ghsv, ghsv_time_unit,
        ghsv_reference_T_K=ghsv_reference_T_K,
        ghsv_reference_P_Pa=ghsv_reference_P_Pa,
    )
    _vprint(
        f"\n【calc_bed_volume】({method}) "
        f"{'参考态 T=' + str(ghsv_reference_T_K) + 'K' if method == 'ghsv' else ''}"
        f" V = {V:.4f} m^3"
    )
    return {"bed_volume_m3": V}


# ═════════════════════════════════════════════
# 工具 10：催化剂质量
# ═════════════════════════════════════════════
def calc_catalyst_mass(
    particle_density_kg_per_m3: float,
    bed_void_fraction: float,
    method: str = "ghsv",
    inlet_molar_flows: dict = None,
    space_velocity_ghsv: float = None,
    ghsv_time_unit: str = "s",
    ghsv_reference_T_K: float = STP_TEMPERATURE_K,
    ghsv_reference_P_Pa: float = STP_PRESSURE_Pa,
    stoichiometric_coefficients: dict = None,
    key_component: str = None,
    key_component_conversion: float = None,
    temperature_K: float = None,
    pressure_Pa: float = None,
    arrhenius_pre_exponential_factor: float = None,
    activation_energy_J_mol: float = None,
    reaction_order: int = 1,
) -> dict:
    """计算 PFR 催化剂装填质量（单位 kg）。

    质量 = 床层体积 × 颗粒密度 × (1−空隙率)。

    参数：
        particle_density_kg_per_m3 — 催化剂颗粒密度 (kg/m³)，须 > 0
        bed_void_fraction          — 床层空隙率 (0, 1)
        ghsv_reference_T_K         — GHSV 参考态温度，默认 STP 273.15K
        ghsv_reference_P_Pa        — GHSV 参考态压力，默认 101325 Pa

    返回:
        {"catalyst_mass_kg": W (kg)}
    """
    if particle_density_kg_per_m3 <= 0:
        raise ValueError(f"颗粒密度须 > 0，得到 {particle_density_kg_per_m3} kg/m³")
    if not (0.0 < bed_void_fraction < 1.0):
        raise ValueError(f"床层空隙率须在 (0, 1) 之间，得到 {bed_void_fraction}")

    V = _compute_bed_volume(
        method, inlet_molar_flows, stoichiometric_coefficients,
        key_component, key_component_conversion, temperature_K, pressure_Pa,
        arrhenius_pre_exponential_factor, activation_energy_J_mol,
        reaction_order, space_velocity_ghsv, ghsv_time_unit,
        ghsv_reference_T_K=ghsv_reference_T_K,
        ghsv_reference_P_Pa=ghsv_reference_P_Pa,
    )
    W = V * particle_density_kg_per_m3 * (1.0 - bed_void_fraction)
    _vprint(f"\n【calc_catalyst_mass】W = {W:.2f} kg")
    return {"catalyst_mass_kg": W}


# ═════════════════════════════════════════════
# 工具 11：单筒内径
# ═════════════════════════════════════════════
def calc_single_tube_diameter(
    inlet_molar_flows: dict,
    temperature_K: float,
    pressure_Pa: float,
    superficial_velocity_m_per_s: float = 1.0,
    stoichiometric_coefficients: dict = None,
    key_component: str = None,
    key_component_conversion: float = None,
) -> dict:
    """计算单筒 PFR 反应器内径（单位 m）。

    由 Q=u×A 反算截面积再求内径，Q 取进出口工况体积流量的较大值（通常为出口）。

    [FIX-9] 原版仅用进口流量。对气相摩尔数增大的反应（如 A→2B），出口流量
    远大于进口，床层出口段气速会超过设计值并产生额外压降，设计应取最大流量。
    若未提供反应参数，则退回到进口流量（与原版行为一致）。

    返回:
        {"single_tube_inner_diameter_m": D (m),
         "single_tube_cross_section_area_m2": A (m^2),
         "design_basis": "inlet" 或 "outlet"（截面设计依据）}
    """
    if (stoichiometric_coefficients is not None and key_component is not None
            and key_component_conversion is not None):
        v_design, v_in, v_out = _compute_max_vol_flow_operating(
            inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
            temperature_K, pressure_Pa
        )
        basis = "outlet" if v_out >= v_in else "inlet"
    else:
        v_design = _compute_inlet_vol_flow_operating(inlet_molar_flows, temperature_K, pressure_Pa)
        basis = "inlet"

    A = v_design / superficial_velocity_m_per_s
    D = math.sqrt(4.0 * A / math.pi)
    _vprint(f"\n【calc_single_tube_diameter】D = {D:.4f} m, A = {A:.4f} m^2, basis={basis}")
    return {
        "single_tube_inner_diameter_m": D,
        "single_tube_cross_section_area_m2": A,
        "design_basis": basis,
    }


# ═════════════════════════════════════════════
# 工具 12：管数
# ═════════════════════════════════════════════
def calc_tube_number(
    inlet_molar_flows: dict,
    temperature_K: float,
    pressure_Pa: float,
    tube_inner_diameter_m: float,
    superficial_velocity_m_per_s: float = 1.0,
    stoichiometric_coefficients: dict = None,
    key_component: str = None,
    key_component_conversion: float = None,
) -> dict:
    """计算列管式 PFR 反应器管数（无量纲，向上取整）。

    [FIX-9] 管数计算基于进出口最大体积流量而非单纯进口流量，
    确保出口段每根管的气速不超过设计值。

    返回:
        {"tube_count": N (int),
         "single_tube_volumetric_flow_m3_per_s": q (m^3/s),
         "design_basis": "inlet" 或 "outlet"}
    """
    if (stoichiometric_coefficients is not None and key_component is not None
            and key_component_conversion is not None):
        v_design, v_in, v_out = _compute_max_vol_flow_operating(
            inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
            temperature_K, pressure_Pa
        )
        basis = "outlet" if v_out >= v_in else "inlet"
    else:
        v_design = _compute_inlet_vol_flow_operating(inlet_molar_flows, temperature_K, pressure_Pa)
        basis = "inlet"

    single_tube_v = (math.pi / 4.0) * tube_inner_diameter_m**2 * superficial_velocity_m_per_s
    tube_count = math.ceil(v_design / single_tube_v)
    _vprint(f"\n【calc_tube_number】N = {tube_count}, q = {single_tube_v:.6f} m^3/s, basis={basis}")
    return {
        "tube_count": tube_count,
        "single_tube_volumetric_flow_m3_per_s": single_tube_v,
        "design_basis": basis,
    }


# ═════════════════════════════════════════════
# 工具 13：列管式床层长度
# ═════════════════════════════════════════════
def calc_bed_length_multitube(
    inlet_molar_flows: dict,
    temperature_K: float,
    pressure_Pa: float,
    tube_inner_diameter_m: float,
    superficial_velocity_m_per_s: float = 1.0,
    method: str = "ghsv",
    space_velocity_ghsv: float = None,
    ghsv_time_unit: str = "s",
    ghsv_reference_T_K: float = STP_TEMPERATURE_K,
    ghsv_reference_P_Pa: float = STP_PRESSURE_Pa,
    stoichiometric_coefficients: dict = None,
    key_component: str = None,
    key_component_conversion: float = None,
    arrhenius_pre_exponential_factor: float = None,
    activation_energy_J_mol: float = None,
    reaction_order: int = 1,
) -> dict:
    """计算列管式 PFR 床层长度（单位 m）。

    床层长度 = 床层体积 / (管数 × 单管截面积)。

    [FIX-9] 管数计算基于进出口最大体积流量（需提供反应参数）。
    [FIX-8] GHSV 参考态可通过 ghsv_reference_T_K / ghsv_reference_P_Pa 指定。

    返回:
        {"bed_length_m": L (m),
         "total_cross_section_area_m2": A_total (m²),
         "tube_count": N}
    """
    # 床层体积
    V = _compute_bed_volume(
        method, inlet_molar_flows, stoichiometric_coefficients,
        key_component, key_component_conversion, temperature_K, pressure_Pa,
        arrhenius_pre_exponential_factor, activation_energy_J_mol,
        reaction_order, space_velocity_ghsv, ghsv_time_unit,
        ghsv_reference_T_K=ghsv_reference_T_K,
        ghsv_reference_P_Pa=ghsv_reference_P_Pa,
    )

    # [FIX-9] 管数：基于进出口最大流量
    if (stoichiometric_coefficients is not None and key_component is not None
            and key_component_conversion is not None):
        v_design, v_in, v_out = _compute_max_vol_flow_operating(
            inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
            temperature_K, pressure_Pa
        )
    else:
        v_design = _compute_inlet_vol_flow_operating(inlet_molar_flows, temperature_K, pressure_Pa)

    single_tube_v = (math.pi / 4.0) * tube_inner_diameter_m**2 * superficial_velocity_m_per_s
    tube_count = math.ceil(v_design / single_tube_v)

    A_total = tube_count * math.pi * tube_inner_diameter_m**2 / 4.0
    L = V / A_total
    _vprint(
        f"\n【calc_bed_length_multitube】tube_count={tube_count}, "
        f"L = {L:.4f} m, A_total = {A_total:.4f} m^2"
    )
    return {
        "bed_length_m": L,
        "total_cross_section_area_m2": A_total,
        "tube_count": tube_count,
    }


# ═════════════════════════════════════════════
# 工具 14：单筒床层长度
# ═════════════════════════════════════════════
def calc_bed_length_single_tube(
    inlet_molar_flows: dict,
    temperature_K: float,
    pressure_Pa: float,
    superficial_velocity_m_per_s: float = 1.0,
    method: str = "ghsv",
    space_velocity_ghsv: float = None,
    ghsv_time_unit: str = "s",
    ghsv_reference_T_K: float = STP_TEMPERATURE_K,
    ghsv_reference_P_Pa: float = STP_PRESSURE_Pa,
    stoichiometric_coefficients: dict = None,
    key_component: str = None,
    key_component_conversion: float = None,
    arrhenius_pre_exponential_factor: float = None,
    activation_energy_J_mol: float = None,
    reaction_order: int = 1,
) -> dict:
    """计算单筒 PFR 床层长度（单位 m）。

    床层长度 = 床层体积 / 截面积。

    [FIX-7] 消除 A→D→A_section 的冗余往返转换（原版多一次 sqrt+平方，引入
    浮点误差），直接复用截面积 A。
    [FIX-9] 截面积基于进出口最大体积流量。
    [FIX-8] 支持 GHSV 参考态参数。

    返回:
        {"bed_length_m": L (m),
         "total_cross_section_area_m2": A (m²)}
    """
    V = _compute_bed_volume(
        method, inlet_molar_flows, stoichiometric_coefficients,
        key_component, key_component_conversion, temperature_K, pressure_Pa,
        arrhenius_pre_exponential_factor, activation_energy_J_mol,
        reaction_order, space_velocity_ghsv, ghsv_time_unit,
        ghsv_reference_T_K=ghsv_reference_T_K,
        ghsv_reference_P_Pa=ghsv_reference_P_Pa,
    )

    # [FIX-9] 截面积基于进出口最大流量
    if (stoichiometric_coefficients is not None and key_component is not None
            and key_component_conversion is not None):
        v_design, v_in, v_out = _compute_max_vol_flow_operating(
            inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
            temperature_K, pressure_Pa
        )
    else:
        v_design = _compute_inlet_vol_flow_operating(inlet_molar_flows, temperature_K, pressure_Pa)

    # [FIX-7] 直接用 A，不再经过 D 反算
    A = v_design / superficial_velocity_m_per_s
    L = V / A
    _vprint(f"\n【calc_bed_length_single_tube】L = {L:.4f} m, A = {A:.4f} m^2")
    return {"bed_length_m": L, "total_cross_section_area_m2": A}


if __name__ == "__main__":
    inlet_molar_flows = {
        'nitrogen': 100_000 / 3600,
        'hydrogen': 130_000 / 3600,
        'ammonia': 0.0,
        'water': 0.0,
    }
    stoichiometric_coefficients = {'nitrogen': -4, 'hydrogen': -5, 'ammonia': 4, 'water': 6}
    key_component = 'nitrogen'
    key_component_conversion = 0.955
    temperature_K = 973.15
    pressure_Pa = 300_000
    space_velocity_ghsv = 4000

    calc_reaction_extent(inlet_molar_flows, stoichiometric_coefficients,
                         key_component, key_component_conversion)
    calc_inlet_total_molar_flow(inlet_molar_flows)
    calc_inlet_volumetric_flow_stp(inlet_molar_flows)
    calc_inlet_volumetric_flow_operating(inlet_molar_flows, temperature_K, pressure_Pa)
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
    calc_bed_volume(method="ghsv", inlet_molar_flows=inlet_molar_flows,
                    space_velocity_ghsv=space_velocity_ghsv, ghsv_time_unit="h")
    calc_catalyst_mass(particle_density_kg_per_m3=3200, bed_void_fraction=0.4,
                       method="ghsv", inlet_molar_flows=inlet_molar_flows,
                       space_velocity_ghsv=space_velocity_ghsv, ghsv_time_unit="h")
    calc_single_tube_diameter(
        inlet_molar_flows, temperature_K, pressure_Pa,
        superficial_velocity_m_per_s=1.0,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    calc_tube_number(
        inlet_molar_flows, temperature_K, pressure_Pa,
        tube_inner_diameter_m=0.05, superficial_velocity_m_per_s=1.0,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    calc_bed_length_multitube(
        inlet_molar_flows, temperature_K, pressure_Pa,
        tube_inner_diameter_m=0.05, superficial_velocity_m_per_s=1.0,
        method="ghsv", space_velocity_ghsv=space_velocity_ghsv, ghsv_time_unit="h",
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    calc_bed_length_single_tube(
        inlet_molar_flows, temperature_K, pressure_Pa,
        superficial_velocity_m_per_s=1.0,
        method="ghsv", space_velocity_ghsv=space_velocity_ghsv, ghsv_time_unit="h",
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )