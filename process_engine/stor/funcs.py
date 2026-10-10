"""化学计量反应器独立工具函数（不依赖其他 public calc_* 函数的输出）。

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
"""


import os
import sys

_VERBOSE = os.getenv("FUNCS_VERBOSE", "0") == "1"


def _vprint(*args, **kwargs):
    if _VERBOSE:
        print(*args, **kwargs)

# 统一物性计算引擎（thermo + CoolProp 双引擎）
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from property.thermo_helper import (
    calc_flash_rachford_rice,
    calc_molar_volume,
    calc_liquid_density,
    get_mixture_properties,
    get_fluid_MW,
)


# ─────────────────────────────────────────────
# 模块级常量
# ─────────────────────────────────────────────
IDEAL_GAS_CONSTANT_J_per_mol_K = 8.314


# ─────────────────────────────────────────────
# 物性辅助函数（私有，可共享，基于 thermo_helper）
# ─────────────────────────────────────────────
def _calc_phase(components, zs, T, P):
    """用 thermo_helper 判断相态，返回 'l' 或 'g'。"""
    result = calc_flash_rachford_rice(components, zs, T, P)
    return 'l' if result.get("beta", 0.0) < 0.5 else 'g'


def _calc_liquid_molar_volume(components, mole_fractions, T, P):
    """液相混合物摩尔体积 (m³/mol)。"""
    # 多组分：Vml = (MW_mix/1000) / rho_L
    if len(components) > 1:
        props = get_mixture_properties(components, mole_fractions, T, P)
        rho_L = props.get("rho_L")
        MW_mix = props.get("MW_mix")
        if rho_L and MW_mix and rho_L > 0:
            return (MW_mix / 1000.0) / rho_L
    # 回退：按单组分加权摩尔体积
    vml = 0.0
    for c, z in zip(components, mole_fractions):
        single = calc_molar_volume(c, T, P, "liquid")
        if not single or single <= 0:
            rho = calc_liquid_density(c, T, P)
            mw = get_fluid_MW(c)
            if rho and mw and rho > 0:
                single = (mw / 1000.0) / rho
        if not single or single <= 0:
            raise ValueError(f"无法计算 '{c}' 的液相摩尔体积")
        vml += z * single
    return vml


# ─────────────────────────────────────────────
# 私有：核心计算（无 print，纯数学，可共享）
# ─────────────────────────────────────────────
def _compute_reaction_extent(inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion):
    inlet_molar_flow_key = inlet_molar_flows[key_component]
    stoich_coeff_key_abs = abs(stoichiometric_coefficients[key_component])
    return (inlet_molar_flow_key * key_component_conversion) / stoich_coeff_key_abs


def _compute_outlet_molar_flows(inlet_molar_flows, stoichiometric_coefficients,
                                 key_component, key_component_conversion):
    reaction_extent = _compute_reaction_extent(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    outlet = {}
    for component_name, inlet_flow in inlet_molar_flows.items():
        stoich_coeff = stoichiometric_coefficients.get(component_name, 0)
        outlet[component_name] = inlet_flow + stoich_coeff * reaction_extent
    return outlet


def _compute_outlet_volumetric_flow(inlet_molar_flows, stoichiometric_coefficients,
                                     key_component, key_component_conversion,
                                     reactor_temp_K, reactor_pressure_Pa):
    outlet_flows = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    # 过滤零流量组分，避免 thermo 在零分率上炸
    active = [(c, f) for c, f in outlet_flows.items() if f > 1e-12]
    if not active:
        raise ValueError("出口所有组分流量均为 0，无法计算体积流量")
    components = [c for c, _ in active]
    flows = [f for _, f in active]
    total = sum(flows)
    mole_fractions = [f / total for f in flows]
    # 统一用 thermo_helper 判相态（单/多组分均可）
    phase = _calc_phase(components, mole_fractions, reactor_temp_K, reactor_pressure_Pa)
    if phase == 'l':
        Vml = _calc_liquid_molar_volume(components, mole_fractions, reactor_temp_K, reactor_pressure_Pa)
        return total * Vml, 'l', Vml
    else:
        v = (total * IDEAL_GAS_CONSTANT_J_per_mol_K * reactor_temp_K) / reactor_pressure_Pa
        return v, 'g', None


# ═════════════════════════════════════════════
# 工具 1：反应进度
# ═════════════════════════════════════════════
def calc_reaction_extent(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
) -> dict:
    """计算化学计量反应器的反应进度 ξ（单位 mol/s）。

    反应进度是描述反应推进程度的基准量，所有组分的变化都遵循统一的计量关系：
    只需知道关键组分转化率，即可通过 ξ 推算所有组分的变化。
    ξ = (关键组分进口流率 × 转化率) / |关键组分计量系数|。
    反应进度越大，说明反应进行得越充分，产物生成越多，反应物消耗越多。

    参数:
        inlet_molar_flows           — {组分名: 进口摩尔流率} (mol/s)
                                       必须包含 stoichiometric_coefficients 中的所有组分，
                                       包括进料流量为 0 的纯产物（如 {"产物": 0.0}）
                                       必须包含 stoichiometric_coefficients 中的所有组分，
                                       包括进料流量为 0 的纯产物（如 {"产物": 0.0}）
        stoichiometric_coefficients — {组分名: 计量系数}，反应物负，产物正
        key_component               — 关键组分名称
        key_component_conversion    — 关键组分转化率 (0~1)

    返回:
        {"reaction_extent_mol_per_s": 反应进度 (mol/s)}
    """
    extent = _compute_reaction_extent(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    _vprint(f"【calc_reaction_extent】ξ = {extent:.2f} mol/s")
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
    """计算化学计量反应器各组分出口摩尔流率（单位 mol/s）。

    根据物料守恒和化学计量关系：出口流率 = 进口流率 + 计量系数 × 反应进度。
    反应物（计量系数为负）出口减少，产物（为正）出口增加，惰性组分（系数为0）
    不变。本工具独立计算，内部完成反应进度推算。

    参数:
        inlet_molar_flows           — {组分名: 进口摩尔流率} (mol/s)
                                       必须包含 stoichiometric_coefficients 中的所有组分，
                                       包括进料流量为 0 的纯产物（如 {"产物": 0.0}）
                                       必须包含 stoichiometric_coefficients 中的所有组分，
                                       包括进料流量为 0 的纯产物（如 {"产物": 0.0}）
        stoichiometric_coefficients — {组分名: 计量系数}
        key_component               — 关键组分名称
        key_component_conversion    — 关键组分转化率 (0~1)

    返回:
        {组分名: 出口摩尔流率 (mol/s)}
    """
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    _vprint("\n【calc_outlet_molar_flows】")
    for name, flow in outlet.items():
        _vprint(f"  F_{name} = {flow:.2f} mol/s")
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
    """计算化学计量反应器出口总摩尔流率（单位 mol/s）。

    由于反应前后分子数可能变化（如 A+B→C，2分子变1分子），出口总流率可能
    与进口不同。总流率是计算体积流量和浓度的基础。本工具独立计算。

    参数:
        inlet_molar_flows           — {组分名: 进口摩尔流率} (mol/s)
                                       必须包含 stoichiometric_coefficients 中的所有组分，
                                       包括进料流量为 0 的纯产物（如 {"产物": 0.0}）
        stoichiometric_coefficients — {组分名: 计量系数}
        key_component               — 关键组分名称
        key_component_conversion    — 关键组分转化率 (0~1)

    返回:
        {"outlet_total_molar_flow_mol_per_s": 出口总摩尔流率 (mol/s)}
    """
    outlet = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    total = sum(outlet.values())
    _vprint(f"\n【calc_outlet_total_molar_flow】F_total = {total:.2f} mol/s")
    return {"outlet_total_molar_flow_mol_per_s": total}


# ═════════════════════════════════════════════
# 工具 4：出口体积流量
# ═════════════════════════════════════════════
def calc_outlet_volumetric_flow(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
    reactor_temp_K: float,
    reactor_pressure_Pa: float,
) -> dict:
    """计算化学计量反应器出口体积流量（单位 m³/s）。

    根据出口组成、温度、压力自动判断相态：气相采用理想气体方程 PV=nRT，
    液相采用 thermo_helper 液体摩尔体积 Vml。本工具独立计算，内部完成出口流量
    与相态判定。

    参数:
        inlet_molar_flows           — {组分名: 进口摩尔流率} (mol/s)
                                       必须包含 stoichiometric_coefficients 中的所有组分，
                                       包括进料流量为 0 的纯产物（如 {"产物": 0.0}）
        stoichiometric_coefficients — {组分名: 计量系数}
        key_component               — 关键组分名称
        key_component_conversion    — 关键组分转化率 (0~1)
        reactor_temp_K              — 反应器操作温度 (K)
        reactor_pressure_Pa         — 反应器操作压力 (Pa)

    返回:
        {"outlet_volumetric_flow_m3_per_s": 出口体积流量 (m³/s)}
    """
    v, phase, Vml = _compute_outlet_volumetric_flow(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
        reactor_temp_K, reactor_pressure_Pa)
    _vprint("\n【calc_outlet_volumetric_flow】")
    if phase == 'l':
        _vprint(f"  相态判定：液相，Vml = {Vml:.6e} m^3/mol")
        _vprint(f"  v = {v:.6f} m^3/s")
    else:
        _vprint(f"  相态判定：气相，理想气体定律")
        _vprint(f"  v = {v:.4f} m^3/s")
    return {"outlet_volumetric_flow_m3_per_s": v}


# ═════════════════════════════════════════════
# 工具 5：出口浓度
# ═════════════════════════════════════════════
def calc_outlet_concentrations(
    inlet_molar_flows: dict,
    stoichiometric_coefficients: dict,
    key_component: str,
    key_component_conversion: float,
    reactor_temp_K: float,
    reactor_pressure_Pa: float,
) -> dict:
    """计算化学计量反应器出口各组分浓度（单位 mol/m³）。

    浓度 = 摩尔流率 / 体积流量，反映反应后各组分在混合物中的密集程度。
    本工具独立计算，内部完成出口流量与体积流量推算。

    参数:
        inlet_molar_flows           — {组分名: 进口摩尔流率} (mol/s)
                                       必须包含 stoichiometric_coefficients 中的所有组分，
                                       包括进料流量为 0 的纯产物（如 {"产物": 0.0}）
        stoichiometric_coefficients — {组分名: 计量系数}
        key_component               — 关键组分名称
        key_component_conversion    — 关键组分转化率 (0~1)
        reactor_temp_K              — 反应器操作温度 (K)
        reactor_pressure_Pa         — 反应器操作压力 (Pa)

    返回:
        {组分名: 出口浓度 (mol/m^3)}
    """
    outlet_flows = _compute_outlet_molar_flows(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion)
    v, _, _ = _compute_outlet_volumetric_flow(
        inlet_molar_flows, stoichiometric_coefficients, key_component, key_component_conversion,
        reactor_temp_K, reactor_pressure_Pa)
    concentrations = {name: flow / v for name, flow in outlet_flows.items()}
    _vprint("\n【calc_outlet_concentrations】")
    for name, c in concentrations.items():
        _vprint(f"  {name} = {c:.2f} mol/m^3")
    return concentrations


if __name__ == "__main__":
    inlet_molar_flows = {
        'Cyclohexene': 100000 / 3600,
        'Acetic Acid': 130000 / 3600,
        'Cyclohexyl Acetate': 0,
    }
    stoichiometric_coefficients = {
        'Cyclohexene': -1, 'Acetic Acid': -1, 'Cyclohexyl Acetate': 1,
    }
    key_component = 'Cyclohexene'
    key_component_conversion = 0.6
    reactor_temp_K = 473.15
    reactor_pressure_Pa = 300000.0

    calc_reaction_extent(inlet_molar_flows, stoichiometric_coefficients,
                         key_component, key_component_conversion)
    calc_outlet_molar_flows(inlet_molar_flows, stoichiometric_coefficients,
                            key_component, key_component_conversion)
    calc_outlet_total_molar_flow(inlet_molar_flows, stoichiometric_coefficients,
                                  key_component, key_component_conversion)
    calc_outlet_volumetric_flow(inlet_molar_flows, stoichiometric_coefficients,
                                 key_component, key_component_conversion,
                                 reactor_temp_K, reactor_pressure_Pa)
    calc_outlet_concentrations(inlet_molar_flows, stoichiometric_coefficients,
                                key_component, key_component_conversion,
                                reactor_temp_K, reactor_pressure_Pa)
