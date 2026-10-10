"""二元/多组分精馏塔独立工具函数。

设计原则：
- 每个 public calc_* 只接收原始题目参数
- public calc_* 之间禁止互相调用
- 私有 _ 辅助函数和常量表可共享
- 每个函数只返回自己的最终结果

热力学模型：thermo_helper 双引擎（thermo + CoolProp）。全程 SI 单位。
"""


import os

_VERBOSE = os.getenv("FUNCS_VERBOSE", "0") == "1"


def _vprint(*args, **kwargs):
    if _VERBOSE:
        print(*args, **kwargs)

import math
import sys

# 统一物性计算引擎（thermo + CoolProp 双引擎）
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from property.thermo_helper import (
    calc_bubble_point_T,
    calc_dew_point_T,
    calc_bubble_point_P,
    calc_dew_point_P,
    calc_vapor_pressure,
    calc_heat_capacity,
    calc_enthalpy_vaporization,
    get_normal_boiling_point,
    get_molecular_weight,
    get_mixture_properties,
)


# ─────────────────────────────────────────────
# 模块级常量
# ─────────────────────────────────────────────
ATMOSPHERIC_PRESSURE_Pa = 101325.0
# 旧固定压降常量已弃用，改用按塔板数估算
PRESSURE_DROP_ATMOSPHERIC_Pa = 50000.0   # kept for backward-compat reference
PRESSURE_DROP_VACUUM_Pa = 5000.0         # kept for backward-compat reference

# 单板压降估计 (Pa/tray)
_DP_PER_TRAY_ATMOSPHERIC_Pa = 700.0      # 0.7 kPa/tray 常压塔
_DP_PER_TRAY_PRESSURIZED_Pa = 500.0      # 0.5 kPa/tray 加压塔
# 用于 operating_conditions 阶段估算全塔压降的典型塔板数
# （此时尚未计算 N_actual，使用工业典型值 32 板）
_N_TRAY_ESTIMATE_FOR_DP = 32

# 特鲁顿规则常数：ΔHvap(J/mol) ≈ TROUTON_CONST * Tb(K)
# 用于 Hvap 缺失时的物理回退
_TROUTON_CONST = 88.0  # J/(mol·K)

COOLANT_OPTIONS = [
    (50.0, "循环水"),
    (25.0, "冷冻水"),
    (10.0, "低温盐水"),
]

_COOLANT_LADDER = [
    (50.0 + 273.15,  "循环水"),
    (20.0 + 273.15,  "低温水"),
    (0.0 + 273.15,   "乙二醇水溶液"),
    (-25.0 + 273.15, "丙烯制冷"),
    (-45.0 + 273.15, "乙烯深冷"),
]
_MAX_COLUMN_PRESSURE_Pa = 5_000_000


# ─────────────────────────────────────────────
# 物性辅助函数（私有，可共享，基于 thermo_helper）
# ─────────────────────────────────────────────
class _ChemInfo:
    """轻量组分信息容器，基于 thermo_helper，属性接口对齐 thermo.Chemical。

    提供 name、Tb(K)、MW(g/mol)、VaporPressure(T)(Pa)；以及按构造时 (T,P) 计算的
    质量基准 Cpl/Cpg(J·kg⁻¹·K⁻¹) 和 Hvap(J·kg⁻¹)，与原 Chemical 字段一致。
    """

    __slots__ = ("name", "_T", "_P", "_Tb", "_MW")

    def __init__(self, name, T=None, P=101325.0):
        self.name = name
        self._T = T
        self._P = P
        self._Tb = None
        self._MW = None

    @property
    def Tb(self):
        if self._Tb is None:
            val = get_normal_boiling_point(self.name)
            if not val:
                raise ValueError(f"无法获取 {self.name} 的正常沸点 Tb")
            self._Tb = val
        return self._Tb

    @property
    def MW(self):
        if self._MW is None:
            val = get_molecular_weight(self.name)
            if not val:
                raise ValueError(f"无法获取 {self.name} 的分子量 MW")
            self._MW = val
        return self._MW

    def _TP(self):
        T = self._T if self._T is not None else 298.15
        return T, self._P

    @property
    def Cpl(self):
        """液相质量热容 J/(kg·K) = 摩尔热容 * 1000 / MW"""
        T, P = self._TP()
        molar = calc_heat_capacity(self.name, T, P, "liquid")
        return (molar * 1000.0 / self.MW) if molar else None

    @property
    def Cpg(self):
        """气相质量热容 J/(kg·K)"""
        T, P = self._TP()
        molar = calc_heat_capacity(self.name, T, P, "gas")
        return (molar * 1000.0 / self.MW) if molar else None

    @property
    def Hvap(self):
        """质量基准汽化焓 J/kg = 摩尔汽化焓 * 1000 / MW"""
        T, _ = self._TP()
        molar = calc_enthalpy_vaporization(self.name, T)
        return (molar * 1000.0 / self.MW) if molar else None

    def VaporPressure(self, temperature_K):
        return calc_vapor_pressure(self.name, temperature_K)


def _to_chemical(component):
    if isinstance(component, _ChemInfo):
        return component
    if isinstance(component, dict):
        for key in ("name", "CAS", "id", "ID"):
            value = component.get(key)
            if value:
                component = value
                break
    try:
        return _ChemInfo(str(component))
    except Exception as exc:
        raise ValueError(f"Invalid component identifier: {component}") from exc


def _normalize_comp_ids(comp_ids):
    if not isinstance(comp_ids, (list, tuple)) or len(comp_ids) < 2:
        raise ValueError("comp_ids must contain at least two component identifiers")
    chem_1 = _to_chemical(comp_ids[0])
    chem_2 = _to_chemical(comp_ids[1])
    if chem_1.Tb > chem_2.Tb:
        chem_1, chem_2 = chem_2, chem_1
    return chem_1, chem_2


def _get_vapor_pressure_Pa(temperature_K, chem_obj):
    val = chem_obj.VaporPressure(temperature_K)
    if not val:
        raise ValueError(f"无法获取 {chem_obj.name} 在 {temperature_K} K 的饱和蒸汽压")
    return val


def _solve_temp_at_pressure(target_pressure_Pa, mole_frac_light_comp, chem_light, chem_heavy, is_dew=True):
    """泡点/露点温度（二元），返回 T(K)。"""
    comps = [chem_light.name, chem_heavy.name]
    zs = [mole_frac_light_comp, 1.0 - mole_frac_light_comp]
    try:
        if is_dew:
            res = calc_dew_point_T(comps, zs, target_pressure_Pa)
        else:
            res = calc_bubble_point_T(comps, zs, target_pressure_Pa)
        if not res or not res.get("converged") or not res.get("T"):
            raise ValueError("泡露点未收敛")
        return res["T"]
    except Exception as exc:
        raise ValueError(
            f"体系 [{chem_light.name}, {chem_heavy.name}] 在 P={target_pressure_Pa:.0f} Pa "
            f"下泡露点计算失败：{type(exc).__name__}: {exc}"
        ) from exc


def _solve_pressure_at_temp(target_temp_K, mole_frac_light_comp, chem_light, chem_heavy, is_dew=True):
    """泡点/露点压力（二元），返回 P(Pa)。"""
    comps = [chem_light.name, chem_heavy.name]
    zs = [mole_frac_light_comp, 1.0 - mole_frac_light_comp]
    try:
        if is_dew:
            res = calc_dew_point_P(comps, zs, target_temp_K)
        else:
            res = calc_bubble_point_P(comps, zs, target_temp_K)
        if not res or not res.get("converged") or not res.get("P"):
            raise ValueError("泡露点压力未收敛")
        return res["P"]
    except Exception as exc:
        raise ValueError(
            f"体系 [{chem_light.name}, {chem_heavy.name}] 在 T={target_temp_K:.2f} K "
            f"下泡露点压力计算失败：{type(exc).__name__}: {exc}"
        ) from exc


def _solve_bubble_T_multicomp(comp_names, zs, P):
    try:
        res = calc_bubble_point_T(comp_names, zs, P)
        if not res or not res.get("converged") or not res.get("T"):
            raise ValueError("未收敛")
        return res["T"]
    except Exception as exc:
        raise ValueError(
            f"体系 {comp_names} 在 P={P:.0f} Pa 下泡点计算失败：{type(exc).__name__}: {exc}。"
            f"可能是体系非理想性强或压力过低，建议提高压力到 0.2 MPa 以上或更换体系。"
        ) from exc


def _solve_dew_T_multicomp(comp_names, zs, P):
    try:
        res = calc_dew_point_T(comp_names, zs, P)
        if not res or not res.get("converged") or not res.get("T"):
            raise ValueError("未收敛")
        return res["T"]
    except Exception as exc:
        raise ValueError(
            f"体系 {comp_names} 在 P={P:.0f} Pa 下露点计算失败：{type(exc).__name__}: {exc}。"
            f"可能是体系非理想性强或压力过低，建议提高压力到 0.2 MPa 以上或更换体系。"
        ) from exc


def _solve_bubble_P_multicomp(comp_names, zs, T):
    try:
        res = calc_bubble_point_P(comp_names, zs, T)
        if not res or not res.get("converged") or not res.get("P"):
            raise ValueError("未收敛")
        return res["P"]
    except Exception as exc:
        raise ValueError(
            f"体系 {comp_names} 在 T={T:.2f} K 下泡点压力计算失败：{type(exc).__name__}: {exc}。"
            f"可能是体系非理想性强或温度超出物性范围。"
        ) from exc


def _calc_relative_volatility(temperature_K, chem_light, chem_heavy):
    p_light = _get_vapor_pressure_Pa(temperature_K, chem_light)
    p_heavy = _get_vapor_pressure_Pa(temperature_K, chem_heavy)
    return p_light / p_heavy


def _estimate_column_pressure_drop_Pa(top_P):
    """按塔板数估算全塔压降 ΔP_total = N * ΔP_per_tray。

    常压塔（top_P ≤ 0.2 MPa）用 0.7 kPa/tray；加压塔用 0.5 kPa/tray。
    operating_conditions 阶段尚未计算 N_actual，使用工业典型估计值 32 板。
    """
    dp_per_tray = _DP_PER_TRAY_ATMOSPHERIC_Pa if top_P <= 200000.0 else _DP_PER_TRAY_PRESSURIZED_Pa
    return _N_TRAY_ESTIMATE_FOR_DP * dp_per_tray


def _hvap_mass_or_trouton(chem_obj):
    """取质量基准汽化焓 J/kg；缺失时用特鲁顿规则回退：ΔHvap_molar ≈ 88*Tb (J/mol)。

    返回 J/kg（质量基准），与 _ChemInfo.Hvap 字段一致。
    """
    hvap = chem_obj.Hvap
    if hvap and hvap > 0:
        return hvap
    # 特鲁顿规则：摩尔汽化焓 ≈ 88 * Tb(K) (J/mol)，再换算为 J/kg
    Tb = chem_obj.Tb
    MW = chem_obj.MW
    return _TROUTON_CONST * Tb * 1000.0 / MW


def _underwood_lhs(theta, alpha_list, feed_mole_frac_list, q):
    return sum(a * z / (a - theta) for a, z in zip(alpha_list, feed_mole_frac_list)) - (1.0 - q)


def _solve_underwood_theta(alpha_list, feed_mole_frac_list, q, tol=1e-10):
    alpha_light, alpha_heavy = alpha_list[0], alpha_list[1]
    eps = 1e-6
    lo = alpha_heavy + eps
    hi = alpha_light - eps
    f_lo = _underwood_lhs(lo, alpha_list, feed_mole_frac_list, q)
    f_hi = _underwood_lhs(hi, alpha_list, feed_mole_frac_list, q)
    if f_lo * f_hi > 0:
        raise RuntimeError("Underwood 方程端点同号")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f_mid = _underwood_lhs(mid, alpha_list, feed_mole_frac_list, q)
        if abs(f_mid) < tol or (hi - lo) < tol:
            return mid
        if f_lo * f_mid < 0:
            hi = mid
        else:
            lo = mid
            f_lo = f_mid
    return 0.5 * (lo + hi)


# ─────────────────────────────────────────────
# 私有：核心计算（无 print，纯数学，可共享）
# ─────────────────────────────────────────────
def _compute_mass_balance(feed_molar_flows, distillate_purity, bottoms_purity,
                          light_key_component, heavy_key_component,
                          light_components, heavy_components):
    """物料衡算核心：返回与 calc_mass_balance 相同的 dict 但无打印。"""
    if light_key_component not in light_components:
        raise ValueError(f"light_key_component '{light_key_component}' 不在 light_components 中")
    if heavy_key_component not in heavy_components:
        raise ValueError(f"heavy_key_component '{heavy_key_component}' 不在 heavy_components 中")
    if set(light_components) & set(heavy_components):
        raise ValueError("light_components 和 heavy_components 有交集")
    all_classified = set(light_components) | set(heavy_components)
    for comp in feed_molar_flows:
        if comp not in all_classified:
            raise ValueError(f"组分 '{comp}' 未分类")

    lnk_flows = {c: feed_molar_flows[c] for c in light_components if c != light_key_component}
    hnk_flows = {c: feed_molar_flows[c] for c in heavy_components if c != heavy_key_component}
    binary_feed = {
        light_key_component: feed_molar_flows[light_key_component],
        heavy_key_component: feed_molar_flows[heavy_key_component],
    }

    chem_light = _ChemInfo(light_key_component)
    chem_heavy = _ChemInfo(heavy_key_component)
    if chem_light.Tb > chem_heavy.Tb:
        # 交换 chem_light/chem_heavy 时必须同步交换 key 名，
        # 否则 binary_feed 仍按原始 key 取值，会与交换后的 MW 错位
        chem_light, chem_heavy = chem_heavy, chem_light
        light_key_component, heavy_key_component = heavy_key_component, light_key_component

    F_light_mol = binary_feed[light_key_component]
    F_heavy_mol = binary_feed[heavy_key_component]
    F_binary_mol = F_light_mol + F_heavy_mol

    MW_light = chem_light.MW / 1000.0
    MW_heavy = chem_heavy.MW / 1000.0

    F_light_kg = F_light_mol * MW_light
    F_heavy_kg = F_heavy_mol * MW_heavy
    F_binary_kg = F_light_kg + F_heavy_kg

    w_F = F_light_kg / F_binary_kg
    w_D = distillate_purity
    w_W = 1.0 - bottoms_purity

    D_kg = F_binary_kg * (w_F - w_W) / (w_D - w_W)
    W_kg = F_binary_kg - D_kg

    D_light_kg = D_kg * w_D
    D_heavy_kg = D_kg - D_light_kg
    W_light_kg = W_kg * w_W
    W_heavy_kg = W_kg - W_light_kg

    D_light_mol = D_light_kg / MW_light
    D_heavy_mol = D_heavy_kg / MW_heavy
    W_light_mol = W_light_kg / MW_light
    W_heavy_mol = W_heavy_kg / MW_heavy

    D_binary_mol = D_light_mol + D_heavy_mol
    W_binary_mol = W_light_mol + W_heavy_mol

    x_F = F_light_mol / F_binary_mol
    x_D = D_light_mol / D_binary_mol
    x_W = W_light_mol / W_binary_mol

    distillate_flows = dict(lnk_flows)
    distillate_flows[light_key_component] = D_light_mol
    distillate_flows[heavy_key_component] = D_heavy_mol

    bottoms_flows = {}
    bottoms_flows[light_key_component] = W_light_mol
    bottoms_flows[heavy_key_component] = W_heavy_mol
    bottoms_flows.update(hnk_flows)

    D_total_mol = sum(distillate_flows.values())
    W_total_mol = sum(bottoms_flows.values())

    F_total = sum(feed_molar_flows.values())
    feed_zs = {comp: flow / F_total for comp, flow in feed_molar_flows.items()}
    distillate_zs = {comp: flow / D_total_mol for comp, flow in distillate_flows.items()}
    bottoms_zs = {comp: flow / W_total_mol for comp, flow in bottoms_flows.items()}

    return {
        "component_ids": [chem_light.name, chem_heavy.name],
        "mole_fractions": {"F": x_F, "D": x_D, "W": x_W},
        "molar_flows_mol_per_s": {"F": F_binary_mol, "D": D_binary_mol, "W": W_binary_mol},
        "mass_flows_kg_per_s": {"F": F_binary_kg, "D": D_kg, "W": W_kg},
        "is_multicomponent": len(feed_molar_flows) > 2,
        "distillate_flows_mol_per_s": distillate_flows,
        "bottoms_flows_mol_per_s": bottoms_flows,
        "distillate_total_mol_per_s": D_total_mol,
        "bottoms_total_mol_per_s": W_total_mol,
        "feed_zs": feed_zs,
        "distillate_zs": distillate_zs,
        "bottoms_zs": bottoms_zs,
    }


def _try_pressurized_topT_with_coolant_ladder(solve_P_func):
    """按冷却介质阶梯依次尝试塔顶温度，返回 (top_T, top_P, coolant_name)。

    跳过条件：flash 不收敛（临界点附近）；或 top_P > 5 MPa（超工业塔上限）。
    全部失败 → ValueError 并提示加闪蒸罐。
    """
    last_exc = None
    tried = []
    for top_T, coolant_name in _COOLANT_LADDER:
        try:
            top_P = solve_P_func(top_T)
            if top_P > _MAX_COLUMN_PRESSURE_Pa:
                tried.append(f"{coolant_name}({top_T-273.15:.0f}℃, P={top_P/1e6:.1f}MPa>5MPa)")
                continue
            return top_T, top_P, coolant_name
        except (ValueError, RuntimeError, UnboundLocalError) as exc:
            last_exc = exc
            tried.append(f"{coolant_name}({top_T-273.15:.0f}℃, flash失败)")
            continue
    raise ValueError(
        f"所有冷却介质都无法让塔顶在 5 MPa 以内液化。已尝试: {', '.join(tried)}。"
        f"最后错误: {type(last_exc).__name__}: {last_exc}。"
        f"建议检查体系是否含不可冷凝组分（如 N2/H2/CH4），上游应使用闪蒸罐先脱除轻气体。"
    )


def _compute_operating_conditions(comp_ids, x_D, x_W, distillate_zs=None, bottoms_zs=None):
    """操作压力与温度核心。"""
    chem_light, chem_heavy = _normalize_comp_ids(comp_ids)
    has_multicomp = distillate_zs is not None and bottoms_zs is not None

    if has_multicomp:
        D_comps = list(distillate_zs.keys())
        D_fracs = list(distillate_zs.values())
        W_comps = list(bottoms_zs.keys())
        W_fracs = list(bottoms_zs.values())
        top_bubble_atm_K = _solve_bubble_T_multicomp(D_comps, D_fracs, ATMOSPHERIC_PRESSURE_Pa)
        if top_bubble_atm_K - 273.15 >= 50.0:
            op_type = "常压精馏 - 循环水冷却"
            coolant_name = "循环水"
            top_P = ATMOSPHERIC_PRESSURE_Pa
            bot_P = top_P + _estimate_column_pressure_drop_Pa(top_P)
            top_T = top_bubble_atm_K
            bot_T = _solve_bubble_T_multicomp(W_comps, W_fracs, bot_P)
        else:
            top_T, top_P, coolant_name = _try_pressurized_topT_with_coolant_ladder(
                lambda T: _solve_bubble_P_multicomp(D_comps, D_fracs, T)
            )
            op_type = f"加压精馏 - {coolant_name}冷却"
            bot_P = top_P + _estimate_column_pressure_drop_Pa(top_P)
            bot_T = _solve_bubble_T_multicomp(W_comps, W_fracs, bot_P)
    else:
        top_bubble_atm_K = _solve_temp_at_pressure(ATMOSPHERIC_PRESSURE_Pa, x_D, chem_light, chem_heavy, is_dew=False)
        if top_bubble_atm_K - 273.15 >= 50.0:
            op_type = "常压精馏 - 循环水冷却"
            coolant_name = "循环水"
            top_P = ATMOSPHERIC_PRESSURE_Pa
            bot_P = top_P + _estimate_column_pressure_drop_Pa(top_P)
            top_T = top_bubble_atm_K
            bot_T = _solve_temp_at_pressure(bot_P, x_W, chem_light, chem_heavy, is_dew=False)
        else:
            top_T, top_P, coolant_name = _try_pressurized_topT_with_coolant_ladder(
                lambda T: _solve_pressure_at_temp(T, x_D, chem_light, chem_heavy, is_dew=False)
            )
            op_type = f"加压精馏 - {coolant_name}冷却"
            bot_P = top_P + _estimate_column_pressure_drop_Pa(top_P)
            bot_T = _solve_temp_at_pressure(bot_P, x_W, chem_light, chem_heavy, is_dew=False)

    return {
        "operation_type": op_type,
        "coolant_type": coolant_name,
        "column_top_pressure_Pa": top_P,
        "column_bottom_pressure_Pa": bot_P,
        "column_top_temperature_K": top_T,
        "column_bottom_temperature_K": bot_T,
    }


def _compute_feed_thermal_condition(feed_temp_K, feed_pressure_Pa, x_F, comp_ids, feed_zs=None):
    """进料热状况核心。"""
    chem_light, chem_heavy = _normalize_comp_ids(comp_ids)
    has_multicomp = feed_zs is not None

    if has_multicomp:
        comps = list(feed_zs.keys())
        fracs = list(feed_zs.values())
        T_bubble = _solve_bubble_T_multicomp(comps, fracs, feed_pressure_Pa)
        T_dew = _solve_dew_T_multicomp(comps, fracs, feed_pressure_Pa)
    else:
        T_bubble = _solve_temp_at_pressure(feed_pressure_Pa, x_F, chem_light, chem_heavy, is_dew=False)
        T_dew = _solve_temp_at_pressure(feed_pressure_Pa, x_F, chem_light, chem_heavy, is_dew=True)

    tol = 0.5
    if feed_temp_K < T_bubble - tol:
        chem_l_obj = _ChemInfo(chem_light.name, T=(feed_temp_K + T_bubble) / 2.0, P=feed_pressure_Pa)
        chem_h_obj = _ChemInfo(chem_heavy.name, T=(feed_temp_K + T_bubble) / 2.0, P=feed_pressure_Pa)
        Cpl = x_F * (chem_l_obj.Cpl or 150.0) + (1 - x_F) * (chem_h_obj.Cpl or 150.0)
        chem_l_bub = _ChemInfo(chem_light.name, T=T_bubble, P=feed_pressure_Pa)
        Hvap = _hvap_mass_or_trouton(chem_l_bub)
        q = 1.0 + Cpl * (T_bubble - feed_temp_K) / Hvap
        state = "过冷液体"
    elif abs(feed_temp_K - T_bubble) <= tol:
        q = 1.0
        state = "饱和液体"
    elif feed_temp_K < T_dew - tol:
        q = (T_dew - feed_temp_K) / (T_dew - T_bubble)
        state = "气液混合"
    elif abs(feed_temp_K - T_dew) <= tol:
        q = 0.0
        state = "饱和蒸汽"
    else:
        chem_l_obj = _ChemInfo(chem_light.name, T=(T_dew + feed_temp_K) / 2.0, P=feed_pressure_Pa)
        chem_h_obj = _ChemInfo(chem_heavy.name, T=(T_dew + feed_temp_K) / 2.0, P=feed_pressure_Pa)
        Cpg = x_F * (chem_l_obj.Cpg or 50.0) + (1 - x_F) * (chem_h_obj.Cpg or 50.0)
        chem_l_dew = _ChemInfo(chem_light.name, T=T_dew, P=feed_pressure_Pa)
        Hvap = _hvap_mass_or_trouton(chem_l_dew)
        q = -Cpg * (feed_temp_K - T_dew) / Hvap
        state = "过热蒸汽"
    return {"q": q, "bubble_point_temperature_K": T_bubble, "dew_point_temperature_K": T_dew, "feed_state_description": state}


def _compute_min_reflux_ratio(comp_ids, x_F, x_D, x_W, top_T, bot_T, feed_T, q,
                              reflux_factor, feed_zs=None, distillate_zs=None):
    """最小回流比核心。

    二元体系采用全塔平均相对挥发度 α_avg = sqrt(α_top * α_bot)，同时用于
    Underwood 第一方程（求 θ）与第二方程（求 R_min），避免进料/塔顶温度下
    α 不一致带来的偏差。多组分体系按传统做法分别用 feed_T 和 top_T 的 α。
    """
    chem_light, chem_heavy = _normalize_comp_ids(comp_ids)
    has_multicomp = feed_zs is not None and distillate_zs is not None

    if has_multicomp:
        feed_comps = list(feed_zs.keys())
        dist_comps = list(distillate_zs.keys())
        Tb_dict = {comp: _ChemInfo(comp).Tb for comp in feed_comps}
        heaviest_comp = max(Tb_dict, key=Tb_dict.get)
        chem_heaviest = _ChemInfo(heaviest_comp)

        P_sat_base_feed = _get_vapor_pressure_Pa(feed_T, chem_heaviest)
        alpha_feed_list = []
        feed_frac_list = []
        for comp in feed_comps:
            chem = _ChemInfo(comp)
            alpha_feed_list.append(_get_vapor_pressure_Pa(feed_T, chem) / P_sat_base_feed)
            feed_frac_list.append(feed_zs[comp])

        P_sat_base_top = _get_vapor_pressure_Pa(top_T, chem_heaviest)
        alpha_top_list = []
        dist_frac_list = []
        for comp in dist_comps:
            chem = _ChemInfo(comp)
            alpha_top_list.append(_get_vapor_pressure_Pa(top_T, chem) / P_sat_base_top)
            dist_frac_list.append(distillate_zs[comp])

        alpha_feed_LK_HK = _calc_relative_volatility(feed_T, chem_light, chem_heavy)

        alpha_LK_feed = _get_vapor_pressure_Pa(feed_T, chem_light) / P_sat_base_feed
        alpha_HK_feed = _get_vapor_pressure_Pa(feed_T, chem_heavy) / P_sat_base_feed
        eps = 1e-6
        theta_lo = min(alpha_LK_feed, alpha_HK_feed) + eps
        theta_hi = max(alpha_LK_feed, alpha_HK_feed) - eps
        f_lo = _underwood_lhs(theta_lo, alpha_feed_list, feed_frac_list, q)
        f_hi = _underwood_lhs(theta_hi, alpha_feed_list, feed_frac_list, q)
        if f_lo * f_hi > 0:
            raise RuntimeError("Underwood 方程端点同号")
        for _ in range(200):
            theta_mid = 0.5 * (theta_lo + theta_hi)
            f_mid = _underwood_lhs(theta_mid, alpha_feed_list, feed_frac_list, q)
            if abs(f_mid) < 1e-10 or (theta_hi - theta_lo) < 1e-10:
                break
            if f_lo * f_mid < 0:
                theta_hi = theta_mid
            else:
                theta_lo = theta_mid
                f_lo = f_mid
        theta = 0.5 * (theta_lo + theta_hi)
        R_min = sum(a * x / (a - theta) for a, x in zip(alpha_top_list, dist_frac_list)) - 1.0
        R_op = round(reflux_factor * R_min, 1)
        return {"R_min": R_min, "R_operating": R_op, "alpha_feed": alpha_feed_LK_HK, "theta": theta}
    else:
        # 全塔平均相对挥发度（塔顶与塔底的几何平均），同时用于 Underwood 两方程
        alpha_top = _calc_relative_volatility(top_T, chem_light, chem_heavy)
        alpha_bot = _calc_relative_volatility(bot_T, chem_light, chem_heavy)
        alpha_avg = math.sqrt(alpha_top * alpha_bot)
        if alpha_avg <= 1.0:
            raise ValueError("全塔平均相对挥发度 <= 1.0")
        theta = _solve_underwood_theta([alpha_avg, 1.0], [x_F, 1.0 - x_F], q)
        R_min = sum(a * x / (a - theta) for a, x in zip([alpha_avg, 1.0], [x_D, 1.0 - x_D])) - 1.0
        R_op = round(reflux_factor * R_min, 1)
        return {"R_min": R_min, "R_operating": R_op, "alpha_feed": alpha_avg, "theta": theta}


def _compute_min_theoretical_stages(x_D, x_W, top_T, bot_T, comp_ids):
    """Fenske 最小理论塔板数核心。"""
    chem_light, chem_heavy = _normalize_comp_ids(comp_ids)
    alpha_top = _calc_relative_volatility(top_T, chem_light, chem_heavy)
    alpha_bot = _calc_relative_volatility(bot_T, chem_light, chem_heavy)
    alpha_geom = math.sqrt(alpha_top * alpha_bot)
    numerator = math.log((x_D / (1.0 - x_D)) * ((1.0 - x_W) / x_W))
    N_min = numerator / math.log(alpha_geom)
    return {"N_min": N_min, "alpha_geometric_avg": alpha_geom}


def _compute_theoretical_stages(R_op, R_min, N_min):
    """Gilliland 理论塔板数核心。"""
    X = (R_op - R_min) / (R_op + 1.0)
    Y = 1.0 - math.exp((1.0 + 54.4 * X) / (11.0 + 117.2 * X) * (X - 1.0) / math.sqrt(X))
    return (Y + N_min) / (1.0 - Y)


def _compute_actual_stages(N_theoretical, alpha_geom, x_F, top_T, bot_T, comp_ids, feed_zs=None):
    """O'Connell 实际塔板数核心。

    N_theoretical（Gilliland）含再沸器 1 级、不含全凝器：
    实际塔板数 = ceil((N_theoretical − 1) / E0)，先扣除再沸器再除以塔板效率，
    得到塔内物理塔板数。
    """
    chem_light, chem_heavy = _normalize_comp_ids(comp_ids)
    avg_T = (top_T + bot_T) / 2.0
    if feed_zs is not None:
        comps = list(feed_zs.keys())
        fracs = list(feed_zs.values())
        bubble_P = _solve_bubble_P_multicomp(comps, fracs, avg_T)
        eval_P = bubble_P + 100.0
        props = get_mixture_properties(comps, fracs, avg_T, eval_P)
        mu_l = props.get("mu_L")
    else:
        bubble_P = _solve_pressure_at_temp(avg_T, x_F, chem_light, chem_heavy, is_dew=False)
        eval_P = bubble_P + 100.0
        props = get_mixture_properties([chem_light.name, chem_heavy.name], [x_F, 1.0 - x_F], avg_T, eval_P)
        mu_l = props.get("mu_L")
    mu_l = mu_l or 0.001
    mu_mPas = mu_l * 1000.0
    E0 = 0.492 * (mu_mPas * alpha_geom) ** (-0.245)
    E0 = max(0.1, min(E0, 1.0))
    # 扣除再沸器 1 级后按效率换算 → 塔内物理塔板数
    N_actual = max(1, math.ceil((N_theoretical - 1.0) / E0))
    return {"N_actual": N_actual, "tray_efficiency": E0}


# ═════════════════════════════════════════════
# 工具 1：物料衡算
# ═════════════════════════════════════════════
def calc_mass_balance(
    feed_molar_flows: dict,
    distillate_purity: float,
    bottoms_purity: float,
    light_key_component: str,
    heavy_key_component: str,
    light_components: list,
    heavy_components: list,
) -> dict:
    """计算精馏塔物料衡算。

    根据进料各组分流量和塔顶/塔底产品纯度，按质量守恒计算塔顶塔底各组分流量与组成。
    支持二元和多组分体系（轻非关键组分全部去塔顶，重非关键组分全部去塔底）。
    返回扁平一层 dict，所有 value 为 number/string/bool 或 {组分名: 数值}。

    返回扁平一层 dict，所有 value 为 number/string/bool 或 {组分名: 数值}。
    """
    mb = _compute_mass_balance(feed_molar_flows, distillate_purity, bottoms_purity,
                                light_key_component, heavy_key_component,
                                light_components, heavy_components)
    result = {
        "is_multicomponent": mb["is_multicomponent"],
        "feed_lk_mole_fraction": mb["mole_fractions"]["F"],
        "distillate_lk_mole_fraction": mb["mole_fractions"]["D"],
        "bottoms_lk_mole_fraction": mb["mole_fractions"]["W"],
        "feed_binary_molar_flow_mol_per_s": mb["molar_flows_mol_per_s"]["F"],
        "distillate_binary_molar_flow_mol_per_s": mb["molar_flows_mol_per_s"]["D"],
        "bottoms_binary_molar_flow_mol_per_s": mb["molar_flows_mol_per_s"]["W"],
        "feed_mass_flow_kg_per_s": mb["mass_flows_kg_per_s"]["F"],
        "distillate_mass_flow_kg_per_s": mb["mass_flows_kg_per_s"]["D"],
        "bottoms_mass_flow_kg_per_s": mb["mass_flows_kg_per_s"]["W"],
        "distillate_flows_mol_per_s": mb["distillate_flows_mol_per_s"],
        "bottoms_flows_mol_per_s": mb["bottoms_flows_mol_per_s"],
        "distillate_total_mol_per_s": mb["distillate_total_mol_per_s"],
        "bottoms_total_mol_per_s": mb["bottoms_total_mol_per_s"],
    }
    _vprint(f"\n【calc_mass_balance】")
    _vprint(f"  组分: {mb['component_ids']}")
    _vprint(f"  二元 x_F={result['feed_lk_mole_fraction']:.4f}, x_D={result['distillate_lk_mole_fraction']:.4f}, x_W={result['bottoms_lk_mole_fraction']:.4f}")
    _vprint(f"  二元 F={result['feed_binary_molar_flow_mol_per_s']:.4f}, D={result['distillate_binary_molar_flow_mol_per_s']:.4f}, W={result['bottoms_binary_molar_flow_mol_per_s']:.4f} mol/s")
    return result


# ═════════════════════════════════════════════
# 工具 2：操作条件
# ═════════════════════════════════════════════
def calc_operating_conditions(
    feed_molar_flows: dict,
    distillate_purity: float,
    bottoms_purity: float,
    light_key_component: str,
    heavy_key_component: str,
    light_components: list,
    heavy_components: list,
) -> dict:
    """判定精馏塔操作压力与塔顶塔底温度。

    根据塔顶组成的露点温度判断是否需要加压操作：若常压下塔顶温度低于冷却水可达
    温度则加压，否则常压。塔底压力 = 塔顶压力 + 压降。返回 operation_type、
    column_top/bottom_pressure_Pa、column_top/bottom_temperature_K。
    本工具独立计算：内部完成物料衡算。

    返回:
        {operation_type, coolant_type, column_top_pressure_Pa, column_bottom_pressure_Pa,
         column_top_temperature_K, column_bottom_temperature_K}
    """
    mb = _compute_mass_balance(feed_molar_flows, distillate_purity, bottoms_purity,
                                light_key_component, heavy_key_component,
                                light_components, heavy_components)
    distillate_zs = mb["distillate_zs"] if mb["is_multicomponent"] else None
    bottoms_zs = mb["bottoms_zs"] if mb["is_multicomponent"] else None
    result = _compute_operating_conditions(
        mb["component_ids"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        distillate_zs=distillate_zs, bottoms_zs=bottoms_zs)
    _vprint(f"\n【calc_operating_conditions】")
    _vprint(f"  {result['operation_type']}")
    _vprint(f"  塔顶 P={result['column_top_pressure_Pa']:.0f} Pa, T={result['column_top_temperature_K']-273.15:.2f}°C")
    _vprint(f"  塔底 P={result['column_bottom_pressure_Pa']:.0f} Pa, T={result['column_bottom_temperature_K']-273.15:.2f}°C")
    return result


# ═════════════════════════════════════════════
# 工具 3：进料热状况
# ═════════════════════════════════════════════
def calc_feed_thermal_condition(
    feed_molar_flows: dict,
    distillate_purity: float,
    bottoms_purity: float,
    light_key_component: str,
    heavy_key_component: str,
    light_components: list,
    heavy_components: list,
    feed_temp_K: float,
    feed_pressure_Pa: float = None,
) -> dict:
    """计算进料热状况参数 q（无量纲）。

    q 表示进料中液相的比例：q=1 为泡点液体，q=0 为露点蒸汽，0<q<1 为汽液两相。
    函数根据进料温度相对于泡点和露点的位置计算 q 值。返回 q、
    bubble_point_temperature_K、dew_point_temperature_K、feed_state_description。

    进料压力若不指定（None），则由塔操作条件内部推导（取塔顶/塔底压力的平均值），
    避免用户输入与塔压不一致造成偏差；该值同时支持减压/加压塔。
    本工具独立计算：内部完成物料衡算、操作条件、泡露点求解。

    返回:
        {q, bubble_point_temperature_K, dew_point_temperature_K, feed_state_description,
         feed_pressure_Pa_used}
    """
    mb = _compute_mass_balance(feed_molar_flows, distillate_purity, bottoms_purity,
                                light_key_component, heavy_key_component,
                                light_components, heavy_components)
    if feed_pressure_Pa is None:
        # 从塔操作条件推导：进料板压力 ≈ (塔顶 + 塔底)/2
        distillate_zs = mb["distillate_zs"] if mb["is_multicomponent"] else None
        bottoms_zs = mb["bottoms_zs"] if mb["is_multicomponent"] else None
        op = _compute_operating_conditions(
            mb["component_ids"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
            distillate_zs=distillate_zs, bottoms_zs=bottoms_zs)
        feed_pressure_Pa = 0.5 * (op["column_top_pressure_Pa"] + op["column_bottom_pressure_Pa"])
    elif feed_pressure_Pa <= 0:
        raise ValueError(f"进料压力必须为正，得到 {feed_pressure_Pa}")
    feed_zs = mb["feed_zs"] if mb["is_multicomponent"] else None
    result = _compute_feed_thermal_condition(
        feed_temp_K, feed_pressure_Pa, mb["mole_fractions"]["F"], mb["component_ids"], feed_zs=feed_zs)
    result["feed_pressure_Pa_used"] = feed_pressure_Pa
    _vprint(f"\n【calc_feed_thermal_condition】")
    _vprint(f"  进料 T={feed_temp_K-273.15:.2f}°C, P={feed_pressure_Pa:.0f} Pa")
    _vprint(f"  泡点={result['bubble_point_temperature_K']-273.15:.2f}°C, 露点={result['dew_point_temperature_K']-273.15:.2f}°C")
    _vprint(f"  状态: {result['feed_state_description']}, q = {result['q']:.4f}")
    return result


# ═════════════════════════════════════════════
# 工具 4：最小回流比
# ═════════════════════════════════════════════
def calc_min_reflux_ratio(
    feed_molar_flows: dict,
    distillate_purity: float,
    bottoms_purity: float,
    light_key_component: str,
    heavy_key_component: str,
    light_components: list,
    heavy_components: list,
    feed_temp_K: float,
    feed_pressure_Pa: float,
    reflux_factor: float = 1.5,
) -> dict:
    """计算最小回流比 R_min 与操作回流比 R_operating（无量纲）。

    采用 Underwood 方程求解，操作回流比 = R_min × reflux_factor（默认 1.5 倍），并取 0.1 精度。
    R_min 是维持指定分离效果所需的最低回流比，实际操作必须大于此值。
    返回 R_min、R_operating、alpha_feed、theta。
    本工具独立计算：内部完成物料衡算、操作条件、进料热状况推算。

    返回:
        {R_min, R_operating, alpha_feed, theta}
    """
    mb = _compute_mass_balance(feed_molar_flows, distillate_purity, bottoms_purity,
                                light_key_component, heavy_key_component,
                                light_components, heavy_components)
    feed_zs = mb["feed_zs"] if mb["is_multicomponent"] else None
    distillate_zs = mb["distillate_zs"] if mb["is_multicomponent"] else None
    bottoms_zs = mb["bottoms_zs"] if mb["is_multicomponent"] else None

    op = _compute_operating_conditions(
        mb["component_ids"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        distillate_zs=distillate_zs, bottoms_zs=bottoms_zs)
    fq = _compute_feed_thermal_condition(
        feed_temp_K, feed_pressure_Pa, mb["mole_fractions"]["F"], mb["component_ids"], feed_zs=feed_zs)
    result = _compute_min_reflux_ratio(
        mb["component_ids"], mb["mole_fractions"]["F"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        op["column_top_temperature_K"], op["column_bottom_temperature_K"],
        feed_temp_K, fq["q"], reflux_factor,
        feed_zs=feed_zs, distillate_zs=distillate_zs)
    _vprint(f"\n【calc_min_reflux_ratio】")
    _vprint(f"  R_min = {result['R_min']:.4f}, R_op = {result['R_operating']:.1f}")
    return result


# ═════════════════════════════════════════════
# 工具 5：最小理论塔板数
# ═════════════════════════════════════════════
def calc_min_theoretical_stages(
    feed_molar_flows: dict,
    distillate_purity: float,
    bottoms_purity: float,
    light_key_component: str,
    heavy_key_component: str,
    light_components: list,
    heavy_components: list,
) -> dict:
    """计算最小理论塔板数 N_min（无量纲）。

    采用 Fenske 方程，在全回流条件下由塔顶塔底组成和相对挥发度求得。
    N_min 是达到指定分离所需的最少理论板数。返回 N_min、alpha_geometric_avg。
    本工具独立计算：内部完成物料衡算与操作条件推算。

    返回:
        {N_min, alpha_geometric_avg}
    """
    mb = _compute_mass_balance(feed_molar_flows, distillate_purity, bottoms_purity,
                                light_key_component, heavy_key_component,
                                light_components, heavy_components)
    distillate_zs = mb["distillate_zs"] if mb["is_multicomponent"] else None
    bottoms_zs = mb["bottoms_zs"] if mb["is_multicomponent"] else None
    op = _compute_operating_conditions(
        mb["component_ids"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        distillate_zs=distillate_zs, bottoms_zs=bottoms_zs)
    result = _compute_min_theoretical_stages(
        mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        op["column_top_temperature_K"], op["column_bottom_temperature_K"], mb["component_ids"])
    _vprint(f"\n【calc_min_theoretical_stages】")
    _vprint(f"  α_geom = {result['alpha_geometric_avg']:.4f}, N_min = {result['N_min']:.2f}")
    return result


# ═════════════════════════════════════════════
# 工具 6：理论塔板数（Gilliland）
# ═════════════════════════════════════════════
def calc_theoretical_stages(
    feed_molar_flows: dict,
    distillate_purity: float,
    bottoms_purity: float,
    light_key_component: str,
    heavy_key_component: str,
    light_components: list,
    heavy_components: list,
    feed_temp_K: float,
    feed_pressure_Pa: float,
    reflux_factor: float = 1.5,
) -> dict:
    """计算理论塔板数 N_theoretical（无量纲）。

    采用 Gilliland 关联式，由操作回流比、最小回流比和最小理论板数求得实际操作
    条件下的理论板数。返回 N_theoretical。
    本工具独立计算：内部完成 mass_balance、operating_conditions、feed_thermal、
    min_reflux、min_stages 全链路推算。

    返回:
        {N_theoretical}
    """
    mb = _compute_mass_balance(feed_molar_flows, distillate_purity, bottoms_purity,
                                light_key_component, heavy_key_component,
                                light_components, heavy_components)
    feed_zs = mb["feed_zs"] if mb["is_multicomponent"] else None
    distillate_zs = mb["distillate_zs"] if mb["is_multicomponent"] else None
    bottoms_zs = mb["bottoms_zs"] if mb["is_multicomponent"] else None

    op = _compute_operating_conditions(
        mb["component_ids"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        distillate_zs=distillate_zs, bottoms_zs=bottoms_zs)
    fq = _compute_feed_thermal_condition(
        feed_temp_K, feed_pressure_Pa, mb["mole_fractions"]["F"], mb["component_ids"], feed_zs=feed_zs)
    rr = _compute_min_reflux_ratio(
        mb["component_ids"], mb["mole_fractions"]["F"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        op["column_top_temperature_K"], op["column_bottom_temperature_K"],
        feed_temp_K, fq["q"], reflux_factor,
        feed_zs=feed_zs, distillate_zs=distillate_zs)
    fs = _compute_min_theoretical_stages(
        mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        op["column_top_temperature_K"], op["column_bottom_temperature_K"], mb["component_ids"])

    N_theoretical = _compute_theoretical_stages(rr["R_operating"], rr["R_min"], fs["N_min"])
    _vprint(f"\n【calc_theoretical_stages】N_theoretical = {N_theoretical:.2f}")
    return {"N_theoretical": N_theoretical}


# ═════════════════════════════════════════════
# 工具 7：实际塔板数
# ═════════════════════════════════════════════
def calc_actual_stages(
    feed_molar_flows: dict,
    distillate_purity: float,
    bottoms_purity: float,
    light_key_component: str,
    heavy_key_component: str,
    light_components: list,
    heavy_components: list,
    feed_temp_K: float,
    feed_pressure_Pa: float,
    reflux_factor: float = 1.5,
) -> dict:
    """计算实际塔板数 N_actual（无量纲）。

    N_actual = N_theoretical / 板效率，板效率采用 O'Connell 关联式由相对挥发度
    和液相粘度估算。返回 N_actual、tray_efficiency。
    本工具独立计算。

    返回:
        {N_actual, tray_efficiency}
    """
    mb = _compute_mass_balance(feed_molar_flows, distillate_purity, bottoms_purity,
                                light_key_component, heavy_key_component,
                                light_components, heavy_components)
    feed_zs = mb["feed_zs"] if mb["is_multicomponent"] else None
    distillate_zs = mb["distillate_zs"] if mb["is_multicomponent"] else None
    bottoms_zs = mb["bottoms_zs"] if mb["is_multicomponent"] else None

    op = _compute_operating_conditions(
        mb["component_ids"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        distillate_zs=distillate_zs, bottoms_zs=bottoms_zs)
    fq = _compute_feed_thermal_condition(
        feed_temp_K, feed_pressure_Pa, mb["mole_fractions"]["F"], mb["component_ids"], feed_zs=feed_zs)
    rr = _compute_min_reflux_ratio(
        mb["component_ids"], mb["mole_fractions"]["F"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        op["column_top_temperature_K"], op["column_bottom_temperature_K"],
        feed_temp_K, fq["q"], reflux_factor,
        feed_zs=feed_zs, distillate_zs=distillate_zs)
    fs = _compute_min_theoretical_stages(
        mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        op["column_top_temperature_K"], op["column_bottom_temperature_K"], mb["component_ids"])
    N_theoretical = _compute_theoretical_stages(rr["R_operating"], rr["R_min"], fs["N_min"])
    result = _compute_actual_stages(
        N_theoretical, fs["alpha_geometric_avg"],
        mb["mole_fractions"]["F"], op["column_top_temperature_K"], op["column_bottom_temperature_K"],
        mb["component_ids"], feed_zs=feed_zs)
    _vprint(f"\n【calc_actual_stages】N_actual = {result['N_actual']}, E0 = {result['tray_efficiency']:.4f}")
    return result


# ═════════════════════════════════════════════
# 工具 8：进料板位置（Kirkbride）
# ═════════════════════════════════════════════
def calc_feed_stage(
    feed_molar_flows: dict,
    distillate_purity: float,
    bottoms_purity: float,
    light_key_component: str,
    heavy_key_component: str,
    light_components: list,
    heavy_components: list,
    feed_temp_K: float,
    feed_pressure_Pa: float,
    reflux_factor: float = 1.5,
) -> dict:
    """计算进料板位置。

    采用 Kirkbride 公式确定精馏段与提馏段板数分配。返回 N_rectifying、
    N_stripping、feed_tray_position_from_top。
    本工具独立计算：内部完成全链路推算。

    返回:
        {N_rectifying, N_stripping, feed_tray_position_from_top}
    """
    mb = _compute_mass_balance(feed_molar_flows, distillate_purity, bottoms_purity,
                                light_key_component, heavy_key_component,
                                light_components, heavy_components)
    feed_zs = mb["feed_zs"] if mb["is_multicomponent"] else None
    distillate_zs = mb["distillate_zs"] if mb["is_multicomponent"] else None
    bottoms_zs = mb["bottoms_zs"] if mb["is_multicomponent"] else None

    op = _compute_operating_conditions(
        mb["component_ids"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        distillate_zs=distillate_zs, bottoms_zs=bottoms_zs)
    fq = _compute_feed_thermal_condition(
        feed_temp_K, feed_pressure_Pa, mb["mole_fractions"]["F"], mb["component_ids"], feed_zs=feed_zs)
    rr = _compute_min_reflux_ratio(
        mb["component_ids"], mb["mole_fractions"]["F"], mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        op["column_top_temperature_K"], op["column_bottom_temperature_K"],
        feed_temp_K, fq["q"], reflux_factor,
        feed_zs=feed_zs, distillate_zs=distillate_zs)
    fs = _compute_min_theoretical_stages(
        mb["mole_fractions"]["D"], mb["mole_fractions"]["W"],
        op["column_top_temperature_K"], op["column_bottom_temperature_K"], mb["component_ids"])
    N_theoretical = _compute_theoretical_stages(rr["R_operating"], rr["R_min"], fs["N_min"])
    actual = _compute_actual_stages(
        N_theoretical, fs["alpha_geometric_avg"],
        mb["mole_fractions"]["F"], op["column_top_temperature_K"], op["column_bottom_temperature_K"],
        mb["component_ids"], feed_zs=feed_zs)

    x_F = mb["mole_fractions"]["F"]
    x_D = mb["mole_fractions"]["D"]
    x_W = mb["mole_fractions"]["W"]
    D_mol_s = mb["molar_flows_mol_per_s"]["D"]
    W_mol_s = mb["molar_flows_mol_per_s"]["W"]
    N_actual = actual["N_actual"]
    E0 = actual["tray_efficiency"]

    K_kirk = (W_mol_s / D_mol_s) * ((1.0 - x_F) / x_F) * (x_W / (1.0 - x_D)) ** 2
    ratio = K_kirk ** 0.206

    # ① 理论级分割（Kirkbride，作用于含再沸器的理论级数）→ 理论进料板编号
    #    （自塔内第一块板起计；提馏段理论级数含再沸器，保底 2 级）
    N_strip_theo = max(2, int(N_theoretical / (1.0 + ratio)))
    N_rect_theo = N_theoretical - N_strip_theo
    feed_stage_theoretical = N_rect_theo + 1.0

    # ② 统一效率近似换算：进料板以上理论板数 / E0 向上取整 → 实际板
    N_rect = max(1, math.ceil((feed_stage_theoretical - 1.0) / E0))
    feed_tray = N_rect + 1
    # 提馏段实际板数 = 总实际板数 − 进料板（均不含进料板，N_rect + 1 + N_strip = N_actual）
    N_strip = N_actual - feed_tray
    if N_strip < 1:
        # 兜底：至少保留 1 块提馏段塔板（压缩精馏段换算裕量）
        N_strip = 1
        feed_tray = N_actual - N_strip
        N_rect = feed_tray - 1

    _vprint(f"\n【calc_feed_stage】理论进料板={feed_stage_theoretical:.2f}, E0={E0:.4f}")
    _vprint(f"【calc_feed_stage】N_rect={N_rect}, N_strip={N_strip}, 进料板={feed_tray} "
            f"(校验: {N_rect}+1+{N_strip}={N_rect + 1 + N_strip} = N_actual={N_actual})")
    return {
        "N_rectifying": N_rect,
        "N_stripping": N_strip,
        "feed_tray_position_from_top": feed_tray,
        "feed_stage_theoretical": feed_stage_theoretical,
        "N_actual": N_actual,
        "tray_efficiency": E0,
    }


if __name__ == "__main__":
    feed_molar_flows = {'N-PENTANE': 5.0, 'N-HEXANE': 15.0, 'N-HEPTANE': 20.0, 'N-OCTANE': 10.0}
    distillate_purity = 0.99
    bottoms_purity = 0.999
    feed_temp_K = 313.15
    feed_pressure_Pa = 200000.0
    light_key_component = 'N-HEXANE'
    heavy_key_component = 'N-HEPTANE'
    light_components = ['N-PENTANE', 'N-HEXANE']
    heavy_components = ['N-HEPTANE', 'N-OCTANE']
    common = dict(
        feed_molar_flows=feed_molar_flows,
        distillate_purity=distillate_purity,
        bottoms_purity=bottoms_purity,
        light_key_component=light_key_component,
        heavy_key_component=heavy_key_component,
        light_components=light_components,
        heavy_components=heavy_components,
    )
    common_with_feed = dict(common, feed_temp_K=feed_temp_K, feed_pressure_Pa=feed_pressure_Pa)

    calc_mass_balance(**common)
    calc_operating_conditions(**common)
    calc_feed_thermal_condition(**common_with_feed)
    calc_min_reflux_ratio(**common_with_feed)
    calc_min_theoretical_stages(**common)
    calc_theoretical_stages(**common_with_feed)
    calc_actual_stages(**common_with_feed)
    calc_feed_stage(**common_with_feed)
