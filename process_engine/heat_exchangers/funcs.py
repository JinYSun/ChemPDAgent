"""换热器独立工具函数（不依赖其他 public calc_* 函数的输出）。

设计原则：
- 每个 public calc_* 只接收原始题目参数
- public calc_* 之间禁止互相调用
- 私有 _ 辅助函数和常量表可共享
- 每个函数只返回自己的最终结果

 
"""

import os
import sys

_VERBOSE = os.getenv("FUNCS_VERBOSE", "0") == "1"


def _vprint(*args, **kwargs):
    if _VERBOSE:
        print(*args, **kwargs)


import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from property.thermo_helper import (
    calc_bubble_point_T,
    calc_dew_point_T,
    calc_flash_rachford_rice,
    calc_enthalpy_vaporization,
    get_mixture_properties,
    get_fluid_MW,
)


# ─────────────────────────────────────────────
# 物性辅助函数（私有）
# ─────────────────────────────────────────────
def _calc_phase(components, zs, T, P):
    """用 thermo_helper 判断相态，返回 'l' 或 'g'。"""
    result = calc_flash_rachford_rice(components, zs, T, P)
    return 'l' if result.get("beta", 0.0) < 0.5 else 'g'


class _MixProps:
    """thermo_helper 物性容器。

    提供 MW(g/mol)、Cp/Cpg/Cpl(J·kg⁻¹·K⁻¹)、Hvapms(J/mol)、mu(Pa·s)。

    [FIX-3] 新增 phase_forced 参数：调用方可强制指定 'l'/'g'，使得显热模式
    在平均温度可能落入两相区时，依然用正确相态的 Cp 和 μ，而不依赖 flash 判断。
    """

    __slots__ = ("components", "zs", "T", "P", "MW", "Cp", "Cpg", "Cpl", "Hvapms", "mu")

    def __init__(self, components, mole_fractions, T, P, phase_forced=None):
        """
        Parameters
        ----------
        phase_forced : str or None
            'l' 强制取液相物性，'g' 强制取气相物性，None 则由 flash 自动判断。
        """
        self.components = components
        self.zs = mole_fractions
        self.T = T
        self.P = P
        MW_mix = get_fluid_MW(components, mole_fractions)  # g/mol
        self.MW = MW_mix
        props = get_mixture_properties(components, mole_fractions, T, P)
        cp_l_molar = props.get("Cp_L")
        cp_v_molar = props.get("Cp_V")
        self.Cpl = (cp_l_molar * 1000.0 / MW_mix) if cp_l_molar else None
        self.Cpg = (cp_v_molar * 1000.0 / MW_mix) if cp_v_molar else None

        # [FIX-3] 相态确定：phase_forced 优先，否则用 flash
        if phase_forced is not None:
            phase = phase_forced
        else:
            phase = _calc_phase(components, mole_fractions, T, P)

        if phase == 'l':
            self.Cp = self.Cpl if self.Cpl else self.Cpg
            self.mu = props.get("mu_L")
        else:
            self.Cp = self.Cpg if self.Cpg else self.Cpl
            self.mu = props.get("mu_V")

        # 各组分纯态汽化潜热 J/mol（T 超过 Tc 时按 0 处理）
        hvap = []
        for c in components:
            try:
                h = calc_enthalpy_vaporization(c, T)
            except Exception:
                h = None
            hvap.append(h if h else 0.0)
        self.Hvapms = hvap


# ─────────────────────────────────────────────
# 永久气体名单（不可凝组分）
# ─────────────────────────────────────────────
_PERMANENT_GAS_NAMES = {
    "hydrogen", "h2",
    "nitrogen", "n2",
    "oxygen", "o2",
    "carbon monoxide", "co",
    "argon", "ar",
    "helium", "he",
    "neon", "ne",
}


def _has_permanent_gas(components, zs, threshold: float = 0.01) -> bool:
    """体系中是否含有摩尔分数 > threshold 的永久气体。"""
    for comp, z in zip(components, zs):
        if z >= threshold and comp.lower().strip() in _PERMANENT_GAS_NAMES:
            return True
    return False


def _calc_bubble_dew_pr(components, zs, P):
    """计算泡点和露点温度，返回 (T_bubble, T_dew)。"""
    if _has_permanent_gas(components, zs):
        _vprint(
            f"  [_calc_bubble_dew_pr] 体系 {components} 含永久气体，"
            f"工况下无液相，按显热气相处理（T_bubble=T_dew=1K）"
        )
        return 1.0, 1.0

    try:
        bub = calc_bubble_point_T(components, zs, P)
        dew = calc_dew_point_T(components, zs, P)
        T_bubble = bub.get("T") if (bub and bub.get("converged")) else None
        T_dew = dew.get("T") if (dew and dew.get("converged")) else None
        if not T_bubble or not T_dew:
            raise ValueError("泡露点未收敛")
    except Exception as exc:
        _vprint(
            f"  [_calc_bubble_dew_pr] 体系 {components} 在 P={P:.0f} Pa 下 "
            f"泡露点计算失败（{type(exc).__name__}: {exc}），fallback 为显热气相"
        )
        return 1.0, 1.0
    return T_bubble, T_dew


def _parse_molar_flows(process_fluid_molar_flows_mol_per_s, T, P, phase_forced=None):
    """解析摩尔流量，返回 (components, mole_fractions, total_molar_flow, mass_flow_kg_per_s, mixture)。

    [FIX-6] 新增返回 total_molar_flow_mol_per_s，避免调用方重复 sum()。
    [FIX-3] 透传 phase_forced 到 _MixProps。
    """
    components = list(process_fluid_molar_flows_mol_per_s.keys())
    flows = list(process_fluid_molar_flows_mol_per_s.values())
    total_molar_flow_mol_per_s = sum(flows)
    mole_fractions = [f / total_molar_flow_mol_per_s for f in flows]
    mixture = _MixProps(components, mole_fractions, T, P, phase_forced=phase_forced)
    mass_flow_kg_per_s = total_molar_flow_mol_per_s * mixture.MW / 1000.0
    return components, mole_fractions, total_molar_flow_mol_per_s, mass_flow_kg_per_s, mixture


# ─────────────────────────────────────────────
# 经验 K 值表（私有）
# ─────────────────────────────────────────────
EMPIRICAL_COOLING_HEAT_TRANSFER_COEFFICIENT_TABLE = {
    ("cooling_water", "light_organic"):  (467.0, 814.0),
    ("cooling_water", "medium_organic"): (290.0, 698.0),
    ("cooling_water", "heavy_organic"):  (116.0, 467.0),
}
EMPIRICAL_CONDENSATION_HEAT_TRANSFER_COEFFICIENT_TABLE = {
    ("cooling_water", "organic_condensation"): (350.0, 1200.0),
}
EMPIRICAL_HEATING_HEAT_TRANSFER_COEFFICIENT_TABLE = {
    ("steam", "light_organic"):  (580.0, 1160.0),
    ("steam", "medium_organic"): (350.0, 930.0),
    ("steam", "heavy_organic"):  (175.0, 700.0),
}
EMPIRICAL_EVAPORATION_HEAT_TRANSFER_COEFFICIENT_TABLE = {
    ("steam", "organic_evaporation"): (580.0, 1750.0),
}

# [FIX-5] 合法 utility_type 白名单，用于友好报错
_VALID_UTILITY_TYPES_COOLING = {"cooling_water"}
_VALID_UTILITY_TYPES_HEATING = {"steam"}
_VALID_UTILITY_TYPES_ALL = _VALID_UTILITY_TYPES_COOLING | _VALID_UTILITY_TYPES_HEATING


# ─────────────────────────────────────────────
# 换热模式判定
# ─────────────────────────────────────────────
def _determine_heat_transfer_mode(T_in, T_out, T_bubble, T_dew):
    """根据进出口温度与泡露点判定换热模式。返回 (mode, mode_cn, is_cooling)。"""
    is_cooling = T_in > T_out
    is_pure = abs(T_dew - T_bubble) < 0.1

    if is_pure:
        T_boil = T_bubble
        if not is_cooling:
            if T_out <= T_boil:
                mode, mode_cn = "sensible", "显热模式（纯液相）"
            elif T_in >= T_boil:
                mode, mode_cn = "sensible", "显热模式（纯气相）"
            else:
                mode, mode_cn = "evaporation", "蒸发模式（纯组分完全汽化）"
        else:
            if T_out >= T_boil:
                mode, mode_cn = "sensible", "显热模式（纯气相）"
            elif T_in <= T_boil:
                mode, mode_cn = "sensible", "显热模式（纯液相）"
            else:
                mode, mode_cn = "condensation", "冷凝模式（纯组分完全冷凝）"
    else:
        if not is_cooling:
            if T_out <= T_bubble:
                mode, mode_cn = "sensible", "显热模式（纯液相）"
            elif T_in <= T_bubble and T_out >= T_dew:
                mode, mode_cn = "evaporation", "蒸发模式（完全汽化）"
            elif T_in <= T_bubble and T_out < T_dew:
                mode, mode_cn = "partial_evaporation", "部分蒸发模式（液相预热+部分汽化）"
            elif T_in >= T_dew:
                mode, mode_cn = "sensible", "显热模式（纯气相）"
            elif T_bubble < T_in < T_dew and T_out >= T_dew:
                mode, mode_cn = "partial_evaporation_no_preheat", "部分蒸发模式（两相区→全汽化+过热）"
            else:
                mode, mode_cn = "two_phase_heating", "两相区加热模式"
        else:
            if T_out >= T_dew:
                mode, mode_cn = "sensible", "显热模式（纯气相）"
            elif T_in >= T_dew and T_out <= T_bubble:
                mode, mode_cn = "condensation", "冷凝模式（完全冷凝）"
            elif T_in >= T_dew and T_out > T_bubble:
                mode, mode_cn = "partial_condensation", "部分冷凝模式（过热冷却+部分冷凝）"
            elif T_in <= T_bubble:
                mode, mode_cn = "sensible", "显热模式（纯液相）"
            elif T_bubble < T_in < T_dew and T_out <= T_bubble:
                mode, mode_cn = "partial_condensation_no_desuperheat", "部分冷凝模式（两相区→全冷凝+过冷）"
            else:
                mode, mode_cn = "two_phase_cooling", "两相区冷却模式"
    return mode, mode_cn, is_cooling


# ─────────────────────────────────────────────
# [FIX-1][FIX-2] 两相区焓差数值积分核心函数
# ─────────────────────────────────────────────
def _hvap_mix_at_T(process_fluid_molar_flows_mol_per_s, T_eval, P, total_molar_flow_mol_per_s):
    """计算混合物在温度 T_eval 处的摩尔汽化焓 ΔHvap_mix (J/mol·进料)。

    [FIX-2] 正确做法：用 flash 计算 T_eval 处的液相组成 xi，以 xi 加权各组分
    纯态汽化潜热，而不是用进料总摩尔分数 zi。

    返回值单位：J / mol进料（可直接乘以总摩尔流量得 W）。
    """
    components = list(process_fluid_molar_flows_mol_per_s.keys())
    flows = list(process_fluid_molar_flows_mol_per_s.values())
    zs = [f / total_molar_flow_mol_per_s for f in flows]

    # flash 获取液相组成
    flash_result = calc_flash_rachford_rice(components, zs, T_eval, P)
    x_liq = flash_result.get("x", zs)   # 液相摩尔分数

    # 各组分纯态汽化潜热
    hvaps = []
    for c in components:
        try:
            h = calc_enthalpy_vaporization(c, T_eval)
        except Exception:
            h = None
        hvaps.append(h if h else 0.0)

    # 用液相组成 xi 加权
    hvap_mix = sum(x_liq[i] * hvaps[i] for i in range(len(components)))
    return hvap_mix   # J/mol（液相基准）


def _integrate_two_phase_enthalpy(
    process_fluid_molar_flows_mol_per_s, T_lo, T_hi, P,
    total_molar_flow_mol_per_s, n_segments=10
):
    """对两相区 [T_lo, T_hi] 进行焓差数值积分，返回热负荷 Q (W)。

    [FIX-1] 取代线性温度插值。在泡露点之间均匀划分 n_segments 个子区间，
    每个子区间用梯形法则：
        Q_seg ≈ Δβ × ΔHvap_mix_avg × F_total_mol
    其中 Δβ 为该子区间内汽化率增量（线性假设仅在子区间内，
    误差与 1/n_segments² 成正比）。

    对宽沸程混合物（泡露点差 > 50K），建议 n_segments ≥ 20。
    """
    if abs(T_hi - T_lo) < 1e-6:
        return 0.0

    T_bubble = T_lo
    T_dew = T_hi
    dT = (T_hi - T_lo) / n_segments
    Q_total = 0.0

    for i in range(n_segments):
        T_a = T_lo + i * dT
        T_b = T_a + dT
        T_mid = (T_a + T_b) / 2.0

        # 梯形法：用区间中点的 ΔHvap_mix 代表该段（Simpson-like midpoint rule）
        hvap_mid = _hvap_mix_at_T(
            process_fluid_molar_flows_mol_per_s, T_mid, P, total_molar_flow_mol_per_s
        )

        # 汽化率增量（线性近似，仅在小子区间内使用）
        d_beta = dT / (T_dew - T_bubble)

        Q_total += total_molar_flow_mol_per_s * hvap_mid * d_beta

    return Q_total


# ─────────────────────────────────────────────
# 核心 Q 计算
# ─────────────────────────────────────────────
def _compute_heat_duty(process_fluid_molar_flows_mol_per_s, T_in, T_out, P, mode,
                       n_segments=10):
    """根据预先判定的 mode 计算热负荷 Q (W)。

    n_segments : 两相区数值积分节点数，默认 10，宽沸程体系建议 ≥ 20。
    """
    # [FIX-6] 统一从 _parse_molar_flows 获取 total_molar_flow，避免重复 sum()
    components, mole_fractions, total_molar_flow_mol_per_s, mass_flow_kg_per_s, _ = \
        _parse_molar_flows(process_fluid_molar_flows_mol_per_s, T_in, P)

    T_bubble, T_dew = _calc_bubble_dew_pr(components, mole_fractions, P)

    # ── 显热：根据温度区间判断实际相态，强制指定 phase_forced ─────────────
    if mode == "sensible":
        # [FIX-3] 判断显热区间对应的实际相态
        T_avg = (T_in + T_out) / 2.0
        is_cooling = T_in > T_out
        if T_avg <= T_bubble:
            phase_forced = 'l'
        elif T_avg >= T_dew:
            phase_forced = 'g'
        else:
            # 平均温度在两相区：按进出口都在同侧判断
            # sensible 模式本身保证了 T_in、T_out 都在同一相区
            if is_cooling:
                # 冷却显热：T_out >= T_dew（纯气）或 T_in <= T_bubble（纯液）
                phase_forced = 'g' if T_out >= T_dew else 'l'
            else:
                phase_forced = 'l' if T_out <= T_bubble else 'g'

        _, _, _, mass_flow_kg_per_s, mixture = _parse_molar_flows(
            process_fluid_molar_flows_mol_per_s, T_avg, P, phase_forced=phase_forced
        )
        return mass_flow_kg_per_s * mixture.Cp * abs(T_in - T_out)

    # ── 冷凝：过热冷却 + 相变潜热（数值积分）+ 子冷 ────────────────────────
    elif mode == "condensation":
        Q_superheat = 0.0
        Q_latent = 0.0
        Q_subcool = 0.0

        # 过热段（纯气相）
        if T_in > T_dew:
            _, _, _, mass_flow_i, mix_sh = _parse_molar_flows(
                process_fluid_molar_flows_mol_per_s,
                (T_in + T_dew) / 2.0, P, phase_forced='g'
            )
            Q_superheat = mass_flow_i * mix_sh.Cpg * (T_in - T_dew)

        # [FIX-1][FIX-2] 相变潜热：数值积分
        Q_latent = _integrate_two_phase_enthalpy(
            process_fluid_molar_flows_mol_per_s,
            T_bubble, T_dew, P, total_molar_flow_mol_per_s, n_segments
        )

        # 子冷段（纯液相）
        if T_out < T_bubble:
            _, _, _, mass_flow_i, mix_sc = _parse_molar_flows(
                process_fluid_molar_flows_mol_per_s,
                (T_bubble + T_out) / 2.0, P, phase_forced='l'
            )
            Q_subcool = mass_flow_i * mix_sc.Cpl * (T_bubble - T_out)

        return Q_superheat + Q_latent + Q_subcool

    # ── 蒸发：预热 + 相变潜热（数值积分）+ 过热 ────────────────────────────
    elif mode == "evaporation":
        Q_preheat = 0.0
        Q_latent = 0.0
        Q_superheat = 0.0

        # 预热段（纯液相）
        if T_in < T_bubble:
            _, _, _, mass_flow_i, mix_ph = _parse_molar_flows(
                process_fluid_molar_flows_mol_per_s,
                (T_in + T_bubble) / 2.0, P, phase_forced='l'
            )
            Q_preheat = mass_flow_i * mix_ph.Cpl * (T_bubble - T_in)

        # [FIX-1][FIX-2] 相变潜热：数值积分
        Q_latent = _integrate_two_phase_enthalpy(
            process_fluid_molar_flows_mol_per_s,
            T_bubble, T_dew, P, total_molar_flow_mol_per_s, n_segments
        )

        # 过热段（纯气相）
        if T_out > T_dew:
            _, _, _, mass_flow_i, mix_sh = _parse_molar_flows(
                process_fluid_molar_flows_mol_per_s,
                (T_dew + T_out) / 2.0, P, phase_forced='g'
            )
            Q_superheat = mass_flow_i * mix_sh.Cpg * (T_out - T_dew)

        return Q_preheat + Q_latent + Q_superheat

    # ── 部分蒸发：预热 + 部分相变（数值积分，T_bubble → T_out） ─────────────
    elif mode == "partial_evaporation":
        # 预热段（液相）
        _, _, _, mass_flow_i, mix_ph = _parse_molar_flows(
            process_fluid_molar_flows_mol_per_s,
            (T_in + T_bubble) / 2.0, P, phase_forced='l'
        )
        Q_preheat = mass_flow_i * mix_ph.Cpl * (T_bubble - T_in)

        # [FIX-1][FIX-2] 部分汽化区间（T_bubble → T_out）数值积分
        Q_partial_vap = _integrate_two_phase_enthalpy(
            process_fluid_molar_flows_mol_per_s,
            T_bubble, T_out, P, total_molar_flow_mol_per_s, n_segments
        )

        return Q_preheat + Q_partial_vap

    # ── 部分蒸发（无预热）：两相区剩余 + 过热 ────────────────────────────────
    elif mode == "partial_evaporation_no_preheat":
        # [FIX-1][FIX-2] T_in（两相区内）→ T_dew 的剩余汽化
        Q_partial_vap = _integrate_two_phase_enthalpy(
            process_fluid_molar_flows_mol_per_s,
            T_in, T_dew, P, total_molar_flow_mol_per_s, n_segments
        )

        # 过热段
        Q_superheat = 0.0
        if T_out > T_dew:
            _, _, _, mass_flow_i, mix_sh = _parse_molar_flows(
                process_fluid_molar_flows_mol_per_s,
                (T_dew + T_out) / 2.0, P, phase_forced='g'
            )
            Q_superheat = mass_flow_i * mix_sh.Cpg * (T_out - T_dew)

        return Q_partial_vap + Q_superheat

    # ── 两相区加热（T_in、T_out 都在泡露点之间）────────────────────────────
    elif mode == "two_phase_heating":
        # [FIX-1][FIX-2] 直接对 [T_in, T_out] 积分
        return _integrate_two_phase_enthalpy(
            process_fluid_molar_flows_mol_per_s,
            T_in, T_out, P, total_molar_flow_mol_per_s, n_segments
        )

    # ── 部分冷凝：过热冷却 + 部分冷凝（T_dew → T_out） ─────────────────────
    elif mode == "partial_condensation":
        Q_desuperheat = 0.0
        if T_in > T_dew:
            _, _, _, mass_flow_i, mix_sh = _parse_molar_flows(
                process_fluid_molar_flows_mol_per_s,
                (T_in + T_dew) / 2.0, P, phase_forced='g'
            )
            Q_desuperheat = mass_flow_i * mix_sh.Cpg * (T_in - T_dew)

        # [FIX-1][FIX-2] T_dew → T_out 部分冷凝数值积分
        Q_partial_cond = _integrate_two_phase_enthalpy(
            process_fluid_molar_flows_mol_per_s,
            T_out, T_dew, P, total_molar_flow_mol_per_s, n_segments
        )

        return Q_desuperheat + Q_partial_cond

    # ── 部分冷凝（无降温）：两相区剩余 + 子冷 ────────────────────────────────
    elif mode == "partial_condensation_no_desuperheat":
        # [FIX-1][FIX-2] T_bubble → T_in（两相区内）的剩余冷凝
        Q_partial_cond = _integrate_two_phase_enthalpy(
            process_fluid_molar_flows_mol_per_s,
            T_bubble, T_in, P, total_molar_flow_mol_per_s, n_segments
        )

        # 子冷段
        Q_subcool = 0.0
        if T_out < T_bubble:
            _, _, _, mass_flow_i, mix_sc = _parse_molar_flows(
                process_fluid_molar_flows_mol_per_s,
                (T_bubble + T_out) / 2.0, P, phase_forced='l'
            )
            Q_subcool = mass_flow_i * mix_sc.Cpl * (T_bubble - T_out)

        return Q_partial_cond + Q_subcool

    # ── 两相区冷却（T_in、T_out 都在泡露点之间）─────────────────────────────
    elif mode == "two_phase_cooling":
        # [FIX-1][FIX-2] 对 [T_out, T_in] 积分（T_out < T_in）
        return _integrate_two_phase_enthalpy(
            process_fluid_molar_flows_mol_per_s,
            T_out, T_in, P, total_molar_flow_mol_per_s, n_segments
        )

    raise ValueError(f"Unknown mode: {mode}")


# ─────────────────────────────────────────────
# 核心 K 计算
# ─────────────────────────────────────────────
def _compute_overall_k(process_fluid_molar_flows_mol_per_s, T_in, T_out, P,
                       mode, is_cooling, utility_type, selection):
    """根据预先判定的 mode 计算总传热系数 K (W/(m²·K))。返回 (k_min, k_max, k_selected)。"""

    # [FIX-5] 校验 utility_type，提前给出友好报错
    if utility_type not in _VALID_UTILITY_TYPES_ALL:
        raise ValueError(
            f"不支持的公用工程类型 '{utility_type}'。"
            f"目前支持：{sorted(_VALID_UTILITY_TYPES_ALL)}。"
            f"如需扩展，请在 EMPIRICAL_*_HEAT_TRANSFER_COEFFICIENT_TABLE 中添加对应条目。"
        )

    if mode in ("condensation", "partial_condensation",
                "partial_condensation_no_desuperheat", "two_phase_cooling"):
        key = ("cooling_water", "organic_condensation")
        if key not in EMPIRICAL_CONDENSATION_HEAT_TRANSFER_COEFFICIENT_TABLE:
            raise ValueError(
                f"冷凝模式查表键 {key} 不存在于 EMPIRICAL_CONDENSATION_HEAT_TRANSFER_COEFFICIENT_TABLE。"
            )
        k_min, k_max = EMPIRICAL_CONDENSATION_HEAT_TRANSFER_COEFFICIENT_TABLE[key]

    elif mode in ("evaporation", "partial_evaporation",
                  "partial_evaporation_no_preheat", "two_phase_heating"):
        key = (utility_type, "organic_evaporation")
        if key not in EMPIRICAL_EVAPORATION_HEAT_TRANSFER_COEFFICIENT_TABLE:
            raise ValueError(
                f"蒸发模式查表键 {key} 不存在于 EMPIRICAL_EVAPORATION_HEAT_TRANSFER_COEFFICIENT_TABLE。"
                f"蒸发模式目前仅支持 utility_type='steam'。"
            )
        k_min, k_max = EMPIRICAL_EVAPORATION_HEAT_TRANSFER_COEFFICIENT_TABLE[key]

    elif is_cooling:
        # [FIX-3] 显热冷却：在平均温度处强制取液相黏度（冷却段一般为液相）
        bulk_temp_K = (T_in + T_out) / 2.0
        # 判断冷却显热区间的相态
        components, mole_fractions, _, _, _ = _parse_molar_flows(
            process_fluid_molar_flows_mol_per_s, bulk_temp_K, P
        )
        T_bubble, T_dew = _calc_bubble_dew_pr(components, mole_fractions, P)
        if bulk_temp_K >= T_dew:
            phase_forced_k = 'g'
        else:
            phase_forced_k = 'l'

        _, _, _, _, mixture = _parse_molar_flows(
            process_fluid_molar_flows_mol_per_s, bulk_temp_K, P,
            phase_forced=phase_forced_k
        )
        mu_mPas = (mixture.mu or 0.0) * 1000.0
        if mu_mPas < 0.5:
            cat = "light_organic"
        elif mu_mPas <= 1.0:
            cat = "medium_organic"
        else:
            cat = "heavy_organic"
        key = ("cooling_water", cat)
        if key not in EMPIRICAL_COOLING_HEAT_TRANSFER_COEFFICIENT_TABLE:
            raise ValueError(
                f"冷却显热查表键 {key} 不存在于 EMPIRICAL_COOLING_HEAT_TRANSFER_COEFFICIENT_TABLE。"
            )
        k_min, k_max = EMPIRICAL_COOLING_HEAT_TRANSFER_COEFFICIENT_TABLE[key]

    else:
        # [FIX-3] 显热加热：在平均温度处强制取液相黏度（加热段一般为液相）
        bulk_temp_K = (T_in + T_out) / 2.0
        components, mole_fractions, _, _, _ = _parse_molar_flows(
            process_fluid_molar_flows_mol_per_s, bulk_temp_K, P
        )
        T_bubble, T_dew = _calc_bubble_dew_pr(components, mole_fractions, P)
        if bulk_temp_K >= T_dew:
            phase_forced_k = 'g'
        else:
            phase_forced_k = 'l'

        _, _, _, _, mixture = _parse_molar_flows(
            process_fluid_molar_flows_mol_per_s, bulk_temp_K, P,
            phase_forced=phase_forced_k
        )
        mu_mPas = (mixture.mu or 0.0) * 1000.0
        if mu_mPas < 0.5:
            cat = "light_organic"
        elif mu_mPas <= 1.0:
            cat = "medium_organic"
        else:
            cat = "heavy_organic"
        key = (utility_type, cat)
        if key not in EMPIRICAL_HEATING_HEAT_TRANSFER_COEFFICIENT_TABLE:
            raise ValueError(
                f"加热显热查表键 {key} 不存在于 EMPIRICAL_HEATING_HEAT_TRANSFER_COEFFICIENT_TABLE。"
                f"加热模式目前仅支持 utility_type='steam'。"
            )
        k_min, k_max = EMPIRICAL_HEATING_HEAT_TRANSFER_COEFFICIENT_TABLE[key]

    sel_map = {"min": k_min, "max": k_max, "mid": (k_min + k_max) / 2.0}
    k_selected = sel_map.get(selection.lower())
    if k_selected is None:
        raise ValueError(f"selection 参数无效：'{selection}'，应为 'min' / 'mid' / 'max'。")
    return k_min, k_max, k_selected


# ─────────────────────────────────────────────
# 核心 LMTD 计算
# ─────────────────────────────────────────────
def _compute_log_mean_temp_difference(T_p_in, T_p_out, T_u_in, T_u_out):
    """计算逆流对数平均温差 LMTD (K)。

    [FIX-4] 补全热端温差（dt1）的合理性校验，防止 dt1 ≤ 0 产生负 LMTD 或 NaN。
    冷却场景逆流定义：
        dt1 = T_p_in  - T_u_out  （热端，工艺入口 vs 公用工程出口）
        dt2 = T_p_out - T_u_in   （冷端，工艺出口 vs 公用工程入口）
    加热场景逆流定义：
        dt1 = T_u_in  - T_p_out  → 用绝对值统一处理
        dt2 = T_u_out - T_p_in   → 用绝对值统一处理
    """
    is_cooling = T_p_in > T_p_out

    # ── 冷端校验（工艺出口 vs 公用工程入口）────────────────────────────────
    if is_cooling and T_p_out < T_u_in:
        raise ValueError(
            f"冷却场景下工艺流体出口温度 ({T_p_out:.2f} K) 低于公用工程入口温度 "
            f"({T_u_in:.2f} K)，违反热力学第二定律。"
            f"冷却介质温度必须始终低于工艺流体。"
        )
    if (not is_cooling) and T_p_out > T_u_in:
        raise ValueError(
            f"加热场景下工艺流体出口温度 ({T_p_out:.2f} K) 高于公用工程入口温度 "
            f"({T_u_in:.2f} K)，违反热力学第二定律。"
            f"加热介质温度必须始终高于工艺流体。"
        )

    # [FIX-4] ── 热端校验（工艺入口 vs 公用工程出口）────────────────────────
    dt1 = abs(T_p_in - T_u_out)
    dt2 = abs(T_p_out - T_u_in)

    if is_cooling and (T_p_in <= T_u_out):
        raise ValueError(
            f"冷却场景下工艺流体入口温度 ({T_p_in:.2f} K) 不高于公用工程出口温度 "
            f"({T_u_out:.2f} K)，热端温差 ΔT1 = {T_p_in - T_u_out:.2f} K ≤ 0，"
            f"违反热力学第二定律，请检查公用工程出口温度设置。"
        )
    if (not is_cooling) and (T_p_in >= T_u_out):
        raise ValueError(
            f"加热场景下工艺流体入口温度 ({T_p_in:.2f} K) 不低于公用工程出口温度 "
            f"({T_u_out:.2f} K)，热端温差 ΔT1 = {T_u_out - T_p_in:.2f} K ≤ 0，"
            f"违反热力学第二定律，请检查公用工程出口温度设置。"
        )

    if dt1 <= 0 or dt2 <= 0:
        raise ValueError(
            f"温差必须大于 0（dt1={dt1:.4f} K, dt2={dt2:.4f} K），"
            f"请检查进出口温度设置。"
        )

    if abs(dt1 - dt2) < 1e-6:
        return dt1
    return (dt1 - dt2) / math.log(dt1 / dt2)


# ═════════════════════════════════════════════
# 工具 1：热负荷
# ═════════════════════════════════════════════
def calc_heat_duty(
    process_fluid_molar_flows_mol_per_s: dict,
    process_fluid_temp_in_K: float,
    process_fluid_temp_out_K: float,
    process_fluid_pressure_Pa: float,
) -> dict:
    """计算换热器热负荷 Q（单位 W）。

    热负荷表示单位时间内需要传递的热量，直接决定换热器规模。函数根据工艺流体
    进出口温度自动判定换热方向（加热/冷却），利用 thermo 库计算泡点露点，按相态
    自动选择计算模式：
      - 纯液相/气相显热：Q = m·Cp·ΔT（强制指定相态，避免两相区物性误判）
      - 完全蒸发/冷凝：分段计算（预热/降温 + 数值积分潜热 + 过热/子冷）
      - 部分蒸发/冷凝/两相区：数值积分（10 节点梯形法）

    修复说明（相比原版）：
      [FIX-1] 两相区潜热改为数值积分（消除宽沸程线性插值误差）
      [FIX-2] 潜热组成加权改用相平衡液相组成 xi
      [FIX-3] 显热段强制指定正确相态取 Cp

    参数:
        process_fluid_molar_flows_mol_per_s — {组分名: 摩尔流量} (mol/s)
        process_fluid_temp_in_K  — 工艺流体初始温度 (K)
        process_fluid_temp_out_K — 工艺流体最终温度 (K)
        process_fluid_pressure_Pa — 工艺流体操作压力 (Pa)

    返回:
        {"heat_duty_W": 换热系统热负荷 (W)}
    """
    T_in = process_fluid_temp_in_K
    T_out = process_fluid_temp_out_K
    P = process_fluid_pressure_Pa

    components, mole_fractions, _, _, _ = _parse_molar_flows(
        process_fluid_molar_flows_mol_per_s, T_in, P)
    T_bubble, T_dew = _calc_bubble_dew_pr(components, mole_fractions, P)
    mode, mode_cn, _ = _determine_heat_transfer_mode(T_in, T_out, T_bubble, T_dew)
    direction = "冷却" if T_in > T_out else "加热"

    _vprint(f"\n【calc_heat_duty】")
    _vprint(f"  自动判定：{direction}方向 - {mode_cn}")
    _vprint(f"  泡点 = {T_bubble - 273.15:.2f}°C, 露点 = {T_dew - 273.15:.2f}°C")

    heat_duty_W = _compute_heat_duty(process_fluid_molar_flows_mol_per_s, T_in, T_out, P, mode)
    _vprint(f"  Q = {heat_duty_W:.2f} W")
    return {"heat_duty_W": heat_duty_W}


# ═════════════════════════════════════════════
# 工具 2：对数平均温差
# ═════════════════════════════════════════════
def calc_log_mean_temp_difference(
    process_fluid_temp_in_K: float,
    process_fluid_temp_out_K: float,
    utility_fluid_temp_in_K: float,
    utility_fluid_temp_out_K: float,
) -> dict:
    """计算换热器对数平均温差 LMTD（单位 K，逆流定义）。

    LMTD 反映换热器两端温差的综合效果，是计算换热面积的关键参数。在逆流换热中，
    ΔT1 = |T工艺进 − T公用出|，ΔT2 = |T工艺出 − T公用进|，
    ΔTm = (ΔT1−ΔT2)/ln(ΔT1/ΔT2)。当两端温差相等时 ΔTm 直接等于该温差。
    函数自动确定热端冷端，支持加热和冷却两个方向。

    修复说明：
      [FIX-4] 补全热端温差 dt1 ≤ 0 的校验，防止产生负 LMTD 或 NaN。

    注意：本函数仅适用于纯逆流换热器。对于 1-2 型管壳式换热器，需在
    结果上乘以 F 校正因子（F < 1），本函数不提供该修正。

    参数:
        process_fluid_temp_in_K  — 工艺流体初始温度 (K)
        process_fluid_temp_out_K — 工艺流体最终温度 (K)
        utility_fluid_temp_in_K  — 公用工程流体初始温度 (K)
        utility_fluid_temp_out_K — 公用工程流体最终温度 (K)

    返回:
        {"log_mean_temp_difference_K": 对数平均温差 (K)}
    """
    lmtd = _compute_log_mean_temp_difference(
        process_fluid_temp_in_K, process_fluid_temp_out_K,
        utility_fluid_temp_in_K, utility_fluid_temp_out_K)
    _vprint(f"\n【calc_log_mean_temp_difference】LMTD = {lmtd:.2f} K")
    return {"log_mean_temp_difference_K": lmtd}


# ═════════════════════════════════════════════
# 工具 3：总传热系数
# ═════════════════════════════════════════════
def calc_overall_heat_transfer_coefficient(
    process_fluid_molar_flows_mol_per_s: dict,
    process_fluid_temp_in_K: float,
    process_fluid_temp_out_K: float,
    process_fluid_pressure_Pa: float,
    utility_type: str = "cooling_water",
    selection: str = "mid",
) -> dict:
    """计算换热器总传热系数 K（单位 W/(m²·K)）。

    总传热系数反映热量从热流体穿过管壁传递到冷流体的综合能力。函数根据工艺
    流体物性自动判定相态和换热模式，结合公用工程类型查询经验传热系数表，
    返回 K 值范围（min/mid/max）。

    修复说明：
      [FIX-3] 显热区间强制指定相态取μ，避免两相区误判
      [FIX-5] utility_type 非法值抛出友好 ValueError（含可用选项列表）

    参数:
        process_fluid_molar_flows_mol_per_s — {组分名: 摩尔流量} (mol/s)
        process_fluid_temp_in_K  — 工艺流体初始温度 (K)
        process_fluid_temp_out_K — 工艺流体最终温度 (K)
        process_fluid_pressure_Pa — 工艺流体操作压力 (Pa)
        utility_type — "cooling_water" 或 "steam"
        selection — "min" / "mid"（默认）/ "max"

    返回:
        {"overall_k_selected_W_per_m2K": 选定传热系数,
         "overall_k_min_W_per_m2K": 最小传热系数,
         "overall_k_max_W_per_m2K": 最大传热系数}  单位 W/(m²·K)
    """
    T_in = process_fluid_temp_in_K
    T_out = process_fluid_temp_out_K
    P = process_fluid_pressure_Pa

    components, mole_fractions, _, _, _ = _parse_molar_flows(
        process_fluid_molar_flows_mol_per_s, T_in, P)
    T_bubble, T_dew = _calc_bubble_dew_pr(components, mole_fractions, P)
    mode, mode_cn, is_cooling = _determine_heat_transfer_mode(T_in, T_out, T_bubble, T_dew)

    k_min, k_max, k_selected = _compute_overall_k(
        process_fluid_molar_flows_mol_per_s, T_in, T_out, P,
        mode, is_cooling, utility_type, selection)

    direction = "冷却" if is_cooling else "加热"
    _vprint(f"\n【calc_overall_heat_transfer_coefficient】")
    _vprint(f"  自动判定：{direction}方向 - {mode_cn}")
    _vprint(f"  公用工程类型：{utility_type}")
    _vprint(f"  K = {k_selected:.2f} W/(m^2*K)")

    return {
        "overall_k_selected_W_per_m2K": k_selected,
        "overall_k_min_W_per_m2K": k_min,
        "overall_k_max_W_per_m2K": k_max,
    }


# ═════════════════════════════════════════════
# 工具 4：换热面积
# ═════════════════════════════════════════════
def calc_heat_transfer_area(
    process_fluid_molar_flows_mol_per_s: dict,
    process_fluid_temp_in_K: float,
    process_fluid_temp_out_K: float,
    process_fluid_pressure_Pa: float,
    utility_fluid_temp_in_K: float,
    utility_fluid_temp_out_K: float,
    utility_type: str = "cooling_water",
    selection: str = "mid",
) -> dict:
    """计算换热器换热面积 A（单位 m²）。

    根据传热基本方程 A = Q / (K × ΔTm) 求得面积。函数内部独立完成
    热负荷 Q、对数平均温差 ΔTm、总传热系数 K 的全部计算后求得面积。

    注意：本函数假设纯逆流换热（F=1）。对于 1-2 型管壳式换热器，
    实际面积 = 本函数结果 / F（F < 1，通常 0.75~0.95），需由调用方修正。

    修复说明：集成 FIX-1 至 FIX-6 的全部修复。

    参数:
        process_fluid_molar_flows_mol_per_s — {组分名: 摩尔流量} (mol/s)
        process_fluid_temp_in_K  — 工艺流体初始温度 (K)
        process_fluid_temp_out_K — 工艺流体最终温度 (K)
        process_fluid_pressure_Pa — 工艺流体操作压力 (Pa)
        utility_fluid_temp_in_K  — 公用工程流体初始温度 (K)
        utility_fluid_temp_out_K — 公用工程流体最终温度 (K)
        utility_type — "cooling_water" 或 "steam"
        selection — "min" / "mid"（默认）/ "max"

    返回:
        {"heat_transfer_area_m2": 换热面积 (m²)}
    """
    T_in = process_fluid_temp_in_K
    T_out = process_fluid_temp_out_K
    P = process_fluid_pressure_Pa

    # 1) 判定换热模式
    components, mole_fractions, _, _, _ = _parse_molar_flows(
        process_fluid_molar_flows_mol_per_s, T_in, P)
    T_bubble, T_dew = _calc_bubble_dew_pr(components, mole_fractions, P)
    mode, _, is_cooling = _determine_heat_transfer_mode(T_in, T_out, T_bubble, T_dew)

    # 2) 计算 Q
    heat_duty_W = _compute_heat_duty(
        process_fluid_molar_flows_mol_per_s, T_in, T_out, P, mode)

    # 3) 计算 LMTD
    lmtd_K = _compute_log_mean_temp_difference(
        T_in, T_out, utility_fluid_temp_in_K, utility_fluid_temp_out_K)

    # 4) 计算 K
    _, _, k_selected = _compute_overall_k(
        process_fluid_molar_flows_mol_per_s, T_in, T_out, P,
        mode, is_cooling, utility_type, selection)

    # 5) 求 A
    heat_transfer_area_m2 = heat_duty_W / (k_selected * lmtd_K)

    _vprint(f"\n【calc_heat_transfer_area】")
    _vprint(f"  内部计算: Q = {heat_duty_W:.2f} W, LMTD = {lmtd_K:.2f} K, "
            f"K = {k_selected:.2f} W/(m^2*K)")
    _vprint(f"  A = {heat_transfer_area_m2:.4f} m^2")
    return {"heat_transfer_area_m2": heat_transfer_area_m2}


# ═════════════════════════════════════════════
# 自测
# ═════════════════════════════════════════════
if __name__ == "__main__":
    process_fluid_molar_flows_mol_per_s = {
        "OLEIC ACID": 100.0, "LINOLEIC ACID": 150.0, "CAPRYLIC ACID": 180.0,
    }
    process_fluid_temp_in_K = 140 + 273.15
    process_fluid_temp_out_K = 40 + 273.15
    process_fluid_pressure_Pa = 300000.0
    utility_fluid_temp_in_K = 30 + 273.15
    utility_fluid_temp_out_K = 40 + 273.15

    print("=== 热负荷 ===")
    r1 = calc_heat_duty(
        process_fluid_molar_flows_mol_per_s=process_fluid_molar_flows_mol_per_s,
        process_fluid_temp_in_K=process_fluid_temp_in_K,
        process_fluid_temp_out_K=process_fluid_temp_out_K,
        process_fluid_pressure_Pa=process_fluid_pressure_Pa,
    )
    print(r1)

    print("\n=== LMTD ===")
    r2 = calc_log_mean_temp_difference(
        process_fluid_temp_in_K=process_fluid_temp_in_K,
        process_fluid_temp_out_K=process_fluid_temp_out_K,
        utility_fluid_temp_in_K=utility_fluid_temp_in_K,
        utility_fluid_temp_out_K=utility_fluid_temp_out_K,
    )
    print(r2)

    print("\n=== 总传热系数 ===")
    r3 = calc_overall_heat_transfer_coefficient(
        process_fluid_molar_flows_mol_per_s=process_fluid_molar_flows_mol_per_s,
        process_fluid_temp_in_K=process_fluid_temp_in_K,
        process_fluid_temp_out_K=process_fluid_temp_out_K,
        process_fluid_pressure_Pa=process_fluid_pressure_Pa,
        utility_type="cooling_water",
    )
    print(r3)

    print("\n=== 换热面积 ===")
    r4 = calc_heat_transfer_area(
        process_fluid_molar_flows_mol_per_s=process_fluid_molar_flows_mol_per_s,
        process_fluid_temp_in_K=process_fluid_temp_in_K,
        process_fluid_temp_out_K=process_fluid_temp_out_K,
        process_fluid_pressure_Pa=process_fluid_pressure_Pa,
        utility_fluid_temp_in_K=utility_fluid_temp_in_K,
        utility_fluid_temp_out_K=utility_fluid_temp_out_K,
        utility_type="cooling_water",
    )
    print(r4)

    print("\n=== 异常测试：非法 utility_type ===")
    try:
        calc_overall_heat_transfer_coefficient(
            process_fluid_molar_flows_mol_per_s=process_fluid_molar_flows_mol_per_s,
            process_fluid_temp_in_K=process_fluid_temp_in_K,
            process_fluid_temp_out_K=process_fluid_temp_out_K,
            process_fluid_pressure_Pa=process_fluid_pressure_Pa,
            utility_type="hot_oil",
        )
    except ValueError as e:
        print(f"✓ 正确捕获: {e}")

    print("\n=== 异常测试：LMTD 热端温差违法 ===")
    try:
        calc_log_mean_temp_difference(
            process_fluid_temp_in_K=process_fluid_temp_in_K,
            process_fluid_temp_out_K=process_fluid_temp_out_K,
            utility_fluid_temp_in_K=utility_fluid_temp_in_K,
            utility_fluid_temp_out_K=160 + 273.15,  # 公用出口高于工艺入口 → 热端温差为负
        )
    except ValueError as e:
        print(f"✓ 正确捕获: {e}")