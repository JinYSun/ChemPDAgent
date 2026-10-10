"""
精馏工具包 — FUG 法精馏塔设计计算
对应研究内容三："精馏塔工具包（包含 Fenske、Underwood、Gilliland 关联）"

深度集成：
- Fenske 最小理论板数（含多组分非关键组分分布）
- Underwood 最小回流比（含 q 因子详细计算，修复 n_LK 降序索引传递 Bug）
- Gilliland 实际板数（采用 Molokanov 连续关联式，修复 Eduljee 边界阶跃不连续问题）
- 进料板位置（传统 Kirkbride 法）
- 换热器热负荷与水力学校验
"""
from __future__ import annotations
import math
from typing import Optional, List, Dict, Tuple

# 导入统一物性引擎 (确保上游文件存在，不改动任何引用)
from .common import _REQUIRED, check_required_params
from .thermo_helper import (
    get_fluid_density, get_fluid_viscosity, get_fluid_Cp,
    get_fluid_vapor_pressure, get_fluid_MW, get_fluid_surface_tension,
    get_chemical_properties as thermo_get_props,
    calc_bubble_point_T as thermo_bubble_T,
    calc_dew_point_T as thermo_dew_T,
    get_mixture_k_values as thermo_K_values,
)


# ============================================================
# 标准塔径系列 (m)
# ============================================================
_STANDARD_DIAMETERS = [
    0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00,
    1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40,2.50, 2.60, 2.70, 2.80, 2.9, 3.0, 3.20, 3.60, 4.00,
    4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 9.00, 10.0
]

# 标准板间距 (m)
_STANDARD_TRAY_SPACING = [0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80]


# ============================================================
# Fenske 最小理论板数（含多组分非关键组分分布）
# ============================================================

def fenske_minimum_stages(
    alpha_LK_HK: float,
    x_LK_dist: float,
    x_HK_dist: float,
    x_LK_btm: float,
    x_HK_btm: float,
    non_key_components: list | None = None,
    alpha_NK_to_HK: list | None = None,
) -> dict:
    """
    Fenske 最小理论板数 N_min（全回流）
    """
    warnings: list[str] = []

    if alpha_LK_HK <= 1.0:
        return {"N_min": 0, "alpha_LK_HK": alpha_LK_HK, "non_key_distribution": None,
                "valid": False, "warnings": [f"α_LK/HK = {alpha_LK_HK:.3f} ≤ 1，分离不可行"]}

    if any(v <= 0 for v in [x_LK_dist, x_HK_dist, x_LK_btm, x_HK_btm]):
        return {"N_min": 0, "alpha_LK_HK": alpha_LK_HK, "non_key_distribution": None,
                "valid": False, "warnings": ["组成不能为零或负"]}

    # Fenske 方程
    numerator = (x_LK_dist / x_HK_dist) * (x_HK_btm / x_LK_btm)
    if numerator <= 0:
        return {"N_min": 0, "alpha_LK_HK": alpha_LK_HK, "non_key_distribution": None,
                "valid": False, "warnings": [f"Fenske 对数项为非正: {numerator:.4e}"]}

    N_min = math.log(numerator) / math.log(alpha_LK_HK)

    if N_min < 2:
        warnings.append(f"N_min = {N_min:.1f} < 2，分离要求可能过低")

    # --- 多组分非关键组分分布 ---
    nk_dist = None
    if non_key_components and alpha_NK_to_HK and len(non_key_components) == len(alpha_NK_to_HK):
        nk_dist = {}
        for name, alpha_NK in zip(non_key_components, alpha_NK_to_HK):
            if alpha_NK <= 0:
                nk_dist[name] = {"x_dist": 0, "x_bottom": 0, "warning": "α_NK ≤ 0"}
                continue
            # 标准非关键组分分布计算
            alpha_N_min = numerator ** (math.log(alpha_NK) / math.log(alpha_LK_HK))
            d_over_b_approx = alpha_N_min * (x_HK_dist / x_HK_btm)
            d_over_f = d_over_b_approx / (1 + d_over_b_approx) if (1 + d_over_b_approx) > 0 else 0
            
            nk_dist[name] = {
                "alpha_to_HK": alpha_NK,
                "d_over_f": d_over_f,
                "b_over_f": 1 - d_over_f,
                "likely_in_distillate": d_over_f > 0.5,
                "likely_in_bottom": d_over_f <= 0.5,
            }

    return {
        "N_min": N_min,
        "alpha_LK_HK": alpha_LK_HK,
        "non_key_distribution": nk_dist,
        "valid": N_min > 0,
        "warnings": warnings,
    }


# ============================================================
# 物料衡算（支持双组分及多组分混合物）
# ============================================================

def material_balance(
    F: float,
    z_i: list[float],
    x_di: list[float] | None = None,
    x_bi: list[float] | None = None,
    alpha_i: list[float] | None = None,
    LK_index: int | None = None,
    HK_index: int | None = None,
    recovery_LK_dist: float | None = None,
    recovery_HK_dist: float | None = None,
    non_key_method: str = "fenske",
    fenske_N_min: float | None = None,
    component_names: list[str] | None = None,
) -> dict:
    """精馏塔物料衡算（支持双组分及多组分混合物）

    三种计算模式:
      模式A — 全组成指定: 给出 x_di 和 x_bi 所有组分 → 直接衡算
      模式B — 关键组分指定: 给出 LK/HK 的回收率或组成 → Fenske 分布非关键组分
      模式C — 回收率指定: 给出 LK 塔顶回收率 + HK 塔釜回收率 → 解析求解

    Args:
        F: 进料总摩尔流量 (kmol/h 或 mol/s)
        z_i: 进料组成 (摩尔分数列表)
        x_di: 塔顶组成 (摩尔分数列表)，模式A时必须完整给出
        x_bi: 塔底组成 (摩尔分数列表)，模式A时必须完整给出
        alpha_i: 各组分相对挥发度 (以HK为基准)，模式B/C时需要
        LK_index: 轻关键组分在 z_i 中的索引 (从0开始)
        HK_index: 重关键组分在 z_i 中的索引 (从0开始)
        recovery_LK_dist: LK 在塔顶的回收率 (0~1)，模式B/C时使用
        recovery_HK_dist: HK 在塔顶的回收率 (0~1)，模式C时使用
                           注: recovery_HK_dist 通常很小（HK主要去塔釜）
        non_key_method: 非关键组分分配方法
            "fenske" — Fenske 全回流方程分配 (需提供 alpha_i 和 fenske_N_min)
            "sharp"  — 清晰分割 (LNK 全去塔顶, HNK 全去塔釜)
        fenske_N_min: Fenske 最小理论板数 (non_key_method="fenske" 时使用)
        component_names: 组分名称列表 (可选，用于结果标注)

    Returns:
        dict:
          D, B: 塔顶/塔底流量
          d_i, b_i: 各组分塔顶/塔底流量 (列表)
          x_di, x_bi: 归一化后的塔顶/塔底组成 (列表)
          recovery_dist_i, recovery_btm_i: 各组分塔顶/塔底回收率
          mode: 实际使用的计算模式
          warnings: 警告信息
    """
    warnings: list[str] = []
    n = len(z_i)

    if component_names is None:
        component_names = [f"C{i+1}" for i in range(n)]

    # --- 模式A: 全组成指定 ---
    if x_di is not None and x_bi is not None:
        return _material_balance_full_spec(F, z_i, x_di, x_bi, component_names, warnings)

    # --- 模式B/C: 关键组分指定 ---
    if alpha_i is None or LK_index is None or HK_index is None:
        return {"valid": False, "error": "模式B/C需要 alpha_i, LK_index, HK_index",
                "warnings": warnings}

    if LK_index == HK_index:
        return {"valid": False, "error": "LK_index 与 HK_index 不能相同", "warnings": warnings}

    alpha_LK = alpha_i[LK_index]
    alpha_HK = alpha_i[HK_index]
    alpha_LK_HK = alpha_LK / alpha_HK

    if alpha_LK_HK <= 1.0:
        return {"valid": False, "error": f"α_LK/HK = {alpha_LK_HK:.3f} ≤ 1, LK应比HK更易挥发",
                "warnings": warnings}

    # 确定 LK 的塔顶回收率
    if recovery_LK_dist is None:
        # 若未给回收率, 尝试从 x_di 中获取 LK 的组成
        if x_di is not None:
            # D = F*z_LK/xD_LK 近似 (假设 LK 几乎全部去塔顶)
            recovery_LK_dist = 0.99  # 默认 99%
            warnings.append("未指定 recovery_LK_dist, 默认取 99%")
        else:
            recovery_LK_dist = 0.99
            warnings.append("未指定 recovery_LK_dist, 默认取 99%")

    # 确定 HK 的塔顶回收率 (HK 主要去塔釜, 塔顶回收率很小)
    if recovery_HK_dist is None:
        recovery_HK_dist = 0.01  # 默认 1%
        warnings.append("未指定 recovery_HK_dist, 默认取 1% (HK 主要去塔釜)")

    # --- 计算各组分在塔顶的回收率 ---
    recovery_dist = [0.0] * n

    for i in range(n):
        if i == LK_index:
            recovery_dist[i] = recovery_LK_dist
        elif i == HK_index:
            recovery_dist[i] = recovery_HK_dist
        else:
            # 非关键组分
            if non_key_method == "sharp":
                # 清晰分割: 比LK更轻的去塔顶, 比HK更重的去塔釜
                if alpha_i[i] > alpha_LK:
                    recovery_dist[i] = 1.0  # LNK 全部去塔顶
                elif alpha_i[i] < alpha_HK:
                    recovery_dist[i] = 0.0  # HNK 全部去塔釜
                else:
                    # 介于 LK 和 HK 之间 (不应出现, 但保守处理)
                    recovery_dist[i] = 0.5
                    warnings.append(f"组分 {component_names[i]} 的挥发度介于 LK/HK 之间, 非清晰分割")
            elif non_key_method == "fenske":
                # Fenske 全回流分配
                alpha_NK = alpha_i[i] / alpha_HK  # 以 HK 为基准
                if fenske_N_min is None:
                    # 未提供 N_min, 用 LK/HK 的回收率反推
                    # (d_LK/b_LK) / (d_HK/b_HK) = alpha_LK_HK^N_min
                    # recovery_LK/(1-recovery_LK) / (recovery_HK/(1-recovery_HK)) = alpha^N_min
                    ratio_LK = recovery_LK_dist / (1 - recovery_LK_dist)
                    ratio_HK = recovery_HK_dist / (1 - recovery_HK_dist)
                    if ratio_HK > 0 and ratio_LK > 0:
                        N_min_calc = math.log(ratio_LK / ratio_HK) / math.log(alpha_LK_HK)
                        fenske_N_min = N_min_calc
                        warnings.append(f"由 LK/HK 回收率反推 Fenske N_min = {N_min_calc:.2f}")
                    else:
                        fenske_N_min = 5.0
                        warnings.append("无法反推 N_min, 默认取 5.0")

                # Fenske 方程: (d_i/b_i) / (d_HK/b_HK) = (alpha_i/alpha_HK)^N_min
                # d_i/b_i = (d_HK/b_HK) * (alpha_i/alpha_HK)^N_min
                ratio_HK = recovery_HK_dist / (1 - recovery_HK_dist)
                alpha_ratio = alpha_i[i] / alpha_HK
                d_over_b_i = ratio_HK * alpha_ratio ** fenske_N_min
                recovery_dist[i] = d_over_b_i / (1 + d_over_b_i)
            else:
                warnings.append(f"未知 non_key_method: {non_key_method}, 使用 sharp")
                recovery_dist[i] = 1.0 if alpha_i[i] > alpha_LK else 0.0

    # --- 由回收率计算 D 和 B ---
    # F*z_i = d_i + b_i, d_i = recovery_dist_i * F * z_i
    d_i = [recovery_dist[i] * F * z_i[i] for i in range(n)]
    b_i = [F * z_i[i] - d_i[i] for i in range(n)]
    D = sum(d_i)
    B = sum(b_i)

    # --- 归一化组成 ---
    x_di_calc = [d / D if D > 0 else 0 for d in d_i]
    x_bi_calc = [b / B if B > 0 else 0 for b in b_i]

    # --- 校验 ---
    checks = []
    sum_z = sum(z_i)
    if abs(sum_z - 1.0) > 1e-6:
        checks.append({"name": "feed_composition_sum", "pass": False, "value": sum_z, "limit": "1.0"})
    else:
        checks.append({"name": "feed_composition_sum", "pass": True, "value": sum_z, "limit": "1.0"})

    sum_xd = sum(x_di_calc)
    sum_xb = sum(x_bi_calc)
    checks.append({"name": "distillate_composition_sum", "pass": abs(sum_xd - 1.0) < 1e-6,
                   "value": sum_xd, "limit": "1.0"})
    checks.append({"name": "bottoms_composition_sum", "pass": abs(sum_xb - 1.0) < 1e-6,
                   "value": sum_xb, "limit": "1.0"})

    fb_check = abs(D + B - F)
    checks.append({"name": "overall_balance", "pass": fb_check < 1e-6,
                   "value": D + B, "limit": f"F={F}"})

    # --- 组分明细表 ---
    component_table = []
    for i in range(n):
        component_table.append({
            "name": component_names[i],
            "z_i": z_i[i],
            "alpha_i": alpha_i[i] if alpha_i else None,
            "f_i": F * z_i[i],
            "d_i": d_i[i],
            "b_i": b_i[i],
            "x_di": x_di_calc[i],
            "x_bi": x_bi_calc[i],
            "recovery_dist": recovery_dist[i],
            "recovery_btm": 1 - recovery_dist[i],
        })

    # 确定实际使用的模式标签
    if recovery_LK_dist is not None and recovery_HK_dist is not None:
        mode = "C (recovery_spec)" if fenske_N_min is not None else "B (key_spec)"
    else:
        mode = "B (key_spec)"

    return {
        "valid": all(c["pass"] for c in checks),
        "mode": mode,
        "F": F,
        "D": D,
        "B": B,
        "D_over_F": D / F if F > 0 else 0,
        "d_i": d_i,
        "b_i": b_i,
        "x_di": x_di_calc,
        "x_bi": x_bi_calc,
        "recovery_dist_i": recovery_dist,
        "recovery_btm_i": [1 - r for r in recovery_dist],
        "component_table": component_table,
        "alpha_LK_HK": alpha_LK_HK,
        "non_key_method": non_key_method,
        "fenske_N_min_used": fenske_N_min,
        "checks": checks,
        "warnings": warnings,
    }


def _material_balance_full_spec(
    F: float,
    z_i: list[float],
    x_di: list[float],
    x_bi: list[float],
    component_names: list[str],
    warnings: list[str],
) -> dict:
    """模式A: 全组成指定 — 直接物料衡算

    F*z_i = D*xD_i + B*xB_i
    对每个组分: D = F*(z_i - xB_i)/(xD_i - xB_i)
    """
    n = len(z_i)

    # 用每个组分计算 D, 取平均值
    D_values = []
    for i in range(n):
        denom = x_di[i] - x_bi[i]
        if abs(denom) > 1e-12:
            D_calc = F * (z_i[i] - x_bi[i]) / denom
            D_values.append(D_calc)

    if not D_values:
        return {"valid": False, "error": "无法计算 D (所有 xD_i ≈ xB_i)", "warnings": warnings}

    D = sum(D_values) / len(D_values)
    B = F - D

    # 检查各组分计算的 D 一致性
    if len(D_values) > 1:
        D_spread = max(D_values) - min(D_values)
        if D_spread / D > 0.05:
            warnings.append(f"各组分计算的 D 不一致 (最大差异 {D_spread/D*100:.1f}%), 组成数据可能不自洽")

    # 计算各组分流率
    d_i = [D * x_di[i] for i in range(n)]
    b_i = [B * x_bi[i] for i in range(n)]

    # 校验: f_i = d_i + b_i
    balance_errors = []
    for i in range(n):
        f_i = F * z_i[i]
        err = abs(f_i - d_i[i] - b_i[i])
        balance_errors.append(err)

    # 回收率
    recovery_dist = []
    for i in range(n):
        f_i = F * z_i[i]
        if f_i > 1e-12:
            recovery_dist.append(d_i[i] / f_i)
        else:
            recovery_dist.append(0.0)

    # 组分明细表
    component_table = []
    for i in range(n):
        component_table.append({
            "name": component_names[i],
            "z_i": z_i[i],
            "f_i": F * z_i[i],
            "d_i": d_i[i],
            "b_i": b_i[i],
            "x_di": x_di[i],
            "x_bi": x_bi[i],
            "recovery_dist": recovery_dist[i],
            "recovery_btm": 1 - recovery_dist[i],
        })

    checks = [
        {"name": "overall_balance", "pass": abs(D + B - F) < 1e-6, "value": D + B, "limit": f"F={F}"},
        {"name": "distillate_composition_sum", "pass": abs(sum(x_di) - 1.0) < 1e-6,
         "value": sum(x_di), "limit": "1.0"},
        {"name": "bottoms_composition_sum", "pass": abs(sum(x_bi) - 1.0) < 1e-6,
         "value": sum(x_bi), "limit": "1.0"},
    ]
    for i, err in enumerate(balance_errors):
        checks.append({"name": f"component_{component_names[i]}_balance",
                       "pass": err < 1e-6, "value": F*z_i[i], "limit": f"d+b={d_i[i]+b_i[i]}"})

    return {
        "valid": all(c["pass"] for c in checks),
        "mode": "A (full_spec)",
        "F": F,
        "D": D,
        "B": B,
        "D_over_F": D / F if F > 0 else 0,
        "d_i": d_i,
        "b_i": b_i,
        "x_di": list(x_di),
        "x_bi": list(x_bi),
        "recovery_dist_i": recovery_dist,
        "recovery_btm_i": [1 - r for r in recovery_dist],
        "component_table": component_table,
        "D_values_per_component": D_values,
        "checks": checks,
        "warnings": warnings,
    }


# ============================================================
# q 因子计算（进料热状态参数，集成到 Underwood 中）
# ============================================================

def calculate_q_factor(
    feed_temperature: float = 0.0,
    bubble_point: float = 0.0,
    dew_point: float = 0.0,
    Cp_feed: float = 0.0,
    latent_heat_feed: float = 0.0,
    quality: float | None = None,
) -> dict:
    """进料热状态参数 q 计算"""
    warnings: list[str] = []

    if quality is not None:
        q = 1.0 - quality
        feed_state = _describe_q_state(q)
        return {"q": q, "feed_state": feed_state, "bubble_point": bubble_point,
                "dew_point": dew_point, "valid": True, "warnings": warnings}

    if bubble_point <= 0 or dew_point <= 0 or feed_temperature <= 0:
        warnings.append("缺少泡露点或进料温度，默认 q=1（饱和液体进料）")
        return {"q": 1.0, "feed_state": "saturated_liquid", "bubble_point": bubble_point,
                "dew_point": dew_point, "valid": False, "warnings": warnings}

    if feed_temperature < bubble_point:
        if Cp_feed > 0 and latent_heat_feed > 0:
            q = 1.0 + Cp_feed * (bubble_point - feed_temperature) / latent_heat_feed
        else:
            q = 1.05  
            warnings.append("缺少 Cp/潜热数据，使用近似 q=1.05（微过冷）")
    elif abs(feed_temperature - bubble_point) < 0.5:
        q = 1.0
    elif feed_temperature < dew_point:
        if dew_point > bubble_point:
            q = (dew_point - feed_temperature) / (dew_point - bubble_point)
            warnings.append("两相区进料采用线性温度插值计算 q 值，可能与非线性闪蒸计算存在微小偏差")
        else:
            q = 0.5
    elif abs(feed_temperature - dew_point) < 0.5:
        q = 0.0
    else:
        if Cp_feed > 0 and latent_heat_feed > 0:
            q = -Cp_feed * (feed_temperature - dew_point) / latent_heat_feed
        else:
            q = -0.1
            warnings.append("缺少 Cp/潜热数据，使用近似 q=-0.1（微过热）")

    return {
        "q": q,
        "feed_state": _describe_q_state(q),
        "bubble_point": bubble_point,
        "dew_point": dew_point,
        "feed_temperature": feed_temperature,
        "valid": True,
        "warnings": warnings,
    }


def _describe_q_state(q: float) -> str:
    if q > 1.0: return "subcooled_liquid"
    elif q == 1.0: return "saturated_liquid"
    elif q > 0: return "vapor_liquid_mixture"
    elif q == 0.0: return "saturated_vapor"
    else: return "superheated_vapor"


# ============================================================
# Underwood 最小回流比
# ============================================================

def underwood_minimum_reflux(
    alpha_i: list[float],
    z_i: list[float],
    x_di: list[float],
    q: float | None = None,
    feed_temperature: float = 0.0,
    bubble_point: float = 0.0,
    dew_point: float = 0.0,
    Cp_feed: float = 0.0,
    latent_heat_feed: float = 0.0,
    quality: float | None = None,
    n_LK: int = 0,
) -> dict:
    """Underwood 最小回流比 R_min"""
    warnings: list[str] = []

    if q is None:
        q_result = calculate_q_factor(
            feed_temperature=feed_temperature, bubble_point=bubble_point,
            dew_point=dew_point, Cp_feed=Cp_feed, latent_heat_feed=latent_heat_feed,
            quality=quality,
        )
        q_val = q_result["q"]
        q_state = q_result["feed_state"]
        warnings.extend(q_result.get("warnings", []))
    else:
        q_val = q
        q_state = _describe_q_state(q_val)

    target = 1 - q_val
    sorted_pairs = sorted(zip(alpha_i, z_i, x_di), key=lambda x: x[0], reverse=True)
    alphas_sorted = [p[0] for p in sorted_pairs]
    z_sorted = [p[1] for p in sorted_pairs]
    x_d_sorted = [p[2] for p in sorted_pairs]

    # 确定 θ 求解的精确上下界（位于降序排列的重关键 HK 和轻关键 LK 挥发度之间）
    # 修正：通过设计函数传入正确的 n_LK 降序位置，避免因 LNK 存在而错配区间
    if n_LK < len(alphas_sorted) - 1:
        theta_lo = alphas_sorted[n_LK + 1] + 1e-5
        theta_hi = alphas_sorted[n_LK] - 1e-5
    else:
        theta_lo = alphas_sorted[-1] * 0.5
        theta_hi = alphas_sorted[-1] - 1e-5

    theta = _solve_theta(alphas_sorted, z_sorted, target, theta_lo, theta_hi)

    if theta is None:
        warnings.append("Underwood θ 求解失败，使用 R_min = 1.0 回退值")
        return {"R_min": 1.0, "theta": 0.0, "q": q_val, "q_state": q_state,
                "valid": False, "warnings": warnings}

    R_min_plus_1 = sum(
        a * xd / (a - theta) for a, xd in zip(alphas_sorted, x_d_sorted) if abs(a - theta) > 1e-10
    )
    R_min = R_min_plus_1 - 1.0

    if R_min < 0:
        warnings.append(f"R_min = {R_min:.3f} < 0，可能非关键组分大量分配到塔顶")
        R_min = 0.0

    return {
        "R_min": R_min,
        "theta": theta,
        "q": q_val,
        "q_state": q_state,
        "n_components": len(alpha_i),
        "valid": R_min >= 0,
        "warnings": warnings,
    }

def _solve_theta(alphas, z_fracs, target, theta_lo, theta_hi, max_iter=200, tol=1e-8):
    def f(theta):
        return sum(a * z / (a - theta) for a, z in zip(alphas, z_fracs) if abs(a - theta) > 1e-10) - target
    try:
        f_lo = f(theta_lo)
        f_hi = f(theta_hi)
    except ZeroDivisionError:
        return None
    if f_lo * f_hi > 0: return None
    for _ in range(max_iter):
        theta_mid = (theta_lo + theta_hi) / 2
        f_mid = f(theta_mid)
        if abs(f_mid) < tol: return theta_mid
        if f_lo * f_mid < 0: theta_hi = theta_mid
        else: theta_lo = theta_mid; f_lo = f_mid
    return (theta_lo + theta_hi) / 2


# ============================================================
# Gilliland 关联（含进料板位置 Kirkbride 法）
# ============================================================

def gilliland_correlation(
    N_min: float,
    R_min: float,
    R_actual: float | None = None,
    R_ratio: float | None = None,
    x_LK_dist: float = 0.0,
    x_HK_dist: float = 0.0,
    x_LK_feed: float = 0.0,
    x_HK_feed: float = 0.0,
    x_LK_btm: float = 0.0,
    D: float = 0.0,
    B: float = 0.0,
    N_stripping: float | None = None,
    feed_stage_method: str = "kirkbride",
) -> dict:
    """Gilliland 关联与 Kirkbride 进料位置

    Kirkbride 公式（标准形式）:
        Nr/Ns = [(B/D) × (zHK/zLK) × (xB_LK/xD_HK)²]^0.206
    需要 D, B（塔顶/塔底流量）和 x_LK_btm（塔底 LK 摩尔分数）。
    若 D/B 未提供，降级为简化公式。
    """
    warnings: list[str] = []

    if R_actual is not None:
        R_act = R_actual
    elif R_ratio is not None:
        R_act = R_ratio * R_min
    else:
        R_act = 1.5 * R_min
        warnings.append(f"未指定回流比，使用 R = 1.5 × R_min = {R_act:.3f}")

    R_ratio_val = R_act / R_min if R_min > 0 else 1.5

    if R_ratio_val <= 1.0:
        return {"N_actual": float("inf"), "R_actual": R_act, "R_ratio": R_ratio_val,
                "X": 0, "Y": 0, "N_rectifying": None, "N_stripping": None,
                "feed_stage": None, "valid": False,
                "warnings": [f"R/R_min = {R_ratio_val:.3f} ≤ 1，无法操作"]}

    X = (R_act - R_min) / (R_act + 1)
    
    # 修正：采用 Molokanov 关联式取代 Eduljee，解决 R -> R_min 处的板数不连续与阶跃性失效问题
    if X <= 1e-7:
        N_actual = float("inf")
        Y_use = 1.0
    else:
        Y_use = 1.0 - math.exp((1.0 + 54.4 * X) / (11.0 + 117.2 * X) * (X - 1.0) / math.sqrt(X))
        if (1 - Y_use) > 1e-10:
            N_actual = (N_min + Y_use) / (1 - Y_use)
        else:
            N_actual = float("inf")

    if N_actual != float("inf"):
        N_actual = max(N_min + 1, math.ceil(N_actual))

    # --- 进料板位置 ---
    N_rect = None
    N_strip = N_stripping
    feed_stage = None

    if N_stripping is not None:
        N_rect = N_actual - N_stripping
        feed_stage = N_rect + 1
    elif all(v > 0 for v in [x_LK_dist, x_HK_dist, x_LK_feed, x_HK_feed]):
        if feed_stage_method == "kirkbride" and N_actual != float("inf"):
            # Kirkbride 方程: Nr/Ns = [(B/D)*(zHK/zLK)*(xB_LK/xD_HK)^2]^0.206
            if D > 0 and B > 0 and x_LK_btm > 0:
                kirkbride_base = (B / D) * (x_HK_feed / x_LK_feed) * (x_LK_btm / x_HK_dist) ** 2
            else:
                # D/B 未提供时降级为简化形式
                kirkbride_base = (x_HK_feed / x_LK_feed) * (x_LK_btm / x_HK_dist if x_LK_btm > 0 else 1.0) ** 2
                if D <= 0 or B <= 0:
                    warnings.append("Kirkbride 公式缺少 D/B 数据，使用简化形式（可能偏差较大）")
            Nr_Ns_ratio = kirkbride_base ** 0.206
            N_strip = max(2, int(N_actual / (1 + Nr_Ns_ratio)))
            N_rect = N_actual - N_strip
            feed_stage = N_rect + 1
        else:
            if N_actual != float("inf"):
                N_strip = max(2, int(N_actual * 0.5))
                N_rect = N_actual - N_strip
                feed_stage = N_rect + 1
    else:
        if N_actual != float("inf"):
            N_rect = max(1, int(N_actual * 0.5))
            N_strip = N_actual - N_rect
            feed_stage = N_rect + 1
            warnings.append("缺少组成数据，进料板位置使用默认中间位置")

    return {
        "N_actual": N_actual,
        "R_actual": R_act,
        "R_ratio": R_ratio_val,
        "X": X,
        "Y": Y_use,
        "N_rectifying": N_rect,
        "N_stripping": N_strip,
        "feed_stage": feed_stage,
        "valid": N_actual > N_min if N_actual != float("inf") else False,
        "warnings": warnings,
    }


# ============================================================
# 塔径计算（Souders-Brown 法）
# ============================================================

# 净面积分率：A_net = A_col × NET_AREA_FRACTION（扣除降液管后的净截面积/全截面积）
NET_AREA_FRACTION = 0.90


def calculate_column_diameter(
    V: float,
    rho_vapor: float,
    rho_liquid: float,
    sigma: float = 0.020,
    tray_spacing: float = 0.60,
    hole_area_fraction: float = 0.10,
    flooding_fraction: float = 0.80,
    L_mass: float | None = None,
    flood_method: str = "simple",
) -> dict:
    """塔径计算（Souders-Brown 法，净面积基准）

    u_max/u_design 以净面积 A_net = 0.90×A_col 为基准（经典 Souders-Brown/Fair 定义，
    扣除降液管占用的截面），塔径由 A_col = A_net/NET_AREA_FRACTION 反推。

    flood_method:
    - "simple": 简化 Csb = 0.002·√(H_t·1000)·(A_h/0.10)^0.2（不含 F_LV 修正）
    - "fair":   Fair (1961) 关联式（Perry's §14 图回归），需提供 L_mass 以计算 F_LV
    """
    warnings: list[str] = []

    F_LV_val = None
    if flood_method == "fair":
        if L_mass is not None and V > 0:
            F_LV_val = (L_mass / V) * math.sqrt(rho_vapor / rho_liquid)
        else:
            warnings.append("fair 法需要 L_mass（液相质量流量）计算 F_LV，回退到简化法")
            flood_method = "simple"

    if flood_method == "fair":
        fair = fair_capacity_factor(tray_spacing, F_LV_val, sigma, hole_area_fraction)
        C_sb = fair["Csb"]
        warnings.extend(fair.get("warnings", []))
        C_base = fair["Csb0"]
        sigma_factor = fair["sigma_factor"]
    else:
        C_base = _calculate_C_base(tray_spacing, hole_area_fraction)
        sigma_factor = (sigma / 0.020) ** 0.2 if sigma > 0 else 1.0
        C_sb = C_base * sigma_factor

    if rho_liquid > rho_vapor and rho_vapor > 0:
        u_max = C_sb * math.sqrt((rho_liquid - rho_vapor) / rho_vapor)
    else:
        warnings.append("液相密度必须大于气相密度")
        return {"diameter": 0, "diameter_calculated": 0, "u_max": 0, "u_design": 0,
                "C_factor": 0, "flooding_fraction": flooding_fraction,
                "valid": False, "warnings": warnings}

    u_design = u_max * flooding_fraction

    Q_vapor = V / rho_vapor if rho_vapor > 0 else 0.0  

    if u_design > 0 and Q_vapor > 0:
        # 液泛气速基准为净面积: A_net = Q_vapor / u_design，A_col = A_net / 净面积分率
        A_net = Q_vapor / u_design
        A_cross = A_net / NET_AREA_FRACTION
        d_calc = math.sqrt(4 * A_cross / math.pi)
    else:
        d_calc = 0.0

    d_standard = _round_to_standard_diameter(d_calc)

    if flooding_fraction > 0.85:
        warnings.append(f"液泛分率 {flooding_fraction:.0%} > 85%，操作弹性不足")
    elif flooding_fraction < 0.60:
        warnings.append(f"液泛分率 {flooding_fraction:.0%} < 60%，塔径偏大")

    return {
        "diameter": d_standard,
        "diameter_calculated": d_calc,
        "u_max": u_max,
        "u_design": u_design,
        "C_factor": C_sb,
        "C_base": C_base,
        "sigma_factor": sigma_factor,
        "F_LV": F_LV_val,
        "flood_method": flood_method,
        "net_area_fraction": NET_AREA_FRACTION,
        "area_basis": "u_max/u_design 基于净面积 A_net = 0.90×A_col；塔径由 A_col = A_net/0.90 反推",
        "flooding_fraction": flooding_fraction,
        "tray_spacing": tray_spacing,
        "valid": d_standard > 0,
        "warnings": warnings,
    }

def _calculate_C_base(tray_spacing: float, hole_area_fraction: float = 0.10) -> float:
    if tray_spacing < 0.30: H_t_eff = 0.30
    elif tray_spacing > 0.90: H_t_eff = 0.90
    else: H_t_eff = tray_spacing

    C_base = 0.002 * (H_t_eff * 1000) ** 0.5 
    if hole_area_fraction > 0:
        C_base *= (hole_area_fraction / 0.10) ** 0.2
    return C_base


def fair_capacity_factor(
    tray_spacing: float,
    F_LV: float,
    sigma: float = 0.020,
    hole_area_fraction: float = 0.10,
) -> dict:
    """Fair (1961) 液泛容量因子 Csb（Perry's §14 通用关联式，图表回归形式）

    Args:
        tray_spacing: 板间距 (m)
        F_LV: 气液两相流参数 (L/V)·√(ρV/ρL)，L/V 为液气质量流量比
        sigma: 表面张力 (N/m)
        hole_area_fraction: 开孔面积/有效面积（图值基于 ≥0.10）

    Csb₀ = 0.0105 + 8.127e-4·T_s^0.775·exp(-1.463·F_LV^0.842)  [T_s: mm, Csb: m/s]
    Csb = Csb₀·(σ/20 mN/m)^0.2
    适用: F_LV 0.01–1.0, 非起泡体系, 开孔面积 ≥ 10% 有效面积, 堰高 ≤ 15% 板间距
    """
    warnings: list[str] = []
    T_s = max(0.15, min(0.90, tray_spacing)) * 1000.0   # mm, Fair 图范围 150–900 mm
    F_LV = max(F_LV, 1e-4)
    Csb0 = 0.0105 + 8.127e-4 * T_s ** 0.775 * math.exp(-1.463 * F_LV ** 0.842)
    sigma_factor = (sigma * 1000.0 / 20.0) ** 0.2 if sigma > 0 else 1.0
    Csb = Csb0 * sigma_factor
    if hole_area_fraction < 0.10:
        warnings.append(f"开孔面积 {hole_area_fraction:.0%} < 10% 有效面积，Fair 图值需另行降低")
    if not (0.01 <= F_LV <= 1.0):
        warnings.append(f"F_LV = {F_LV:.4f} 超出 Fair 关联式适用范围 0.01–1.0，结果仅作参考")
    return {
        "Csb": Csb, "Csb0": Csb0, "F_LV": F_LV, "sigma_factor": sigma_factor,
        "correlation": "Fair (1961), Perry's §14 chart regression",
        "warnings": warnings,
    }

def _round_to_standard_diameter(d: float) -> float:
    for sd in _STANDARD_DIAMETERS:
        if sd >= d: return sd
    return math.ceil(d * 2) / 2


# ============================================================
# 塔板效率估算（O'Connell 关联）
# ============================================================

def estimate_tray_efficiency(
    alpha_LK_HK: float,
    mu_feed: float,
    temperature: float = 298.15,
) -> dict:
    """塔板效率估算（O'Connell 关联）
    
    Args:
        alpha_LK_HK: 轻关键组分对重关键组分的相对挥发度
        mu_feed: 进料液相黏度 (Pa·s)，内部自动转换为 cP
        temperature: 温度 (K)
    
    O'Connell 关联式: E_o = 0.492 × (α×μ)^(-0.245)
    其中 μ 单位为 cP (厘泊), 1 Pa·s = 1000 cP
    """
    mu_cP = mu_feed * 1000.0
    alpha_mu = alpha_LK_HK * mu_cP

    if alpha_mu <= 0:
        return {"efficiency": 0.5, "efficiency_percent": 50.0, "alpha_mu_product": alpha_mu,
                "correlation": "O'Connell (fallback)", "valid": False,
                "warning": "α×μ ≤ 0，使用默认效率 50%"}

    E_o = 0.492 * alpha_mu ** (-0.245)
    E_o = max(0.15, min(0.95, E_o)) 

    warning = None
    if alpha_mu < 0.1: warning = f"α×μ = {alpha_mu:.4f} < 0.1，超出关联范围"
    elif alpha_mu > 100: warning = f"α×μ = {alpha_mu:.1f} > 100，超出关联范围"

    return {
        "efficiency": E_o, "efficiency_percent": E_o * 100, "alpha_mu_product": alpha_mu,
        "alpha_LK_HK": alpha_LK_HK, "mu_feed_cP": mu_feed * 1000,
        "correlation": "O'Connell 1946", "valid": 0.1 <= alpha_mu <= 100, "warning": warning,
    }


# ============================================================
# 塔板水力学校验
# ============================================================

def calculate_tray_hydraulics(
    V_mass: float,
    L_mass: float,
    rho_vapor: float,
    rho_liquid: float,
    sigma: Optional[float] = None,
    mu_liquid: Optional[float] = None,
    D_column: float = 1.0,
    tray_spacing: float = 0.60,
    hole_diameter: float = 0.005,
    hole_area_fraction: float = 0.10,
    weir_height: float = 0.050,
    weir_length: float | None = None,
    flood_method: str = "simple",
) -> dict:
    """塔板水力学校验
    
    sigma / mu_liquid 不提供时:
    - C_sb 计算跳过表面张力修正
    - 干板压降仍可计算（仅依赖气相密度和孔速）

    flood_method: "simple"（简化式）/ "fair"（Fair 1961 关联式，含 F_LV 修正）
    """
    warnings: list[str] = []

    if weir_length is None:
        weir_length = 0.7 * D_column

    A_column = math.pi / 4 * D_column ** 2
    A_net = A_column * NET_AREA_FRACTION   # 净面积（扣除降液管）
    A_holes = A_column * hole_area_fraction

    # 液泛校核气速以净面积为基准（与 Souders-Brown/Fair 关联式定义一致）
    u_vapor = V_mass / (rho_vapor * A_net) if (rho_vapor > 0 and A_net > 0) else 0.0
    u_hole = V_mass / (rho_vapor * A_holes) if (rho_vapor > 0 and A_holes > 0) else 0.0

    F_factor = u_vapor * math.sqrt(rho_vapor) if rho_vapor > 0 else 0.0
    F_hole = u_hole * math.sqrt(rho_vapor) if rho_vapor > 0 else 0.0

    # C_sb 计算
    F_LV_val = None
    if flood_method == "fair" and V_mass > 0 and L_mass is not None:
        F_LV_val = (L_mass / V_mass) * math.sqrt(rho_vapor / rho_liquid)
        fair = fair_capacity_factor(tray_spacing, F_LV_val, sigma if sigma else 0.020, hole_area_fraction)
        C_sb = fair["Csb"]
        warnings.extend(fair.get("warnings", []))
    else:
        # 简化法：仅在提供 sigma 时施加表面张力修正
        C_base_raw = _calculate_C_base(tray_spacing, hole_area_fraction)
        if sigma is not None and sigma > 0:
            sigma_factor = (sigma / 0.020) ** 0.2
            C_sb = C_base_raw * sigma_factor
        else:
            C_sb = C_base_raw  # 无表面张力修正

    u_flood = C_sb * math.sqrt((rho_liquid - rho_vapor) / rho_vapor) if rho_vapor > 0 else 0.0
    flooding_pct = (u_vapor / u_flood * 100) if u_flood > 0 else 0.0

    C_d = 0.73  
    if rho_vapor > 0 and A_holes > 0:
        dP_dry = (1.0 / C_d ** 2) * (rho_vapor / 2) * u_hole ** 2
    else:
        dP_dry = 0.0

    L_per_weir = L_mass / (rho_liquid * weir_length) if (rho_liquid > 0 and weir_length > 0) else 0.0
    h_ow = 0.666 * (L_per_weir) ** (2 / 3) if L_per_weir > 0 else 0.0  
    
    # 单板水力学考虑筛板上的充气系数 (Aeration Factor beta，典型值 0.6)
    beta = 0.6  
    h_cl = beta * (weir_height + h_ow)  

    dP_wet = rho_liquid * 9.81 * h_cl  
    dP_total = dP_dry + dP_wet

    F_hole_min = 5.0  
    weeping = F_hole < F_hole_min

    weeping_check = {
        "F_hole": F_hole, "F_hole_min": F_hole_min, "weeping": weeping,
        "status": "OK" if not weeping else "LEAKING - reduce hole area or increase vapor load",
    }

    if flooding_pct > 85: warnings.append(f"液泛分率 {flooding_pct:.1f}% > 85%，需增大塔径")
    elif flooding_pct < 50: warnings.append(f"液泛分率 {flooding_pct:.1f}% < 50%，塔径偏大")
    if weeping: warnings.append(f"筛孔 F={F_hole:.1f} < {F_hole_min}，漏液风险")
    if dP_total > 1500: warnings.append(f"单板压降 {dP_total:.0f} Pa 偏高")

    return {
        "flooding_percent": flooding_pct, "pressure_drop_per_tray": dP_total,
        "pressure_drop_dry": dP_dry, "pressure_drop_wet": dP_wet,
        "clear_liquid_height": h_cl, "h_ow": h_ow,
        "weeping_check": weeping_check,
        "vapor_velocity": u_vapor, "hole_velocity": u_hole,
        "F_factor": F_factor, "F_hole": F_hole,
        "C_sb": C_sb, "u_flood": u_flood, "F_LV": F_LV_val, "flood_method": flood_method,
        "A_column": A_column, "A_net": A_net, "A_holes": A_holes,
        "area_basis": f"u_vapor/u_flood/FF 基于净面积 A_net = {NET_AREA_FRACTION:.2f}×A_col；u_hole 基于开孔面积",
        "valid": 50 <= flooding_pct <= 85, "warnings": warnings,
    }


# ============================================================
# 筛板塔综合校核（液泛 + 夹带 + 压降 + 漏液）
# ============================================================

def check_sieve_tray(
    D_column: float,
    tray_spacing: float,
    hole_diameter: float,
    hole_area_fraction: float,
    rho_vapor: float,
    rho_liquid: float,
    Q_vapor: float,
    sigma: Optional[float] = None,
    mu_liquid: Optional[float] = None,
    L_mass: Optional[float] = None,
    weir_height: float = 0.050,
) -> dict:
    """筛板精馏塔综合校核

    对已有塔进行液泛、雾沫夹带、压降、漏液四项校核。
    sigma / mu_liquid / L_mass 未提供时自动降级，不使用默认值。

    Args:
        D_column: 塔径 (m)
        tray_spacing: 板间距 (m)
        hole_diameter: 孔径 (m)
        hole_area_fraction: 开孔率 (-)
        rho_vapor: 气相密度 (kg/m3)
        rho_liquid: 液相密度 (kg/m3)
        Q_vapor: 气相体积流量 (m3/s)
        sigma: 表面张力 (N/m)，可选；不提供则夹带校核回退 Souders-Brown 法
        mu_liquid: 液相黏度 (Pa.s)，可选
        L_mass: 液相质量流量 (kg/s)，可选；不提供时跳过堰上液层高度计算
        weir_height: 堰高 (m)

    Returns:
        综合校核 dict，包含 flooding / entrainment / pressure_drop / weeping 子结果
    """
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    D_column = safe_float(D_column, name='D_column')
    tray_spacing = safe_float(tray_spacing, name='tray_spacing')
    hole_diameter = safe_float(hole_diameter, name='hole_diameter')
    hole_area_fraction = safe_float(hole_area_fraction, name='hole_area_fraction')
    rho_vapor = safe_float(rho_vapor, name='rho_vapor')
    rho_liquid = safe_float(rho_liquid, name='rho_liquid')
    Q_vapor = safe_float(Q_vapor, name='Q_vapor')
    sigma = safe_float(sigma, name='sigma')
    mu_liquid = safe_float(mu_liquid, name='mu_liquid')
    L_mass = safe_float(L_mass, name='L_mass')
    weir_height = safe_float(weir_height, default=0.050, name='weir_height')
    
    from .flash_drum import entrainment_fraction as _entrainment_fraction

    warnings: list[str] = []

    # ---- 基本气速参数 ----
    A_cross = math.pi / 4 * D_column ** 2
    A_net = A_cross * NET_AREA_FRACTION   # 净面积（扣除降液管）
    A_holes = A_cross * hole_area_fraction
    V_mass = Q_vapor * rho_vapor  # kg/s

    u_superficial = Q_vapor / A_net   # 液泛校核气速以净面积为基准
    u_hole = Q_vapor / A_holes if A_holes > 0 else 0.0
    F_factor = u_superficial * math.sqrt(rho_vapor)
    F_hole = u_hole * math.sqrt(rho_vapor)

    velocity_params = {
        "A_cross": round(A_cross, 4),
        "A_net": round(A_net, 4),
        "A_holes": round(A_holes, 4),
        "u_superficial": round(u_superficial, 4),
        "u_hole": round(u_hole, 2),
        "F_factor": round(F_factor, 3),
        "F_hole": round(F_hole, 2),
    }

    # ---- 液泛校核 (Souders-Brown 法) ----
    C_base_raw = _calculate_C_base(tray_spacing, hole_area_fraction)
    if sigma is not None and sigma > 0:
        sigma_factor = (sigma / 0.020) ** 0.2
        C_sb = C_base_raw * sigma_factor
    else:
        C_sb = C_base_raw
        if sigma is None:
            warnings.append("未提供表面张力，C_sb 未经表面张力修正")

    u_flood = C_sb * math.sqrt((rho_liquid - rho_vapor) / rho_vapor) if rho_vapor > 0 else 0.0
    flood_ratio = u_superficial / u_flood if u_flood > 0 else 0.0

    flooding = {
        "C_sb": round(C_sb, 4),
        "u_flood": round(u_flood, 4),
        "u_superficial": round(u_superficial, 4),
        "flood_ratio": round(flood_ratio, 4),
        "flood_percent": round(flood_ratio * 100, 1),
        "pass": flood_ratio <= 0.85,
        "note": (f"液泛分率 {flood_ratio*100:.1f}%"
                 + (" <= 85%，合格" if flood_ratio <= 0.85
                    else " >> 85%，已液泛!" if flood_ratio > 1.0
                    else "，超标")),
    }
    if flood_ratio > 0.85:
        warnings.append(f"液泛分率 {flood_ratio*100:.1f}% > 85%，塔已液泛或接近液泛")

    # ---- 雾沫夹带校核 ----
    ent_result = _entrainment_fraction(
        vapor_velocity=u_superficial,
        surface_tension=sigma,
        liquid_viscosity=mu_liquid,
        gas_density=rho_vapor,
        liquid_density=rho_liquid,
        separator_diameter=D_column,
    )
    ent_method = ent_result.get("method", "unknown")
    ent_frac = ent_result.get("entrainment_fraction")
    entrainment = {
        "method": ent_method,
        "entrainment_fraction": ent_frac,
        "velocity_ratio": ent_result.get("velocity_ratio"),
        "critical_droplet_diameter_um": ent_result.get("critical_droplet_diameter_um"),
        "pass": ent_result.get("velocity_ratio", 999) < 0.8 if ent_method == "Souders-Brown" else None,
        "note": ent_result.get("note", ""),
    }
    if ent_frac is not None and ent_frac >= 0.01:
        warnings.append(f"夹带分率 {ent_frac*100:.2f}% >= 1%，超标")

    # ---- 压降校核 ----
    C_d = 0.73
    dP_dry = (1.0 / C_d ** 2) * (rho_vapor / 2) * u_hole ** 2 if A_holes > 0 else 0.0

    # 堰上液层高度（需 L_mass）
    h_ow = None
    h_cl = None
    dP_wet = None
    if L_mass is not None and L_mass > 0:
        weir_length = 0.7 * D_column
        L_per_weir = L_mass / (rho_liquid * weir_length) if (rho_liquid > 0 and weir_length > 0) else 0.0
        h_ow = 0.666 * L_per_weir ** (2 / 3) if L_per_weir > 0 else 0.0
        beta = 0.6
        h_cl = beta * (weir_height + h_ow)
        dP_wet = rho_liquid * 9.81 * h_cl
    else:
        warnings.append("未提供 L_mass，跳过堰上液层高度和湿板压降计算")

    dP_total = dP_dry + (dP_wet if dP_wet is not None else 0.0)
    pressure_drop = {
        "dP_dry": round(dP_dry, 1),
        "dP_wet": round(dP_wet, 1) if dP_wet is not None else None,
        "dP_total": round(dP_total, 1),
        "h_ow": round(h_ow * 1000, 1) if h_ow is not None else None,
        "h_cl": round(h_cl * 1000, 1) if h_cl is not None else None,
        "pass": dP_total < 1500,
        "note": f"单板压降 {dP_total:.1f} Pa" + (" < 1500 Pa，合格" if dP_total < 1500 else " >= 1500 Pa，偏高"),
    }
    if dP_total >= 1500:
        warnings.append(f"单板压降 {dP_total:.0f} Pa >= 1500 Pa")

    # ---- 漏液校核 ----
    F_hole_min = 5.0
    weeping = {
        "F_hole": round(F_hole, 2),
        "F_hole_min": F_hole_min,
        "pass": F_hole >= F_hole_min,
        "note": f"F_hole = {F_hole:.2f}" + (" >= 5，无漏液" if F_hole >= F_hole_min else " < 5，漏液风险"),
    }
    if F_hole < F_hole_min:
        warnings.append(f"筛孔 F_factor = {F_hole:.1f} < {F_hole_min}，漏液风险")

    # ---- 综合判定 ----
    items = {
        "flooding": flooding["pass"],
        "entrainment": entrainment["pass"],
        "pressure_drop": pressure_drop["pass"],
        "weeping": weeping["pass"],
    }
    n_pass = sum(1 for v in items.values() if v)
    n_total = len(items)
    overall_pass = n_pass == n_total

    return {
        "velocity_params": velocity_params,
        "flooding": flooding,
        "entrainment": entrainment,
        "pressure_drop": pressure_drop,
        "weeping": weeping,
        "overall_pass": overall_pass,
        "n_pass": n_pass,
        "n_total": n_total,
        "warnings": warnings,
    }


# ============================================================
# 完整 FUG 精馏塔设计函数
# ============================================================

def design_distillation_column(
    alpha_LK_HK: float,
    x_LK_dist: float,
    x_HK_dist: float,
    x_LK_btm: float,
    x_HK_btm: float,
    z_i: list[float] | None = None,
    alpha_i: list[float] | None = None,
    x_di: list[float] | None = None,
    q: float | None = None,
    feed_temperature: float = 0.0,
    bubble_point: float = 0.0,
    dew_point: float = 0.0,
    Cp_feed: float = 0.0,
    latent_heat_feed: float = 0.0,
    quality: float | None = None,
    R_actual: float | None = None,
    R_ratio: float  = 1.2,
    non_key_components: list | None = None,
    alpha_NK_to_HK: list | None = None,
    x_LK_feed: float = 0.0,
    x_HK_feed: float = 0.0,
    V_mass: float = _REQUIRED,
    L_mass: float = _REQUIRED,
    rho_vapor: float = _REQUIRED,
    rho_liquid: float = _REQUIRED,
    sigma: float = 0.020,
    mu_feed: float = 0.001,
    tray_spacing: float = _REQUIRED,
    hole_area_fraction: float = 0.10,
    flooding_fraction: float = _REQUIRED,
    weir_height: float = 0.050,
    hole_diameter: float = 0.005,
    tray_efficiency: float | None = None,
    P_top: float | None = None,
    flood_method: str = "simple",
) -> dict:
    """完整 FUG 精馏塔设计计算

    板数口径说明：
    - N_theoretical（Gilliland 理论级数）含再沸器 1 级、不含全凝器；
    - N_actual_trays 为塔内物理塔板数 = ceil((N_theoretical - 1) / 塔板效率)；
    - feed_stage 为理论进料板编号（自塔内第一块板起计），
      feed_stage_actual = ceil((feed_stage - 1) / 塔板效率) + 1（统一效率近似换算）；
    - N_rectifying_actual / N_stripping_actual 均不含进料板，
      满足 N_rectifying_actual + 1 + N_stripping_actual = N_actual_trays。
    """
    missing = check_required_params(locals(), {
        "V_mass": "塔内气相质量流量 (kg/s)，由设备层计算提供",
        "L_mass": "塔内液相质量流量 (kg/s)，由设备层计算提供",
        "rho_vapor": "气相密度 (kg/m³)，影响塔径计算",
        "rho_liquid": "液相密度 (kg/m³)，影响塔径和液泛校核",
        "tray_spacing": "塔板间距 (m)，影响塔径和液泛校核",
        "flooding_fraction": "目标液泛分率 (0.5~0.85)，影响塔径校核",
    })
    if missing:
        return missing

    all_warnings: list[str] = []

    fenske = fenske_minimum_stages(
        alpha_LK_HK=alpha_LK_HK, x_LK_dist=x_LK_dist, x_HK_dist=x_HK_dist,
        x_LK_btm=x_LK_btm, x_HK_btm=x_HK_btm,
        non_key_components=non_key_components, alpha_NK_to_HK=alpha_NK_to_HK,
    )
    all_warnings.extend(fenske.get("warnings", []))
    if not fenske["valid"]:
        return {"valid": False, "error": "Fenske 计算失败", "details": fenske, "warnings": all_warnings}

    n_LK_val = 0
    if z_i and alpha_i and x_di:
        uw_alpha, uw_z, uw_xd = alpha_i, z_i, x_di
        sorted_alphas_temp = sorted(alpha_i, reverse=True)
        n_LK_val = min(range(len(sorted_alphas_temp)), key=lambda idx: abs(sorted_alphas_temp[idx] - alpha_LK_HK))
    else:
        uw_alpha = [alpha_LK_HK, 1.0]
        uw_z = [x_LK_feed or 0.5, x_HK_feed or 0.5]
        uw_xd = [x_LK_dist, x_HK_dist]
        n_LK_val = 0

    underwood = underwood_minimum_reflux(
        alpha_i=uw_alpha, z_i=uw_z, x_di=uw_xd, q=q, 
        feed_temperature=feed_temperature, bubble_point=bubble_point, dew_point=dew_point, 
        Cp_feed=Cp_feed, latent_heat_feed=latent_heat_feed, quality=quality, n_LK=n_LK_val,
    )
    all_warnings.extend(underwood.get("warnings", []))

    # 由物料衡算精确计算 D 和 B，用于 Kirkbride 进料板位置
    D_val = 0.0
    B_val = 0.0
    if z_i and len(z_i) > 0:
        # 确定 LK 在 z_i 中的索引（α_i 降序排列，LK 位于 n_LK_val）
        lk_idx = n_LK_val if n_LK_val < len(z_i) else 0
        z_LK = z_i[lk_idx]
        # 由 LK 组分物料衡算: F*z_LK = D*x_LK_dist + B*x_LK_btm, F = D + B
        # → D = F*(z_LK - x_LK_btm) / (x_LK_dist - x_LK_btm)
        F_est = 1.0  # 归一化基准
        denom = x_LK_dist - x_LK_btm
        if abs(denom) > 1e-10:
            D_val = F_est * (z_LK - x_LK_btm) / denom
            B_val = F_est - D_val
            D_val = max(0.0, min(D_val, F_est))  # 限制在 [0, F]
            B_val = F_est - D_val
        else:
            # x_LK_dist ≈ x_LK_btm 时无法由 LK 衡算，降级为估算
            for i in range(len(z_i)):
                if alpha_i and alpha_i[i] >= alpha_LK_HK:
                    D_val += F_est * z_i[i] * 0.99
                elif alpha_i and alpha_i[i] <= 1.0:
                    D_val += F_est * z_i[i] * 0.01
                else:
                    D_val += F_est * z_i[i] * 0.5
            B_val = F_est - D_val

    gilliland = gilliland_correlation(
        N_min=fenske["N_min"], R_min=underwood["R_min"], R_actual=R_actual, R_ratio=R_ratio,
        x_LK_dist=x_LK_dist, x_HK_dist=x_HK_dist,
        x_LK_feed=x_LK_feed or (z_i[0] if z_i else 0.0),
        x_HK_feed=x_HK_feed or (z_i[-1] if z_i and len(z_i) > 1 else 0.0),
        x_LK_btm=x_LK_btm,
        D=D_val, B=B_val,
    )
    all_warnings.extend(gilliland.get("warnings", []))

    if tray_efficiency is not None and 0 < tray_efficiency <= 1:
        efficiency = {"efficiency": tray_efficiency, "efficiency_percent": tray_efficiency * 100, "method": "user_specified"}
    else:
        efficiency = estimate_tray_efficiency(alpha_LK_HK, mu_feed)
        efficiency["method"] = "oconnell"
    N_theoretical = gilliland["N_actual"]
    
    if N_theoretical == float("inf"):
        N_actual_trays = float("inf")
        feed_stage_actual = None
        N_rectifying_actual = None
        N_stripping_actual = None
    else:
        # N_theoretical 含再沸器 1 级：先扣除再沸器，再除以塔板效率 → 塔内物理塔板数
        N_actual_trays = max(1, int(math.ceil((N_theoretical - 1) / efficiency["efficiency"])))
        # 实际进料板位置：理论进料板编号自塔内第一块板起计，
        # 按统一效率换算进料板以上的板数（近似）
        feed_stage_theo = gilliland.get("feed_stage")
        if feed_stage_theo:
            N_rectifying_actual = max(1, int(math.ceil((feed_stage_theo - 1) / efficiency["efficiency"])))
            feed_stage_actual = N_rectifying_actual + 1
            N_stripping_actual = N_actual_trays - feed_stage_actual
        else:
            feed_stage_actual = None
            N_rectifying_actual = None
            N_stripping_actual = None

    diameter = calculate_column_diameter(
        V=V_mass, rho_vapor=rho_vapor, rho_liquid=rho_liquid, sigma=sigma,
        tray_spacing=tray_spacing, hole_area_fraction=hole_area_fraction, flooding_fraction=flooding_fraction,
        L_mass=L_mass, flood_method=flood_method,
    )
    all_warnings.extend(diameter.get("warnings", []))

    hydraulics = calculate_tray_hydraulics(
        V_mass=V_mass, L_mass=L_mass, rho_vapor=rho_vapor, rho_liquid=rho_liquid,
        sigma=sigma, mu_liquid=mu_feed, D_column=diameter["diameter"],
        tray_spacing=tray_spacing, hole_diameter=hole_diameter,
        hole_area_fraction=hole_area_fraction, weir_height=weir_height,
        flood_method=flood_method,
    )
    all_warnings.extend(hydraulics.get("warnings", []))

    if N_actual_trays != float("inf"):
        H_column = (N_actual_trays * tray_spacing) + 2 * tray_spacing + 3.0  
        H_D_ratio = H_column / diameter["diameter"] if diameter["diameter"] > 0 else 0
    else:
        H_column = float("inf")
        H_D_ratio = float("inf")

    # 总压降与塔底压力：按修正后的物理塔板数累加（统一单板压降时 ΔP_total = N_actual_trays × ΔP_tray）
    dP_total = None
    P_bottom = None
    if N_actual_trays != float("inf"):
        dP_total = N_actual_trays * hydraulics["pressure_drop_per_tray"]
        if P_top is not None:
            P_bottom = P_top + dP_total

    constraints = _check_distillation_constraints(fenske, underwood, gilliland, efficiency, diameter, hydraulics, H_D_ratio)

    return {
        "valid": True if N_theoretical != float("inf") else False,
        "N_min": fenske["N_min"], "R_min": underwood["R_min"], "R_actual": gilliland["R_actual"],
        "N_theoretical": N_theoretical, "N_actual_trays": N_actual_trays,
        "feed_stage": gilliland.get("feed_stage"), "N_rectifying": gilliland.get("N_rectifying"), "N_stripping": gilliland.get("N_stripping"),
        "feed_stage_actual": feed_stage_actual,
        "N_rectifying_actual": N_rectifying_actual, "N_stripping_actual": N_stripping_actual,
        "tray_efficiency": efficiency["efficiency"], "tray_efficiency_percent": efficiency["efficiency_percent"],
        "D_column": diameter["diameter"], "D_column_calculated": diameter["diameter_calculated"],
        "H_column": H_column, "H_D_ratio": H_D_ratio,
        "q": underwood["q"], "q_state": underwood.get("q_state", ""),
        "flooding_percent": hydraulics["flooding_percent"],
        "pressure_drop_per_tray": hydraulics["pressure_drop_per_tray"],
        "pressure_drop_total": dP_total,
        "P_top": P_top, "P_bottom": P_bottom,
        "F_factor": hydraulics["F_factor"],
        "non_key_distribution": fenske.get("non_key_distribution"),
        "constraint_checks": constraints, "warnings": all_warnings,
        "fenske": fenske, "underwood": underwood, "gilliland": gilliland,
        "efficiency_result": efficiency, "diameter_result": diameter, "hydraulics_result": hydraulics,
    }

def _check_distillation_constraints(fenske_r, underwood_r, gilliland_r, efficiency_r, diameter_r, hydraulics_r, H_D_ratio) -> dict:
    checks = {}
    checks["N_min_positive"] = {"value": fenske_r["N_min"], "pass": fenske_r["N_min"] > 0, "criterion": "N_min > 0"}
    R_ratio = gilliland_r.get("R_ratio", 0)
    checks["R_ratio"] = {"value": R_ratio, "pass": 1.1 <= R_ratio <= 3.0, "criterion": "1.1 ≤ R/R_min ≤ 3.0"}
    checks["efficiency"] = {"value": efficiency_r["efficiency_percent"], "pass": 30 <= efficiency_r["efficiency_percent"] <= 90, "criterion": "E_o = 30-90%"}
    checks["flooding"] = {"value": hydraulics_r["flooding_percent"], "pass": 50 <= hydraulics_r["flooding_percent"] <= 85, "criterion": "液泛分率 50-85%"}
    checks["pressure_drop"] = {"value": hydraulics_r["pressure_drop_per_tray"], "pass": hydraulics_r["pressure_drop_per_tray"] <= 1500, "criterion": "单板压降 ≤ 1500 Pa"}
    checks["weeping"] = {"value": hydraulics_r["weeping_check"]["F_hole"], "pass": not hydraulics_r["weeping_check"]["weeping"], "criterion": "F_hole ≥ 5 (无漏液)"}
    checks["H_D_ratio"] = {"value": H_D_ratio, "pass": 2 <= H_D_ratio <= 30, "criterion": "H/D = 2-30"}
    checks["column_diameter"] = {"value": diameter_r["diameter"], "pass": diameter_r["diameter"] >= 0.3, "criterion": "D ≥ 0.3 m"}

    n_pass = sum(1 for c in checks.values() if c["pass"])
    return {"checks": checks, "n_pass": n_pass, "n_total": len(checks), "pass_rate": n_pass / len(checks) if len(checks) > 0 else 0}


# ============================================================
# 动态物性代理与辅助函数 (彻底剔除静态本地库，代理至标准统一物性引擎)
# ============================================================

def get_component_properties(name: str) -> dict:
    """由物性引擎 thermo_get_props 获取组分参考物性"""
    try:
        props = thermo_get_props(name)
        if props:
            return {
                "Tb": props.get("Tb", 350.0),
                "MW": props.get("MW", 50.0),
                "Pc": props.get("Pc", 40.0),
                "omega": props.get("omega", 0.2),
            }
    except Exception:
        pass
    return {"Tb": 350.0, "MW": 50.0, "Pc": 40.0, "omega": 0.2, "note": "unknown component"}


def antoine_equation(T: float, A: float, B: float, C: float, Tmin: float = 0, Tmax: float = 1000, unit_P: str = "kPa") -> float:
    T_C = T - 273.15
    try:
        if T_C <= -C: T_C = -C + 0.01
        P_mmHg = 10 ** (A - B / (T_C + C))
        P = P_mmHg * 0.133322 if unit_P == "kPa" else P_mmHg
    except:
        P = 101.325
    return float(P)


def get_saturation_pressure(component: str, T: float) -> float:
    """由物性引擎 get_fluid_vapor_pressure 计算饱和蒸气压 (kPa)"""
    try:
        p_pa = get_fluid_vapor_pressure(component, T)
        if p_pa is not None and p_pa > 0:
            return p_pa / 1000.0
    except Exception:
        pass
    Tb = get_component_properties(component).get("Tb", 350.0)
    if T > 0 and Tb > 0 and 0.5 < T / Tb < 1.0:
        return float(101.325 * math.exp(40000/8.314 * (1/Tb - 1/T)))
    return 101.325


def wilson_k_values(T: float, P: float, components: list) -> dict:
    return {c: float(get_saturation_pressure(c, T) * 1000 / P) if P > 0 else 1.0 for c in components}


def bubble_point_temperature(P: float, components: list, z: list, T_guess: float = 350.0) -> dict:
    """泡点温度 — 优先走 thermo_helper UNIFAC 活度系数模型（非理想体系），
    失败时回退到 Raoult 定律迭代（理想溶液，γ=1）。"""
    # 优先：UNIFAC 活度系数路径
    try:
        r = thermo_bubble_T(components, z, P)
        if r.get("converged") and r.get("T"):
            return {"Tb": float(r["T"]), "K_values": dict(zip(components, r["K"])),
                    "iterations": r.get("iterations", 0), "method": r.get("method", "unifac")}
    except Exception:
        pass
    # 回退：Raoult 迭代（原实现）
    T = T_guess
    for _ in range(50):
        K = [wilson_k_values(T, P, components)[c] for c in components]
        f = sum(z[i] * (K[i] - 1) for i in range(len(z)))
        if abs(f) < 1e-6: break
        K2 = [wilson_k_values(T + 0.1, P, components)[c] for c in components]
        df_dT = (sum(z[i] * (K2[i] - 1) for i in range(len(z))) - f) / 0.1
        if abs(df_dT) < 1e-10: break
        T -= f / df_dT
    return {"Tb": float(T), "K_values": wilson_k_values(T, P, components),
            "iterations": _ + 1, "method": "raoult"}


def dew_point_temperature(P: float, components: list, y: list, T_guess: float = 350.0) -> dict:
    """露点温度 — 优先走 thermo_helper UNIFAC 活度系数模型（非理想体系），
    失败时回退到 Raoult 定律迭代（理想溶液，γ=1）。"""
    # 优先：UNIFAC 活度系数路径
    try:
        r = thermo_dew_T(components, y, P)
        if r.get("converged") and r.get("T"):
            return {"Td": float(r["T"]), "K_values": dict(zip(components, r["K"])),
                    "iterations": r.get("iterations", 0), "method": r.get("method", "unifac")}
    except Exception:
        pass
    # 回退：Raoult 迭代（原实现）
    T = T_guess
    for _ in range(50):
        K = [wilson_k_values(T, P, components)[c] for c in components]
        f = sum(y[i] / K[i] for i in range(len(y))) - 1.0
        if abs(f) < 1e-6: break
        K2 = [wilson_k_values(T + 0.1, P, components)[c] for c in components]
        df_dT = (sum(y[i] / K2[i] for i in range(len(y))) - 1.0 - f) / 0.1
        if abs(df_dT) < 1e-10: break
        T -= f / df_dT
    return {"Td": float(T), "K_values": wilson_k_values(T, P, components),
            "iterations": _ + 1, "method": "raoult"}


# ============================================================
# 填料塔压降计算（修正为湿床层 Leva 公式）
# ============================================================

def calculate_packed_column_pressure_drop(
    Z: float, G: float, rho_g: float, rho_l: float, mu_g: float,
    epsilon: float = 0.7, packing_type: str = "raschig_ring", dp: float = 0.025,
) -> dict:
    """填料塔压降计算 — 使用经典的 Leva 湿床层关联式 (Irrigated Packed Bed)，
    由于原始接口未直接传入液相流量，此处的液相负荷按精馏常态 L/G = 0.8 估算。
    这克服了干床 Ergun 方程在气液两相逆流中严重低估压降的根本缺陷。
    """
    warnings = []
    # G 在原设计中为气相质量通量 (kg/(m²·s))
    G_flux = G
    # 典型液气比 L/G ≈ 0.80
    L_flux = 0.8 * G_flux
    
    # 英制换算
    G_prime = G_flux * 0.204816
    L_prime = L_flux * 0.204816
    rho_g_lb = rho_g * 0.062428
    
    if packing_type == "raschig_ring":
        alpha, beta = 0.52, 0.16
    elif packing_type == "pall_ring":
        alpha, beta = 0.21, 0.06
    else:
        alpha, beta = 0.15, 0.04
        
    dp_inch_ft = alpha * (10 ** (beta * L_prime)) * (G_prime ** 2) / rho_g_lb if rho_g_lb > 0 else 0
    dp_Pa_m = dp_inch_ft * 816.6
    delta_P = dp_Pa_m * Z
    
    dp_flood_m = 1800.0
    flooding_fraction = min(dp_Pa_m / dp_flood_m, 1.0)
    
    if dp_Pa_m > 1200.0:
        warnings.append(f"填料塔压降梯度 {dp_Pa_m:.1f} Pa/m 偏高，可能已接近载液或泛点区间")
        
    return {
        "delta_P": float(delta_P),
        "pressure_gradient": float(dp_Pa_m),
        "flooding_fraction": float(flooding_fraction),
        "velocity": G / rho_g if rho_g > 0 else 0,
        "Re": rho_g * (G / rho_g) * dp / (mu_g * (1 - epsilon)) if (mu_g * (1 - epsilon)) > 0 else 0,
        "packing_type": packing_type,
        "Z": Z,
        "warnings": warnings,
    }


# ============================================================
# 工具层兼容别名
# ============================================================

def column_diameter(
    V: float, rho_vapor: float, rho_liquid: float,
    sigma: float = 0.020, tray_spacing: float = 0.60,
    hole_area_fraction: float = 0.10, flooding_fraction: float = 0.80,
) -> dict:
    return calculate_column_diameter(
        V=V, rho_vapor=rho_vapor, rho_liquid=rho_liquid,
        sigma=sigma, tray_spacing=tray_spacing,
        hole_area_fraction=hole_area_fraction, flooding_fraction=flooding_fraction,
    )


def calculate_tray_pressure_drop(
    V_mass: float, L_mass: float, rho_vapor: float, rho_liquid: float,
    sigma: float = 0.020, mu_liquid: float = 0.001, D_column: float = 1.0,
    tray_spacing: float = 0.60, hole_diameter: float = 0.005,
    hole_area_fraction: float = 0.10, weir_height: float = 0.050,
) -> dict:
    hyd = calculate_tray_hydraulics(
        V_mass=V_mass, L_mass=L_mass, rho_vapor=rho_vapor, rho_liquid=rho_liquid,
        sigma=sigma, mu_liquid=mu_liquid, D_column=D_column,
        tray_spacing=tray_spacing, hole_diameter=hole_diameter,
        hole_area_fraction=hole_area_fraction, weir_height=weir_height,
    )
    return {
        "delta_P": hyd["pressure_drop_per_tray"],
        "delta_P_dry": hyd["pressure_drop_dry"],
        "delta_P_wet": hyd["pressure_drop_wet"],
        "flooding_percent": hyd["flooding_percent"],
        "vapor_velocity": hyd["vapor_velocity"],
        "valid": hyd["valid"],
        "warnings": hyd.get("warnings", []),
    }


def tray_efficiency(
    alpha_LK_HK: float, mu_feed: float, temperature: float = 298.15,
) -> dict:
    return estimate_tray_efficiency(
        alpha_LK_HK=alpha_LK_HK, mu_feed=mu_feed, temperature=temperature,
    )


def tray_efficiency_oconnell(
    alpha_LK_HK: float, mu_feed: float, temperature: float = 298.15,
) -> dict:
    return estimate_tray_efficiency(
        alpha_LK_HK=alpha_LK_HK, mu_feed=mu_feed, temperature=temperature,
    )


def predict_saturation_pressure(component: str, T: float) -> dict:
    Psat = get_saturation_pressure(component, T)
    return {
        "component": component,
        "T": T,
        "T_C": T - 273.15,
        "Psat_kPa": Psat,
        "Psat_Pa": Psat * 1000.0,
        "Psat_bar": Psat / 100.0,
        "valid": Psat > 0,
    }


def relative_volatility_from_pressure(
    component_A: str, component_B: str, T: float,
) -> dict:
    Psat_A = get_saturation_pressure(component_A, T)
    Psat_B = get_saturation_pressure(component_B, T)
    if Psat_B <= 0:
        return {
            "alpha": 1.0, "component_A": component_A, "component_B": component_B,
            "T": T, "Psat_A_kPa": Psat_A, "Psat_B_kPa": Psat_B,
            "valid": False, "warning": "Psat_B ≤ 0, 默认 α=1",
        }
    alpha = Psat_A / Psat_B
    return {
        "alpha": alpha, "component_A": component_A, "component_B": component_B,
        "T": T, "T_C": T - 273.15,
        "Psat_A_kPa": Psat_A, "Psat_B_kPa": Psat_B,
        "valid": alpha > 0,
    }


def predict_mixture_viscosity(
    components: list, z: list, T: float, P: float = 101325.0,
) -> dict:
    n = len(components)
    mu_vals_Pas = []
    for name in components:
        mu_Pas = get_fluid_viscosity(name, T, P, phase="liquid")
        mu_vals_Pas.append(max(mu_Pas, 1e-10))
    log_mu = sum(z[i] * math.log(mu_vals_Pas[i]) for i in range(n))
    mu_mix_Pas = math.exp(log_mu)
    mu_mix_cP = mu_mix_Pas * 1000.0
    return {
        "mu_l": mu_mix_cP,
        "mu_l_Pas": mu_mix_Pas,
        "T": T,
        "components": components,
        "valid": mu_mix_Pas > 0,
    }


def predict_mixture_molecular_weight(
    components: list, z: list,
) -> dict:
    MW_total = 0.0
    for i, name in enumerate(components):
        MW_i = get_fluid_MW(name)
        MW_total += z[i] * MW_i
    return {
        "MW_mix": MW_total,
        "MW": MW_total,
        "components": components,
        "valid": MW_total > 0,
    }


def feed_stage_estimate(
    N_actual: int, x_LK_dist: float, x_HK_dist: float,
    x_LK_feed: float, x_HK_feed: float,
    x_LK_btm: float = 0.0,
    D: float = 0.0,
    B: float = 0.0,
    feed_stage_method: str = "kirkbride",
) -> dict:
    """进料板位置估算（Kirkbride 法）

    Kirkbride 方程: Nr/Ns = [(B/D)*(zHK/zLK)*(xB_LK/xD_HK)^2]^0.206
    """
    if all(v > 0 for v in [x_LK_dist, x_HK_dist, x_LK_feed, x_HK_feed]):
        if feed_stage_method == "kirkbride":
            # 标准 Kirkbride 方程
            if D > 0 and B > 0 and x_LK_btm > 0:
                kirkbride_base = (B / D) * (x_HK_feed / x_LK_feed) * (x_LK_btm / x_HK_dist) ** 2
            else:
                kirkbride_base = (x_HK_feed / x_LK_feed) * (x_LK_btm / x_HK_dist if x_LK_btm > 0 else 1.0) ** 2
            Nr_Ns_ratio = kirkbride_base ** 0.206
            N_strip = max(2, int(N_actual / (1 + Nr_Ns_ratio)))
            N_rect = N_actual - N_strip
            feed_stage = N_rect + 1
        else:
            N_strip = max(2, int(N_actual * 0.5))
            N_rect = N_actual - N_strip
            feed_stage = N_rect + 1
    else:
        N_strip = max(2, int(N_actual * 0.5))
        N_rect = N_actual - N_strip
        feed_stage = N_rect + 1
    return {
        "feed_stage": feed_stage,
        "N_rectifying": N_rect,
        "N_stripping": N_strip,
        "N_actual": N_actual,
        "method": feed_stage_method,
        "valid": feed_stage > 1,
    }


# ============================================================
# 精馏塔热负荷计算（冷凝器、再沸器）
# ============================================================

def calculate_condenser_duty(
    D_mol_s: float,
    R_actual: float,
    components: list,
    x_D: list,
    T_top: float,
    P: float = 101325.0,
    Hvap_values: list | None = None,
) -> dict:
    warnings: list[str] = []
    V_top = D_mol_s * (R_actual + 1.0)
    L_top = D_mol_s * R_actual
    
    if Hvap_values is None or len(Hvap_values) != len(components):
        try:
            from .thermo_helper import calc_enthalpy_vaporization
            Hvap_values = []
            for i, comp in enumerate(components):
                try:
                    hvap = calc_enthalpy_vaporization(comp, T_top)
                    if hvap is not None and hvap > 0:
                        Hvap_values.append(hvap)
                    else:
                        props = get_component_properties(comp)
                        Tb = props.get("Tb", T_top)
                        Hvap_values.append(88.0 * Tb)  
                        warnings.append(f"组分 {comp} 汽化焓使用 Trouton 规则估算")
                except Exception as e:
                    props = get_component_properties(comp)
                    Tb = props.get("Tb", T_top)
                    Hvap_values.append(88.0 * Tb)
                    warnings.append(f"组分 {comp} 汽化焓计算失败，使用 Trouton 规则: {str(e)}")
        except ImportError:
            Hvap_values = []
            for comp in components:
                props = get_component_properties(comp)
                Tb = props.get("Tb", T_top)
                Hvap_values.append(88.0 * Tb)
            warnings.append("使用 Trouton 规则估算汽化焓")
    
    Hvap_mix = sum(x_D[i] * Hvap_values[i] for i in range(len(components)))
    Q_condenser = V_top * Hvap_mix
    Q_condenser_kW = Q_condenser / 1000.0
    
    return {
        "Q_condenser_W": Q_condenser,
        "Q_condenser_kW": Q_condenser_kW,
        "V_top_mol_s": V_top,
        "L_top_mol_s": L_top,
        "Hvap_mix_J_mol": Hvap_mix,
        "Hvap_components_J_mol": Hvap_values,
        "components": components,
        "x_D": x_D,
        "T_top_K": T_top,
        "valid": Q_condenser > 0,
        "warnings": warnings,
    }


def calculate_reboiler_duty(
    B_mol_s: float,
    R_actual: float,
    components: list,
    x_B: list,
    T_bottom: float,
    q: float = 1.0,
    z_feed: list | None = None,
    P: float = 101325.0,
    Hvap_values: list | None = None,
    D_mol_s: float | None = None,  # 保持兼容性，放在末尾作为可选参数
) -> dict:
    """计算塔底再沸器热负荷"""
    warnings: list[str] = []
    
    if Hvap_values is None or len(Hvap_values) != len(components):
        try:
            from .thermo_helper import calc_enthalpy_vaporization
            Hvap_values = []
            for i, comp in enumerate(components):
                try:
                    hvap = calc_enthalpy_vaporization(comp, T_bottom)
                    if hvap is not None and hvap > 0:
                        Hvap_values.append(hvap)
                    else:
                        props = get_component_properties(comp)
                        Tb = props.get("Tb", T_bottom)
                        Hvap_values.append(88.0 * Tb)
                        warnings.append(f"组分 {comp} 汽化焓使用 Trouton 规则估算")
                except Exception as e:
                    props = get_component_properties(comp)
                    Tb = props.get("Tb", T_bottom)
                    Hvap_values.append(88.0 * Tb)
                    warnings.append(f"组分 {comp} 汽化焓计算失败，使用 Trouton 规则: {str(e)}")
        except ImportError:
            Hvap_values = []
            for comp in components:
                props = get_component_properties(comp)
                Tb = props.get("Tb", T_bottom)
                Hvap_values.append(88.0 * Tb)
            warnings.append("使用 Trouton 规则估算汽化焓")
    
    Hvap_mix = sum(x_B[i] * Hvap_values[i] for i in range(len(components)))
    
    # 严格采用恒摩尔流 (CMO) 的物料衡算 V' = D*(R+1) - (1-q)*F 计算提馏段气相流量
    if D_mol_s is not None:
        V_strip = D_mol_s * (R_actual + 1.0) - (1.0 - q) * (D_mol_s + B_mol_s)
    else:
        V_strip = B_mol_s * (R_actual + 1.0) - (1.0 - q) * (2.0 * B_mol_s)
        warnings.append("由于未提供 D_mol_s，提馏段流量已降级为 D_mol_s ≈ B_mol_s 的等摩尔假定进行估算")
        
    Q_reboiler = V_strip * Hvap_mix
    Q_reboiler_kW = Q_reboiler / 1000.0
    L_strip = V_strip + B_mol_s
    
    return {
        "Q_reboiler_W": Q_reboiler,
        "Q_reboiler_kW": Q_reboiler_kW,
        "V_strip_mol_s": V_strip,
        "L_strip_mol_s": L_strip,
        "Hvap_mix_J_mol": Hvap_mix,
        "Hvap_components_J_mol": Hvap_values,
        "components": components,
        "x_B": x_B,
        "T_bottom_K": T_bottom,
        "q": q,
        "valid": Q_reboiler > 0,
        "warnings": warnings,
    }


def calculate_column_heat_duty(
    D_mol_s: float,
    B_mol_s: float,
    R_actual: float,
    components: list,
    x_D: list,
    x_B: list,
    T_top: float,
    T_bottom: float,
    q: float = 1.0,
    P: float = 101325.0,
) -> dict:
    """计算精馏塔总热负荷（冷凝器 + 再沸器）"""
    condenser_result = calculate_condenser_duty(
        D_mol_s=D_mol_s,
        R_actual=R_actual,
        components=components,
        x_D=x_D,
        T_top=T_top,
        P=P,
    )
    
    reboiler_result = calculate_reboiler_duty(
        B_mol_s=B_mol_s,
        R_actual=R_actual,
        components=components,
        x_B=x_B,
        T_bottom=T_bottom,
        q=q,
        P=P,
        D_mol_s=D_mol_s,  # 正确传递
    )
    
    Q_total = condenser_result["Q_condenser_W"] + reboiler_result["Q_reboiler_W"]
    warnings = condenser_result.get("warnings", []) + reboiler_result.get("warnings", [])
    
    return {
        "Q_condenser_W": condenser_result["Q_condenser_W"],
        "Q_condenser_kW": condenser_result["Q_condenser_kW"],
        "Q_reboiler_W": reboiler_result["Q_reboiler_W"],
        "Q_reboiler_kW": reboiler_result["Q_reboiler_kW"],
        "Q_total_W": Q_total,
        "Q_total_kW": Q_total / 1000.0,
        "condenser_result": condenser_result,
        "reboiler_result": reboiler_result,
        "valid": condenser_result["valid"] and reboiler_result["valid"],
        "warnings": warnings,
    }