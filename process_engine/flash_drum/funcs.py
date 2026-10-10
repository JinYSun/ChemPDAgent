"""闪蒸罐工具函数。

设计哲学：
- 闪蒸罐是单自由度气液分离设备：给定 P + 进料组成 → 只剩 1 个自由度（T）
- 指定任一关键组分在某相的回收率 → 唯一确定 T → 全部相平衡分配确定
- 工程上推荐用挥发度最低的待回收组分作关键组分（即"分离瓶颈"），反算 T 最稳定

热力学模型：thermo_helper 双引擎（thermo + CoolProp）。全程 SI 单位。

修复记录（化工角度）：
  [FIX-1] 搜索范围约束在两相区内（泡点~露点），避免在单相区做无意义 flash；
           同时补充温度收敛判据（|T_high - T_low| < 0.01 K），防止回收率
           平台区无限迭代。
  [FIX-2] 二分法单调性验证：在进入迭代前用两端点实测回收率验证单调方向，
           当出现非单调（mul体系、近沸点体系）时给出明确错误而非返回错误结果。
  [FIX-3] beta 边界时 y/x 取值保护：全液（beta≈0）时气相回收率强制为 0，
           全气（beta≈1）时液相回收率强制为 0，不再使用物理意义不明的 y/x。
  [FIX-4] 可达性判断的饱和情况改为明确 ValueError，不静默改变语义；
           同时将"饱和"温度限制在物理合理的两相区边界而非固定 100K/800K。
  [FIX-5] 输入校验：组分名大小写规范化匹配 + 各流量非负检查。
  [FIX-6] _build_result 统一用 beta 边界阈值（EPS_VF）保护 y/x 取值，
           确保全液/全气时输出结构一致且物理正确。
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
    calc_bubble_point_T,
    calc_dew_point_T,
)

# beta 边界判定阈值：beta < EPS_VF 视为全液，beta > 1 - EPS_VF 视为全气
_EPS_VF = 1e-9

# 二分法温度收敛阈值 (K)
_T_TOL = 0.01

# 回收率收敛阈值（绝对误差）
_REC_TOL = 1e-4

# 二分法最大迭代次数
_MAX_ITER = 120


# ─────────────────────────────────────────────
# 闪蒸器封装（私有）
# ─────────────────────────────────────────────
def _safe_flash_TP(components, T, P, zs):
    """带友好错误的等温闪蒸，返回 {beta, x, y, K, phase, method}。"""
    try:
        result = calc_flash_rachford_rice(components, zs, T, P)
        if not result or result.get("beta") is None:
            raise ValueError("Rachford-Rice 未收敛")
        return result
    except Exception as exc:
        raise ValueError(
            f"闪蒸失败 T={T:.2f} K, P={P:.0f} Pa: {type(exc).__name__}: {exc}。"
            f"可能是体系非理想性强或操作条件超出物性范围。"
        ) from exc


def _component_recovery(result, key_idx, F_total, key_inlet_flow, recovery_phase):
    """从 flash 结果计算关键组分在指定相的回收率。

    [FIX-3] 当 beta 处于边界（全液/全气）时，直接返回 0 或 1，
    不使用物理意义不明确的 y/x 值（thermo 在单相时 y≈z 或 x≈z）。

    result 含:
      beta — 汽化率（气相摩尔分率）
      y    — 气相摩尔分数列表
      x    — 液相摩尔分数列表
    """
    beta = result.get("beta", 0.0)

    if recovery_phase == "vapor":
        # [FIX-3] 全液时气相回收率严格为 0
        if beta <= _EPS_VF:
            return 0.0
        # [FIX-3] 全气时气相回收率严格为 1
        if beta >= 1.0 - _EPS_VF:
            return 1.0
        y = result.get("y", [])
        y_key = y[key_idx] if key_idx < len(y) else 0.0
        vapor_total = F_total * beta
        return (y_key * vapor_total) / key_inlet_flow
    else:
        # [FIX-3] 全气时液相回收率严格为 0
        if beta >= 1.0 - _EPS_VF:
            return 0.0
        # [FIX-3] 全液时液相回收率严格为 1
        if beta <= _EPS_VF:
            return 1.0
        x = result.get("x", [])
        x_key = x[key_idx] if key_idx < len(x) else 0.0
        liquid_total = F_total * (1.0 - beta)
        return (x_key * liquid_total) / key_inlet_flow


def _build_result(T_final, final, F_total, components,
                  flash_pressure_Pa, key_component, target, recovery_phase):
    """根据 flash 结果构建最终输出 dict（气相/液相各组分流量）。

    [FIX-6] 用 _EPS_VF 统一保护 y/x 取值，全液时 vapor_flows 全为 0 且
    liquid_flows 用 x（全液 flash 返回的 x 等于 z，物理正确）；
    全气时反之。
    """
    VF = final.get("beta", 0.0)
    ys = final.get("y", [])
    xs = final.get("x", [])
    vapor_total = F_total * VF
    liquid_total = F_total * (1.0 - VF)

    # [FIX-6] 全液时 y 无意义，用零流量；全气时 x 无意义，同理
    if VF <= _EPS_VF:
        vapor_flows = {comp: 0.0 for comp in components}
        liquid_flows = {comp: xs[i] * liquid_total if i < len(xs) else 0.0
                        for i, comp in enumerate(components)}
    elif VF >= 1.0 - _EPS_VF:
        vapor_flows = {comp: ys[i] * vapor_total if i < len(ys) else 0.0
                       for i, comp in enumerate(components)}
        liquid_flows = {comp: 0.0 for comp in components}
    else:
        vapor_flows = {comp: ys[i] * vapor_total if i < len(ys) else 0.0
                       for i, comp in enumerate(components)}
        liquid_flows = {comp: xs[i] * liquid_total if i < len(xs) else 0.0
                        for i, comp in enumerate(components)}

    _vprint(f"\n【calc_flash_drum】")
    _vprint(f"  P = {flash_pressure_Pa:.0f} Pa ({flash_pressure_Pa/1e6:.3f} MPa)")
    _vprint(f"  关键组分 {key_component} {recovery_phase} 回收率目标: {target:.4f}")
    _vprint(f"  反算闪蒸温度 T = {T_final:.2f} K ({T_final-273.15:.2f}°C)")
    _vprint(f"  汽化率 VF = {VF:.4f}")
    _vprint(f"  气相总流量: {vapor_total:.2f} mol/s")
    _vprint(f"  液相总流量: {liquid_total:.2f} mol/s")

    return {
        "flash_temperature_K": T_final,
        "vapor_molar_flows_mol_per_s": vapor_flows,
        "liquid_molar_flows_mol_per_s": liquid_flows,
    }


# ─────────────────────────────────────────────
# [FIX-1] 私有：获取体系两相区温度边界（泡点、露点）
# ─────────────────────────────────────────────
def _get_two_phase_bounds(components, zs, P):
    """计算体系在压力 P 下的泡点和露点温度，返回 (T_bubble, T_dew)。

    [FIX-1] 二分法搜索范围应限定在两相区 [T_bubble, T_dew] 内。在此区间外
    flash 的 beta=0（T < T_bubble）或 beta=1（T > T_dew）是边界值，
    不能用于回收率的单调性判断，会让二分法方向判断失效。

    失败时返回 (None, None)，由调用方 fallback 到宽搜索范围并给出警告。
    """
    try:
        bub = calc_bubble_point_T(components, zs, P)
        dew = calc_dew_point_T(components, zs, P)
        T_bubble = bub.get("T") if (bub and bub.get("converged")) else None
        T_dew = dew.get("T") if (dew and dew.get("converged")) else None
        if T_bubble and T_dew and T_dew > T_bubble:
            return T_bubble, T_dew
    except Exception:
        pass
    return None, None


# ═════════════════════════════════════════════
# 工具：闪蒸罐设计计算
# ═════════════════════════════════════════════
def calc_flash_drum(
    flash_pressure_Pa: float,
    inlet_molar_flows_mol_per_s: dict,
    key_component: str,
    key_component_recovery: float,
    recovery_phase: str = "vapor",
) -> dict:
    """计算闪蒸罐反算操作温度并给出气液出口流量分配。

    场景：等压冷却闪蒸（永久气体 + 有机物分离）。给定上游工艺压力 P 和进料组成，
    指定某个关键组分在指定相的目标回收率，二分法反算所需闪蒸温度 T，
    再做一次完整 flash 得到气液各组分摩尔流量。

    单自由度系统：给定 P 和进料后，一个回收率约束唯一确定 T，所有组分的相平衡
    分配随之确定。

    工程提示：关键组分应选挥发度最低的待回收组分（分离瓶颈），使回收率对 T
    的敏感度最高，二分收敛最稳定。

    参数:
        flash_pressure_Pa             — 闪蒸压力 (Pa)，由上游工艺约束决定
        inlet_molar_flows_mol_per_s   — {组分名: 进料摩尔流量} (mol/s)
        key_component                 — 关键组分名（必须在 inlet_molar_flows_mol_per_s 的键中）
        key_component_recovery        — 关键组分目标回收率 (0~1，开区间)
        recovery_phase                — "vapor" 或 "liquid"，关键组分从哪个相回收

    返回:
        {"flash_temperature_K": 闪蒸温度 T (K),
         "vapor_molar_flows_mol_per_s":  {组分: 气相摩尔流量 mol/s},
         "liquid_molar_flows_mol_per_s": {组分: 液相摩尔流量 mol/s}}
    """
    # ── 输入校验 ─────────────────────────────────────────────────────────────
    if recovery_phase not in ("vapor", "liquid"):
        raise ValueError(
            f"recovery_phase 必须是 'vapor' 或 'liquid'，得到 '{recovery_phase}'"
        )
    if not (0.0 < key_component_recovery < 1.0):
        raise ValueError(
            f"key_component_recovery 必须在开区间 (0, 1) 内，得到 {key_component_recovery}"
        )
    if len(inlet_molar_flows_mol_per_s) < 2:
        raise ValueError(
            f"闪蒸分离需要至少 2 个组分（不同挥发度），"
            f"进料只有 {len(inlet_molar_flows_mol_per_s)} 个组分。"
        )

    # [FIX-5] 流量非负检查
    for comp, flow in inlet_molar_flows_mol_per_s.items():
        if flow < 0:
            raise ValueError(
                f"组分 '{comp}' 的摩尔流量为负值 ({flow} mol/s)，请检查输入。"
            )

    components = list(inlet_molar_flows_mol_per_s.keys())
    flows = list(inlet_molar_flows_mol_per_s.values())
    F_total = sum(flows)
    if F_total <= 0:
        raise ValueError("进料总流量必须大于 0。")

    # [FIX-5] 组分名大小写规范化匹配
    # key_component 可能与字典键大小写不同（thermo 不区分大小写但 dict 区分）
    key_component_matched = None
    for comp in components:
        if comp.lower().strip() == key_component.lower().strip():
            key_component_matched = comp
            break
    if key_component_matched is None:
        raise ValueError(
            f"关键组分 '{key_component}' 不在进料组分中（大小写不敏感匹配）: "
            f"{components}"
        )
    key_component = key_component_matched  # 使用字典中的实际键

    zs = [f / F_total for f in flows]
    key_idx = components.index(key_component)
    key_inlet_flow = flows[key_idx]
    if key_inlet_flow <= 0:
        raise ValueError(
            f"关键组分 '{key_component}' 的进料流量为 {key_inlet_flow} mol/s，"
            f"无法计算回收率。"
        )

    # ── [FIX-1] 确定搜索温度范围（限定在两相区） ────────────────────────────
    T_bubble, T_dew = _get_two_phase_bounds(components, zs, flash_pressure_Pa)

    if T_bubble is not None and T_dew is not None:
        # 在两相区边界内留 0.1 K 余量，避免边界 flash 不收敛
        T_low = T_bubble + 0.1
        T_high = T_dew - 0.1
        _vprint(
            f"  两相区范围: T_bubble={T_bubble:.2f}K, T_dew={T_dew:.2f}K → "
            f"搜索 [{T_low:.2f}, {T_high:.2f}] K"
        )
        if T_high <= T_low:
            # 泡露点过于接近（纯组分或极窄沸程）：用单点 flash
            T_mid_single = (T_bubble + T_dew) / 2.0
            res_single = _safe_flash_TP(components, T_mid_single, flash_pressure_Pa, zs)
            rec_single = _component_recovery(
                res_single, key_idx, F_total, key_inlet_flow, recovery_phase
            )
            _vprint(
                f"  注意：泡露点差 < 0.2K，近似纯组分体系，"
                f"T={T_mid_single:.2f}K 时回收率={rec_single:.4f}"
            )
            return _build_result(
                T_mid_single, res_single, F_total, components,
                flash_pressure_Pa, key_component,
                key_component_recovery, recovery_phase
            )
    else:
        # 泡露点计算失败（如含永久气体体系）：fallback 到宽范围并给出警告
        T_low, T_high = 150.0, 700.0
        _vprint(
            f"  警告：泡露点计算失败，使用宽温度范围 [{T_low}, {T_high}] K。"
            f"结果可能在单相区，请核查。"
        )

    target = key_component_recovery

    # ── [FIX-2] 两端点回收率 + 单调性验证 ───────────────────────────────────
    res_low = _safe_flash_TP(components, T_low, flash_pressure_Pa, zs)
    res_high = _safe_flash_TP(components, T_high, flash_pressure_Pa, zs)
    rec_low = _component_recovery(res_low, key_idx, F_total, key_inlet_flow, recovery_phase)
    rec_high = _component_recovery(res_high, key_idx, F_total, key_inlet_flow, recovery_phase)

    _vprint(
        f"  端点回收率: T_low={T_low:.1f}K → rec={rec_low:.4f}, "
        f"T_high={T_high:.1f}K → rec={rec_high:.4f}"
    )

    # [FIX-2] 验证单调性（气相回收率应随 T 单调递增；液相单调递减）
    if recovery_phase == "vapor":
        is_monotone = rec_high >= rec_low
    else:
        is_monotone = rec_high <= rec_low

    if not is_monotone:
        raise ValueError(
            f"关键组分 '{key_component}' 在 [{T_low:.1f}K, {T_high:.1f}K] 范围内"
            f"{recovery_phase} 回收率不单调（"
            f"T_low 时 {rec_low:.4f}，T_high 时 {rec_high:.4f}）。"
            f"可能原因：体系存在强非理想性、液液分层、或关键组分选取不当。"
            f"建议更换关键组分或缩小温度搜索范围。"
        )

    # ── [FIX-4] 可达性检验（改为明确 ValueError，不静默改变语义）────────────
    if recovery_phase == "vapor":
        # 气相回收率单调递增于 T
        if rec_high < target:
            raise ValueError(
                f"目标气相回收率 {target:.4f} 在 P={flash_pressure_Pa:.0f} Pa 下"
                f"不可达：两相区上限 T={T_high:.2f}K 时关键组分气相回收率仅 "
                f"{rec_high:.4f}。建议降低压力、降低目标回收率，"
                f"或改用其他分离手段。"
            )
        if rec_low >= target:
            raise ValueError(
                f"目标气相回收率 {target:.4f} 过低：即使在两相区下限 "
                f"T={T_low:.2f}K（接近泡点）时气相回收率已达 {rec_low:.4f}。"
                f"此时两相区内任何温度均满足目标，温度无唯一解。"
                f"建议提高目标回收率，或改用挥发度更低的组分作关键组分。"
            )
    else:
        # 液相回收率单调递减于 T
        if rec_low < target:
            raise ValueError(
                f"目标液相回收率 {target:.4f} 在 P={flash_pressure_Pa:.0f} Pa 下"
                f"不可达：两相区下限 T={T_low:.2f}K（接近泡点）时关键组分液相"
                f"回收率仅 {rec_low:.4f}。建议提高压力、降低目标回收率，"
                f"或改用其他分离手段。"
            )
        if rec_high >= target:
            raise ValueError(
                f"目标液相回收率 {target:.4f} 过低：即使在两相区上限 "
                f"T={T_high:.2f}K（接近露点）时液相回收率已达 {rec_high:.4f}。"
                f"此时两相区内任何温度均满足目标，温度无唯一解。"
                f"建议提高目标回收率，或改用挥发度更高的组分作关键组分。"
            )

    # ── 二分法反算 T ──────────────────────────────────────────────────────────
    # [FIX-1] 双重收敛判据：回收率误差 < _REC_TOL  OR  温度区间 < _T_TOL
    T_final = None
    for _ in range(_MAX_ITER):
        T_mid = 0.5 * (T_low + T_high)

        # [FIX-1] 温度收敛优先（防止回收率平台区无限迭代）
        if (T_high - T_low) < _T_TOL:
            T_final = T_mid
            break

        res_mid = _safe_flash_TP(components, T_mid, flash_pressure_Pa, zs)
        rec_mid = _component_recovery(res_mid, key_idx, F_total, key_inlet_flow, recovery_phase)

        if abs(rec_mid - target) < _REC_TOL:
            T_final = T_mid
            break

        if recovery_phase == "vapor":
            # 气相回收率单调递增于 T
            if rec_mid < target:
                T_low = T_mid
            else:
                T_high = T_mid
        else:
            # 液相回收率单调递减于 T
            if rec_mid < target:
                T_high = T_mid
            else:
                T_low = T_mid

    if T_final is None:
        T_final = 0.5 * (T_low + T_high)

    # ── 用收敛 T 做最终 flash，输出气液流量 ─────────────────────────────────
    final = _safe_flash_TP(components, T_final, flash_pressure_Pa, zs)
    return _build_result(
        T_final, final, F_total, components,
        flash_pressure_Pa, key_component, target, recovery_phase
    )


# ═════════════════════════════════════════════
# 工具：等温闪蒸（TP flash，给定 T、P 求气液分配）
# ═════════════════════════════════════════════
def calc_flash_drum_TP(
    flash_temperature_K: float,
    flash_pressure_Pa: float,
    inlet_molar_flows_mol_per_s: dict,
) -> dict:
    """给定闪蒸温度 T 和压力 P，直接求气液分配（等温闪蒸，不反算 T）。

    场景：上游温度、压力均已知（如换热器把物流定温后进闪蒸罐），需要
    在该 (T, P) 下把进料分成气液两相。这与 calc_flash_drum 相反——后者已知
    回收率反算 T，本工具已知 T 正算分配。

    单自由度系统：进料组成固定后，锁定 (T, P) 即唯一确定汽化率 VF 和各组分
    的气液分配。内部做一次 Rachford-Rice 等温闪蒸，按 VF 与相组成 y/x
    折算各组分气液摩尔流量。全液（VF≈0）/全气（VF≈1）时用零流量保护，
    与 calc_flash_drum 的边界处理一致。

    参数:
        flash_temperature_K           — 闪蒸温度 T (K)
        flash_pressure_Pa             — 闪蒸压力 P (Pa)
        inlet_molar_flows_mol_per_s   — {组分名: 进料摩尔流量} (mol/s)

    返回:
        {"flash_temperature_K": T (K),
         "vapor_fraction": 汽化率 VF (0~1，摩尔基准),
         "vapor_molar_flows_mol_per_s":  {组分: 气相摩尔流量 mol/s},
         "liquid_molar_flows_mol_per_s": {组分: 液相摩尔流量 mol/s}}
    """
    # ── 输入校验 ─────────────────────────────────────────────────────────────
    if flash_temperature_K <= 0:
        raise ValueError(f"闪蒸温度必须为正 (K)，得到 {flash_temperature_K}")
    if flash_pressure_Pa <= 0:
        raise ValueError(f"闪蒸压力必须为正 (Pa)，得到 {flash_pressure_Pa}")
    if len(inlet_molar_flows_mol_per_s) < 1:
        raise ValueError("进料至少需要 1 个组分。")

    # [FIX-5] 流量非负检查
    for comp, flow in inlet_molar_flows_mol_per_s.items():
        if flow < 0:
            raise ValueError(
                f"组分 '{comp}' 的摩尔流量为负值 ({flow} mol/s)，请检查输入。"
            )

    components = list(inlet_molar_flows_mol_per_s.keys())
    flows = list(inlet_molar_flows_mol_per_s.values())
    F_total = sum(flows)
    if F_total <= 0:
        raise ValueError("进料总流量必须大于 0。")

    zs = [f / F_total for f in flows]

    # ── 单次等温闪蒸（复用私有底层） ─────────────────────────────────────────
    result = _safe_flash_TP(components, flash_temperature_K, flash_pressure_Pa, zs)

    # ── 按 VF 与相组成折算气液流量（复用 _build_result 的边界保护逻辑）──────
    VF = result.get("beta", 0.0)
    ys = result.get("y", [])
    xs = result.get("x", [])
    vapor_total = F_total * VF
    liquid_total = F_total * (1.0 - VF)

    if VF <= _EPS_VF:
        # 全液：y 无物理意义，气相零流量
        vapor_flows = {comp: 0.0 for comp in components}
        liquid_flows = {comp: xs[i] * liquid_total if i < len(xs) else 0.0
                        for i, comp in enumerate(components)}
    elif VF >= 1.0 - _EPS_VF:
        # 全气：x 无物理意义，液相零流量
        vapor_flows = {comp: ys[i] * vapor_total if i < len(ys) else 0.0
                       for i, comp in enumerate(components)}
        liquid_flows = {comp: 0.0 for comp in components}
    else:
        vapor_flows = {comp: ys[i] * vapor_total if i < len(ys) else 0.0
                       for i, comp in enumerate(components)}
        liquid_flows = {comp: xs[i] * liquid_total if i < len(xs) else 0.0
                        for i, comp in enumerate(components)}

    _vprint(f"\n【calc_flash_drum_TP】")
    _vprint(f"  T = {flash_temperature_K:.2f} K ({flash_temperature_K-273.15:.2f}°C)")
    _vprint(f"  P = {flash_pressure_Pa:.0f} Pa ({flash_pressure_Pa/1e6:.4f} MPa)")
    _vprint(f"  汽化率 VF = {VF:.4f}")
    _vprint(f"  气相总流量: {vapor_total:.4f} mol/s")
    _vprint(f"  液相总流量: {liquid_total:.4f} mol/s")

    return {
        "flash_temperature_K": flash_temperature_K,
        "vapor_fraction": VF,
        "vapor_molar_flows_mol_per_s": vapor_flows,
        "liquid_molar_flows_mol_per_s": liquid_flows,
    }


if __name__ == "__main__":
    # 示例 1：回收率反算 T（原有工具）
    print("=== calc_flash_drum（回收率反算 T）===")
    inlet = {"Propene": 10.0, "Benzene": 5.0, "Cumene": 5.0}
    r1 = calc_flash_drum(
        flash_pressure_Pa=99300.0,
        inlet_molar_flows_mol_per_s=inlet,
        key_component="Cumene",
        key_component_recovery=0.9,
        recovery_phase="liquid",
    )
    print(r1)

    # 示例 2：等温闪蒸（新工具，给定 T、P）
    print("\n=== calc_flash_drum_TP（给定 T=54°C, P=99.3kPa）===")
    r2 = calc_flash_drum_TP(
        flash_temperature_K=54 + 273.15,
        flash_pressure_Pa=99300.0,
        inlet_molar_flows_mol_per_s=inlet,
    )
    print(f"VF = {r2['vapor_fraction']:.4f}")
    print(f"气相: {r2['vapor_molar_flows_mol_per_s']}")
    print(f"液相: {r2['liquid_molar_flows_mol_per_s']}")

    # 示例 3：异常测试（非正温度）
    print("\n=== 异常测试：非正温度 ===")
    try:
        calc_flash_drum_TP(-5.0, 99300.0, inlet)
    except ValueError as e:
        print(f"[OK] 正确捕获: {e}")