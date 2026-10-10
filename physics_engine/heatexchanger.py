"""
换热器工具包 — 对数平均温差、热负荷平衡、换热面积、水力学设计
对应研究内容二："换热器工具包（包含对数平均温差计算、热负荷平衡校验）"

包含：
- LMTD 法：对数平均温差、热负荷、换热面积、F 因子修正（多壳程配置）
- ε-NTU 法：效能 - 传热单元数
- 传热系数估算（含表 2-9/2-10 动态污垢热阻、混合物计算）
- 压降计算（管程 + 壳程）
- 换热器类型选择
- 约束校验与优化建议
"""
from __future__ import annotations
import numpy as np
import math
from typing import Optional, List, Dict, Tuple, Union

from .common import _REQUIRED, check_required_params


# ============================================================
# 污垢热阻数据库 (基于 教材表 2-9 和 表 2-10)
# ============================================================

def _lookup_single_fouling(spec: str | tuple, T_heat_C: float = 100.0, v_flow: float = 1.0) -> float:
    """
    针对单一流体查表计算污垢热阻 (m²·K/W)。
    依据：教材表 2-9（水类）、表 2-10（其他流体）。
    """
    if spec in ("default", "none", "clean", ""):
        spec = "tap"

    category = None
    sub = None
    if isinstance(spec, tuple) and len(spec) == 2:
        category, sub = spec
    elif isinstance(spec, str):
        if "/" in spec:
            parts = spec.split("/")
            category, sub = parts[0].strip(), parts[1].strip()
        else:
            sub = spec.strip().lower()

    # ── 1. 水类流体 (表 2-9 动态逻辑) ─────────────────────────
    _WATER_ALIASES = {
        "sea": "sea", "seawater": "sea",
        "tap": "tap", "well": "tap", "boiler_soft": "tap", "boiler_feed": "tap",
        "cooling_water": "tap", "chilled": "tap", "water": "tap",
        "distilled": "distilled", "demineralized": "distilled",
        "hard": "hard", "hard_water": "hard",
        "river": "river", "brine": "brine"
    }
    
    is_water = False
    if category == "water" or sub in _WATER_ALIASES:
        is_water = True
        sub = _WATER_ALIASES.get(sub, "tap")

    if is_water and sub != "brine":  # brine(盐水)在表2-10中
        is_high_temp = T_heat_C >= 115.0
        is_fast = v_flow > 1.0
        
        if sub == "sea":  # 海水
            return 1.7197e-4 if is_high_temp else 0.8598e-4
        elif sub == "distilled":  # 蒸馏水
            return 0.8598e-4
        elif sub == "hard":  # 硬水
            return 8.5980e-4 if is_high_temp else 5.1590e-4
        elif sub == "river":  # 河水
            if is_high_temp:
                return 5.1590e-4 if is_fast else 6.8788e-4
            else:
                return 3.4394e-4 if is_fast else 5.1590e-4
        else:  # 自来水、井水、锅炉软水、循环水 (默认水类档位)
            return 3.4394e-4 if is_high_temp else 1.7197e-4

    # ── 2. 其他流体 (表 2-10 定值映射) ────────────────────────
    _TABLE_2_10 = {
        "organic_vapor": 0.8598e-4, "solvent_vapor": 1.7197e-4, 
        "natural_gas": 1.7197e-4, "coke_oven_gas": 1.7197e-4, 
        "steam": 0.8598e-4, "water_vapor": 0.8598e-4,
        "air": 3.4394e-4, "flue_gas": 3.4394e-4,
        "organic": 1.7197e-4, "brine": 1.7197e-4, "molten_salt": 0.8598e-4,
        "vegetable_oil": 5.1590e-4,
        "crude_oil": 7.7687e-4,  # 原油中值
        "diesel": 4.2992e-4,     # 柴油中值
        "naphtha": 1.7197e-4, "kerosene": 1.7197e-4, "gasoline": 1.7197e-4,
        "heavy_oil": 8.5980e-4, "bitumen": 1.7197e-3, "asphalt": 1.7197e-3,
        "methanol": 1.7197e-4, "ethanol": 1.7197e-4, "benzene": 1.7197e-4,
        "toluene": 1.7197e-4, "acetone": 1.7197e-4,
        "r22": 0.8598e-4, "r134a": 0.8598e-4, "ammonia": 0.8598e-4,
        "propylene": 0.8598e-4, "ethylene": 0.8598e-4,
    }
    return _TABLE_2_10.get(sub, 1.7197e-4)


def calculate_fouling_resistance(
    spec: Union[str, tuple, list, float],
    z: list[float] | None = None,
    T_heat_C: float = 100.0,
    v_flow: float = 1.0,
    mixture_rule: str = "weighted"
) -> float:
    """
    计算流体的污垢热阻 (m²·K/W)，支持单一流体及多物质混合物。
    """
    if isinstance(spec, (int, float)) and spec > 0:
        return float(spec)

    if isinstance(spec, list):
        if z is None:
            z = [1.0 / len(spec)] * len(spec)
        if len(spec) != len(z):
            raise ValueError("流体组分列表与组成(z)列表长度必须一致")
        
        rfs = [_lookup_single_fouling(comp, T_heat_C, v_flow) for comp in spec]
        
        if mixture_rule == "weighted":
            return sum(frac * rf for frac, rf in zip(z, rfs))
        elif mixture_rule == "max":
            return max(rfs)
            
    return _lookup_single_fouling(spec, T_heat_C, v_flow)


# ============================================================
# 典型总传热系数范围 (W/m²·K)
# ============================================================
_TYPICAL_U: dict[tuple, tuple] = {
    ("water", "water"):                     (582,  1163),
    ("water_low_v", "water_low_v"):         (582,   698),
    ("water_high_v", "water_high_v"):       (814,  1163),
    ("water", "light_organic"):             (467,   814),
    ("light_organic", "water"):             (467,   814),
    ("water", "medium_organic"):            (290,   698),
    ("medium_organic", "water"):            (290,   698),
    ("water", "heavy_organic"):             (116,   467),
    ("heavy_organic", "water"):             (116,   467),
    ("brine", "light_organic"):             (233,   582),
    ("light_organic", "brine"):             (233,   582),
    ("organic_solvent", "organic_solvent"): (198,   233),
    ("light_organic", "light_organic"):     (233,   465),
    ("medium_organic", "medium_organic"):   (116,   349),
    ("heavy_organic", "heavy_organic"):      (58,   233),
    ("organic", "organic"):                 (116,   465),
    ("organic", "water"):                   (233,   698),
    ("water", "condensing_steam_pressure"): (2326, 4652),
    ("condensing_steam_pressure", "water"): (2326, 4652),
    ("water", "condensing_steam"):          (1745, 3489),
    ("condensing_steam", "water"):          (1745, 3489),
    ("water_solution_light", "condensing_steam"): (1071, 1163),
    ("water_solution_heavy", "condensing_steam"): (582,  2908),
    ("light_organic", "condensing_steam"):  (582,  1193),
    ("condensing_steam", "light_organic"):  (582,  1193),
    ("medium_organic", "condensing_steam"): (291,   582),
    ("condensing_steam", "medium_organic"): (291,   582),
    ("heavy_organic", "condensing_steam"):  (114,   349),
    ("condensing_steam", "heavy_organic"):  (114,   349),
    ("water", "condensing_mixed"):          (582,  1163),
    ("water", "condensing_heavy_organic"):  (116,   349),
    ("condensing_steam", "organic"):        (500,  1500),
    ("boiling_water", "water"):             (1000, 3000),
    ("boiling_organic", "water"):           (300,  1000),
    ("gas", "gas"):                         (10,    50),
    ("gas", "water"):                       (20,   300),
    ("condensing_organic", "water"):        (116,   814),
}

_STANDARD_TUBE_OD = [0.010, 0.012, 0.015, 0.019, 0.025, 0.032, 0.038, 0.050, 0.057]
_STANDARD_TUBE_LENGTH = [0.5, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0, 9.0]
_STANDARD_SHELL_DIAMETER = [
    0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
    0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60,
    1.80, 2.00, 2.40, 2.80, 3.00
]


# ============================================================
# 默认冷却介质参数库
# ============================================================
_DEFAULT_COOLANTS: dict[str, dict] = {
    "cooling_water": {
        "fluid_name":  "water",
        "T_in_K":      298.15,
        "T_out_K":     313.15,
        "delta_T_K":   10.0,
        "P_Pa":        300000.0,
        "fouling_key": ("water", "cooling_water"),
        "note":        "循环冷却水：供水入口 25°C，出口 小于80°C，优先选择，冷却25-80°C的物流。",
    },
    "air_cooler": {
        "fluid_name":  "air",
        "T_in_K":      298.15,
        "T_out_K":     323.15,
        "delta_T_K":   15.0,
        "P_Pa":        101325.0,
        "fouling_key": ("gas", "air"),
        "note":        "空气冷却器：入口 25°C（夏季），出口 50~60°C。适用于缺水地区或高温物料，除非特别要求，否则不选择。",
    },
    "chilled_water_7": {
        "fluid_name":  "water",
        "T_in_K":      280.15,
        "T_out_K":     290.15,
        "delta_T_K":   10.0,
        "P_Pa":        200000.0,
        "fouling_key": ("water", "chilled"),
        "note":        "7°C 冷冻水：供水 7°C，回水 17°C。适用于 20~50°C 冷却需求，普通冷水机组即可提供。",
    },
    "chilled_water_5": {
        "fluid_name":  "water",
        "T_in_K":      278.15,
        "T_out_K":     288.15,
        "delta_T_K":   10.0,
        "P_Pa":        200000.0,
        "fouling_key": ("water", "chilled"),
        "note":        "冷冻水，适用于冷却（0-25°C），轻度制冷需求。",
    },
    "brine_minus5": {
        "fluid_name":  "water",
        "T_in_K":      268.15,
        "T_out_K":     273.15,
        "delta_T_K":   5.0,
        "P_Pa":        200000.0,
        "fouling_key": ("water", "brine"),
        "note":        "-5°C 盐水：供液 -5°C，回液 0°C。适用于 -5~5°C 冷却需求。",
    },
    "brine_minus15": {
        "fluid_name":  "water",
        "T_in_K":      258.15,
        "T_out_K":     263.15,
        "delta_T_K":   5.0,
        "P_Pa":        200000.0,
        "fouling_key": ("water", "brine"),
        "note":        "-15°C 盐水：供液 -15°C，回液 -10°C。适用于 -15~0°C 冷却需求。",
    },
    "brine_minus30": {
        "fluid_name":  "water",
        "T_in_K":      243.15,
        "T_out_K":     248.15,
        "delta_T_K":   5.0,
        "P_Pa":        200000.0,
        "fouling_key": ("water", "brine"),
        "note":        "-30°C 盐水：供液 -30°C，回液 -25°C。适用于 -30~-10°C 冷却需求。",
    },
    "propylene_refrigerant": {
        "fluid_name":  "water",
        "T_in_K":      233.15,
        "T_out_K":     233.15,
        "delta_T_K":   0.0,
        "P_Pa":        101325.0,
        "fouling_key": ("refrigerant", "propylene"),
        "note":        "丙烯制冷：蒸发温度约 -40°C（恒温），适用于轻烃分离。",
    },
    "ethylene_refrigerant": {
        "fluid_name":  "water",
        "T_in_K":      173.15,
        "T_out_K":     173.15,
        "delta_T_K":   0.0,
        "P_Pa":        101325.0,
        "fouling_key": ("refrigerant", "ethylene"),
        "note":        "乙烯制冷：蒸发温度约 -100°C（恒温），适用于深冷分离。",
    },
}

_COOLANT_SELECTION_LADDER: list[tuple[float, str]] = [
    (323.15,  "cooling_water"),
    (293.15,  "chilled_water_7"),
    (278.15,  "chilled_water_5"),
    (268.15,  "brine_minus5"),
    (258.15,  "brine_minus15"),
    (243.15,  "brine_minus30"),
    (233.15,  "propylene_refrigerant"),
    (0.0,     "ethylene_refrigerant"),
]

def get_default_coolant(T_hot_out_K: float, prefer_air: bool = False) -> dict:
    if prefer_air and T_hot_out_K >= 323.15:
        key = "air_cooler"
    else:
        key = "ethylene_refrigerant"
        for temp_threshold, coolant_key in sorted(
            _COOLANT_SELECTION_LADDER, key=lambda x: x[0], reverse=True
        ):
            if T_hot_out_K >= temp_threshold:
                key = coolant_key
                break

    cfg = _DEFAULT_COOLANTS[key].copy()
    return {
        "coolant_type":  key,
        "fluid_name":    cfg["fluid_name"],
        "T_cold_in_K":   cfg["T_in_K"],
        "T_cold_out_K":  cfg["T_out_K"],
        "delta_T_K":     cfg["delta_T_K"],
        "P_cold_Pa":     cfg["P_Pa"],
        "fouling_key":   cfg["fouling_key"],
        "note":          cfg["note"],
        "auto_selected": True,
    }

def _round_to_standard(value: float, series: list) -> float:
    if not series:
        return value
    for s in sorted(series):
        if s >= value:
            return s
    return series[-1]


# ============================================================
# LMTD & F 因子
# ============================================================

def lmtd(
    T_hot_in: float,
    T_hot_out: float,
    T_cold_in: float,
    T_cold_out: float,
    flow_arrangement: str = "counter",
    configuration: str = "1_shell_2_tube",
    apply_correction: bool = True,
) -> dict:
    warnings: list[str] = []

    if flow_arrangement == "counter":
        dT1 = T_hot_in - T_cold_out
        dT2 = T_hot_out - T_cold_in
    else:
        dT1 = T_hot_in - T_cold_in
        dT2 = T_hot_out - T_cold_out

    if dT1 <= 0 or dT2 <= 0:
        return {
            "LMTD": 0.0, "LMTD_raw": 0.0, "F": None, "F_quality": None,
            "dT1": dT1, "dT2": dT2, "P": None, "R": None,
            "approach_temp_min": min(dT1, dT2),
            "valid": False,
            "warnings": [f"温度交叉！端差 dT1={dT1:.1f} K, dT2={dT2:.1f} K, LMTD 无意义"],
        }

    if abs(dT1 - dT2) < 0.01:
        LMTD_raw = (dT1 + dT2) / 2.0
    else:
        LMTD_raw = (dT1 - dT2) / math.log(dT1 / dT2)

    F = None
    F_quality = None
    P_val = None
    R_val = None

    if apply_correction:
        F_result = _compute_F_factor(T_hot_in, T_hot_out, T_cold_in, T_cold_out, configuration)
        F = F_result["F"]
        F_quality = F_result.get("quality")
        P_val = F_result.get("P")
        R_val = F_result.get("R")

        if F_result["valid"]:
            if F < 0.5: warnings.append(f"F={F:.3f} 较差，建议增加壳程数或调整温度参数")
            elif F < 0.75: warnings.append(f"F={F:.3f} 可接受但偏低，建议增加壳程数")
        else:
            warnings.append(f"F 因子计算无效: {F_result['message']}")
            F = 1.0

    LMTD = LMTD_raw * F if F is not None else LMTD_raw
    approach_min = min(dT1, dT2)
    if approach_min < 5.0:
        warnings.append(f"最小接近温度 {approach_min:.1f} K < 5 K，换热面积可能过大")

    return {
        "LMTD": LMTD, "LMTD_raw": LMTD_raw, "F": F, "F_quality": F_quality,
        "dT1": dT1, "dT2": dT2, "P": P_val, "R": R_val,
        "approach_temp_min": approach_min, "valid": dT1 > 0 and dT2 > 0 and LMTD > 0,
        "warnings": warnings,
    }

def _compute_F_factor(T_hot_in: float, T_hot_out: float, T_cold_in: float, T_cold_out: float, configuration: str = "1_shell_2_tube") -> dict:
    delta_T_max = T_hot_in - T_cold_in
    if delta_T_max <= 0:
        return {"F": 0.0, "P": 0.0, "R": 0.0, "valid": False, "message": "温度交叉"}

    delta_T_cold = T_cold_out - T_cold_in
    P = delta_T_cold / delta_T_max
    R = (T_hot_in - T_hot_out) / delta_T_cold if delta_T_cold > 1e-6 else 1.0

    if P < 0 or P > 1: return {"F": 0.0, "P": P, "R": R, "valid": False, "message": f"P超出范围"}

    if configuration.startswith("1_shell"): F = _F_bowman_1_shell(R, P)
    elif configuration.startswith("2_shell"): F = _F_bowman_n_shell(R, P, n_shells=2)
    elif configuration.startswith("3_shell"): F = _F_bowman_n_shell(R, P, n_shells=3)
    elif configuration.startswith("4_shell"): F = _F_bowman_n_shell(R, P, n_shells=4)
    elif configuration == "cross_flow": F = _F_cross_flow(R, P)
    else: F = _F_bowman_1_shell(R, P)

    if F < 0 or F > 1.01: return {"F": F, "P": P, "R": R, "valid": False, "message": f"F={F:.4f}超出范围"}
    F = min(1.0, max(0.0, F))

    if F >= 0.90: quality = "excellent"
    elif F >= 0.75: quality = "good"
    elif F >= 0.50: quality = "acceptable"
    else: quality = "poor"
    return {"F": F, "P": P, "R": R, "valid": (F > 0), "quality": quality, "configuration": configuration}

def _F_bowman_1_shell(R: float, P: float) -> float:
    try:
        if P <= 0 or P >= 1: return 0.0
        sqrt_r2p1 = math.sqrt(R ** 2 + 1)
        if abs(R - 1.0) < 1e-6:
            A = 2.0 / P - 2.0 + sqrt_r2p1
            B = 2.0 / P - 2.0 - sqrt_r2p1
            # 修正：严格拦截无效边界（B <= 0 意味着存在严重的温度交叉，在 1-2 型管壳中不可行），不允许使用 abs 取绝对值掩盖温度交叉错误
            if A <= 0 or B <= 0: return 0.0
            arg = A / B
            if arg < 1e-12 or abs(math.log(arg)) < 1e-12: return 0.0
            F = (sqrt_r2p1 * P) / ((1 - P) * math.log(arg))
        else:
            arg_num = (1 - P) / (1 - P * R) if (1 - P * R) > 0 else -1.0
            if arg_num <= 0: return 0.0
            A = 2 / P - (1 + R) + sqrt_r2p1
            B = 2 / P - (1 + R) - sqrt_r2p1
            # 修正：严格拦截无效边界（A <= 0 或 B <= 0 意味着发生了温度交叉现象），杜绝强取绝对值计算出虚假的良好设计
            if A <= 0 or B <= 0: return 0.0
            arg_den = A / B
            F = (sqrt_r2p1 / (R - 1)) * math.log(arg_num) / math.log(arg_den)
        return max(0.0, min(1.0, F))
    except: return 0.0

def _F_bowman_n_shell(R: float, P: float, n_shells: int = 2) -> float:
    try:
        if abs(R - 1.0) < 1e-6:
            P_1 = P / (n_shells - (n_shells - 1) * P)
        else:
            ratio = ((1 - P * R) / (1 - P)) ** (1.0 / n_shells)
            if ratio <= 0: return 0.0
            P_1 = (ratio - 1) / (ratio - R)
        return _F_bowman_1_shell(R, P_1)
    except: return 0.0

def _F_cross_flow(R: float, P: float) -> float:
    try:
        if P <= 0 or P >= 1: return 0.0
        NTU = -math.log(1 - P)
        if NTU <= 0: return 1.0
        if abs(R - 1.0) < 1e-6:
            F = NTU / (2 * math.log(1 + NTU / 2))
        else:
            F = (1 / (R - 1)) * math.log((1 - P) / (1 - P * R)) / NTU
        return max(0.0, min(1.0, F))
    except: return 0.0

def suggest_configuration(T_hot_in: float, T_hot_out: float, T_cold_in: float, T_cold_out: float, min_F: float = 0.75) -> dict:
    configs = [("1_shell_2_tube", 1), ("2_shell_4_tube", 2), ("3_shell_6_tube", 3), ("4_shell_8_tube", 4)]
    results = []
    for config, n in configs:
        r = _compute_F_factor(T_hot_in, T_hot_out, T_cold_in, T_cold_out, config)
        if r["valid"]: results.append({"config": config, "n_shells": n, "F": r["F"], "quality": r.get("quality", "")})
    
    if not results:
        return {"recommended_config": "custom", "F_factor": 0.0, "n_shells": 0, "alternatives": [], "reasoning": "无法找到合适配置"}

    for r in results:
        if r["F"] >= min_F:
            alts = [x["config"] for x in results if x["F"] >= min_F]
            return {"recommended_config": r["config"], "F_factor": r["F"], "n_shells": r["n_shells"], "alternatives": alts, "reasoning": f"满足 F>={min_F}"}
    
    best = max(results, key=lambda x: x["F"])
    return {"recommended_config": best["config"], "F_factor": best["F"], "n_shells": best["n_shells"], "alternatives": [x["config"] for x in results], "reasoning": "全部低于下限"}


# ============================================================
# 相变检测与潜热计算
# ============================================================

def detect_phase_change(
    fluid_name: str,
    T_in: float,
    T_out: float,
    P: float,
    side: str = "hot",
) -> dict:
    """检测流体是否发生相变

    根据入口/出口温度与饱和温度的关系判断相变类型。

    Args:
        fluid_name: 流体名称 (如 "water", "methanol")
        T_in: 入口温度 (K)
        T_out: 出口温度 (K)
        P: 操作压力 (Pa)
        side: "hot" 或 "cold"

    Returns:
        dict: {
            "has_phase_change": bool,      # 是否发生相变
            "phase_type": str,             # "condensing" / "evaporating" / "single"
            "T_sat": float,                # 饱和温度 (K)
            "latent_heat": float,          # 潜热 (J/kg)
            "superheat": float,            # 过热度 (K)
            "subcool": float,              # 过冷度 (K)
            "phase_fraction": float,       # 相变分率 (0-1)
        }
    """
    from physics_engine.thermo_helper import calc_vapor_pressure, calc_enthalpy_vaporization, get_molecular_weight

    result = {
        "has_phase_change": False,
        "phase_type": "single",
        "T_sat": None,
        "latent_heat": 0.0,
        "superheat": 0.0,
        "subcool": 0.0,
        "phase_fraction": 0.0,
    }

    # 查找饱和温度 (通过蒸气压 = 操作压力反推)
    # 使用二分法求解 T_sat
    T_min, T_max = 200.0, 700.0
    T_sat = None

    for _ in range(50):
        T_mid = (T_min + T_max) / 2
        try:
            P_sat = calc_vapor_pressure(fluid_name, T_mid)
            if P_sat is None or P_sat <= 0:
                break
            if abs(P_sat - P) < 1.0:  # 收敛
                T_sat = T_mid
                break
            elif P_sat < P:
                T_min = T_mid
            else:
                T_max = T_mid
        except Exception:
            break

    if T_sat is None:
        return result

    result["T_sat"] = T_sat

    # 获取汽化潜热
    try:
        dH_vap_molar = calc_enthalpy_vaporization(fluid_name, T_sat)  # J/mol
        if dH_vap_molar and dH_vap_molar > 0:
            MW = get_molecular_weight(fluid_name)  # g/mol
            if MW and MW > 0:
                result["latent_heat"] = dH_vap_molar / (MW / 1000.0)  # J/kg
    except Exception:
        pass

    # 判断相变类型
    T_high = max(T_in, T_out)
    T_low = min(T_in, T_out)

    if side == "hot":
        # 热侧：温度下降
        if T_high > T_sat > T_low:
            # 发生冷凝
            result["has_phase_change"] = True
            result["phase_type"] = "condensing"
            # 计算过热度 (气相入口) 和过冷度 (液相出口)
            if T_in > T_sat:
                result["superheat"] = T_in - T_sat
            if T_out < T_sat:
                result["subcool"] = T_sat - T_out
            # 相变分率：简化为完全冷凝 (1.0)
            # 实际工程中可通过热平衡精确计算
            result["phase_fraction"] = 1.0
        elif T_low >= T_sat:
            # 全气相，无相变
            pass
        elif T_high <= T_sat:
            # 全液相，无相变
            pass
    else:
        # 冷侧：温度上升
        if T_low < T_sat < T_high:
            # 发生蒸发/沸腾
            result["has_phase_change"] = True
            result["phase_type"] = "evaporating"
            if T_in < T_sat:
                result["subcool"] = T_sat - T_in
            if T_out > T_sat:
                result["superheat"] = T_out - T_sat
            # 相变分率：简化为完全汽化 (1.0)
            result["phase_fraction"] = 1.0

    return result


# ============================================================
# 热负荷计算
# ============================================================

def heat_duty(
    m_hot: float = 0.0, Cp_hot: float = 0.0, T_hot_in: float = 0.0, T_hot_out: float = 0.0,
    m_cold: float = 0.0, Cp_cold: float = 0.0, T_cold_in: float = 0.0, T_cold_out: float = 0.0,
    latent_heat_hot: float = 0.0, latent_heat_cold: float = 0.0,
    phase_hot: str = "single", phase_cold: str = "single",
    hot_subcool: float = 0.0, cold_superheat: float = 0.0,
) -> dict:
    warnings: list[str] = []

    Q_hot_sensible = m_hot * Cp_hot * abs(T_hot_in - T_hot_out)
    Q_hot_phase = m_hot * latent_heat_hot if latent_heat_hot > 0 else 0.0
    
    # 修正：相变总热负荷为相变焓与过热段降温/冷凝液过冷段降温等显热负荷的严密物理之和，不应粗暴舍弃
    if phase_hot in ("condensing", "evaporating"):
        Q_hot = Q_hot_phase + Q_hot_sensible
    else: 
        Q_hot = Q_hot_sensible

    Q_cold_sensible = m_cold * Cp_cold * abs(T_cold_out - T_cold_in)
    Q_cold_phase = m_cold * latent_heat_cold if latent_heat_cold > 0 else 0.0

    # 修正：冷流体相变总热负荷为液体预热/沸腾相变/气体过热段的物理之和，不应粗暴舍弃预热显热
    if phase_cold in ("evaporating", "condensing"):
        Q_cold = Q_cold_phase + Q_cold_sensible
    else: 
        Q_cold = Q_cold_sensible

    Q_avg = (Q_hot + Q_cold) / 2.0 if (Q_hot + Q_cold) > 0 else 0.0
    balance_error = abs(Q_hot - Q_cold) / max(Q_hot, Q_cold, 1e-6) * 100 if max(Q_hot, Q_cold) > 0 else 0.0

    if balance_error > 10: warnings.append(f"热平衡偏差 {balance_error:.1f}% > 10%")
    elif balance_error > 5: warnings.append(f"热平衡偏差 {balance_error:.1f}% > 5%")

    return {
        "Q_hot": Q_hot, "Q_cold": Q_cold, "Q_avg": Q_avg,
        "Q_hot_sensible": Q_hot_sensible, "Q_hot_phase": Q_hot_phase,
        "Q_cold_sensible": Q_cold_sensible, "Q_cold_phase": Q_cold_phase,
        "balance_error": balance_error, "phase_hot": phase_hot, "phase_cold": phase_cold,
        "valid": Q_avg > 0, "warnings": warnings,
    }


# ============================================================
# 总传热系数估算（支持混合物与动态污垢）
# ============================================================

def overall_htc_estimate(
    h_tube: float,
    h_shell: float,
    d_o: float = 0.025,
    d_i: float = 0.021,
    k_wall: float = 50.0,
    fouling_tube: str | tuple | list = "default",
    fouling_shell: str | tuple | list = "default",
    z_tube: list[float] | None = None,
    z_shell: list[float] | None = None,
    R_f_tube: float | None = None,
    R_f_shell: float | None = None,
    T_heat_C: float = 100.0,
    v_tube: float = 1.0,
    v_shell: float = 1.0,
) -> dict:
    """
    总传热系数估算，动态计算并加入污垢热阻。
    """
    Rf_t = R_f_tube if R_f_tube is not None else calculate_fouling_resistance(fouling_tube, z_tube, T_heat_C, v_tube)
    Rf_s = R_f_shell if R_f_shell is not None else calculate_fouling_resistance(fouling_shell, z_shell, T_heat_C, v_shell)

    R_conv_tube = (1.0 / h_tube) * (d_o / d_i) if (h_tube > 0 and d_i > 0) else 0.0
    R_conv_shell = 1.0 / h_shell if h_shell > 0 else 0.0
    R_wall = (d_o * math.log(d_o / d_i)) / (2 * k_wall) if (d_o > 0 and d_i > 0 and k_wall > 0) else 0.0
    
    R_foul_t = Rf_t * (d_o / d_i) if d_i > 0 else Rf_t
    R_foul_s = Rf_s

    R_total = R_conv_tube + R_foul_t + R_wall + R_foul_s + R_conv_shell
    U = 1.0 / R_total if R_total > 0 else 0.0

    fouling_fraction = (R_foul_t + R_foul_s) / R_total * 100 if R_total > 0 else 0.0

    return {
        "U": U,
        "R_total": R_total,
        "R_tube": R_conv_tube,
        "R_shell": R_conv_shell,
        "R_wall": R_wall,
        "R_fouling_tube": R_foul_t,
        "R_fouling_shell": R_foul_s,
        "R_fouling_total": R_foul_t + R_foul_s,
        "fouling_fraction": fouling_fraction,
        "valid": U > 0,
    }


def typical_U_range(category: str = "water/water") -> dict:
    parts  = category.split("/")
    key    = tuple(p.strip() for p in parts)
    key_rev= (key[1], key[0]) if len(key) == 2 else None

    for k in (key, key_rev):
        if k and k in _TYPICAL_U:
            lo, hi = _TYPICAL_U[k]
            return {"category": category, "U_min": lo, "U_max": hi, "U_typical": round((lo + hi) / 2), "valid": True}

    _CATEGORY_FALLBACK = {
        "organic": "medium_organic", "light organic": "light_organic",
        "heavy organic": "heavy_organic", "brine": "brine",
        "condensing_steam": "condensing_steam", "steam": "condensing_steam",
    }
    for i, part in enumerate(key):
        fallback = _CATEGORY_FALLBACK.get(part)
        if fallback:
            new_key = list(key); new_key[i] = fallback; new_key = tuple(new_key)
            new_rev = (new_key[1], new_key[0])
            for k in (new_key, new_rev):
                if k in _TYPICAL_U:
                    lo, hi = _TYPICAL_U[k]
                    return {"category": category, "U_min": lo, "U_max": hi, "U_typical": round((lo + hi) / 2), "valid": True, "note": f"降级匹配到 {k}"}

    return {"category": category, "U_min": 116, "U_max": 698, "U_typical": 350, "valid": False, "message": "使用通用范围"}


# ============================================================
# 换热面积计算
# ============================================================

def area_required(
    Q: float, U: float, LMTD: float,
    overspec_percent: float = 15.0,
    d_o: float = 0.025, tube_length: float = 6.0,
    tube_passes: int = 2, tube_layout: str = "triangular",
    # ── 流速约束参数 ──
    m_flow: float = 0.0,
    rho: float = 0.0,
    d_i: float = 0.0,
    V_min: float = 0.5,
    # ── 强制参数标志 ──
    force_tube_passes: bool = False,
    force_baffle_spacing: float = None,
) -> dict:
    warnings: list[str] = []
    if U <= 0 or LMTD <= 0:
        return {"A_required": 0, "A_design": 0, "overspec_percent": 0, "N_tubes": 0, "D_shell_est": 0, "length_diameter_ratio": 0, "tube_velocity_estimate": None, "valid": False, "warnings": ["U 或 LMTD 为零"]}

    A_required = Q / (U * LMTD)
    A_design = A_required * (1 + overspec_percent / 100.0)

    single_tube_area = math.pi * d_o * tube_length if d_o > 0 and tube_length > 0 else 0.01
    N_tubes = max(1, int(math.ceil(A_design / single_tube_area)))

    # ── 流速约束校核与自动调整 ──
    tube_velocity = None
    velocity_constrained = False
    if m_flow > 0 and rho > 0 and d_i > 0:
        A_single_tube = math.pi / 4 * d_i ** 2

        def _calc_velocity(n_tubes, passes):
            n_per_pass = n_tubes // passes
            if n_per_pass < 1:
                return 0.0
            return m_flow / (rho * A_single_tube * n_per_pass)

        tube_velocity = _calc_velocity(N_tubes, tube_passes)

        if tube_velocity < V_min and not force_tube_passes:
            velocity_constrained = True
            # 策略1: 增加管程数 → 减少每程管数 → 提高流速
            for new_passes in [4, 6, 8]:
                n_per_pass = N_tubes // new_passes
                if n_per_pass < 1:
                    continue
                v_new = _calc_velocity(N_tubes, new_passes)
                if v_new >= V_min:
                    tube_passes = new_passes
                    tube_velocity = v_new
                    break
            else:
                # 策略2: 减少管数 + 增加管长补偿面积
                n_per_pass_min = max(1, math.ceil(m_flow / (rho * A_single_tube * V_min * tube_passes)))
                N_tubes_min = n_per_pass_min * tube_passes
                if N_tubes_min < N_tubes:
                    A_total_min = N_tubes_min * single_tube_area
                    if A_total_min < A_design:
                        tube_length_needed = A_design / (N_tubes_min * math.pi * d_o)
                        for L_std in [0.5, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0, 9.0]:
                            if L_std >= tube_length_needed:
                                tube_length = L_std
                                single_tube_area = math.pi * d_o * tube_length
                                break
                    N_tubes = N_tubes_min
                    tube_velocity = _calc_velocity(N_tubes, tube_passes)

    pitch_ratio = 1.25
    k_layout = 0.65 if tube_layout == "triangular" else 0.80
    D_shell_est = d_o * pitch_ratio * math.sqrt(N_tubes / k_layout) if N_tubes > 0 else 0.0
    D_shell_std = _round_to_standard(D_shell_est, _STANDARD_SHELL_DIAMETER)

    L_D = tube_length / D_shell_std if D_shell_std > 0 else 0
    if L_D < 3: warnings.append(f"长径比 {L_D:.1f} < 3，换热器偏短")
    elif L_D > 15: warnings.append(f"长径比 {L_D:.1f} > 15，换热器偏长")

    # ── 长径比约束自动调整 ──
    LD_CONSTRAINED = False
    _STANDARD_TUBE_LENGTHS = [1.0, 1.5, 2.0, 3.0, 4.5, 6.0, 9.0]

    if L_D > 15 and D_shell_std > 0:
        # L/D 过大 → 减小管长（增大壳径方向调整）
        LD_CONSTRAINED = True
        for new_L in _STANDARD_TUBE_LENGTHS:
            if new_L < tube_length:
                # 重新计算: 更短的管 → 更多的管数 → 更大的壳径
                new_single_tube_area = math.pi * d_o * new_L
                new_N_tubes = max(1, int(math.ceil(A_design / new_single_tube_area)))
                new_D_shell_est = d_o * pitch_ratio * math.sqrt(new_N_tubes / k_layout) if new_N_tubes > 0 else 0.0
                new_D_shell_std = _round_to_standard(new_D_shell_est, _STANDARD_SHELL_DIAMETER)
                new_L_D = new_L / new_D_shell_std if new_D_shell_std > 0 else 99
                if new_L_D <= 15:
                    warnings.append(f"[L/D约束] 管长 {tube_length:.1f}→{new_L:.1f}m, N_tubes {N_tubes}→{new_N_tubes}, D_shell {D_shell_std:.2f}→{new_D_shell_std:.2f}m, L/D {L_D:.1f}→{new_L_D:.1f}")
                    tube_length = new_L
                    N_tubes = new_N_tubes
                    D_shell_est = new_D_shell_est
                    D_shell_std = new_D_shell_std
                    single_tube_area = new_single_tube_area
                    L_D = new_L_D
                    break
        else:
            # 所有标准管长都不满足，用最短的
            new_L = _STANDARD_TUBE_LENGTHS[0]
            new_single_tube_area = math.pi * d_o * new_L
            new_N_tubes = max(1, int(math.ceil(A_design / new_single_tube_area)))
            new_D_shell_est = d_o * pitch_ratio * math.sqrt(new_N_tubes / k_layout) if new_N_tubes > 0 else 0.0
            new_D_shell_std = _round_to_standard(new_D_shell_est, _STANDARD_SHELL_DIAMETER)
            new_L_D = new_L / new_D_shell_std if new_D_shell_std > 0 else 99
            warnings.append(f"[L/D约束] 管长缩至最小 {tube_length:.1f}→{new_L:.1f}m, L/D {L_D:.1f}→{new_L_D:.1f} (仍超限)")
            tube_length = new_L
            N_tubes = new_N_tubes
            D_shell_est = new_D_shell_est
            D_shell_std = new_D_shell_std
            single_tube_area = new_single_tube_area
            L_D = new_L_D

    elif L_D < 3 and D_shell_std > 0:
        # L/D 过小 → 增加管长（减少管数→减小壳径）
        LD_CONSTRAINED = True
        for new_L in _STANDARD_TUBE_LENGTHS:
            if new_L > tube_length:
                new_single_tube_area = math.pi * d_o * new_L
                new_N_tubes = max(1, int(math.ceil(A_design / new_single_tube_area)))
                new_D_shell_est = d_o * pitch_ratio * math.sqrt(new_N_tubes / k_layout) if new_N_tubes > 0 else 0.0
                new_D_shell_std = _round_to_standard(new_D_shell_est, _STANDARD_SHELL_DIAMETER)
                new_L_D = new_L / new_D_shell_std if new_D_shell_std > 0 else 99
                if new_L_D >= 3:
                    warnings.append(f"[L/D约束] 管长 {tube_length:.1f}→{new_L:.1f}m, N_tubes {N_tubes}→{new_N_tubes}, D_shell {D_shell_std:.2f}→{new_D_shell_std:.2f}m, L/D {L_D:.1f}→{new_L_D:.1f}")
                    tube_length = new_L
                    N_tubes = new_N_tubes
                    D_shell_est = new_D_shell_est
                    D_shell_std = new_D_shell_std
                    single_tube_area = new_single_tube_area
                    L_D = new_L_D
                    break
        else:
            # 所有标准管长都不满足，用最长的
            new_L = _STANDARD_TUBE_LENGTHS[-1]
            new_single_tube_area = math.pi * d_o * new_L
            new_N_tubes = max(1, int(math.ceil(A_design / new_single_tube_area)))
            new_D_shell_est = d_o * pitch_ratio * math.sqrt(new_N_tubes / k_layout) if new_N_tubes > 0 else 0.0
            new_D_shell_std = _round_to_standard(new_D_shell_est, _STANDARD_SHELL_DIAMETER)
            new_L_D = new_L / new_D_shell_std if new_D_shell_std > 0 else 99
            warnings.append(f"[L/D约束] 管长增至最大 {tube_length:.1f}→{new_L:.1f}m, L/D {L_D:.1f}→{new_L_D:.1f} (仍超限)")
            tube_length = new_L
            N_tubes = new_N_tubes
            D_shell_est = new_D_shell_est
            D_shell_std = new_D_shell_std
            single_tube_area = new_single_tube_area
            L_D = new_L_D
    if overspec_percent < 10: warnings.append(f"面积裕度 {overspec_percent:.0f}% < 10%，偏小")
    elif overspec_percent > 30: warnings.append(f"面积裕度 {overspec_percent:.0f}% > 30%，偏大")
    if velocity_constrained:
        warnings.append(f"流速约束 V≥{V_min}m/s 已触发，管程数调整为 {tube_passes}")

    return {
        "A_required": A_required, "A_design": A_design, "overspec_percent": overspec_percent,
        "N_tubes": N_tubes, "D_shell_est": D_shell_est, "D_shell_standard": D_shell_std,
        "tube_length": tube_length, "tube_passes": tube_passes, "length_diameter_ratio": L_D,
        "tube_velocity": tube_velocity, "V_min": V_min if (m_flow > 0 and rho > 0 and d_i > 0) else None,
        "valid": A_required > 0, "warnings": warnings,
    }


# ============================================================
# 压降计算（管程 + 壳程）
# ============================================================

def pressure_drop(
    side: str, m: float, rho: float, mu: float, d: float, L: float,
    n_passes: int = 1, roughness: float = 4.5e-5, D_shell: float = 0.0,
    B: float = 0.0, N_baffles: int = 0, d_o: float = 0.025, pitch: float = 0.03125,
    n_tubes: int = 1, layout: str = "triangular",
) -> dict:
    warnings: list[str] = []
    if side == "tube":
        return _pressure_drop_tube(m, rho, mu, d, L, n_passes, roughness, warnings, n_tubes)
    else:
        return _pressure_drop_shell(m, rho, mu, d, D_shell, L, B, N_baffles, d_o, pitch, layout, warnings)

def _pressure_drop_tube(m: float, rho: float, mu: float, d_i: float, L: float, n_passes: int, roughness: float, warnings: list, n_tubes: int = 1) -> dict:
    # 每程并行管数
    n_per_pass = max(1, n_tubes // max(1, n_passes))
    A_flow = n_per_pass * math.pi / 4 * d_i ** 2
    u = m / (rho * A_flow) if (rho > 0 and A_flow > 0) else 0.0
    Re = rho * u * d_i / mu if mu > 0 else 0.0
    f = _churchill_friction(Re, roughness / d_i if d_i > 0 else 0)
    dP_friction = f * (L / d_i) * (rho * u ** 2 / 2) * n_passes if d_i > 0 else 0.0
    dP_return = 4 * (rho * u ** 2 / 2) * (n_passes - 1)
    dP_nozzle = 1.5 * (rho * u ** 2 / 2)
    dP_total = dP_friction + dP_return + dP_nozzle

    if dP_total > 100000: warnings.append(f"管程压降 {dP_total/1000:.1f} kPa > 100 kPa")
    return {"delta_P": dP_total, "delta_P_friction": dP_friction, "delta_P_return": dP_return, "delta_P_nozzle": dP_nozzle, "friction_factor": f, "velocity": u, "Re": Re, "valid": Re > 0, "warnings": warnings}

def calculate_shell_side_pressure_drop(
    m: float,
    rho: float,
    mu: float,
    D_shell: float,
    d_o: float,
    pitch: float,
    layout: str = "triangular",
    baffle_spacing: float = None,
    tube_length: float = None,
    n_baffles: int = None,
) -> dict:
    """壳侧压降计算（Kern 法，折流板分布）
    """
    from physics_engine.common import safe_float, safe_int
    
    # ── 防御性类型转换 ──
    m = safe_float(m, name='m')
    rho = safe_float(rho, name='rho')
    mu = safe_float(mu, name='mu')
    D_shell = safe_float(D_shell, name='D_shell')
    d_o = safe_float(d_o, name='d_o')
    pitch = safe_float(pitch, name='pitch')
    baffle_spacing = safe_float(baffle_spacing, name='baffle_spacing')
    tube_length = safe_float(tube_length, name='tube_length')
    n_baffles = safe_int(n_baffles, name='n_baffles')
    
    warnings = []
    
    # 自动计算折流板参数
    if baffle_spacing is None and tube_length is not None and n_baffles is not None:
        B = tube_length / (n_baffles + 1)
    elif baffle_spacing is not None:
        B = baffle_spacing
    else:
        B = D_shell  # 默认
    
    if n_baffles is None and tube_length is not None and B > 0:
        n_baffles = int(tube_length / B) - 1
    
    # 当量直径
    if layout == "triangular":
        d_e = 4 * (0.866 * pitch ** 2 - math.pi * d_o ** 2 / 4) / (math.pi * d_o)
    else:
        d_e = 4 * (pitch ** 2 - math.pi * d_o ** 2 / 4) / (math.pi * d_o)
    
    L = tube_length if tube_length is not None else D_shell * 3
    
    result = _pressure_drop_shell(m, rho, mu, d_e, D_shell, L, B, n_baffles, d_o, pitch, layout, warnings)
    result["delta_P_Pa"] = result.pop("delta_P")
    result["velocity_m_s"] = result.pop("velocity")
    return result


def _pressure_drop_shell(m: float, rho: float, mu: float, d_e: float, D_shell: float, L: float, B: float, N_baffles: int, d_o: float, pitch: float, layout: str, warnings: list) -> dict:
    if D_shell <= 0 or B <= 0: return {"delta_P": 0, "friction_factor": 0, "velocity": 0, "Re": 0, "valid": False, "warnings": ["缺少壳程几何参数"]}
    
    # 如果 d_e 不合理（<5mm），根据管排列方式重算
    if d_e < 0.005 and pitch > 0 and d_o > 0:
        if layout == "triangular":
            d_e = 4 * (0.866 * pitch ** 2 - math.pi * d_o ** 2 / 4) / (math.pi * d_o)
        else:
            d_e = 4 * (pitch ** 2 - math.pi * d_o ** 2 / 4) / (math.pi * d_o)

    C = pitch - d_o
    A_shell = (D_shell * C * B / pitch) if C > 0 else (D_shell * B * 0.25)
    u_s = m / (rho * A_shell) if (rho > 0 and A_shell > 0) else 0.0
    Re_s = rho * u_s * d_e / mu if (mu > 0 and d_e > 0) else 0.0

    f_s = (1.728 * Re_s ** (-0.188) if Re_s > 100 else 100.0 / Re_s) if Re_s > 0 else 0.0
    N_cross = N_baffles + 1 if N_baffles > 0 else int(L / B) + 1
    dP_shell = f_s * (D_shell / d_e) * N_cross * (rho * u_s ** 2 / 2) if d_e > 0 else 0.0

    if dP_shell > 70000: warnings.append(f"壳程压降 {dP_shell/1000:.1f} kPa > 70 kPa")
    return {"delta_P": dP_shell, "friction_factor": f_s, "velocity": u_s, "Re": Re_s, "d_equivalent": d_e, "N_cross_passes": N_cross, "N_baffles": N_baffles, "valid": Re_s > 0, "warnings": warnings}

def _churchill_friction(Re: float, rel_rough: float = 0.0) -> float:
    if Re < 1: return 64.0
    if Re < 2300: return 64.0 / Re
    if Re > 4000:
        if rel_rough > 0: return 0.25 / (math.log10(rel_rough / 3.7 + 5.74 / Re ** 0.9)) ** 2
        else: return 0.316 / Re ** 0.25
    return (64.0 / 2300) * (1 - (Re - 2300) / 1700) + (0.316 / 4000 ** 0.25) * ((Re - 2300) / 1700)


# ============================================================
# ε-NTU 法
# ============================================================

def epsilon_ntu(
    C_hot: float, C_cold: float, UA: float,
    T_hot_in: float = None, T_cold_in: float = None,
    flow_arrangement: str = "counter",
) -> dict:
    """ε-NTU 法完整计算：效能 ε、传热量 Q、出口温度。
    """
    C_min = min(C_hot, C_cold)
    C_max = max(C_hot, C_cold)
    C_r = C_min / C_max if C_max > 0 else 0.0
    NTU = UA / C_min if C_min > 0 else 0.0

    if flow_arrangement == "counter":
        if C_r < 1e-6: epsilon = 1 - math.exp(-NTU)
        elif abs(C_r - 1.0) < 1e-6: epsilon = NTU / (1 + NTU)
        else: epsilon = (1 - math.exp(-NTU * (1 - C_r))) / (1 - C_r * math.exp(-NTU * (1 - C_r)))
    elif flow_arrangement == "parallel":
        epsilon = (1 - math.exp(-NTU * (1 + C_r))) / (1 + C_r)
    elif flow_arrangement == "cross_flow":
        if C_r < 1e-6: epsilon = 1 - math.exp(-NTU)
        else: epsilon = 1 - math.exp((math.exp(-NTU * C_r) - 1) / C_r)
    else: epsilon = 0.0

    epsilon = max(0.0, min(1.0, epsilon))

    # 最大可能传热量和实际传热量
    T_hot_out = None
    T_cold_out = None
    if T_hot_in is not None and T_cold_in is not None and C_min > 0:
        delta_T_max = T_hot_in - T_cold_in
        Q_max = C_min * delta_T_max
        Q_actual = epsilon * Q_max
        T_hot_out = T_hot_in - Q_actual / C_hot if C_hot > 0 else T_hot_in
        T_cold_out = T_cold_in + Q_actual / C_cold if C_cold > 0 else T_cold_in
    else:
        Q_max = C_min * 1.0  
        Q_actual = epsilon * Q_max

    return {
        "epsilon": epsilon, "NTU": NTU, "C_r": C_r,
        "C_min": C_min, "C_max": C_max,
        "Q_max": Q_max, "Q_actual": Q_actual,
        "T_hot_out": T_hot_out, "T_cold_out": T_cold_out,
        "valid": NTU > 0,
    }


# ============================================================
# 对流换热系数
# ============================================================

def tube_side_htc(
    m: float, rho: float, mu: float, k: float, Cp: float, d_i: float,
    n_tubes: int = 1, n_passes: int = 1, tube_length: float = 6.0, mu_wall: float | None = None,
) -> dict:
    A_flow = n_tubes / n_passes * math.pi / 4 * d_i ** 2
    u = m / (rho * A_flow) if (rho > 0 and A_flow > 0) else 0.0
    Re = rho * u * d_i / mu if mu > 0 else 0.0
    Pr = mu * Cp / k if k > 0 else 0.0

    if Re < 2300:
        L_eff = max(tube_length, d_i)
        Nu = 1.86 * (Re * Pr * d_i / L_eff) ** (1 / 3)
        if mu_wall and mu > 0: Nu *= (mu / mu_wall) ** 0.14
    elif Re < 4000:
        Nu_lam = 1.86 * (2300 * Pr) ** (1 / 3)
        Nu_turb = 0.023 * Re ** 0.8 * Pr ** 0.4
        Nu = Nu_lam * (1 - (Re - 2300) / 1700) + Nu_turb * ((Re - 2300) / 1700)
    else:
        Nu = 0.023 * Re ** 0.8 * Pr ** 0.4
        if mu_wall and mu > 0: Nu *= (mu / mu_wall) ** 0.14 

    h = Nu * k / d_i if d_i > 0 else 0.0
    return {"h": h, "Nu": Nu, "Re": Re, "Pr": Pr, "velocity": u, "flow_regime": "laminar" if Re < 2300 else ("transition" if Re < 4000 else "turbulent"), "valid": h > 0}

# ============================================================
# 折流板设计
# ============================================================

def design_baffles(D_shell: float, L_tube: float, baffle_cut: float = 0.25) -> dict:
    """折流板设计：计算间距、数量

    TEMA 推荐折流板间距 B = 0.2~1.0 × D_shell，常用 0.4 × D_shell。
    弓形折流板切口一般取 20~35%，默认 25%。
    """
    if D_shell <= 0 or L_tube <= 0:
        return {"baffle_spacing_m": 0.0, "N_baffles": 0, "N_cross_passes": 1, "baffle_cut": baffle_cut, "valid": False}

    B = 0.4 * D_shell
    B = max(0.15, min(B, 1.0))  # 典型工业限制
    N_baffles = max(1, int(L_tube / B) - 1)
    
    return {
        "baffle_spacing_m": B,
        "N_baffles": N_baffles,
        "N_cross_passes": N_baffles + 1,
        "baffle_cut": baffle_cut,
        "valid": True
    }


# ============================================================
# 传热系数计算辅助
# ============================================================

def viscosity_correction(mu_wall: float, mu_bulk: float, n: float = 0.14) -> float:
    """黏度修正因子
    修正：Sieder-Tate 黏度修正分子分母位置颠倒，应为 (mu_bulk / mu_wall) ** n，符合传统教材
    """
    if mu_wall <= 0: return 1.0
    return (mu_bulk / mu_wall) ** n


# ============================================================
# 壳程对流传热系数
# ============================================================

def shell_side_htc(
    m: float, rho: float, mu: float, k: float, Cp: float, D_shell: float,
    B: float, d_o: float, pitch: float = 0.03125, layout: str = "triangular",
) -> dict:
    C = pitch - d_o
    A_s = D_shell * C * B / pitch if pitch > 0 else 0.0
    u_s = m / (rho * A_s) if (rho > 0 and A_s > 0) else 0.0
    
    # 壳程当量直径
    if layout == "triangular":
        d_e = 4 * (0.866 * pitch ** 2 - 0.5 * math.pi * d_o ** 2 / 4) / (0.5 * math.pi * d_o)
    else:
        d_e = 4 * (pitch ** 2 - math.pi * d_o ** 2 / 4) / (math.pi * d_o)
        
    Re_s = rho * u_s * d_e / mu if (mu > 0 and d_e > 0) else 0.0
    Pr_s = mu * Cp / k if k > 0 else 0.0

    if Re_s > 0:
        j_H = 0.36 * Re_s ** (-0.39) if Re_s > 200 else 0.5 * Re_s ** (-0.5)
        Nu_s = j_H * Re_s * Pr_s ** (1 / 3)
    else: Nu_s = 0.0

    h_s = Nu_s * k / d_e if (d_e > 0 and k > 0) else 0.0
    return {"h": h_s, "Nu": Nu_s, "Re": Re_s, "Pr": Pr_s, "velocity": u_s, "d_equivalent": d_e, "flow_regime": "laminar" if Re_s < 200 else "turbulent", "valid": h_s > 0}


# ============================================================
# 换热器类型与相变支持
# ============================================================

def condensate_heat_transfer(T_wall: float, T_sat: float, rho_l: float, rho_v: float, mu_l: float, k_l: float, h_fg: float, L: float = 1.0, angle: float = 90.0, Cp_l: float = 4180.0) -> dict:
    delta_T = T_sat - T_wall
    if delta_T <= 0: return {"h": 0.0, "valid": False, "error": "Delta T <= 0"}
    h_fg_mod = h_fg + 0.68 * Cp_l * delta_T
    C = 0.725 if abs(angle - 90) < 1 else 0.943
    try: h = C * ((rho_l * (rho_l - rho_v) * 9.81 * h_fg_mod * k_l**3) / (mu_l * delta_T * L)) ** 0.25
    except: h = 0.0
    return {"h": float(h), "valid": h > 0, "delta_T": delta_T}

def boiling_heat_transfer(delta_T: float, C_sf: float = 0.0132, n: float = 1.0, mu_l: float = 2.82e-4, h_fg: float = 2257e3, rho_v: float = 0.598, rho_l: float = 958, sigma: float = 0.0589, C_p_l: float = 4217, Pr_l: float = 1.76) -> dict:
    if delta_T <= 0: return {"h": 0.0, "q_flux": 0.0, "valid": False}
    try:
        q_flux = mu_l * h_fg * math.sqrt(9.81 * (rho_l - rho_v) / sigma) * ((C_p_l * delta_T) / (C_sf * h_fg * Pr_l**n)) ** 3
        h = q_flux / delta_T
    except: h = 0.0; q_flux = 0.0
    return {"h": float(h), "q_flux": float(q_flux), "delta_T": float(delta_T), "valid": h > 0}

def sieder_tate_nu(Re: float, Pr: float, mu_wall: float, mu_bulk: float, d_i: float = 0.021, L: float = 3.0) -> float:
    # 修正：本处黏度修正一并修正为 (mu_bulk / mu_wall) ** 0.14
    if Re < 2300: Nu = 1.86 * (Re * Pr * (d_i / L)) ** (1/3) * (mu_bulk / mu_wall) ** 0.14
    elif Re < 10000: Nu = 0.116 * (Re**0.67 - 125) * Pr**0.43 * (1 + (d_i / L)**(2/3)) * (mu_bulk / mu_wall)**0.14
    else: Nu = 0.023 * Re**0.8 * Pr**(1/3) * (mu_bulk / mu_wall)**0.14
    return float(Nu)

def fouling_impact(U_clean: float, R_f_tube: float = 0.00018, R_f_shell: float = 0.00018, A: float = 1.0, delta_t: float = 8000) -> dict:
    R_f_total = R_f_tube + R_f_shell
    U_fouled = 1.0 / (1.0 / U_clean + R_f_total)
    efficiency = U_fouled / U_clean if U_clean > 0 else 1.0
    fouling_rate = R_f_total / delta_t 
    if efficiency < 0.8: action = "Clean immediately"
    elif efficiency < 0.9: action = "Schedule cleaning"
    else: action = "Normal"
    return {"U_fouled": float(U_fouled), "U_clean": float(U_clean), "efficiency": float(efficiency), "R_f_total": float(R_f_total), "fouling_rate": float(fouling_rate), "action": action}

def heat_transfer_with_fouling(h_tube: float, h_shell: float, d_o: float = 0.025, d_i: float = 0.021, k_wall: float = 50.0, R_f_tube: float = 0.00018, R_f_shell: float = 0.00018) -> dict:
    R_conv_tube = (1.0 / h_tube) * (d_o / d_i) if (h_tube > 0 and d_i > 0) else 0.0
    R_conv_shell = 1.0 / h_shell if h_shell > 0 else 0.0
    R_wall = (d_o * math.log(d_o / d_i)) / (2 * k_wall) if (k_wall > 0 and d_i > 0) else 0.0
    R_foul_t = R_f_tube * (d_o / d_i) if d_i > 0 else R_f_tube
    U_clean = 1.0 / (R_conv_tube + R_conv_shell + R_wall) if (R_conv_tube + R_conv_shell + R_wall) > 0 else 0.0
    U_fouled = 1.0 / (R_conv_tube + R_conv_shell + R_wall + R_foul_t + R_f_shell) if (R_conv_tube + R_conv_shell + R_wall + R_foul_t + R_f_shell) > 0 else 0.0
    efficiency = U_fouled / U_clean if U_clean > 0 else 1.0
    return {"U": float(U_fouled), "U_clean": float(U_clean), "efficiency": float(efficiency), "fouling_resistance": R_f_tube + R_f_shell}

def select_exchanger_type(hot_fluid: str = "", cold_fluid: str = "", T_hot: float = 298.15, T_cold: float = 298.15, P_hot: float = 101325, P_cold: float = 101325, Q: float = 0, phase_hot: str = "single", phase_cold: str = "single", corrosion_risk: str = "low") -> dict:
    reasons: dict[str, int] = {}
    types = ["shell_and_tube_BEM", "shell_and_tube_AES", "shell_and_tube_AKM", "shell_and_tube_BEU", "double_pipe", "plate", "air_cooled", "spiral", "plate_fin", "spiral_tube"]
    for t in types: reasons[t] = 0

    if abs(T_hot - T_cold) > 50:
        reasons["shell_and_tube_BEM"] += 3; reasons["shell_and_tube_AES"] += 3
    if P_hot > 1e6 or P_cold > 1e6:
        reasons["shell_and_tube_BEM"] += 4; reasons["double_pipe"] += 2; reasons["plate"] -= 3
    if phase_hot == "condensing":
        reasons["shell_and_tube_AES"] += 3; reasons["shell_and_tube_AKM"] += 2; reasons["air_cooled"] += 2
    if phase_cold == "evaporating":
        reasons["shell_and_tube_AKM"] += 3; reasons["shell_and_tube_BEU"] += 2
    if corrosion_risk in ("high", "medium"):
        reasons["shell_and_tube_AES"] += 2; reasons["plate"] += 1
    if 0 < Q < 50000:
        reasons["double_pipe"] += 3; reasons["plate"] += 2
    if "air" in cold_fluid.lower() or "air" in hot_fluid.lower():
        reasons["air_cooled"] += 4

    sorted_types = sorted(reasons.items(), key=lambda x: x[1], reverse=True)
    return {"recommended": sorted_types[0][0], "score": sorted_types[0][1], "alternatives": [t[0] for t in sorted_types[1:3]] if len(sorted_types) > 1 else [], "reasoning": "; ".join([f"{r}({c})" for r, c in sorted_types]), "all_factors": reasons}


# ============================================================
# 完整换热器设计函数 (基于动态更新流速与动态折流板间距)
# ============================================================

def design_heat_exchanger(
    T_hot_in: float, T_hot_out: float,
    T_cold_in: float, T_cold_out: float,
    m_hot: float, m_cold: float,
    Cp_hot: float = _REQUIRED, Cp_cold: float = _REQUIRED,
    rho_hot: float = _REQUIRED, rho_cold: float = _REQUIRED,
    mu_hot: float = _REQUIRED, mu_cold: float = _REQUIRED,
    k_hot: float = _REQUIRED, k_cold: float = _REQUIRED,
    latent_heat_hot: float = 0.0, latent_heat_cold: float = 0.0,
    phase_hot: str = "single", phase_cold: str = "single",
    configuration: str = "1_shell_2_tube",
    d_o: float = 0.025, d_i: float = 0.021,
    tube_length: float = 6.0, tube_passes: int = 2,
    tube_layout: str = "triangular",
    overspec_percent: float = 15.0,
    k_wall: float = 50.0,
    fouling_tube: str | tuple | list = "default",
    fouling_shell: str | tuple | list = "default",
    z_tube: list[float] | None = None,
    z_shell: list[float] | None = None,
    baffle_spacing: float = 0.6,  # 修正：将 tray_spacing 纠正为 baffle_spacing（折流板间距）
    # ── 强制参数标志（用于约束迭代设计）──
    force_tube_passes: bool = False,
    force_baffle_spacing: float = None,
) -> dict:
    """
    完整换热器设计计算
    """
    # 检查必填参数
    missing = check_required_params(locals(), {
        "Cp_hot": "热流体比热容 (J/kg·K)，由设备层从物性引擎获取",
        "Cp_cold": "冷流体比热容 (J/kg·K)，由设备层从物性引擎获取",
        "rho_hot": "热流体密度 (kg/m³)，由设备层从物性引擎获取",
        "rho_cold": "冷流体密度 (kg/m³)，由设备层从物性引擎获取",
        "mu_hot": "热流体动力粘度 (Pa·s)，由设备层从物性引擎获取",
        "mu_cold": "冷流体动力粘度 (Pa·s)，由设备层从物性引擎获取",
        "k_hot": "热流体导热系数 (W/m·K)，由设备层从物性引擎获取",
        "k_cold": "冷流体导热系数 (W/m·K)，由设备层从物性引擎获取",
    })
    if missing:
        return missing

    all_warnings: list[str] = []

    # 1. 热负荷
    duty = heat_duty(
        m_hot=m_hot, Cp_hot=Cp_hot, T_hot_in=T_hot_in, T_hot_out=T_hot_out,
        m_cold=m_cold, Cp_cold=Cp_cold, T_cold_in=T_cold_in, T_cold_out=T_cold_out,
        latent_heat_hot=latent_heat_hot, latent_heat_cold=latent_heat_cold,
        phase_hot=phase_hot, phase_cold=phase_cold,
    )
    all_warnings.extend(duty.get("warnings", []))

    # 2. LMTD + F 因子
    lmtd_result = lmtd(
        T_hot_in, T_hot_out, T_cold_in, T_cold_out,
        flow_arrangement="counter",
        configuration=configuration,
        apply_correction=True,
    )
    all_warnings.extend(lmtd_result.get("warnings", []))

    if not lmtd_result["valid"]:
        return {"valid": False, "error": "LMTD 计算无效", "details": lmtd_result, "warnings": all_warnings}

    # 用于动态查表的热侧平均摄氏度
    T_hot_avg_C = ((T_hot_in + T_hot_out) / 2.0) - 273.15

    # 3. 传热系数 — 迭代自洽
    n_tubes_guess = 100
    h_tube = 3000.0
    h_shell = 5000.0
    u_tube = 1.0
    u_shell = 1.0
    
    # 动态管间距设置
    pitch = 1.25 * d_o

    for iteration in range(5):
        # 计算管程换热与流速
        h_tube_r = tube_side_htc(m_hot, rho_hot, mu_hot, k_hot, Cp_hot, d_i,
                                 n_tubes=n_tubes_guess, n_passes=tube_passes,
                                 tube_length=tube_length)
        h_tube = h_tube_r["h"] if h_tube_r["valid"] else 3000.0
        u_tube = h_tube_r.get("velocity", 1.0)

        # 估算壳径以计算壳程换热
        k_layout = 0.65 if tube_layout == "triangular" else 0.80
        D_shell_est = d_o * 1.25 * math.sqrt(n_tubes_guess / k_layout) if n_tubes_guess > 0 else 0.5
        
        # 修正：将折流板设计放入迭代自洽循环内，而不是使用外部硬编码的 tray_spacing = 0.6 计算 h_shell。
        # 否则循环中估算的 h_shell（基于 0.6m）与最终折流设计后的 h_shell（基于 0.15m 等真实折流距）差别极大，导致迭代失真。
        if force_baffle_spacing is not None and force_baffle_spacing > 0:
            B_temp = force_baffle_spacing
        else:
            baffle_temp = design_baffles(D_shell_est, tube_length)
            B_temp = baffle_temp["baffle_spacing_m"] if baffle_temp["valid"] else baffle_spacing
        
        # 计算壳程换热与流速（采用迭代对应的动态 B_temp）
        h_shell_r = shell_side_htc(m_cold, rho_cold, mu_cold, k_cold, Cp_cold,
                                    D_shell_est, B_temp, d_o, pitch=pitch, layout=tube_layout)
        h_shell = h_shell_r["h"] if h_shell_r["valid"] else 5000.0
        u_shell = h_shell_r.get("velocity", 1.0)

        # 动态更新污垢与整体 U 值
        u_result = overall_htc_estimate(
            h_tube=h_tube, h_shell=h_shell,
            d_o=d_o, d_i=d_i, k_wall=k_wall,
            fouling_tube=fouling_tube, fouling_shell=fouling_shell,
            z_tube=z_tube, z_shell=z_shell,
            T_heat_C=T_hot_avg_C, v_tube=u_tube, v_shell=u_shell
        )

        area = area_required(
            Q=duty["Q_avg"], U=u_result["U"], LMTD=lmtd_result["LMTD"],
            overspec_percent=overspec_percent,
            d_o=d_o, tube_length=tube_length, tube_passes=tube_passes, tube_layout=tube_layout,
            m_flow=m_hot, rho=rho_hot, d_i=d_i, V_min=0.5,
            force_tube_passes=force_tube_passes,
        )
        
        n_tubes_actual = area.get("N_tubes", 100)
        # 收敛判断：管数变化 < 5%
        if abs(n_tubes_actual - n_tubes_guess) / max(n_tubes_actual, 1) < 0.05:
            break
        n_tubes_guess = n_tubes_actual
        
    all_warnings.extend(area.get("warnings", []))

    # 使用流速约束调整后的管程数
    tube_passes = area.get("tube_passes", tube_passes)
    # 使用长径比约束调整后的管长
    tube_length = area.get("tube_length", tube_length)

    # 4b. 确定最终折流板设计
    D_shell_final = area.get("D_shell_standard", area.get("D_shell_est", 0.5))
    if force_baffle_spacing is not None and force_baffle_spacing > 0:
        B = force_baffle_spacing
        # 仍需从 design_baffles 获取 N_baffles 等其他参数
        baffle = design_baffles(D_shell_final, tube_length)
        baffle["baffle_spacing_m"] = B  # 覆盖间距为强制值
    else:
        baffle = design_baffles(D_shell_final, tube_length)
        B = baffle["baffle_spacing_m"]
    N_baffles = baffle["N_baffles"]

    # 5. 压降计算
    d_e_shell = h_shell_r.get("d_equivalent", 0.02)  
    dp_tube = pressure_drop("tube", m_hot, rho_hot, mu_hot, d_i, tube_length, tube_passes,
                            n_tubes=area.get("N_tubes", n_tubes_guess))
    dp_shell = pressure_drop("shell", m_cold, rho_cold, mu_cold, d=d_e_shell, L=tube_length,
                              D_shell=D_shell_final, B=B, N_baffles=N_baffles, d_o=d_o,
                              pitch=pitch, layout=tube_layout)

    # 6. 约束校验
    constraint_checks = _check_constraints(lmtd_result, u_result, area, dp_tube, dp_shell)

    return {
        "valid": True,
        "Q": duty["Q_avg"],
        "LMTD": lmtd_result["LMTD"],
        "LMTD_raw": lmtd_result["LMTD_raw"],
        "F": lmtd_result["F"],
        "F_quality": lmtd_result.get("F_quality"),
        "U": u_result["U"],
        "A_required": area["A_required"],
        "A_design": area["A_design"],
        "N_tubes": area["N_tubes"],
        "tube_passes": area.get("tube_passes", tube_passes),
        "D_shell": area.get("D_shell_standard", area.get("D_shell_est", 0)),
        "overspec_percent": overspec_percent,
        "length_diameter_ratio": area.get("length_diameter_ratio", 0),
        "delta_P_tube": dp_tube.get("delta_P", 0),
        "delta_P_shell": dp_shell.get("delta_P", 0),
        "h_tube": h_tube,
        "h_shell": h_shell,
        "fouling_fraction": u_result.get("fouling_fraction", 0),
        "tube_layout": tube_layout,
        "baffle_spacing_m": B,
        "N_baffles": N_baffles,
        "N_cross_passes": baffle["N_cross_passes"],
        "baffle_cut": baffle["baffle_cut"],
        "D_shell_est": area.get("D_shell_est", 0),
        "d_e_shell": d_e_shell,
        "constraint_checks": constraint_checks,
        "warnings": all_warnings,
        "duty": duty,
        "lmtd_result": lmtd_result,
        "u_result": u_result,
        "area_result": area,
        "dp_tube_result": dp_tube,
        "dp_shell_result": dp_shell,
        "h_tube_detail": h_tube_r,
        "h_shell_detail": h_shell_r,
    }


def _check_constraints(lmtd_r, u_r, area_r, dp_t, dp_s) -> dict:
    checks = {}
    checks["LMTD_positive"] = {"value": lmtd_r.get("LMTD", 0), "pass": lmtd_r.get("LMTD", 0) > 0, "criterion": "LMTD > 0"}
    checks["F_factor"] = {"value": lmtd_r.get("F"), "pass": (lmtd_r.get("F") or 0) >= 0.5, "criterion": "F >= 0.50"}
    checks["approach_temp"] = {"value": lmtd_r.get("approach_temp_min", 0), "pass": lmtd_r.get("approach_temp_min", 0) >= 5.0, "criterion": "ΔT_min >= 5 K"}
    checks["overspec"] = {"value": area_r.get("overspec_percent", 0), "pass": 10 <= area_r.get("overspec_percent", 0) <= 30, "criterion": "面积裕度 10-30%"}
    
    L_D = area_r.get("length_diameter_ratio", 0)
    checks["L_D_ratio"] = {"value": L_D, "pass": 3 <= L_D <= 15, "criterion": "L/D = 3-15"}
    checks["delta_P_tube"] = {"value": dp_t.get("delta_P", 0), "pass": dp_t.get("delta_P", 0) <= 100000, "criterion": "ΔP_tube <= 100 kPa"}
    checks["delta_P_shell"] = {"value": dp_s.get("delta_P", 0), "pass": dp_s.get("delta_P", 0) <= 70000, "criterion": "ΔP_shell <= 70 kPa"}

    n_pass = sum(1 for c in checks.values() if c["pass"])
    n_total = len(checks)
    return {"checks": checks, "n_pass": n_pass, "n_total": n_total, "pass_rate": n_pass / n_total if n_total > 0 else 0}


# ============================================================
# 工具层兼容别名
# ============================================================

def lmtd_with_correction(T_hot_in: float, T_hot_out: float, T_cold_in: float, T_cold_out: float, flow_arrangement: str = "counter", configuration: str = "1_shell_2_tube") -> dict:
    return lmtd(T_hot_in=T_hot_in, T_hot_out=T_hot_out, T_cold_in=T_cold_in, T_cold_out=T_cold_out, flow_arrangement=flow_arrangement, configuration=configuration, apply_correction=True)

def lmtd_correction_factor(R: float, P: float, configuration: str = "1_shell_2_tube") -> dict:
    F = _F_bowman_1_shell(R, P)
    return {"F": F, "R": R, "P": P, "configuration": configuration, "acceptable": F >= 0.75, "valid": F > 0}

def heat_duty_from_flow(m: float, Cp: float, T_in: float, T_out: float, latent_heat: float = 0.0, phase: str = "single") -> dict:
    Q_sensible = m * Cp * (T_out - T_in)
    Q_latent = m * latent_heat if phase in ("condensing", "evaporating") else 0.0
    Q = Q_sensible + Q_latent
    return {"Q": Q, "Q_sensible": Q_sensible, "Q_latent": Q_latent, "m": m, "T_in": T_in, "T_out": T_out, "valid": abs(Q) > 1e-6}

def heat_exchanger_area(Q: float, U: float, LMTD: float, overspec_percent: float = 15.0) -> dict:
    return area_required(Q=Q, U=U, LMTD=LMTD, overspec_percent=overspec_percent)

def overall_heat_transfer_coefficient(
    h_tube: float, h_shell: float, d_o: float = 0.025, d_i: float = 0.021, k_wall: float = 50.0,
    fouling_tube: str | tuple | list = "default", fouling_shell: str | tuple | list = "default",
    z_tube: list[float] | None = None, z_shell: list[float] | None = None,
    R_f_tube: float | None = None, R_f_shell: float | None = None,
) -> dict:
    return overall_htc_estimate(
        h_tube=h_tube, h_shell=h_shell, d_o=d_o, d_i=d_i, k_wall=k_wall,
        fouling_tube=fouling_tube, fouling_shell=fouling_shell, z_tube=z_tube, z_shell=z_shell,
        R_f_tube=R_f_tube, R_f_shell=R_f_shell
    )

def fouling_resistance(fluid_category: str = "water", fluid_type: str = "tap", z: list | None = None, T_heat_C: float = 100.0, v_flow: float = 1.0) -> dict:
    spec = (fluid_category, fluid_type) if not isinstance(fluid_category, list) else fluid_category
    Rf = calculate_fouling_resistance(spec, z, T_heat_C, v_flow)
    return {"R_f": Rf, "fluid_category": fluid_category, "fluid_type": fluid_type, "valid": True}

def effectiveness_ntu(C_hot: float, C_cold: float, UA: float, T_hot_in: float = None, T_cold_in: float = None, flow_arrangement: str = "counter") -> dict:
    return epsilon_ntu(C_hot=C_hot, C_cold=C_cold, UA=UA, T_hot_in=T_hot_in, T_cold_in=T_cold_in, flow_arrangement=flow_arrangement)

def estimate_U(category: str = "water/water") -> dict:
    u_range = typical_U_range(category)
    return {"U_min": u_range.get("U_min", 0), "U_max": u_range.get("U_max", 0), "U_typical": u_range.get("U_typical", 0), "category": category, "valid": True}

def heat_exchanger_design(
    T_hot_in: float, T_hot_out: float, T_cold_in: float, T_cold_out: float, m_hot: float, m_cold: float,
    Cp_hot: float = _REQUIRED, Cp_cold: float = _REQUIRED, rho_hot: float = _REQUIRED, rho_cold: float = _REQUIRED,
    mu_hot: float = _REQUIRED, mu_cold: float = _REQUIRED, k_hot: float = _REQUIRED, k_cold: float = _REQUIRED,
    latent_heat_hot: float = 0.0, latent_heat_cold: float = 0.0, phase_hot: str = "single", phase_cold: str = "single",
    configuration: str = "1_shell_2_tube", d_o: float = 0.025, d_i: float = 0.021, tube_length: float = 6.0,
    tube_passes: int = 2, tube_layout: str = "triangular", overspec_percent: float = 15.0, k_wall: float = 50.0,
    fouling_tube: str | tuple | list = "default", fouling_shell: str | tuple | list = "default",
    z_tube: list[float] | None = None, z_shell: list[float] | None = None, baffle_spacing: float = 0.6,
) -> dict:
    return design_heat_exchanger(
        T_hot_in=T_hot_in, T_hot_out=T_hot_out, T_cold_in=T_cold_in, T_cold_out=T_cold_out,
        m_hot=m_hot, m_cold=m_cold, Cp_hot=Cp_hot, Cp_cold=Cp_cold, rho_hot=rho_hot, rho_cold=rho_cold,
        mu_hot=mu_hot, mu_cold=mu_cold, k_hot=k_hot, k_cold=k_cold, latent_heat_hot=latent_heat_hot, latent_heat_cold=latent_heat_cold,
        phase_hot=phase_hot, phase_cold=phase_cold, configuration=configuration, d_o=d_o, d_i=d_i,
        tube_length=tube_length, tube_passes=tube_passes, tube_layout=tube_layout,
        overspec_percent=overspec_percent, k_wall=k_wall, fouling_tube=fouling_tube, fouling_shell=fouling_shell,
        z_tube=z_tube, z_shell=z_shell, baffle_spacing=baffle_spacing
    )