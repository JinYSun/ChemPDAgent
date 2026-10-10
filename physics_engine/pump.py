"""
泵选型工具包 — 离心泵/容积泵选型、NPSH计算、管路特性
包含：
- 泵功率计算（水力功率、轴功率、电机功率）
- NPSH 计算（NPSHa, NPSHr）
- 泵类型选择（离心泵/容积泵）
- 管路特性曲线
- 泵性能曲线匹配
- 效率估算
- 气蚀余量校核
- 标准泵系列选型
"""
from __future__ import annotations
import numpy as np
import math
from typing import Optional, List, Dict, Tuple

# Import thermo helper
from .common import _REQUIRED, check_required_params


# ============================================================
# 泵类型与标准系列
# ============================================================

# 泵类型
PUMP_CENTRIFUGAL = "centrifugal"             # IS/IR 清水离心泵
PUMP_CHEMICAL = "chemical"                    # IH 化工离心泵 (GB/T 5656)
PUMP_DOUBLE_SUCTION = "double_suction"         # S/SH 双吸离心泵
PUMP_MULTISTAGE = "multistage"                 # D/DG 多级离心泵 (GB/T 5657)
PUMP_MAGNETIC = "magnetic"                    # CQ 磁力驱动泵
PUMP_POSITIVE_DISPLACEMENT = "positive_displacement"
PUMP_GEAR = "gear"
PUMP_SCREW = "screw"
PUMP_DIAPHRAGM = "diaphragm"
PUMP_MIXED_AXIAL = "mixed_flow_axial"
PUMP_VORTEX = "vortex"
PUMP_METERING = "metering"
PUMP_SUBMERSIBLE = "submersible"  # deprecated, 用 PUMP_SEWAGE 代替
PUMP_SEWAGE = "sewage"                        # QW/WQ 潜污泵
PUMP_SELF_PRIMING = "self_priming"              # ZX 自吸离心泵
PUMP_FLUOROPLASTIC = "fluoroplastic"            # IHF 氟塑料衬里离心泵

# ── 离心泵标准系列 ── GB/T 5660/5661, IS 系列清水泵（流量 m³/h, 扬程 m）
_CENTRIFUGAL_PUMP_SERIES = [
    # 2900 rpm 小流量低扬程
    {"model": "IS50-32-125",  "Q_min": 3,   "Q_max": 12,  "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.55},
    {"model": "IS50-32-160",  "Q_min": 3,   "Q_max": 12.5,"H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.49},
    {"model": "IS50-32-200",  "Q_min": 3,   "Q_max": 12.5,"H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.42},
    {"model": "IS50-32-250",  "Q_min": 3,   "Q_max": 12.5,"H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.33},
    {"model": "IS65-50-125",  "Q_min": 7.5, "Q_max": 30,  "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.63},
    {"model": "IS65-50-160",  "Q_min": 7.5, "Q_max": 30,  "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.59},
    {"model": "IS65-40-200",  "Q_min": 7.5, "Q_max": 30,  "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.53},
    {"model": "IS65-40-250",  "Q_min": 7.5, "Q_max": 30,  "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.43},
    {"model": "IS65-40-315",  "Q_min": 7.5, "Q_max": 30,  "H_min": 32, "H_max": 125, "n": 2900, "eff": 0.35},
    {"model": "IS80-65-125",  "Q_min": 15,  "Q_max": 60,  "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.71},
    {"model": "IS80-65-160",  "Q_min": 15,  "Q_max": 60,  "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.69},
    {"model": "IS80-50-200",  "Q_min": 15,  "Q_max": 60,  "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.65},
    {"model": "IS80-50-250",  "Q_min": 15,  "Q_max": 60,  "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.57},
    {"model": "IS80-50-315",  "Q_min": 15,  "Q_max": 60,  "H_min": 32, "H_max": 125, "n": 2900, "eff": 0.49},
    {"model": "IS100-80-125", "Q_min": 30,  "Q_max": 120, "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.75},
    {"model": "IS100-80-160", "Q_min": 30,  "Q_max": 120, "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.74},
    {"model": "IS100-65-200", "Q_min": 30,  "Q_max": 120, "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.72},
    {"model": "IS100-65-250", "Q_min": 30,  "Q_max": 120, "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.66},
    {"model": "IS100-65-315", "Q_min": 30,  "Q_max": 120, "H_min": 32, "H_max": 125, "n": 2900, "eff": 0.60},
    {"model": "IS125-100-200","Q_min": 60,  "Q_max": 240, "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.78},
    {"model": "IS125-100-250","Q_min": 60,  "Q_max": 240, "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.75},
    {"model": "IS125-100-315","Q_min": 60,  "Q_max": 240, "H_min": 32, "H_max": 125, "n": 2900, "eff": 0.71},
    # 1480 rpm 大流量中高扬程
    {"model": "IS150-125-250","Q_min": 120, "Q_max": 360, "H_min": 20, "H_max": 80,  "n": 1480, "eff": 0.78},
    {"model": "IS150-125-315","Q_min": 120, "Q_max": 360, "H_min": 32, "H_max": 125, "n": 1480, "eff": 0.75},
    {"model": "IS150-125-400","Q_min": 120, "Q_max": 360, "H_min": 50, "H_max": 125, "n": 1480, "eff": 0.70},
    {"model": "IS200-150-250","Q_min": 240, "Q_max": 600, "H_min": 20, "H_max": 80,  "n": 1480, "eff": 0.82},
    {"model": "IS200-150-315","Q_min": 240, "Q_max": 600, "H_min": 32, "H_max": 125, "n": 1480, "eff": 0.80},
    {"model": "IS200-150-400","Q_min": 240, "Q_max": 600, "H_min": 50, "H_max": 125, "n": 1480, "eff": 0.78},
]

# ── 计量泵标准系列 ── GB/T 7782, J/JZ/JD/JY 系列
# 流量 L/h → m³/h（Q_min/Q_max 为 m³/h），压力 MPa → 扬程 m（1 MPa ≈ 102 m 水柱）
_METERING_PUMP_SERIES = [
    # 柱塞式计量泵 J 系列 — 高压精密计量
    {"model": "J-X 0.5/50",  "Q_min": 0.0005, "Q_max": 0.0005, "H_min": 5100, "H_max": 5100, "n": 100, "eff": 0.85, "type": "plunger",    "P_max_MPa": 50,  "accuracy": 0.01},
    {"model": "J-Z 2/10",    "Q_min": 0.002,  "Q_max": 0.002,  "H_min": 1020, "H_max": 1020, "n": 100, "eff": 0.85, "type": "plunger",    "P_max_MPa": 10,  "accuracy": 0.01},
    {"model": "J-Z 10/5",    "Q_min": 0.010,  "Q_max": 0.010,  "H_min": 510,  "H_max": 510,  "n": 100, "eff": 0.85, "type": "plunger",    "P_max_MPa": 5,   "accuracy": 0.01},
    {"model": "J-D 50/2.5",  "Q_min": 0.050,  "Q_max": 0.050,  "H_min": 255,  "H_max": 255,  "n": 80,  "eff": 0.85, "type": "plunger",    "P_max_MPa": 2.5, "accuracy": 0.01},
    {"model": "J-D 200/1.0", "Q_min": 0.200,  "Q_max": 0.200,  "H_min": 102,  "H_max": 102,  "n": 80,  "eff": 0.85, "type": "plunger",    "P_max_MPa": 1.0, "accuracy": 0.01},
    {"model": "J-D 500/0.5", "Q_min": 0.500,  "Q_max": 0.500,  "H_min": 51,   "H_max": 51,   "n": 80,  "eff": 0.85, "type": "plunger",    "P_max_MPa": 0.5, "accuracy": 0.01},
    # 液压隔膜计量泵 JY 系列 — 耐腐蚀、无泄漏
    {"model": "JY 10/0.5",   "Q_min": 0.010,  "Q_max": 0.010,  "H_min": 51,   "H_max": 51,   "n": 100, "eff": 0.80, "type": "diaphragm",  "P_max_MPa": 0.5, "accuracy": 0.01},
    {"model": "JY 50/0.3",   "Q_min": 0.050,  "Q_max": 0.050,  "H_min": 30.6, "H_max": 30.6, "n": 80,  "eff": 0.80, "type": "diaphragm",  "P_max_MPa": 0.3, "accuracy": 0.01},
    {"model": "JY 200/0.3",  "Q_min": 0.200,  "Q_max": 0.200,  "H_min": 30.6, "H_max": 30.6, "n": 80,  "eff": 0.80, "type": "diaphragm",  "P_max_MPa": 0.3, "accuracy": 0.01},
    {"model": "JY 500/0.2",  "Q_min": 0.500,  "Q_max": 0.500,  "H_min": 20.4, "H_max": 20.4, "n": 80,  "eff": 0.80, "type": "diaphragm",  "P_max_MPa": 0.2, "accuracy": 0.01},
    # 机械隔膜计量泵 JD 系列 — 中低压、结构简单
    {"model": "JD 50/0.5",   "Q_min": 0.050,  "Q_max": 0.050,  "H_min": 51,   "H_max": 51,   "n": 100, "eff": 0.80, "type": "mech_diaph", "P_max_MPa": 0.5, "accuracy": 0.02},
    {"model": "JD 200/0.3",  "Q_min": 0.200,  "Q_max": 0.200,  "H_min": 30.6, "H_max": 30.6, "n": 80,  "eff": 0.80, "type": "mech_diaph", "P_max_MPa": 0.3, "accuracy": 0.02},
    {"model": "JD 500/0.2",  "Q_min": 0.500,  "Q_max": 0.500,  "H_min": 20.4, "H_max": 20.4, "n": 80,  "eff": 0.80, "type": "mech_diaph", "P_max_MPa": 0.2, "accuracy": 0.02},
    {"model": "JD 1000/0.1", "Q_min": 1.000,  "Q_max": 1.000,  "H_min": 10.2, "H_max": 10.2, "n": 60,  "eff": 0.80, "type": "mech_diaph", "P_max_MPa": 0.1, "accuracy": 0.02},
]

# ── 漩涡泵标准系列 ── GB/T 29531, W/WB 系列
# 小流量高扬程，自吸性好，适合含气介质
_VORTEX_PUMP_SERIES = [
    # W 型开式漩涡泵
    {"model": "20W-20",   "Q_min": 0.2,  "Q_max": 0.5,  "H_min": 15, "H_max": 25, "n": 2900, "eff": 0.25},
    {"model": "25W-25",   "Q_min": 0.4,  "Q_max": 1.0,  "H_min": 20, "H_max": 30, "n": 2900, "eff": 0.28},
    {"model": "32W-30",   "Q_min": 0.8,  "Q_max": 2.0,  "H_min": 25, "H_max": 40, "n": 2900, "eff": 0.30},
    {"model": "32W-75",   "Q_min": 0.8,  "Q_max": 2.0,  "H_min": 60, "H_max": 90, "n": 2900, "eff": 0.22},
    {"model": "40W-40",   "Q_min": 1.5,  "Q_max": 3.5,  "H_min": 30, "H_max": 50, "n": 2900, "eff": 0.32},
    {"model": "40W-90",   "Q_min": 1.5,  "Q_max": 3.5,  "H_min": 70, "H_max": 110,"n": 2900, "eff": 0.22},
    {"model": "50W-45",   "Q_min": 3.0,  "Q_max": 7.0,  "H_min": 35, "H_max": 55, "n": 2900, "eff": 0.35},
    # WB 型闭式漩涡泵（自吸性能更好）
    {"model": "1.5WB-120", "Q_min": 0.1, "Q_max": 0.3,  "H_min": 100,"H_max": 140,"n": 2900, "eff": 0.18},
    {"model": "25WB-70",  "Q_min": 0.4,  "Q_max": 1.0,  "H_min": 55, "H_max": 85, "n": 2900, "eff": 0.22},
    {"model": "32WB-30",  "Q_min": 0.8,  "Q_max": 2.0,  "H_min": 22, "H_max": 38, "n": 2900, "eff": 0.30},
    {"model": "40WB-55",  "Q_min": 1.5,  "Q_max": 3.5,  "H_min": 40, "H_max": 65, "n": 2900, "eff": 0.28},
    {"model": "50WB-45",  "Q_min": 3.0,  "Q_max": 7.0,  "H_min": 35, "H_max": 55, "n": 2900, "eff": 0.32},
]

# ── 齿轮泵标准系列 ── GB/T 25631, KCB/2CY/CB 系列
# 高黏度、含气介质、自吸性好
_GEAR_PUMP_SERIES = [
    # KCB 型渐开线齿轮泵 — 低压通用
    {"model": "KCB-18.3",   "Q_min": 0.5,  "Q_max": 1.2,  "H_min": 100, "H_max": 148, "n": 1400, "eff": 0.45, "P_MPa": 1.45},
    {"model": "KCB-33.3",   "Q_min": 1.2,  "Q_max": 2.5,  "H_min": 100, "H_max": 148, "n": 1420, "eff": 0.50, "P_MPa": 1.45},
    {"model": "KCB-55",     "Q_min": 2.0,  "Q_max": 4.0,  "H_min": 25,  "H_max": 40,  "n": 1400, "eff": 0.55, "P_MPa": 0.33},
    {"model": "KCB-83.3",   "Q_min": 3.5,  "Q_max": 6.5,  "H_min": 25,  "H_max": 40,  "n": 1420, "eff": 0.58, "P_MPa": 0.33},
    {"model": "KCB-200",    "Q_min": 8,    "Q_max": 16,   "H_min": 25,  "H_max": 40,  "n": 1450, "eff": 0.62, "P_MPa": 0.33},
    {"model": "KCB-300",    "Q_min": 12,   "Q_max": 22,   "H_min": 25,  "H_max": 45,  "n": 970,  "eff": 0.65, "P_MPa": 0.36},
    {"model": "KCB-483.3",  "Q_min": 20,   "Q_max": 36,   "H_min": 25,  "H_max": 45,  "n": 970,  "eff": 0.65, "P_MPa": 0.36},
    {"model": "KCB-633",    "Q_min": 28,   "Q_max": 48,   "H_min": 20,  "H_max": 35,  "n": 970,  "eff": 0.65, "P_MPa": 0.28},
    {"model": "KCB-960",    "Q_min": 40,   "Q_max": 70,   "H_min": 20,  "H_max": 35,  "n": 970,  "eff": 0.65, "P_MPa": 0.28},
    # 2CY 型内啮合齿轮泵 — 高压
    {"model": "2CY-1.08/2.5","Q_min": 0.5, "Q_max": 1.2,  "H_min": 200, "H_max": 255, "n": 1420, "eff": 0.40, "P_MPa": 2.5},
    {"model": "2CY-2.1/2.5","Q_min": 1.2,  "Q_max": 2.5,  "H_min": 200, "H_max": 255, "n": 1420, "eff": 0.45, "P_MPa": 2.5},
    {"model": "2CY-3/2.5",  "Q_min": 2.0,  "Q_max": 4.0,  "H_min": 200, "H_max": 255, "n": 1440, "eff": 0.50, "P_MPa": 2.5},
    {"model": "2CY-7.5/2.5","Q_min": 5.0,  "Q_max": 10,   "H_min": 200, "H_max": 255, "n": 960,  "eff": 0.55, "P_MPa": 2.5},
    {"model": "2CY-12/2.5", "Q_min": 8,    "Q_max": 16,   "H_min": 200, "H_max": 255, "n": 970,  "eff": 0.55, "P_MPa": 2.5},
    {"model": "2CY-18/2.5", "Q_min": 12,   "Q_max": 22,   "H_min": 200, "H_max": 255, "n": 970,  "eff": 0.55, "P_MPa": 2.5},
]

# ── 单螺杆泵标准系列 ── GB/T 25695, G/GF/GN 系列
# 高黏度、含固体颗粒、自吸性好
_SCREW_PUMP_SERIES = [
    # G 型单螺杆泵通用系列
    {"model": "G20-1",   "Q_min": 0.2,  "Q_max": 0.8,  "H_min": 30,  "H_max": 60,  "n": 960,  "eff": 0.50, "P_MPa": 0.6, "solids_pct": 10},
    {"model": "G25-1",   "Q_min": 0.5,  "Q_max": 1.5,  "H_min": 30,  "H_max": 60,  "n": 960,  "eff": 0.55, "P_MPa": 0.6, "solids_pct": 15},
    {"model": "G30-1",   "Q_min": 1.0,  "Q_max": 3.0,  "H_min": 30,  "H_max": 60,  "n": 960,  "eff": 0.58, "P_MPa": 0.6, "solids_pct": 20},
    {"model": "G35-1",   "Q_min": 2.0,  "Q_max": 5.0,  "H_min": 30,  "H_max": 60,  "n": 960,  "eff": 0.60, "P_MPa": 0.6, "solids_pct": 20},
    {"model": "G40-1",   "Q_min": 3.0,  "Q_max": 8.0,  "H_min": 30,  "H_max": 60,  "n": 720,  "eff": 0.62, "P_MPa": 0.6, "solids_pct": 25},
    {"model": "G50-1",   "Q_min": 5.0,  "Q_max": 15,   "H_min": 30,  "H_max": 60,  "n": 720,  "eff": 0.62, "P_MPa": 0.6, "solids_pct": 25},
    {"model": "G70-1",   "Q_min": 10,   "Q_max": 30,   "H_min": 30,  "H_max": 60,  "n": 720,  "eff": 0.65, "P_MPa": 0.6, "solids_pct": 30},
    {"model": "G85-1",   "Q_min": 20,   "Q_max": 50,   "H_min": 30,  "H_max": 60,  "n": 480,  "eff": 0.65, "P_MPa": 0.6, "solids_pct": 30},
    {"model": "G105-1",  "Q_min": 30,   "Q_max": 80,   "H_min": 30,  "H_max": 60,  "n": 380,  "eff": 0.65, "P_MPa": 0.6, "solids_pct": 30},
    # 高压二级系列
    {"model": "G25-2",   "Q_min": 0.3,  "Q_max": 1.0,  "H_min": 60,  "H_max": 120, "n": 960,  "eff": 0.45, "P_MPa": 1.2, "solids_pct": 10},
    {"model": "G35-2",   "Q_min": 1.5,  "Q_max": 4.0,  "H_min": 60,  "H_max": 120, "n": 960,  "eff": 0.50, "P_MPa": 1.2, "solids_pct": 15},
    {"model": "G50-2",   "Q_min": 3.0,  "Q_max": 10,   "H_min": 60,  "H_max": 120, "n": 720,  "eff": 0.55, "P_MPa": 1.2, "solids_pct": 20},
    {"model": "G70-2",   "Q_min": 8,    "Q_max": 22,   "H_min": 60,  "H_max": 120, "n": 720,  "eff": 0.55, "P_MPa": 1.2, "solids_pct": 25},
]

# ── 混流泵 / 轴流泵标准系列 ── GB/T 13008, HW/HB/ZLB 系列
# 大流量低扬程
_MIXED_AXIAL_PUMP_SERIES = [
    # HW 型卧式混流泵（中小型）
    {"model": "HW-150",  "Q_min": 50,   "Q_max": 120,  "H_min": 5,  "H_max": 12, "n": 1450, "eff": 0.70, "type": "mixed"},
    {"model": "HW-200",  "Q_min": 100,  "Q_max": 250,  "H_min": 5,  "H_max": 14, "n": 1450, "eff": 0.75, "type": "mixed"},
    {"model": "HW-250",  "Q_min": 180,  "Q_max": 400,  "H_min": 6,  "H_max": 16, "n": 1450, "eff": 0.78, "type": "mixed"},
    {"model": "HW-300",  "Q_min": 300,  "Q_max": 650,  "H_min": 6,  "H_max": 18, "n": 980,  "eff": 0.80, "type": "mixed"},
    {"model": "HW-350",  "Q_min": 450,  "Q_max": 1000, "H_min": 6,  "H_max": 18, "n": 980,  "eff": 0.80, "type": "mixed"},
    {"model": "HW-400",  "Q_min": 600,  "Q_max": 1500, "H_min": 5,  "H_max": 15, "n": 730,  "eff": 0.82, "type": "mixed"},
    # HB 型大型立式混流泵
    {"model": "HB-500",  "Q_min": 1000, "Q_max": 2500, "H_min": 5,  "H_max": 20, "n": 590,  "eff": 0.83, "type": "mixed"},
    {"model": "HB-600",  "Q_min": 1500, "Q_max": 4000, "H_min": 5,  "H_max": 22, "n": 490,  "eff": 0.84, "type": "mixed"},
    {"model": "HB-800",  "Q_min": 3000, "Q_max": 8000, "H_min": 5,  "H_max": 25, "n": 370,  "eff": 0.85, "type": "mixed"},
    {"model": "HB-1000", "Q_min": 5000, "Q_max": 15000,"H_min": 5,  "H_max": 25, "n": 300,  "eff": 0.85, "type": "mixed"},
    # ZLB 型立式轴流泵
    {"model": "ZLB-350", "Q_min": 300,  "Q_max": 800,  "H_min": 2,  "H_max": 6,  "n": 1450, "eff": 0.78, "type": "axial"},
    {"model": "ZLB-500", "Q_min": 700,  "Q_max": 2000, "H_min": 2,  "H_max": 7,  "n": 980,  "eff": 0.80, "type": "axial"},
    {"model": "ZLB-600", "Q_min": 1200, "Q_max": 3500, "H_min": 2,  "H_max": 8,  "n": 730,  "eff": 0.82, "type": "axial"},
    {"model": "ZLB-800", "Q_min": 2500, "Q_max": 7000, "H_min": 2,  "H_max": 8,  "n": 590,  "eff": 0.83, "type": "axial"},
    {"model": "ZLB-1000","Q_min": 5000, "Q_max": 15000,"H_min": 1,  "H_max": 6,  "n": 490,  "eff": 0.82, "type": "axial"},
]

# ── 隔膜泵标准系列 ── QBY/DBY 系列气动/电动隔膜泵
# 耐腐蚀、含固体颗粒、含气介质
_DIAPHRAGM_PUMP_SERIES = [
    # QBY 型气动隔膜泵 — 无电机，适合防爆场合
    {"model": "QBY-10",  "Q_min": 0.2,  "Q_max": 0.8,  "H_min": 30,  "H_max": 50,  "n": 0, "eff": 0.30, "drive": "pneumatic", "solids_mm": 1},
    {"model": "QBY-15",  "Q_min": 0.5,  "Q_max": 1.5,  "H_min": 30,  "H_max": 50,  "n": 0, "eff": 0.35, "drive": "pneumatic", "solids_mm": 1},
    {"model": "QBY-25",  "Q_min": 1.0,  "Q_max": 3.5,  "H_min": 30,  "H_max": 50,  "n": 0, "eff": 0.40, "drive": "pneumatic", "solids_mm": 2.5},
    {"model": "QBY-40",  "Q_min": 3.0,  "Q_max": 8.0,  "H_min": 30,  "H_max": 50,  "n": 0, "eff": 0.45, "drive": "pneumatic", "solids_mm": 2.5},
    {"model": "QBY-50",  "Q_min": 5.0,  "Q_max": 15,   "H_min": 30,  "H_max": 50,  "n": 0, "eff": 0.48, "drive": "pneumatic", "solids_mm": 4},
    {"model": "QBY-65",  "Q_min": 10,   "Q_max": 25,   "H_min": 30,  "H_max": 50,  "n": 0, "eff": 0.50, "drive": "pneumatic", "solids_mm": 4},
    {"model": "QBY-80",  "Q_min": 15,   "Q_max": 40,   "H_min": 30,  "H_max": 50,  "n": 0, "eff": 0.50, "drive": "pneumatic", "solids_mm": 5},
    {"model": "QBY-100", "Q_min": 25,   "Q_max": 60,   "H_min": 30,  "H_max": 50,  "n": 0, "eff": 0.50, "drive": "pneumatic", "solids_mm": 6},
    # DBY 型电动隔膜泵
    {"model": "DBY-10",  "Q_min": 0.2,  "Q_max": 0.8,  "H_min": 30,  "H_max": 50,  "n": 1400, "eff": 0.35, "drive": "electric", "solids_mm": 1},
    {"model": "DBY-15",  "Q_min": 0.5,  "Q_max": 1.5,  "H_min": 30,  "H_max": 50,  "n": 1400, "eff": 0.40, "drive": "electric", "solids_mm": 1},
    {"model": "DBY-25",  "Q_min": 1.0,  "Q_max": 3.5,  "H_min": 30,  "H_max": 50,  "n": 1400, "eff": 0.45, "drive": "electric", "solids_mm": 2.5},
    {"model": "DBY-40",  "Q_min": 3.0,  "Q_max": 8.0,  "H_min": 30,  "H_max": 50,  "n": 1400, "eff": 0.48, "drive": "electric", "solids_mm": 2.5},
    {"model": "DBY-50",  "Q_min": 5.0,  "Q_max": 15,   "H_min": 30,  "H_max": 50,  "n": 1400, "eff": 0.50, "drive": "electric", "solids_mm": 4},
    {"model": "DBY-65",  "Q_min": 10,   "Q_max": 25,   "H_min": 30,  "H_max": 50,  "n": 960,  "eff": 0.50, "drive": "electric", "solids_mm": 4},
    {"model": "DBY-80",  "Q_min": 15,   "Q_max": 40,   "H_min": 30,  "H_max": 50,  "n": 960,  "eff": 0.50, "drive": "electric", "solids_mm": 5},
]

# ── 化工离心泵标准系列 ── GB/T 5656, IH 系列不锈钢耐腐蚀泵
# 流量 3.4–460 m³/h, 扬程 3.6–132 m, 2900/1450 rpm
_IH_CHEMICAL_PUMP_SERIES = [
    # 2900 rpm
    {"model": "IH50-32-125",  "Q_min": 6.3,  "Q_max": 12.5, "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.51, "material": "304/316L"},
    {"model": "IH50-32-160",  "Q_min": 6.3,  "Q_max": 12.5, "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.46, "material": "304/316L"},
    {"model": "IH50-32-200",  "Q_min": 6.3,  "Q_max": 12.5, "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.39, "material": "304/316L"},
    {"model": "IH50-32-250",  "Q_min": 6.3,  "Q_max": 12.5, "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.31, "material": "304/316L"},
    {"model": "IH65-50-125",  "Q_min": 12.5, "Q_max": 25,  "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.60, "material": "304/316L"},
    {"model": "IH65-50-160",  "Q_min": 12.5, "Q_max": 25,  "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.57, "material": "304/316L"},
    {"model": "IH65-40-200",  "Q_min": 12.5, "Q_max": 25,  "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.50, "material": "304/316L"},
    {"model": "IH65-40-250",  "Q_min": 12.5, "Q_max": 25,  "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.41, "material": "304/316L"},
    {"model": "IH65-40-315",  "Q_min": 12.5, "Q_max": 25,  "H_min": 32, "H_max": 125, "n": 2900, "eff": 0.33, "material": "304/316L"},
    {"model": "IH80-65-125",  "Q_min": 25,   "Q_max": 50,  "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.68, "material": "304/316L"},
    {"model": "IH80-65-160",  "Q_min": 25,   "Q_max": 50,  "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.67, "material": "304/316L"},
    {"model": "IH80-50-200",  "Q_min": 25,   "Q_max": 50,  "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.63, "material": "304/316L"},
    {"model": "IH80-50-250",  "Q_min": 25,   "Q_max": 50,  "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.55, "material": "304/316L"},
    {"model": "IH80-50-315",  "Q_min": 25,   "Q_max": 50,  "H_min": 32, "H_max": 125, "n": 2900, "eff": 0.47, "material": "304/316L"},
    {"model": "IH100-80-125", "Q_min": 50,   "Q_max": 100, "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.73, "material": "304/316L"},
    {"model": "IH100-80-160", "Q_min": 50,   "Q_max": 100, "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.72, "material": "304/316L"},
    {"model": "IH100-65-200", "Q_min": 50,   "Q_max": 100, "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.70, "material": "304/316L"},
    {"model": "IH100-65-250", "Q_min": 50,   "Q_max": 100, "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.64, "material": "304/316L"},
    {"model": "IH100-65-315", "Q_min": 50,   "Q_max": 100, "H_min": 32, "H_max": 125, "n": 2900, "eff": 0.58, "material": "304/316L"},
    {"model": "IH125-100-200","Q_min": 100,  "Q_max": 200, "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.75, "material": "304/316L"},
    {"model": "IH125-100-250","Q_min": 100,  "Q_max": 200, "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.73, "material": "304/316L"},
    {"model": "IH125-100-315","Q_min": 100,  "Q_max": 200, "H_min": 32, "H_max": 125, "n": 2900, "eff": 0.68, "material": "304/316L"},
    {"model": "IH125-100-400","Q_min": 100,  "Q_max": 200, "H_min": 50, "H_max": 125, "n": 1450, "eff": 0.62, "material": "304/316L"},
    # 1450 rpm
    {"model": "IH150-125-250","Q_min": 200,  "Q_max": 400, "H_min": 20, "H_max": 80,  "n": 1450, "eff": 0.76, "material": "304/316L"},
    {"model": "IH150-125-315","Q_min": 200,  "Q_max": 400, "H_min": 32, "H_max": 125, "n": 1450, "eff": 0.73, "material": "304/316L"},
    {"model": "IH150-125-400","Q_min": 200,  "Q_max": 400, "H_min": 50, "H_max": 125, "n": 1450, "eff": 0.68, "material": "304/316L"},
    {"model": "IH200-150-315","Q_min": 400,  "Q_max": 600, "H_min": 32, "H_max": 125, "n": 1450, "eff": 0.78, "material": "304/316L"},
    {"model": "IH200-150-400","Q_min": 400,  "Q_max": 600, "H_min": 50, "H_max": 125, "n": 1450, "eff": 0.75, "material": "304/316L"},
]

# ── 双吸离心泵标准系列 ── GB/T 5657, S/SH 系列中开式双吸泵
# 流量 72–3170 m³/h, 扬程 10–125 m
_S_DOUBLE_SUCTION_SERIES = [
    {"model": "100S90",   "Q_min": 72,  "Q_max": 80,   "H_min": 75, "H_max": 95,  "n": 2950, "eff": 0.65},
    {"model": "150S100",  "Q_min": 140, "Q_max": 180,  "H_min": 85, "H_max": 110, "n": 2950, "eff": 0.73},
    {"model": "150S78",   "Q_min": 140, "Q_max": 180,  "H_min": 65, "H_max": 85,  "n": 2950, "eff": 0.74},
    {"model": "150S50",   "Q_min": 140, "Q_max": 180,  "H_min": 40, "H_max": 55,  "n": 2950, "eff": 0.80},
    {"model": "200S95",   "Q_min": 250, "Q_max": 310,  "H_min": 80, "H_max": 100, "n": 2950, "eff": 0.79},
    {"model": "200S63",   "Q_min": 250, "Q_max": 310,  "H_min": 50, "H_max": 70,  "n": 2950, "eff": 0.83},
    {"model": "200S42",   "Q_min": 250, "Q_max": 310,  "H_min": 35, "H_max": 48,  "n": 2950, "eff": 0.84},
    {"model": "250S65",   "Q_min": 420, "Q_max": 530,  "H_min": 55, "H_max": 70,  "n": 1450, "eff": 0.79},
    {"model": "250S39",   "Q_min": 420, "Q_max": 530,  "H_min": 33, "H_max": 43,  "n": 1450, "eff": 0.84},
    {"model": "250S24",   "Q_min": 420, "Q_max": 530,  "H_min": 20, "H_max": 27,  "n": 1450, "eff": 0.86},
    {"model": "250S14",   "Q_min": 420, "Q_max": 530,  "H_min": 11, "H_max": 16,  "n": 1450, "eff": 0.86},
    {"model": "300S110",  "Q_min": 700, "Q_max": 870,  "H_min": 95, "H_max": 120, "n": 1450, "eff": 0.80},
    {"model": "300S90",   "Q_min": 700, "Q_max": 870,  "H_min": 78, "H_max": 98,  "n": 1450, "eff": 0.80},
    {"model": "300S58",   "Q_min": 700, "Q_max": 870,  "H_min": 49, "H_max": 64,  "n": 1450, "eff": 0.84},
    {"model": "300S32",   "Q_min": 700, "Q_max": 870,  "H_min": 26, "H_max": 35,  "n": 1450, "eff": 0.87},
    {"model": "300S19",   "Q_min": 700, "Q_max": 870,  "H_min": 16, "H_max": 22,  "n": 1450, "eff": 0.87},
    {"model": "300S12",   "Q_min": 700, "Q_max": 870,  "H_min": 10, "H_max": 14,  "n": 1450, "eff": 0.85},
    {"model": "350S125",  "Q_min": 1100,"Q_max": 1400, "H_min": 105,"H_max": 135, "n": 1450, "eff": 0.81},
    {"model": "350S75",   "Q_min": 1100,"Q_max": 1400, "H_min": 65, "H_max": 82,  "n": 1450, "eff": 0.85},
    {"model": "350S44",   "Q_min": 1100,"Q_max": 1400, "H_min": 36, "H_max": 48,  "n": 1450, "eff": 0.88},
    {"model": "350S26",   "Q_min": 1100,"Q_max": 1400, "H_min": 21, "H_max": 28,  "n": 1450, "eff": 0.88},
    {"model": "350S16",   "Q_min": 1100,"Q_max": 1400, "H_min": 13, "H_max": 18,  "n": 1450, "eff": 0.85},
    {"model": "500S98",   "Q_min": 1800,"Q_max": 2200, "H_min": 83, "H_max": 105, "n": 970,  "eff": 0.80},
    {"model": "500S59",   "Q_min": 1800,"Q_max": 2200, "H_min": 49, "H_max": 65,  "n": 970,  "eff": 0.83},
    {"model": "500S35",   "Q_min": 1800,"Q_max": 2200, "H_min": 27, "H_max": 38,  "n": 970,  "eff": 0.88},
    {"model": "500S22",   "Q_min": 1800,"Q_max": 2200, "H_min": 17, "H_max": 24,  "n": 970,  "eff": 0.84},
    {"model": "500S13",   "Q_min": 1800,"Q_max": 2200, "H_min": 11, "H_max": 15,  "n": 970,  "eff": 0.83},
    {"model": "600S75",   "Q_min": 2800,"Q_max": 3400, "H_min": 65, "H_max": 82,  "n": 970,  "eff": 0.89},
    {"model": "600S47",   "Q_min": 2800,"Q_max": 3400, "H_min": 40, "H_max": 52,  "n": 970,  "eff": 0.88},
    {"model": "600S32",   "Q_min": 2800,"Q_max": 3400, "H_min": 27, "H_max": 35,  "n": 970,  "eff": 0.89},
    {"model": "600S22",   "Q_min": 2800,"Q_max": 3400, "H_min": 18, "H_max": 25,  "n": 970,  "eff": 0.85},
]

# ── 多级离心泵标准系列 ── GB/T 5657, D/DG 系列卧式多级泵
# 每级扬程 25–50 m, 级数 3–12
_D_MULTISTAGE_SERIES = [
    # 格式: 型号, Q_min/Q_max (单级流量), H_min/H_max (总扬程, 3-12级), 转速, 效率
    {"model": "D6-25",   "Q_min": 3.75, "Q_max": 7.5,  "H_min": 76,  "H_max": 300,  "n": 2950, "eff": 0.45, "stages": "3-12"},
    {"model": "D12-25",  "Q_min": 7.5,  "Q_max": 15,   "H_min": 76,  "H_max": 300,  "n": 2950, "eff": 0.54, "stages": "3-12"},
    {"model": "D12-50",  "Q_min": 7.5,  "Q_max": 15,   "H_min": 150, "H_max": 600,  "n": 2950, "eff": 0.44, "stages": "3-12"},
    {"model": "D25-30",  "Q_min": 15,   "Q_max": 30,   "H_min": 90,  "H_max": 360,  "n": 2950, "eff": 0.60, "stages": "3-12"},
    {"model": "D25-50",  "Q_min": 15,   "Q_max": 30,   "H_min": 150, "H_max": 600,  "n": 2950, "eff": 0.55, "stages": "3-12"},
    {"model": "D46-30",  "Q_min": 30,   "Q_max": 55,   "H_min": 90,  "H_max": 360,  "n": 2950, "eff": 0.68, "stages": "3-12"},
    {"model": "D46-50",  "Q_min": 30,   "Q_max": 55,   "H_min": 150, "H_max": 600,  "n": 2950, "eff": 0.64, "stages": "3-12"},
    {"model": "D85-45",  "Q_min": 55,   "Q_max": 100,  "H_min": 135, "H_max": 540,  "n": 2950, "eff": 0.72, "stages": "3-12"},
    {"model": "D85-67",  "Q_min": 55,   "Q_max": 100,  "H_min": 200, "H_max": 800,  "n": 2950, "eff": 0.68, "stages": "3-12"},
    {"model": "D155-30", "Q_min": 100,  "Q_max": 190,  "H_min": 90,  "H_max": 360,  "n": 1480, "eff": 0.77, "stages": "3-12"},
    {"model": "D155-67", "Q_min": 100,  "Q_max": 190,  "H_min": 200, "H_max": 800,  "n": 2950, "eff": 0.74, "stages": "3-12"},
    {"model": "D280-43", "Q_min": 190,  "Q_max": 350,  "H_min": 129, "H_max": 516,  "n": 1480, "eff": 0.80, "stages": "3-12"},
    {"model": "D280-65", "Q_min": 190,  "Q_max": 350,  "H_min": 195, "H_max": 780,  "n": 1480, "eff": 0.78, "stages": "3-12"},
    {"model": "D450-60", "Q_min": 350,  "Q_max": 550,  "H_min": 180, "H_max": 720,  "n": 1480, "eff": 0.80, "stages": "3-12"},
]

# ── 磁力驱动泵标准系列 ── CQ 型不锈钢磁力泵（无泄漏）
# 流量 0.3–100 m³/h, 扬程 3–50 m, 2900 rpm
_CQ_MAGNETIC_SERIES = [
    {"model": "CQ16-12-50",  "Q_min": 0.3,  "Q_max": 0.8,  "H_min": 3,  "H_max": 8,   "n": 2900, "eff": 0.20, "material": "304/316L"},
    {"model": "CQ20-14-90",  "Q_min": 0.6,  "Q_max": 1.5,  "H_min": 8,  "H_max": 16,  "n": 2900, "eff": 0.25, "material": "304/316L"},
    {"model": "CQ25-20-100", "Q_min": 1.5,  "Q_max": 3.5,  "H_min": 8,  "H_max": 16,  "n": 2900, "eff": 0.32, "material": "304/316L"},
    {"model": "CQ32-20-110", "Q_min": 3.0,  "Q_max": 6.5,  "H_min": 10, "H_max": 20,  "n": 2900, "eff": 0.40, "material": "304/316L"},
    {"model": "CQ40-25-120", "Q_min": 6.0,  "Q_max": 12,   "H_min": 15, "H_max": 28,  "n": 2900, "eff": 0.48, "material": "304/316L"},
    {"model": "CQ40-25-160", "Q_min": 6.0,  "Q_max": 12,   "H_min": 22, "H_max": 38,  "n": 2900, "eff": 0.42, "material": "304/316L"},
    {"model": "CQ50-32-125", "Q_min": 10,   "Q_max": 20,   "H_min": 15, "H_max": 25,  "n": 2900, "eff": 0.55, "material": "304/316L"},
    {"model": "CQ50-32-160", "Q_min": 10,   "Q_max": 20,   "H_min": 25, "H_max": 38,  "n": 2900, "eff": 0.50, "material": "304/316L"},
    {"model": "CQ50-32-200", "Q_min": 10,   "Q_max": 20,   "H_min": 40, "H_max": 55,  "n": 2900, "eff": 0.44, "material": "304/316L"},
    {"model": "CQ65-50-125", "Q_min": 20,   "Q_max": 35,   "H_min": 15, "H_max": 25,  "n": 2900, "eff": 0.62, "material": "304/316L"},
    {"model": "CQ65-50-160", "Q_min": 20,   "Q_max": 35,   "H_min": 25, "H_max": 38,  "n": 2900, "eff": 0.58, "material": "304/316L"},
    {"model": "CQ80-65-125", "Q_min": 35,   "Q_max": 55,   "H_min": 15, "H_max": 25,  "n": 2900, "eff": 0.66, "material": "304/316L"},
    {"model": "CQ80-65-160", "Q_min": 35,   "Q_max": 55,   "H_min": 25, "H_max": 38,  "n": 2900, "eff": 0.63, "material": "304/316L"},
    {"model": "CQ100-80-125","Q_min": 55,   "Q_max": 100,  "H_min": 15, "H_max": 25,  "n": 2900, "eff": 0.70, "material": "304/316L"},
    {"model": "CQ100-80-160","Q_min": 55,   "Q_max": 100,  "H_min": 25, "H_max": 38,  "n": 2900, "eff": 0.68, "material": "304/316L"},
]

# ── 潜污泵标准系列 ── GB/T 24674, QW/WQ 系列潜水排污泵
# 流量 5–500 m³/h, 扬程 5–40 m
_QW_SEWAGE_SERIES = [
    {"model": "50QW10-10-0.75", "Q_min": 5,   "Q_max": 12,  "H_min": 7,  "H_max": 12,  "n": 2900, "eff": 0.45, "P_kW": 0.75},
    {"model": "50QW15-15-1.5",  "Q_min": 10,  "Q_max": 20,  "H_min": 10, "H_max": 18,  "n": 2900, "eff": 0.50, "P_kW": 1.5},
    {"model": "50QW25-10-1.5",  "Q_min": 15,  "Q_max": 32,  "H_min": 7,  "H_max": 13,  "n": 2900, "eff": 0.55, "P_kW": 1.5},
    {"model": "65QW25-15-2.2",  "Q_min": 18,  "Q_max": 35,  "H_min": 10, "H_max": 20,  "n": 2900, "eff": 0.55, "P_kW": 2.2},
    {"model": "65QW37-13-3",    "Q_min": 25,  "Q_max": 45,  "H_min": 8,  "H_max": 16,  "n": 2900, "eff": 0.60, "P_kW": 3.0},
    {"model": "80QW43-13-3",    "Q_min": 30,  "Q_max": 55,  "H_min": 8,  "H_max": 16,  "n": 2900, "eff": 0.63, "P_kW": 3.0},
    {"model": "80QW50-20-5.5",  "Q_min": 35,  "Q_max": 65,  "H_min": 14, "H_max": 25,  "n": 2900, "eff": 0.60, "P_kW": 5.5},
    {"model": "100QW65-15-5.5", "Q_min": 45,  "Q_max": 80,  "H_min": 10, "H_max": 20,  "n": 2900, "eff": 0.65, "P_kW": 5.5},
    {"model": "100QW70-22-7.5", "Q_min": 50,  "Q_max": 90,  "H_min": 15, "H_max": 28,  "n": 2900, "eff": 0.62, "P_kW": 7.5},
    {"model": "100QW100-15-7.5","Q_min": 70,  "Q_max": 130, "H_min": 10, "H_max": 20,  "n": 2900, "eff": 0.67, "P_kW": 7.5},
    {"model": "150QW100-25-11", "Q_min": 75,  "Q_max": 140, "H_min": 18, "H_max": 32,  "n": 1460, "eff": 0.65, "P_kW": 11},
    {"model": "150QW145-10-7.5","Q_min": 100, "Q_max": 180, "H_min": 7,  "H_max": 14,  "n": 1460, "eff": 0.70, "P_kW": 7.5},
    {"model": "150QW200-10-11", "Q_min": 140, "Q_max": 250, "H_min": 7,  "H_max": 14,  "n": 1460, "eff": 0.72, "P_kW": 11},
    {"model": "200QW250-15-18.5","Q_min":180,"Q_max": 320, "H_min": 10, "H_max": 20,  "n": 1460, "eff": 0.72, "P_kW": 18.5},
    {"model": "200QW300-10-15",  "Q_min": 220,"Q_max": 380, "H_min": 7,  "H_max": 14,  "n": 1460, "eff": 0.74, "P_kW": 15},
    {"model": "250QW400-15-30",  "Q_min": 280,"Q_max": 500, "H_min": 10, "H_max": 20,  "n": 980,  "eff": 0.75, "P_kW": 30},
]

# ── 自吸离心泵标准系列 ── ZX 系列卧式自吸泵（自吸高度 5-6.5 m）
# 流量 3–400 m³/h, 扬程 5–132 m, 2900/1450 rpm
_ZX_SELF_PRIMING_SERIES = [
    # 2900 rpm
    {"model": "ZX25-32-125",  "Q_min": 3,   "Q_max": 6.5,  "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.43, "self_prime_m": 5.0},
    {"model": "ZX25-32-160",  "Q_min": 3,   "Q_max": 6.5,  "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.37, "self_prime_m": 5.0},
    {"model": "ZX32-25-125",  "Q_min": 4,   "Q_max": 8.5,  "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.47, "self_prime_m": 5.0},
    {"model": "ZX32-25-160",  "Q_min": 4,   "Q_max": 8.5,  "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.42, "self_prime_m": 5.0},
    {"model": "ZX40-32-125",  "Q_min": 5,   "Q_max": 10,   "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.51, "self_prime_m": 5.5},
    {"model": "ZX40-32-160",  "Q_min": 5,   "Q_max": 10,   "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.47, "self_prime_m": 5.5},
    {"model": "ZX50-32-125",  "Q_min": 8,   "Q_max": 16,   "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.57, "self_prime_m": 5.5},
    {"model": "ZX50-32-160",  "Q_min": 8,   "Q_max": 16,   "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.54, "self_prime_m": 5.5},
    {"model": "ZX50-32-200",  "Q_min": 8,   "Q_max": 16,   "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.47, "self_prime_m": 5.5},
    {"model": "ZX65-50-125",  "Q_min": 15,  "Q_max": 30,   "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.62, "self_prime_m": 6.0},
    {"model": "ZX65-50-160",  "Q_min": 15,  "Q_max": 30,   "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.58, "self_prime_m": 6.0},
    {"model": "ZX65-40-200",  "Q_min": 15,  "Q_max": 30,   "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.52, "self_prime_m": 6.0},
    {"model": "ZX65-40-250",  "Q_min": 15,  "Q_max": 30,   "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.42, "self_prime_m": 6.0},
    {"model": "ZX80-65-125",  "Q_min": 30,  "Q_max": 60,   "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.68, "self_prime_m": 6.0},
    {"model": "ZX80-65-160",  "Q_min": 30,  "Q_max": 60,   "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.66, "self_prime_m": 6.0},
    {"model": "ZX80-50-200",  "Q_min": 30,  "Q_max": 60,   "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.62, "self_prime_m": 6.0},
    {"model": "ZX80-50-250",  "Q_min": 30,  "Q_max": 60,   "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.54, "self_prime_m": 6.0},
    {"model": "ZX100-80-125", "Q_min": 60,  "Q_max": 120,  "H_min": 5,  "H_max": 20,  "n": 2900, "eff": 0.73, "self_prime_m": 6.5},
    {"model": "ZX100-80-160", "Q_min": 60,  "Q_max": 120,  "H_min": 8,  "H_max": 32,  "n": 2900, "eff": 0.71, "self_prime_m": 6.5},
    {"model": "ZX100-65-200", "Q_min": 60,  "Q_max": 120,  "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.69, "self_prime_m": 6.5},
    {"model": "ZX100-65-250", "Q_min": 60,  "Q_max": 120,  "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.63, "self_prime_m": 6.5},
    {"model": "ZX125-100-200","Q_min": 120, "Q_max": 240,  "H_min": 12, "H_max": 50,  "n": 2900, "eff": 0.76, "self_prime_m": 6.5},
    {"model": "ZX125-100-250","Q_min": 120, "Q_max": 240,  "H_min": 20, "H_max": 80,  "n": 2900, "eff": 0.73, "self_prime_m": 6.5},
    {"model": "ZX125-100-315","Q_min": 120, "Q_max": 240,  "H_min": 32, "H_max": 125, "n": 1450, "eff": 0.69, "self_prime_m": 6.0},
    # 1450 rpm 大流量
    {"model": "ZX150-125-250","Q_min": 200, "Q_max": 400,  "H_min": 20, "H_max": 80,  "n": 1450, "eff": 0.76, "self_prime_m": 5.5},
]

# ── 氟塑料衬里离心泵标准系列 ── IHF 系列 (F46/FEP 衬里, 耐强酸强碱)
# 流量 3.5–130 m³/h, 扬程 4.5–87 m, 2900/1450 rpm
_IHF_FLUOROPLASTIC_SERIES = [
    # 2900 rpm
    {"model": "IHF32-25-125",  "Q_min": 3.5,  "Q_max": 6.5,  "H_min": 18, "H_max": 21,  "n": 2900, "eff": 0.40, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF32-25-160",  "Q_min": 3.5,  "Q_max": 6.5,  "H_min": 30, "H_max": 33,  "n": 2900, "eff": 0.34, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF40-32-125",  "Q_min": 4.4,  "Q_max": 8.3,  "H_min": 18, "H_max": 21,  "n": 2900, "eff": 0.40, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF40-32-160",  "Q_min": 4.4,  "Q_max": 8.3,  "H_min": 30, "H_max": 33,  "n": 2900, "eff": 0.34, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF50-32-125",  "Q_min": 8.8,  "Q_max": 16.3, "H_min": 17.5,"H_max": 21.5,"n": 2900, "eff": 0.45, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF50-32-160",  "Q_min": 8.8,  "Q_max": 16.3, "H_min": 30, "H_max": 33,  "n": 2900, "eff": 0.41, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF50-32-200",  "Q_min": 8.8,  "Q_max": 16.3, "H_min": 48, "H_max": 52,  "n": 2900, "eff": 0.34, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF50-32-250",  "Q_min": 8.8,  "Q_max": 16.3, "H_min": 76, "H_max": 82,  "n": 2900, "eff": 0.27, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF65-50-125",  "Q_min": 17.5, "Q_max": 32,   "H_min": 17.5,"H_max": 21.5,"n": 2900, "eff": 0.56, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF65-50-160",  "Q_min": 17.5, "Q_max": 32,   "H_min": 27.5,"H_max": 33,  "n": 2900, "eff": 0.50, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF65-40-200",  "Q_min": 17.5, "Q_max": 32,   "H_min": 45.5,"H_max": 52,  "n": 2900, "eff": 0.45, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF65-40-250",  "Q_min": 17.5, "Q_max": 32,   "H_min": 76, "H_max": 82,  "n": 2900, "eff": 0.39, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF80-65-125",  "Q_min": 35,   "Q_max": 65,   "H_min": 17, "H_max": 21.5,"n": 2900, "eff": 0.64, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF80-65-160",  "Q_min": 35,   "Q_max": 65,   "H_min": 27.5,"H_max": 33,  "n": 2900, "eff": 0.60, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF80-50-200",  "Q_min": 35,   "Q_max": 65,   "H_min": 45.5,"H_max": 52,  "n": 2900, "eff": 0.52, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF80-50-250",  "Q_min": 35,   "Q_max": 65,   "H_min": 72, "H_max": 82,  "n": 2900, "eff": 0.40, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF100-80-125", "Q_min": 70,   "Q_max": 130,  "H_min": 14, "H_max": 23,  "n": 2900, "eff": 0.65, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF100-80-160", "Q_min": 70,   "Q_max": 130,  "H_min": 24, "H_max": 34,  "n": 2900, "eff": 0.65, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF100-65-200", "Q_min": 70,   "Q_max": 130,  "H_min": 42, "H_max": 52,  "n": 2900, "eff": 0.64, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF100-65-250", "Q_min": 70,   "Q_max": 130,  "H_min": 68, "H_max": 87,  "n": 2900, "eff": 0.62, "material": "F46/FEP", "temp_max": 150},
    # 1450 rpm
    {"model": "IHF125-100-200","Q_min": 130,  "Q_max": 240,  "H_min": 10, "H_max": 13,  "n": 1450, "eff": 0.63, "material": "F46/FEP", "temp_max": 150},
    {"model": "IHF125-100-250","Q_min": 130,  "Q_max": 240,  "H_min": 18, "H_max": 22,  "n": 1450, "eff": 0.63, "material": "F46/FEP", "temp_max": 150},
]

# 液体物性数据库（典型值）
_LIQUID_PROPERTIES = {
    "water": {"rho": 1000, "mu": 0.001, "vapor_pressure_20C": 2338},
    "diesel": {"rho": 850, "mu": 0.003, "vapor_pressure_20C": 500},
    "gasoline": {"rho": 740, "mu": 0.0006, "vapor_pressure_20C": 60000},
    "methanol": {"rho": 792, "mu": 0.0006, "vapor_pressure_20C": 12800},
    "ethanol": {"rho": 789, "mu": 0.0012, "vapor_pressure_20C": 5900},
    "benzene": {"rho": 879, "mu": 0.00065, "vapor_pressure_20C": 10000},
}


# ============================================================
# 泵功率计算
# ============================================================

def pump_hydraulic_power(Q: float, H: float, rho: float = 1000.0) -> float:
    """计算水力功率
    
    Args:
        Q: 体积流量 (m³/s)
        H: 扬程 (m)
        rho: 液体密度 (kg/m³)
        
    Returns:
        水力功率 (W)
    """
    g = 9.81
    P_hydraulic = rho * g * Q * H
    return P_hydraulic


def pump_shaft_power(P_hydraulic: float, efficiency: float) -> float:
    """计算轴功率
    
    Args:
        P_hydraulic: 水力功率 (W)
        efficiency: 泵效率
        
    Returns:
        轴功率 (W)
    """
    if efficiency <= 0 or efficiency > 1:
        raise ValueError("效率必须在 (0, 1] 范围内")
    return P_hydraulic / efficiency


def motor_power(P_shaft: float, motor_efficiency: float = 0.9,
                service_factor: float = 1.1) -> float:
    """计算电机功率（含安全系数）
    
    Args:
        P_shaft: 轴功率 (W)
        motor_efficiency: 电机效率
        service_factor: 服务系数（安全系数）
        
    Returns:
        电机功率 (W)
    """
    P_motor = P_shaft / motor_efficiency * service_factor
    return P_motor


def pump_power_full(Q: float, H: float, rho: float = 1000.0,
                    efficiency: float = 0.7, motor_efficiency: float = 0.9,
                    service_factor: float = 1.1) -> dict:
    """完整泵功率计算
    
    Returns:
        {
            "P_hydraulic": 水力功率 (W),
            "P_shaft": 轴功率 (W),
            "P_motor": 电机功率 (W),
            "efficiency": 泵效率
        }
    """
    P_hyd = pump_hydraulic_power(Q, H, rho)
    P_shaft = pump_shaft_power(P_hyd, efficiency)
    P_motor = motor_power(P_shaft, motor_efficiency, service_factor)
    
    return {
        "P_hydraulic": P_hyd,
        "P_shaft": P_shaft,
        "P_motor": P_motor,
        "efficiency": efficiency,
        "motor_efficiency": motor_efficiency,
        "service_factor": service_factor,
    }


# ============================================================
# 系统总扬程计算
# ============================================================

def calculate_system_head(
    delta_z: float = 0.0,
    delta_P: float = 0.0,
    h_friction: float = 0.0,
    rho: float = None,
    # ── 可选：用户给定流体信息，自动查询物性 ──
    fluid: str = None,
    T: float = None,
    P: float = 101325.0,
) -> dict:
    """系统总扬程计算
    
    H_total = Δz + ΔP/(ρg) + h_friction
    
    Args:
        delta_z: 排出液面与吸入液面高差 (m)，排出高于吸入为正
        delta_P: 排出与吸入空间压差 (Pa)，排出压力高为正
        h_friction: 管路总阻力损失 (m液柱)，含直管摩擦+局部阻力
        rho: 流体密度 (kg/m³)，给定后跳过查询
        fluid: 流体名称，用于自动查询密度
        T: 温度 (K)
        P: 压力 (Pa)，默认 101325
    
    Returns:
        系统总扬程及各分量
    """
    from physics_engine.thermo_helper import get_fluid_density
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    delta_z = safe_float(delta_z, default=0.0, name='delta_z')
    delta_P = safe_float(delta_P, default=0.0, name='delta_P')
    h_friction = safe_float(h_friction, default=0.0, name='h_friction')
    rho = safe_float(rho, name='rho')
    T = safe_float(T, name='T')
    P = safe_float(P, default=101325.0, name='P')
    
    conversion_notes = []
    g = 9.81
    
    # ── 密度：用户给定值优先 ──
    if rho is not None:
        rho_val = float(rho)
        conversion_notes.append(f"[用户给定] 密度 rho = {rho_val} kg/m³")
    elif fluid is not None and T is not None:
        try:
            rho_val = get_fluid_density(fluid, T, P, phase="liquid")
            conversion_notes.append(f"[自动查询] {fluid} 在 {T}K 的密度 = {rho_val:.2f} kg/m³")
        except Exception as e:
            return {"error": f"密度查询失败：{e}。请手动提供 rho 参数。"}
    else:
        rho_val = 1000.0
        conversion_notes.append("[默认] 使用水密度 1000 kg/m³")
    
    # 静扬程
    H_static = delta_z
    
    # 压力扬程
    H_pressure = delta_P / (rho_val * g) if rho_val > 0 else 0.0
    
    # 阻力损失扬程
    H_friction = h_friction
    
    # 总扬程
    H_total = H_static + H_pressure + H_friction
    
    return {
        "H_total_m": round(H_total, 3),
        "H_static_m": round(H_static, 3),
        "H_pressure_m": round(H_pressure, 3),
        "H_friction_m": round(H_friction, 3),
        "rho_kg_m3": round(rho_val, 2),
        "conversion_notes": conversion_notes,
    }


# ============================================================
# 离心泵功率估算（便捷接口）
# ============================================================

def size_centrifugal_pump_power(
    Q: float,
    H: float,
    rho: float = None,
    delta_P: Optional[float] = None,
    efficiency: Optional[float] = None,
    # ── 可选：用户给定流体信息，自动查询物性 ──
    fluid: str = None,
    T: float = None,
    P: float = 101325.0,
) -> dict:
    """离心泵功率估算
    
    由流量、扬程、压差计算轴功率。
    密度支持用户给定优先，缺失时自动查询或默认1000。
    
    Args:
        Q: 体积流量 (m³/s)
        H: 扬程 (m)
        rho: 流体密度 (kg/m³)，给定后跳过查询
        delta_P: 压差 (Pa)，可选，若提供则 H = delta_P / (rho·g)
        efficiency: 泵效率，不指定则自动估算
        fluid: 流体名称，用于自动查询密度
        T: 温度 (K)
        P: 压力 (Pa)，默认 101325
        
    Returns:
        {
            "P_hydraulic_W": 水力功率 (W),
            "P_shaft_W": 轴功率 (W),
            "P_motor_W": 电机功率 (W),
            "efficiency": 泵效率,
            "H_m": 扬程 (m),
        }
    """
    from physics_engine.thermo_helper import get_fluid_density
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    Q = safe_float(Q, name='Q')
    H = safe_float(H, name='H')
    rho = safe_float(rho, name='rho')
    delta_P = safe_float(delta_P, name='delta_P')
    efficiency = safe_float(efficiency, name='efficiency')
    T = safe_float(T, name='T')
    P = safe_float(P, default=101325.0, name='P')
    
    conversion_notes = []
    g = 9.81
    
    # ── 密度：用户给定值优先 ──
    if rho is not None:
        rho_val = float(rho)
        conversion_notes.append(f"[用户给定] 密度 rho = {rho_val} kg/m³")
    elif fluid is not None and T is not None:
        try:
            rho_val = get_fluid_density(fluid, T, P, phase="liquid")
            conversion_notes.append(f"[自动查询] {fluid} 在 {T}K 的密度 = {rho_val:.2f} kg/m³")
        except Exception:
            rho_val = 1000.0
            conversion_notes.append(f"[默认] {fluid} 密度查询失败，使用 1000 kg/m³")
    else:
        rho_val = 1000.0
        conversion_notes.append("[默认] 使用水密度 1000 kg/m³")
    
    if delta_P is not None:
        H = delta_P / (rho_val * g)
    
    if efficiency is None:
        efficiency = estimate_pump_efficiency(Q, H)
    
    P_hyd = pump_hydraulic_power(Q, H, rho_val)
    P_shaft = pump_shaft_power(P_hyd, efficiency)
    P_motor = motor_power(P_shaft)
    
    return {
        "P_hydraulic_W": P_hyd,
        "P_shaft_W": P_shaft,
        "P_motor_W": P_motor,
        "efficiency": efficiency,
        "H_m": H,
        "Q_m3_s": Q,
        "rho_kg_m3": rho_val,
        "conversion_notes": conversion_notes,
    }


# ============================================================
# 泵 NPSH 校核（便捷接口）
# ============================================================

def calculate_pump_npsh(
    P_suction: float,
    P_vapor: float = None,
    rho: float = None,
    suction_head: float = 0.0,
    friction_loss: float = 0.0,
    Q: float = 0.0,
    n: float = 2900,
    margin: float = 0.5,
    # ── 可选：用户给定流体信息，自动查询物性 ──
    fluid: str = None,
    T: float = None,
    P: float = 101325.0,
) -> dict:
    """泵 NPSH 校核
    
    计算可用 NPSH_a，判断是否汽蚀。
    密度和蒸气压支持用户给定优先，缺失时自动查询。
    
    Args:
        P_suction: 泵入口压力 (Pa, 绝对)
        P_vapor: 液体蒸气压 (Pa, 绝对)，给定后跳过查询
        rho: 液体密度 (kg/m³)，给定后跳过查询
        suction_head: 吸入高度 (m)，正=倒灌
        friction_loss: 吸入管路摩擦损失 (m)
        Q: 流量 (m³/s)
        n: 转速 (rpm)
        margin: 安全余量 (m)
        fluid: 流体名称，用于自动查询
        T: 温度 (K)
        P: 压力 (Pa)，默认 101325
        
    Returns:
        {
            "NPSHa_m": 可用汽蚀余量 (m),
            "NPSHr_m": 必需汽蚀余量 (m),
            "margin_m": 实际余量 (m),
            "safe": 是否安全,
            "cavitation_risk": 汽蚀风险,
        }
    """
    from physics_engine.thermo_helper import get_fluid_density, get_fluid_vapor_pressure
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    P_suction = safe_float(P_suction, name='P_suction')
    P_vapor = safe_float(P_vapor, name='P_vapor')
    rho = safe_float(rho, name='rho')
    suction_head = safe_float(suction_head, default=0.0, name='suction_head')
    friction_loss = safe_float(friction_loss, default=0.0, name='friction_loss')
    Q = safe_float(Q, default=0.0, name='Q')
    n = safe_float(n, default=2900, name='n')
    margin = safe_float(margin, default=0.5, name='margin')
    NPSH_required = safe_float(NPSH_required, name='NPSH_required')
    T = safe_float(T, name='T')
    P = safe_float(P, default=101325.0, name='P')
    
    conversion_notes = []
    
    # ── 密度：用户给定值优先 ──
    if rho is not None:
        rho_val = float(rho)
        conversion_notes.append(f"[用户给定] 密度 rho = {rho_val} kg/m³")
    elif fluid is not None and T is not None:
        try:
            rho_val = get_fluid_density(fluid, T, P, phase="liquid")
            conversion_notes.append(f"[自动查询] {fluid} 在 {T}K 的密度 = {rho_val:.2f} kg/m³")
        except Exception:
            rho_val = 1000.0
            conversion_notes.append(f"[默认] {fluid} 密度查询失败，使用 1000 kg/m³")
    else:
        rho_val = 1000.0
        conversion_notes.append("[默认] 使用水密度 1000 kg/m³")
    
    # ── 蒸气压：用户给定值优先 ──
    if P_vapor is not None:
        Pv_val = float(P_vapor)
        conversion_notes.append(f"[用户给定] 蒸气压 Pv = {Pv_val} Pa")
    elif fluid is not None and T is not None:
        try:
            Pv_val = get_fluid_vapor_pressure(fluid, T)
            conversion_notes.append(f"[自动查询] {fluid} 在 {T}K 的蒸气压 = {Pv_val:.1f} Pa")
        except Exception:
            return {"error": f"蒸气压查询失败。请手动提供 P_vapor 参数。"}
    else:
        return {"error": "必须提供 P_vapor 或 (fluid + T) 参数"}
    
    npsha_val = npsha(P_suction, Pv_val, rho_val, suction_head, friction_loss)
    npshr_val = npshr_estimate(Q, n)
    check_result = npsh_check(npsha_val, npshr_val, margin)
    
    if check_result["safe"]:
        risk = "low"
    elif npsha_val > npshr_val:
        risk = "moderate"
    else:
        risk = "high"
    
    return {
        "NPSHa_m": npsha_val,
        "NPSHr_m": npshr_val,
        "margin_m": check_result["margin_actual"],
        "safe": check_result["safe"],
        "cavitation_risk": risk,
        "conversion_notes": conversion_notes,
    }


# ============================================================
# 泵型推荐（按比转速）
# ============================================================

def select_pump_type(
    Q: float,
    H: float,
    n: float = 2900,
) -> dict:
    """按比转速 n_s 推荐泵型
    
    比转速 n_s = n·√Q / H^(3/4)（n in rpm, Q in m³/s, H in m）
    
    Args:
        Q: 流量 (m³/s)
        H: 扬程 (m)
        n: 转速 (rpm)
        
    Returns:
        {
            "recommended_type": 推荐泵型,
            "specific_speed": 比转速,
            "reason": 推荐理由,
        }
    """
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    Q = safe_float(Q, name='Q')
    H = safe_float(H, name='H')
    n = safe_float(n, default=2900, name='n')
    
    if Q <= 0 or H <= 0:
        return {"error": "流量和扬程必须为正数"}
    
    # 比转速（无量纲形式）
    n_s = n * math.sqrt(Q) / (H ** 0.75)
    
    if n_s < 500:
        pump_type = "centrifugal_radial"
        reason = "低比转速，适合径向叶轮离心泵（高扬程小流量）"
    elif n_s < 1500:
        pump_type = "centrifugal"
        reason = "中低比转速，适合标准离心泵"
    elif n_s < 3000:
        pump_type = "centrifugal_mixed"
        reason = "中高比转速，适合混流离心泵"
    elif n_s < 6000:
        pump_type = "mixed_flow"
        reason = "高比转速，适合混流泵"
    else:
        pump_type = "axial_flow"
        reason = "极高比转速，适合轴流泵（大流量低扬程）"
    
    return {
        "recommended_type": pump_type,
        "specific_speed": n_s,
        "reason": reason,
    }


# ============================================================
# NPSH 计算
# ============================================================

def npsha(P_suction: float, P_vapor: float, rho: float,
           suction_head: float, friction_loss: float = 0.0) -> float:
    """计算有效汽蚀余量 NPSHa
    
    NPSHa = (P_suction - P_vapor) / (ρg) + suction_head - friction_loss
    
    Args:
        P_suction: 泵入口压力 (Pa)
        P_vapor: 液体蒸气压 (Pa)
        rho: 液体密度 (kg/m³)
        suction_head: 吸入高度（正=倒灌，负=吸上）(m)
        friction_loss: 吸入管路摩擦损失 (m)
        
    Returns:
        NPSHa (m)
    """
    g = 9.81
    npsha_value = (P_suction - P_vapor) / (rho * g) + suction_head - friction_loss
    return npsha_value


def npshr_estimate(Q: float, n: float, NPSHR_type: str = "centrifugal") -> float:
    """估算必需汽蚀余量 NPSHr（简化关联式）
    
    Args:
        Q: 流量 (m³/s)
        n: 泵转速 (rpm)
        NPSHR_type: 泵类型
        
    Returns:
        NPSHr (m)
    """
    if NPSHR_type == "centrifugal":
        # 简化公式：NPSHr ≈ k * Q^0.5 / n^1.5
        # k 为常数，取决于泵设计
        k = 0.5
        n_rps = n / 60  # 转每秒
        NPSHr = k * math.sqrt(Q) / (n_rps ** 1.5)
        return max(NPSHr, 2.0)  # 最小 2m
    else:
        # 容积泵通常有更低的 NPSHr
        return 1.0


def npsh_check(npsha_value: float, npshr_value: float,
               margin: float = 0.5) -> dict:
    """NPSH 校核
    
    Args:
        npsha_value: NPSHa (m)
        npshr_value: NPSHr (m)
        margin: 安全余量 (m)
        
    Returns:
        {
            "safe": 是否安全,
            "margin_actual": 实际余量,
            "margin_required": 要求余量,
            "npsha": NPSHa,
            "npshr": NPSHr
        }
    """
    actual_margin = npsha_value - npshr_value
    safe = actual_margin >= margin
    
    return {
        "safe": safe,
        "margin_actual": actual_margin,
        "margin_required": margin,
        "npsha": npsha_value,
        "npshr": npshr_value,
    }


# ============================================================
# 泵效率估算
# ============================================================

def estimate_pump_efficiency(Q: float, H: float, n: float = 2900) -> float:
    """估算泵效率（基于经验关联式）
    
    Args:
        Q: 流量 (m³/s)
        H: 扬程 (m)
        n: 转速 (rpm)
        
    Returns:
        效率 (0-1)
    """
    # 比转速
    n_rps = n / 60
    N_s = n_rps * math.sqrt(Q) / (H ** 0.75)  # 无量纲比转速（简化）
    
    # 效率与比转速的关系（离心泵）
    if N_s < 500:  # 低比转速
        eff = 0.65
    elif N_s < 1500:  # 中比转速
        eff = 0.75
    elif N_s < 3000:  # 高比转速
        eff = 0.80
    else:  # 混流泵/轴流泵
        eff = 0.85
    
    # 流量修正（小流量泵效率低）
    if Q < 0.005:  # < 18 m³/h
        eff *= 0.9
    elif Q > 0.1:  # > 360 m³/h
        eff *= 1.02
    
    return min(eff, 0.90)


# ============================================================
# 泵选型
# ============================================================

def select_centrifugal_pump(Q_required: float, H_required: float,
                             n: float = 2900) -> dict:
    """选择离心泵型号
    
    Args:
        Q_required: 所需流量 (m³/s)
        H_required: 所需扬程 (m)
        n: 转速 (rpm)
        
    Returns:
        {
            "model": 型号,
            "Q_rated": 额定流量 (m³/s),
            "H_rated": 额定扬程 (m),
            "efficiency": 效率,
            "n": 转速 (rpm),
            "suitable": 是否合适
        }
    """
    Q_m3h = Q_required * 3600  # 转 m³/h
    
    best_match = None
    best_score = float('inf')
    
    for pump in _CENTRIFUGAL_PUMP_SERIES:
        # 检查流量范围
        if pump["Q_min"] <= Q_m3h <= pump["Q_max"]:
            # 检查扬程范围
            if pump["H_min"] <= H_required <= pump["H_max"]:
                # 计算匹配度（距离）
                Q_center = (pump["Q_min"] + pump["Q_max"]) / 2
                H_center = (pump["H_min"] + pump["H_max"]) / 2
                score = math.sqrt(((Q_m3h - Q_center) / max(Q_center, 1))**2 + 
                               ((H_required - H_center) / max(H_center, 1))**2)
                
                if score < best_score:
                    best_score = score
                    best_match = pump.copy()
                    best_match["suitable"] = True
                    best_match["Q_rated"] = Q_m3h
                    best_match["H_rated"] = H_required
        
        # 如果没找到完全匹配的，找最接近的
        if best_match is None:
            for pump in _CENTRIFUGAL_PUMP_SERIES:
                if Q_m3h <= pump["Q_max"] and H_required <= pump["H_max"]:
                    best_match = pump.copy()
                    best_match["suitable"] = True
                    best_match["Q_rated"] = Q_m3h
                    best_match["H_rated"] = H_required
                    break
        
        if best_match is None:
            # 返回最大泵
            pump = _CENTRIFUGAL_PUMP_SERIES[-1]
            best_match = pump.copy()
            best_match["suitable"] = False
            best_match["Q_rated"] = Q_m3h
            best_match["H_rated"] = H_required
            best_match["warning"] = "超出标准系列范围"
        
        # 转换单位
        best_match["Q_rated_m3s"] = best_match["Q_rated"] / 3600
        best_match["n"] = pump.get("n", n)
        best_match["efficiency"] = pump.get("eff", 0.7)
        
        return best_match


def _select_pump_from_series(series: list, Q_m3h: float, H_required: float,
                              pump_type_name: str) -> dict:
    """通用泵型选型 — 从系列目录中匹配最佳型号

    Args:
        series: 泵系列数据列表
        Q_m3h: 所需流量 (m³/h)
        H_required: 所需扬程 (m)
        pump_type_name: 泵类型名称（用于返回结果）

    Returns:
        匹配结果字典
    """
    best_match = None
    best_score = float('inf')

    # 精确匹配：流量和扬程都在范围内
    for pump in series:
        if pump["Q_min"] <= Q_m3h <= pump["Q_max"] and pump["H_min"] <= H_required <= pump["H_max"]:
            Q_center = (pump["Q_min"] + pump["Q_max"]) / 2
            H_center = (pump["H_min"] + pump["H_max"]) / 2
            score = math.sqrt(
                ((Q_m3h - Q_center) / max(Q_center, 0.1)) ** 2 +
                ((H_required - H_center) / max(H_center, 0.1)) ** 2
            )
            if score < best_score:
                best_score = score
                best_match = pump.copy()
                best_match["suitable"] = True
                best_match["Q_rated"] = Q_m3h
                best_match["H_rated"] = H_required

    # 宽松匹配：流量在范围内，扬程不超过上限
    if best_match is None:
        for pump in series:
            if pump["Q_min"] <= Q_m3h <= pump["Q_max"] and H_required <= pump["H_max"]:
                best_match = pump.copy()
                best_match["suitable"] = True
                best_match["Q_rated"] = Q_m3h
                best_match["H_rated"] = H_required
                break

    # 超出范围：返回最大型号
    if best_match is None:
        pump = series[-1]
        best_match = pump.copy()
        best_match["suitable"] = False
        best_match["Q_rated"] = Q_m3h
        best_match["H_rated"] = H_required
        best_match["warning"] = f"超出{pump_type_name}标准系列范围"

    best_match["Q_rated_m3s"] = best_match["Q_rated"] / 3600
    best_match["n"] = best_match.get("n", 0)
    best_match["efficiency"] = best_match.get("eff", 0.5)
    return best_match


def select_metering_pump(Q_required: float, H_required: float) -> dict:
    """选择计量泵型号 — GB/T 7782, J/JZ/JD/JY 系列"""
    return _select_pump_from_series(_METERING_PUMP_SERIES, Q_required * 3600, H_required, "计量泵")


def select_vortex_pump(Q_required: float, H_required: float) -> dict:
    """选择漩涡泵型号 — GB/T 29531, W/WB 系列"""
    return _select_pump_from_series(_VORTEX_PUMP_SERIES, Q_required * 3600, H_required, "漩涡泵")


def select_gear_pump(Q_required: float, H_required: float) -> dict:
    """选择齿轮泵型号 — GB/T 25631, KCB/2CY 系列"""
    return _select_pump_from_series(_GEAR_PUMP_SERIES, Q_required * 3600, H_required, "齿轮泵")


def select_screw_pump(Q_required: float, H_required: float) -> dict:
    """选择螺杆泵型号 — GB/T 25695, G/GF 系列"""
    return _select_pump_from_series(_SCREW_PUMP_SERIES, Q_required * 3600, H_required, "螺杆泵")


def select_mixed_axial_pump(Q_required: float, H_required: float) -> dict:
    """选择混流泵/轴流泵型号 — GB/T 13008, HW/HB/ZLB 系列"""
    return _select_pump_from_series(_MIXED_AXIAL_PUMP_SERIES, Q_required * 3600, H_required, "混流泵/轴流泵")


def select_diaphragm_pump(Q_required: float, H_required: float) -> dict:
    """选择隔膜泵型号 — QBY/DBY 系列"""
    return _select_pump_from_series(_DIAPHRAGM_PUMP_SERIES, Q_required * 3600, H_required, "隔膜泵")


def select_chemical_pump(Q_required: float, H_required: float) -> dict:
    """选择化工离心泵型号 — GB/T 5656, IH 系列"""
    return _select_pump_from_series(_IH_CHEMICAL_PUMP_SERIES, Q_required * 3600, H_required, "化工离心泵")


def select_double_suction_pump(Q_required: float, H_required: float) -> dict:
    """选择双吸离心泵型号 — S/SH 系列"""
    return _select_pump_from_series(_S_DOUBLE_SUCTION_SERIES, Q_required * 3600, H_required, "双吸离心泵")


def select_multistage_pump(Q_required: float, H_required: float) -> dict:
    """选择多级离心泵型号 — GB/T 5657, D 系列

    多级泵扬程范围较宽 (3-12级), 匹配时以总扬程为准。
    """
    return _select_pump_from_series(_D_MULTISTAGE_SERIES, Q_required * 3600, H_required, "多级离心泵")


def select_magnetic_pump(Q_required: float, H_required: float) -> dict:
    """选择磁力驱动泵型号 — CQ 系列"""
    return _select_pump_from_series(_CQ_MAGNETIC_SERIES, Q_required * 3600, H_required, "磁力驱动泵")


def select_sewage_pump(Q_required: float, H_required: float) -> dict:
    """选择潜污泵型号 — GB/T 24674, QW/WQ 系列"""
    return _select_pump_from_series(_QW_SEWAGE_SERIES, Q_required * 3600, H_required, "潜污泵")


def select_self_priming_pump(Q_required: float, H_required: float) -> dict:
    """选择自吸离心泵型号 — ZX 系列 (自吸高度 5-6.5 m)"""
    return _select_pump_from_series(_ZX_SELF_PRIMING_SERIES, Q_required * 3600, H_required, "自吸离心泵")


def select_fluoroplastic_pump(Q_required: float, H_required: float) -> dict:
    """选择氟塑料衬里离心泵型号 — IHF 系列 (耐强酸强碱)"""
    return _select_pump_from_series(_IHF_FLUOROPLASTIC_SERIES, Q_required * 3600, H_required, "氟塑料衬里离心泵")


def select_pump_by_type(pump_type: str, Q_required: float, H_required: float) -> dict:
    """根据泵类型调用对应的选型函数

    Args:
        pump_type: 泵类型常量
        Q_required: 所需流量 (m³/s)
        H_required: 所需扬程 (m)

    Returns:
        选型结果，无匹配时返回 None
    """
    selectors = {
        PUMP_CENTRIFUGAL: select_centrifugal_pump,
        PUMP_CHEMICAL: select_chemical_pump,
        PUMP_DOUBLE_SUCTION: select_double_suction_pump,
        PUMP_MULTISTAGE: select_multistage_pump,
        PUMP_MAGNETIC: select_magnetic_pump,
        PUMP_SEWAGE: select_sewage_pump,
        PUMP_METERING: select_metering_pump,
        PUMP_VORTEX: select_vortex_pump,
        PUMP_GEAR: select_gear_pump,
        PUMP_SCREW: select_screw_pump,
        PUMP_POSITIVE_DISPLACEMENT: select_gear_pump,  # 容积式默认齿轮泵
        PUMP_MIXED_AXIAL: select_mixed_axial_pump,
        PUMP_DIAPHRAGM: select_diaphragm_pump,
        PUMP_SELF_PRIMING: select_self_priming_pump,
        PUMP_FLUOROPLASTIC: select_fluoroplastic_pump,
    }
    selector = selectors.get(pump_type)
    if selector:
        return selector(Q_required, H_required)
    return None


def suggest_pump_type(Q: float, H: float, viscosity: float = 0.001,
                      rho: float = 1000.0,
                      self_priming: bool = False, metering_required: bool = False,
                      gas_content: float = 0.0) -> dict:
    """基于流程图的泵自动选型算法

    决策流程：
    1. 计量要求 → 计量泵
    2. 黏度 > 650 mm²/s → 容积式泵 (Q≤600) 或超出范围
    3. 扬程与比转数 ns 判断：
       - H ≤ 20: 计算 ns，ns≥300 → 混流泵/轴流泵
       - H > 20: 计算 ns，ns≥300 → 混流泵/轴流泵，ns<300 → 进入含气量判断
    4. 含气量 > 5% → 根据 H、Q、黏度选择漩涡泵或容积式泵
    5. 含气量 ≤ 5% → 离心泵

    泵类型选择依据：
    - 计量泵: 需要精确计量输送流量的工况（如加药、配料）
    - 容积式泵: 高黏度介质(>650 mm²/s)、含气量高(>5%)、需自吸能力
    - 漩涡泵: 含气量>5%、低扬程(H≤150)、小流量(Q≤10 m³/h)、低黏度(≤37.4 mm²/s)
    - 离心泵: 常规工况，含气量≤5%，黏度适中，流量范围广
    - 混流泵/轴流泵: 比转数 ns≥300，大流量低扬程工况

    Args:
        Q: 流量 (m³/s)
        H: 扬程 (m)
        viscosity: 动力粘度 (Pa·s)
        rho: 液体密度 (kg/m³), 用于计算运动粘度
        self_priming: 是否需要自吸
        metering_required: 是否有计量要求
        gas_content: 含气量百分比 (例如 6% 输入 6)

    Returns:
        {
            "pump_type": 建议泵类型,
            "reason": 选择依据,
            "alternative": 备选类型,
            "ns": 比转数（如果计算）
        }
    """
    Q_m3h = Q * 3600
    # 运动粘度 mm²/s = 动力粘度 Pa·s / 密度 kg/m³ × 1e6
    viscosity_mm2s = viscosity / max(rho, 0.1) * 1e6

    # ── 1. 计量要求判断 ──
    # 适用场景：化工加药、食品配料、水处理等需要精确计量的工况
    if metering_required:
        return {
            "pump_type": PUMP_METERING,
            "reason": "有计量要求，选择计量泵（精度高，可调节流量）",
            "alternative": PUMP_DIAPHRAGM,
            "ns": None,
        }

    # ── 2. 黏度超高判断 ──
    # 黏度 > 650 mm²/s 时，离心泵效率急剧下降，必须用容积式泵
    if viscosity_mm2s > 650:
        if Q_m3h <= 600:
            return {
                "pump_type": PUMP_POSITIVE_DISPLACEMENT,
                "reason": f"黏度超高 ({viscosity_mm2s:.0f} mm²/s > 650)，离心泵效率过低，"
                          f"选择容积式泵（齿轮泵/螺杆泵）",
                "alternative": PUMP_SCREW,
                "ns": None,
            }
        else:
            return {
                "pump_type": "超出范围",
                "reason": f"黏度 {viscosity_mm2s:.0f} mm²/s 且流量 {Q_m3h:.0f} m³/h > 600，"
                          f"超出常规泵选型范围，需特殊设计",
                "alternative": None,
                "ns": None,
            }

    # ── 3. 扬程与比转数判断 ──
    # 比转数 ns = 3.65*n*sqrt(Q_m³/s) / H^(3/4)
    # 径向离心泵: ns=10-80, 混流泵: ns=80-160, 轴流泵: ns=160-400
    if H <= 0:
        raise ValueError("扬程(H)必须大于0")

    n = 1450  # 固定转速 r/min
    ns = (3.65 * n * math.sqrt(Q)) / (H ** 0.75)  # Q in m³/s

    # 高扬程 (>125m) → 多级泵
    if H > 125:
        return {
            "pump_type": PUMP_MULTISTAGE,
            "reason": f"扬程 H={H:.1f} m > 125 m，单级离心泵难以满足，"
                      f"选择多级离心泵 (D/DG 系列)",
            "alternative": PUMP_POSITIVE_DISPLACEMENT,
            "ns": ns,
        }

    # 比转数判断: ns >= 160 → 轴流泵; ns >= 80 → 混流泵
    if ns >= 160:
        return {
            "pump_type": PUMP_MIXED_AXIAL,
            "reason": f"比转数 ns={ns:.1f} ≥ 160，大流量低扬程工况，"
                      f"选择混流泵/轴流泵 (HW/ZLB 系列)",
            "alternative": PUMP_DOUBLE_SUCTION,
            "ns": ns,
        }

    if ns >= 80:
        # 中比转数: 双吸泵优先 (大流量)
        if Q_m3h >= 70:
            return {
                "pump_type": PUMP_DOUBLE_SUCTION,
                "reason": f"比转数 ns={ns:.1f} ∈ [80,160)，流量 Q={Q_m3h:.1f} m³/h ≥ 70，"
                          f"推荐双吸离心泵 (S 系列) 以平衡轴向力",
                "alternative": PUMP_CENTRIFUGAL,
                "ns": ns,
            }
        return {
            "pump_type": PUMP_CENTRIFUGAL,
            "reason": f"比转数 ns={ns:.1f} ∈ [80,160)，中比转数工况，"
                      f"选择离心泵 (IS 系列)",
            "alternative": PUMP_DOUBLE_SUCTION,
            "ns": ns,
        }

    # ns < 80: 径向离心泵适用范围
    # ── 4. 含气量判断 ──
    # 含气量 > 5% 时，离心泵容易气蚀，需选用耐气蚀泵型
    if gas_content > 5:
        if H <= 150:
            if Q_m3h <= 10:
                if viscosity_mm2s <= 37.4:
                    return {
                        "pump_type": PUMP_VORTEX,
                        "reason": f"含气量 {gas_content}% > 5%，小流量(Q={Q_m3h:.1f} m³/h ≤ 10)、"
                                  f"低黏度({viscosity_mm2s:.1f} ≤ 37.4 mm²/s)，"
                                  f"漩涡泵自吸性好，适合含气介质",
                        "alternative": PUMP_DIAPHRAGM,
                        "ns": ns,
                    }
                else:
                    return {
                        "pump_type": PUMP_POSITIVE_DISPLACEMENT,
                        "reason": f"含气量 {gas_content}% > 5%，黏度 {viscosity_mm2s:.1f} mm²/s > 37.4，"
                                  f"选择容积式泵（螺杆泵耐含气、耐高黏）",
                        "alternative": PUMP_SCREW,
                        "ns": ns,
                    }
            else:
                return {
                    "pump_type": PUMP_POSITIVE_DISPLACEMENT,
                    "reason": f"含气量 {gas_content}% > 5%，流量 Q={Q_m3h:.1f} m³/h > 10，"
                              f"选择容积式泵（隔膜泵/螺杆泵耐含气）",
                    "alternative": PUMP_DIAPHRAGM,
                    "ns": ns,
                }
        else:
            return {
                "pump_type": PUMP_POSITIVE_DISPLACEMENT,
                "reason": f"含气量 {gas_content}% > 5%，扬程 H={H:.1f} m > 150，"
                          f"高扬程含气工况选择容积式泵",
                "alternative": PUMP_GEAR,
                "ns": ns,
            }

    # ── 5. 含气量 ≤ 5%：离心泵 ──
    # 大流量 (>600 m³/h) → 双吸泵
    if Q_m3h >= 600:
        return {
            "pump_type": PUMP_DOUBLE_SUCTION,
            "reason": f"流量 Q={Q_m3h:.0f} m³/h ≥ 600，大流量工况推荐双吸离心泵 (S 系列)",
            "alternative": PUMP_CENTRIFUGAL,
            "ns": ns,
        }

    return {
        "pump_type": PUMP_CENTRIFUGAL,
        "reason": f"比转数 ns={ns:.1f}，含气量 {gas_content}% ≤ 5%，"
                  f"黏度 {viscosity_mm2s:.1f} mm²/s ≤ 650，"
                  f"常规工况选择离心泵 (IS 系列，效率高、维护简单)",
        "alternative": PUMP_MULTISTAGE,
        "ns": ns,
    }


# ============================================================
# 管路特性曲线
# ============================================================

def pipe_system_curve(Q: np.ndarray, pipe_diameter: float, pipe_length: float,
                      elevation_diff: float, viscosity: float = 0.001,
                      roughness: float = 0.0001, g: float = 9.81) -> np.ndarray:
    """计算管路特性曲线 H = H_static + K * Q²
    
    Args:
        Q: 流量数组 (m³/s)
        pipe_diameter: 管径 (m)
        pipe_length: 管长 (m)
        elevation_diff: 高程差 (m)
        viscosity: 粘度 (Pa·s)
        roughness: 管壁粗糙度 (m)
        g: 重力加速度
        
    Returns:
        扬程数组 (m)
    """
    rho = 1000  # 假设水
    H_static = elevation_diff
    
    # 计算摩擦系数 K
    # 流速
    A = math.pi * (pipe_diameter / 2)**2
    # 假设一个典型流速计算 Re
    Q_typical = np.mean(Q) if len(Q) > 0 else 0.01
    v_typical = Q_typical / A if A > 0 else 1.0
    Re = rho * v_typical * pipe_diameter / viscosity
    
    # Colebrook-White 方程（简化：取 Swamee-Jain 近似）
    if Re < 2000:
        f = 64 / Re
    else:
        # Swamee-Jain
        f = 0.25 / (math.log10(roughness / (3.7 * pipe_diameter) + 5.74 / (Re**0.9)))**2
    
    # K = f * L / (2 * g * D * A²)
    K = f * pipe_length / (2 * g * pipe_diameter * A**2)
    
    # 管路特性
    H = H_static + K * Q**2
    
    return H


# ============================================================
# 泵-系统曲线工作点匹配
# ============================================================

def match_pump_system_curve(
    pump_Q: np.ndarray,
    pump_H: np.ndarray,
    H_static: float,
    K_system: float,
    n_points: int = 100,
) -> dict:
    """泵特性曲线与系统曲线工作点匹配
    
    求泵 H-Q 曲线与系统曲线 H_sys = H_static + K·Q² 的交点。
    
    Args:
        pump_Q: 泵流量数组 (m³/s)
        pump_H: 泵扬程数组 (m)，与 pump_Q 对应
        H_static: 系统静扬程 (m)
        K_system: 系统阻力系数 (s²/m⁵)
        n_points: 插值点数
        
    Returns:
        {
            "operating_Q_m3_s": 工作点流量 (m³/s),
            "operating_Q_m3_h": 工作点流量 (m³/h),
            "operating_H_m": 工作点扬程 (m),
            "pump_curve": 泵曲线数据点列表 [{"Q_m3_s": ..., "H_m": ...}, ...],
            "system_curve": 系统曲线数据点列表,
            "efficiency_at_operating": 工作点效率估计,
        }
    """
    # 拟合泵曲线为二次多项式 H = a - b·Q - c·Q²
    coeffs = np.polyfit(pump_Q, pump_H, 2)
    a, b, c = coeffs[2], coeffs[1], coeffs[0]
    
    # 系统曲线 H_sys = H_static + K·Q²
    # 求交点: a + b·Q + c·Q² = H_static + K·Q²
    # (c - K)·Q² + b·Q + (a - H_static) = 0
    
    A_eq = c - K_system
    B_eq = b
    C_eq = a - H_static
    
    if abs(A_eq) > 1e-12:
        discriminant = B_eq ** 2 - 4 * A_eq * C_eq
        if discriminant < 0:
            return {
                "error": "无交点：泵曲线与系统曲线不相交",
                "operating_Q_m3_s": None,
                "operating_H_m": None,
            }
        Q1 = (-B_eq + math.sqrt(discriminant)) / (2 * A_eq)
        Q2 = (-B_eq - math.sqrt(discriminant)) / (2 * A_eq)
        
        # 选择正流量解
        Q_op = Q1 if Q1 > 0 else Q2
    else:
        # 线性情况
        Q_op = -C_eq / B_eq if abs(B_eq) > 1e-12 else 0
    
    if Q_op < 0:
        return {
            "error": "无有效工作点（流量为负）",
            "operating_Q_m3_s": None,
            "operating_H_m": None,
        }
    
    H_op = H_static + K_system * Q_op ** 2
    
    # 生成曲线数据
    Q_range = np.linspace(0, max(pump_Q) * 1.1, n_points)
    H_pump_curve = np.polyval(coeffs, Q_range)
    H_sys_curve = H_static + K_system * Q_range ** 2
    
    pump_curve_data = [
        {"Q_m3_s": float(q), "H_m": float(h)}
        for q, h in zip(Q_range, H_pump_curve) if h > 0
    ]
    system_curve_data = [
        {"Q_m3_s": float(q), "H_m": float(h)}
        for q, h in zip(Q_range, H_sys_curve)
    ]
    
    # 工作点效率估计
    eff_op = estimate_pump_efficiency(Q_op, H_op)
    
    return {
        "operating_Q_m3_s": Q_op,
        "operating_Q_m3_h": Q_op * 3600,
        "operating_H_m": H_op,
        "pump_curve": pump_curve_data,
        "system_curve": system_curve_data,
        "efficiency_at_operating": eff_op,
        "H_static": H_static,
        "K_system": K_system,
    }


def calculate_pump_system_curve(
    Q_points: np.ndarray,
    H_static: float,
    pipe_diameter: float,
    pipe_length: float,
    elevation_diff: float = 0.0,
    viscosity: float = 0.001,
    roughness: float = 0.0001,
    fitting_K: float = 0.0,
) -> dict:
    """泵系统曲线计算 H_sys(Q) = H_static + K·Q²，多流量点
    
    Args:
        Q_points: 流量点数组 (m³/s)
        H_static: 静扬程 (m)
        pipe_diameter: 管道内径 (m)
        pipe_length: 管道总长 (m)
        elevation_diff: 高程差 (m)
        viscosity: 流体黏度 (Pa·s)
        roughness: 管道粗糙度 (m)
        fitting_K: 管件总阻力系数
        
    Returns:
        {
            "H_static_total": 总静扬程 (m),
            "K": 系统阻力系数 (s²/m⁵),
            "curve_points": [{"Q_m3_s": ..., "H_m": ...}, ...],
            "friction_drop_at_mean_Q": 平均流量下摩擦压降 (m),
        }
    """
    g = 9.81
    rho = 1000
    
    # 总静扬程
    H_static_total = H_static + elevation_diff
    
    # 管道截面积
    A = math.pi * (pipe_diameter / 2) ** 2
    
    # 计算典型雷诺数和摩擦因子
    Q_mean = float(np.mean(Q_points))
    v_mean = Q_mean / A if A > 0 else 0
    Re = rho * v_mean * pipe_diameter / viscosity if viscosity > 0 else 1e5
    
    if Re < 2300:
        f = 64 / Re
    else:
        rel_rough = roughness / pipe_diameter if pipe_diameter > 0 else 0
        f = 0.25 / (math.log10(rel_rough / 3.7 + 5.74 / Re ** 0.9)) ** 2
    
    # 系统阻力系数
    # H_friction = f·L/D·v²/(2g) + ΣK·v²/(2g)
    # v = Q/A → v² = Q²/A²
    # H_friction = (f·L/D + ΣK) / (2g·A²) · Q²
    K_pipe = f * pipe_length / pipe_diameter if pipe_diameter > 0 else 0
    K_total = (K_pipe + fitting_K) / (2 * g * A ** 2)
    
    # 计算各流量点扬程
    curve_points = []
    for Q in Q_points:
        H = H_static_total + K_total * Q ** 2
        curve_points.append({"Q_m3_s": float(Q), "H_m": float(H)})
    
    # 平均流量下摩擦压降
    H_fric_mean = K_total * Q_mean ** 2
    
    return {
        "H_static_total": H_static_total,
        "K": K_total,
        "curve_points": curve_points,
        "friction_drop_at_mean_Q": H_fric_mean,
        "f": f,
        "Re": Re,
    }


# ============================================================
# 泵选型完整流程
# ============================================================

def pump_selection_full(Q_required: float, H_required: float,
                        fluid: str = _REQUIRED, T: float = _REQUIRED,
                        P_suction: float = 101325, suction_head: float = 0,
                        n: float = 2900,
                        metering_required: bool = False,
                        gas_content: float = 0.0,
                        pump_type: str = None) -> dict:
    """完整泵选型流程

    Args:
        Q_required: 所需流量 (m³/s)
        H_required: 所需扬程 (m)
        fluid: 流体类型（必填，如 'water'）
        T: 温度 (K)（必填）
        P_suction: 泵入口压力 (Pa)，默认大气压 101325 Pa
        suction_head: 吸入高度 (m)
        n: 转速 (rpm)
        metering_required: 是否有计量要求
        gas_content: 含气量百分比 (例如 6 表示 6%)
        pump_type: 泵类型（None 则自动推荐）

    Returns:
        选型结果字典
    """
    # 检查必填参数
    missing = check_required_params(locals(), {
        "fluid": "流体名称（如 'water', 'methanol'），用于查询物性",
        "T": "流体温度 (K)，用于蒸气压修正",
    })
    if missing:
        return missing
    # 1. 获取流体物性
    props = _LIQUID_PROPERTIES.get(fluid)
    if props is None:
        return {
            "needs_human": True,
            "human_prompt": f"流体 '{fluid}' 不在内置物性数据库中。请提供以下物性参数：\n"
                           f"  • rho：密度 (kg/m³)\n"
                           f"  • mu：动力粘度 (Pa·s)\n"
                           f"  • P_vapor：蒸气压 (Pa，@20°C)\n"
                           f"当前支持的流体：{list(_LIQUID_PROPERTIES.keys())}",
            "error": "unknown_fluid",
        }
    rho = props["rho"]
    mu = props["mu"]
    P_vapor = props["vapor_pressure_20C"] * math.exp(-5000 * (1/T - 1/293.15))  # 简化蒸气压温度修正

    # 2. 建议泵类型 (允许外部覆盖)
    if pump_type:
        pump_suggestion = {"pump_type": pump_type, "reason": "外部指定泵类型 (反思修正)", "alternative": None, "ns": None}
    else:
        pump_suggestion = suggest_pump_type(
            Q_required, H_required, mu,
            rho=rho,
            metering_required=metering_required,
            gas_content=gas_content,
        )
    pump_type = pump_suggestion["pump_type"]
    
    # 3. 效率估算
    efficiency = estimate_pump_efficiency(Q_required, H_required, n)
    
    # 4. 功率计算
    power = pump_power_full(Q_required, H_required, rho, efficiency)
    
    # 5. NPSH 计算
    npsha_value = npsha(P_suction, P_vapor, rho, suction_head, friction_loss=0.5)
    npshr_value = npshr_estimate(Q_required, n)
    npsh_result = npsh_check(npsha_value, npshr_value)
    
    # 6. 泵型号匹配 (所有泵类型)
    pump_model = select_pump_by_type(pump_type, Q_required, H_required)
    
    result = {
        "pump_type": pump_type,
        "pump_suggestion": pump_suggestion,
        "efficiency": efficiency,
        "power": power,
        "npsh": {
            "npsha": npsha_value,
            "npshr": npshr_value,
            "check": npsh_result,
        },
        "fluid_properties": {
            "rho": rho,
            "mu": mu,
            "P_vapor": P_vapor,
        },
    }
    
    if pump_model:
        result["pump_model"] = pump_model
        result["model_suitable"] = pump_model.get("suitable", False)
    
    return result


# ============================================================
# 设计校验
# ============================================================

def validate_pump_selection(selection: dict) -> dict:
    """校验泵选型
    
    Args:
        selection: 选型结果字典
        
    Returns:
        {
            "valid": 是否有效,
            "checks": [检查项],
            "warnings": [警告]
        }
    """
    checks = []
    warnings = []
    
    # 1. NPSH 校核
    npsh_check_result = selection.get("npsh", {}).get("check", {})
    if npsh_check_result.get("safe", False):
        checks.append({"name": "npsh", "pass": True, "value": npsh_check_result.get("margin_actual", 0), "limit": ">0.5m"})
    else:
        checks.append({"name": "npsh", "pass": False, "value": npsh_check_result.get("margin_actual", 0), "limit": ">0.5m"})
        warnings.append("NPSH 余量不足，可能发生气蚀")
    
    # 2. 效率检查
    eff = selection.get("efficiency", 0)
    if eff >= 0.6:
        checks.append({"name": "efficiency", "pass": True, "value": eff, "limit": ">0.6"})
    else:
        checks.append({"name": "efficiency", "pass": False, "value": eff, "limit": ">0.6"})
        warnings.append(f"泵效率 {eff:.2f} 较低")
    
    # 3. 电机功率安全系数
    power = selection.get("power", {})
    P_motor = power.get("P_motor", 0)
    P_shaft = power.get("P_shaft", 0)
    if P_motor > 0 and P_shaft > 0:
        actual_sf = P_motor / P_shaft
        if actual_sf >= 1.1:
            checks.append({"name": "service_factor", "pass": True, "value": actual_sf, "limit": ">=1.1"})
        else:
            checks.append({"name": "service_factor", "pass": False, "value": actual_sf, "limit": ">=1.1"})
            warnings.append("电机功率安全系数不足")
    
    # 4. 泵型号匹配
    if "pump_model" in selection:
        suitable = selection.get("model_suitable", False)
        if suitable:
            checks.append({"name": "model_match", "pass": True, "value": True, "limit": "True"})
        else:
            checks.append({"name": "model_match", "pass": False, "value": False, "limit": "True"})
            warnings.append("泵型号超出标准系列范围")
    
    n_pass = sum(1 for c in checks if c["pass"])
    n_total = len(checks)
    pass_rate = n_pass / n_total if n_total > 0 else 0
    
    valid = pass_rate >= 0.8
    
    return {
        "valid": valid,
        "pass_rate": pass_rate,
        "n_pass": n_pass,
        "n_total": n_total,
        "checks": checks,
        "warnings": warnings,
    }
