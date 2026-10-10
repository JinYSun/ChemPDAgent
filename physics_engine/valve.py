"""
阀门计算工具包 — 控制阀、安全阀
包含：
- 液体控制阀 Cv/Kv 计算（IEC 60534-2-1 标准，含空化判别，精准区分公英制系数）
- 气体安全阀泄放面积（API 520 Part I 标准，深度修正公制常数与压力单位）
"""
from __future__ import annotations
import math
from typing import Optional, Dict


# ============================================================
# 液体控制阀 Cv 计算（IEC 60534-2-1）
# ============================================================

def size_control_valve_liquid(
    Q: float,
    rho: float = None,
    P1: float = None,
    P2: float = None,
    Pv: float = None,
    Pc: Optional[float] = None,
    valve_type: str = "globe",
    xT: Optional[float] = None,
    FL: Optional[float] = None,
    NPSHa: Optional[float] = None,
    # ── 可选：用户给定流体信息，自动查询物性 ──
    fluid: str = None,
    T: float = None,
    P: float = 101325.0,
) -> dict:
    """液体控制阀流量系数计算（IEC 60534-2-1 标准）
    
    含空化/闪蒸判别，严格区分并计算公制 Kv 与美制 Cv。
    
    物性参数支持用户给定优先：
    - 若提供 rho/Pv，直接使用
    - 若未提供但提供 fluid/T，自动从数据库查询
    """
    from physics_engine.thermo_helper import get_fluid_density, get_fluid_vapor_pressure
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    Q = safe_float(Q, name='Q')
    rho = safe_float(rho, name='rho')
    P1 = safe_float(P1, name='P1')
    P2 = safe_float(P2, name='P2')
    Pv = safe_float(Pv, name='Pv')
    Pc = safe_float(Pc, name='Pc')
    xT = safe_float(xT, name='xT')
    FL = safe_float(FL, name='FL')
    NPSHa = safe_float(NPSHa, name='NPSHa')
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
        except Exception as e:
            return {"error": f"密度查询失败：{e}。请手动提供 rho 参数。"}
    else:
        return {"error": "必须提供 rho 或 (fluid + T) 参数"}
    
    # ── 饱和蒸气压：用户给定值优先 ──
    if Pv is not None:
        Pv_val = float(Pv)
        conversion_notes.append(f"[用户给定] 饱和蒸气压 Pv = {Pv_val} Pa")
    elif fluid is not None and T is not None:
        try:
            Pv_val = get_fluid_vapor_pressure(fluid, T)
            conversion_notes.append(f"[自动查询] {fluid} 在 {T}K 的饱和蒸气压 = {Pv_val:.1f} Pa")
        except Exception as e:
            return {"error": f"饱和蒸气压查询失败：{e}。请手动提供 Pv 参数。"}
    else:
        return {"error": "必须提供 Pv 或 (fluid + T) 参数"}
    
    # ── 压力参数检查 ──
    if P1 is None or P2 is None:
        return {"error": "必须提供 P1（阀前压力）和 P2（阀后压力）"}
    P1_val = float(P1)
    P2_val = float(P2)
    
    g = 9.81
    
    # 比重 (以 4°C 水为基准)
    SG = rho_val / 1000.0
    
    delta_P = P1_val - P2_val
    if delta_P <= 0:
        raise ValueError("阀前压力 P1 必须大于阀后压力 P2")
    
    # 阀门压力恢复系数默认值 (IEC 60534)
    FL_defaults = {
        "globe": 0.9,
        "ball": 0.65,
        "butterfly": 0.6,
        "plug": 0.8,
        "gate": 0.9,
    }
    
    xT_defaults = {
        "globe": 0.72,
        "ball": 0.65,
        "butterfly": 0.35,
        "plug": 0.65,
        "gate": 0.9,
    }
    
    FL_val = FL if FL is not None else FL_defaults.get(valve_type, 0.9)
    xT_val = xT if xT is not None else xT_defaults.get(valve_type, 0.72)
    
    # ── IEC 60534-2-1 标准三级判别 ──
    
    # 液体临界压力比因子 FF (IEC 60534-2-1)
    if Pc is not None and Pc > 0:
        FF = 0.96 - 0.28 * math.sqrt(Pv_val / Pc)
    else:
        FF = 0.96
    
    # ============================================================
    # 判别 1：闪蒸 (Flashing)
    # 判据：P2 < Pv → 阀后液体持续汽化
    # ============================================================
    flashing = P2_val < Pv_val
    
    # ============================================================
    # 判别 2：阻塞流 (Choked Flow) — IEC 60534 含 FF 修正
    # ΔP_choked = FL² × (P1 - FF × Pv)
    # ============================================================
    delta_P_choked = FL_val ** 2 * (P1_val - FF * Pv_val)
    choked = delta_P > delta_P_choked
    
    # 阻塞流裕度比 = ΔP / ΔP_choked，越接近 1.0 越危险
    choked_margin = delta_P / delta_P_choked if delta_P_choked > 0 else float('inf')
    
    # ============================================================
    # 判别 3：空化 (Cavitation) — FL-based 严格判据
    # 缩流断面压力 P_vc = P1 - ΔP / FL²
    # 空化条件：P_vc < Pv
    # ============================================================
    P_vc = P1_val - delta_P / (FL_val ** 2)
    cavitation = P_vc < Pv_val
    
    # 空化安全裕度 σ_FL = (P_vc - Pv) / ΔP
    # σ_FL > 0 安全, σ_FL < 0 空化, σ_FL < 0.05 临界
    sigma_FL = (P_vc - Pv_val) / delta_P if delta_P > 0 else float('inf')
    
    # 辅助参考指标（非主判据）
    sigma_actual = (P2_val - Pv_val) / delta_P if delta_P > 0 else float('inf')
    
    # ============================================================
    # 综合风险等级判定
    # ============================================================
    if flashing:
        cavitation_risk = "severe"       # 闪蒸 — 最严重
        flow_condition = "flashing"
    elif cavitation:
        if choked:
            cavitation_risk = "choked_cavitation"  # 阻塞流 + 空化
        else:
            cavitation_risk = "high"     # 空化但未阻塞
        flow_condition = "choked" if choked else "non_choked"
    elif choked:
        cavitation_risk = "choked_no_cavitation"  # 阻塞流但缩流面未汽化（罕见）
        flow_condition = "choked"
    elif sigma_FL < 0.05:
        cavitation_risk = "marginal"     # 临界状态，工况波动可能引发空化
        flow_condition = "non_choked"
    else:
        cavitation_risk = "low"          # 安全
        flow_condition = "non_choked"
    
    # 修正：公制下 Q (m3/h) 与 P (bar) 算出的流量系数定义为 Kv
    delta_P_eff = delta_P_choked if choked else delta_P
    delta_P_bar = delta_P_eff / 1e5  # Pa → bar
    Kv = Q * math.sqrt(SG / delta_P_bar)
    
    # 修正：通过标准系数 1.156 严格由 Kv 换算得到美制流量系数 Cv
    Cv = 1.156 * Kv
    
    # 安全裕度系数 (通常最大设计流量预留 30% 裕量)
    Kv_required = Kv * 1.3
    Cv_required = Cv * 1.3
    
    # NPSHa 汽蚀余量校核
    if NPSHa is not None:
        NPSH_required = (P1_val - Pv_val) / (rho_val * g) - delta_P_eff / (rho_val * g)
        if NPSHa < NPSH_required:
            cavitation_risk = "severe"
    
    return {
        "Kv": Kv,
        "Kv_required": Kv_required,
        "Cv": Cv,
        "Cv_required": Cv_required,
        "valve_type": valve_type,
        "flow_condition": flow_condition,
        "delta_P_Pa": delta_P,
        "delta_P_choked_Pa": delta_P_choked,
        "choked_margin": round(choked_margin, 4),
        "P_vc_Pa": round(P_vc, 1),
        "sigma_FL": round(sigma_FL, 4),
        "sigma_actual": round(sigma_actual, 4),
        "cavitation_risk": cavitation_risk,
        "flashing": flashing,
        "choked": choked,
        "cavitation": cavitation,
        "recommended_Cv": Cv_required,
        "recommended_Kv": Kv_required,
        "conversion_notes": conversion_notes,
        "FL": FL_val,
        "xT": xT_val,
        "FF": round(FF, 4),
    }


# ============================================================
# 气体安全阀泄放面积（API 520 Part I）
# ============================================================

def size_safety_valve_gas(
    W: float,
    P_set: float,
    P_back: float,
    T: float,
    MW: float,
    k: float = 1.4,
    Z: float = 1.0,
    C: Optional[float] = None,
    Kd: float = 0.975,
    Kb: float = 1.0,
    Kc: float = 1.0,
) -> dict:
    """气体/蒸汽安全阀泄放面积计算（API 520 Part I 标准公制法）
    
    A = W · √(T·Z) / (C · K · P1 · Kb · Kc · √(MW))
    
    Args:
        W: 泄放量 (kg/h)
        P_set: 设定压力 (Pa, 表压)
        P_back: 背压 (Pa, 表压)
        T: 泄放温度 (K)
        MW: 分子量 (kg/kmol，即 g/mol)
        k: 绝热指数 Cp/Cv
        Z: 气体压缩因子
        C: 气体常数，不指定则根据 k 自动采用 API 520 公制标准公式计算
        Kd: 额定排放系数（气体默认 0.975）
        Kb: 背压修正系数
        Kc: 破裂盘修正系数（无破裂盘为 1.0）
        
    Returns:
        安全阀泄放设计计算结果字典
    """
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    W = safe_float(W, name='W')
    P_set = safe_float(P_set, name='P_set')
    P_back = safe_float(P_back, name='P_back')
    T = safe_float(T, name='T')
    MW = safe_float(MW, name='MW')
    k = safe_float(k, default=1.4, name='k')
    Z = safe_float(Z, default=1.0, name='Z')
    C = safe_float(C, name='C')
    Kd = safe_float(Kd, default=0.975, name='Kd')
    Kb = safe_float(Kb, default=1.0, name='Kb')
    Kc = safe_float(Kc, default=1.0, name='Kc')
    
    P_atm = 101325.0  # 标准大气压 Pa
    
    # 1. 泄放压力（绝对压力，Pa）= 设定压力 × 1.1（10%超压）+ 大气压
    P1 = P_set * 1.1 + P_atm
    
    # 修正：API 520 SI 公制面积公式要求上游泄放压力 P1 的物理单位必须为 kPa（绝对压力）
    P1_kPa = P1 / 1000.0
    
    # 2. 修正：气体常数 C 值计算采用标准的 API 520 公制前置常数 0.03948，
    # 彻底杜绝原先误用英制前置常数 520 导致的安全阀泄放面积被缩小一千万倍的灾难性 Bug
    if C is None:
        C = 0.03948 * math.sqrt(k * (2.0 / (k + 1.0)) ** ((k + 1.0) / (k - 1.0)))
    
    # 3. 临界流动状态判别
    P_critical_ratio = (2.0 / (k + 1.0)) ** (k / (k - 1.0))
    P_critical = P_critical_ratio * P1
    
    # 确定实际泄放背压（绝对压力）
    P2_abs = P_back + P_atm
    
    if P2_abs <= P_critical:
        flow_type = "critical"
        Kb_actual = Kb
    else:
        flow_type = "subcritical"
        # 亚临界流背压修正系数估算（基于临界流面积向亚临界外推）
        r = P2_abs / P1
        Kb_actual = math.sqrt(1.0 - r) / math.sqrt(1.0 - P_critical_ratio)
    
    # 4. 泄放面积计算（API 520 标准公制公式）
    # A (mm²) = W * sqrt(T * Z) / (C * Kd * Kb * Kc * P1_kPa * sqrt(MW))
    A_mm2 = W * math.sqrt(T * Z) / (C * Kd * Kb_actual * Kc * P1_kPa * math.sqrt(MW))
    A_in2 = A_mm2 / 645.16  # mm² → in²
    
    # 5. 标准孔口代号匹配（API 526）
    standard_orifices = [
        ("D", 0.110),   # in²
        ("E", 0.196),
        ("F", 0.307),
        ("G", 0.503),
        ("H", 0.785),
        ("J", 1.287),
        ("K", 1.838),
        ("L", 2.853),
        ("M", 3.600),
        ("N", 4.340),
        ("P", 6.380),
        ("Q", 11.050),
        ("R", 16.000),
        ("T", 26.000),
    ]
    
    recommended = None
    for code, area_in2 in standard_orifices:
        if area_in2 >= A_in2:
            recommended = code
            break
    
    if recommended is None:
        recommended = "custom"
    
    return {
        "required_area_mm2": A_mm2,
        "required_area_in2": A_in2,
        "P1_Pa": P1,
        "P1_kPa": P1_kPa,
        "flow_type": flow_type,
        "C": C,
        "Kb_actual": Kb_actual,
        "recommended_orifice": recommended,
        "critical_pressure_ratio": P_critical_ratio,
        "critical_pressure_Pa": P_critical,
    }