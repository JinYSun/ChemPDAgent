"""
管道计算工具包 — 管道压降、管径选择、管件局部阻力
包含：
- 单相流体直管段压降（Darcy-Weisbach + Churchill 摩擦因子，修复极端雷诺数下溢 Bug）
- 气液两相流直管段压降（重构为严密、精准的 Friedel 经典状态方程关联式）
- 管件局部阻力压降（弯头、阀门、三通、进出口等，带空安全防护）
- 管径选择（符合化工工艺标准流速限制的管径精确设计法）
"""
from __future__ import annotations
import math
from typing import Optional, List, Dict, Tuple

from .thermo_helper import (
    get_fluid_density,
    get_fluid_viscosity,
)


# ============================================================
# 标准管径系列 (GB/T 17395, 外径 × 壁厚, mm)
# ============================================================

_STANDARD_PIPE_SIZES = [
    # (DN, OD_mm, schedule, ID_mm)
    (15,  21.3, "SCH40", 15.8),
    (20,  26.7, "SCH40", 20.7),
    (25,  33.4, "SCH40", 26.6),
    (32,  42.2, "SCH40", 35.0),
    (40,  48.3, "SCH40", 40.9),
    (50,  60.3, "SCH40", 52.5),
    (65,  73.0, "SCH40", 62.7),
    (80,  88.9, "SCH40", 77.9),
    (100, 114.3, "SCH40", 102.3),
    (125, 139.7, "SCH40", 128.2),
    (150, 168.3, "SCH40", 154.1),
    (200, 219.1, "SCH40", 202.7),
    (250, 273.0, "SCH40", 254.5),
    (300, 323.9, "SCH40", 303.2),
    (350, 355.6, "SCH40", 333.4),
    (400, 406.4, "SCH40", 381.0),
    (450, 457.0, "SCH40", 429.0),
    (500, 508.0, "SCH40", 477.8),
    (600, 610.0, "SCH40", 575.6),
]

# 推荐流速范围 (m/s)
_RECOMMENDED_VELOCITY = {
    "liquid": (0.5, 3.0),       # 一般液体
    "viscous_liquid": (0.1, 1.0), # 高黏度液体
    "gas_low": (5, 15),          # 低压气体
    "gas_medium": (15, 25),      # 中压气体
    "steam": (25, 40),           # 蒸汽
}

# ============================================================
# 管件局部阻力系数 K 值
# ============================================================

FITTING_K_VALUES = {
    # 阀门（全开）
    "gate_valve": 0.15,          # 闸阀
    "globe_valve": 10.0,         # 截止阀
    "ball_valve": 0.05,          # 球阀
    "butterfly_valve": 0.3,      # 蝶阀
    "check_valve": 2.0,          # 止回阀
    "plug_valve": 0.4,           # 旋塞阀
    "diaphragm_valve": 2.5,      # 隔膜阀
    # 弯头
    "elbow_90_long": 0.3,        # 90° 长半径弯头
    "elbow_90_short": 0.5,       # 90° 短半径弯头
    "elbow_45": 0.2,             # 45° 弯头
    "elbow_180_close": 0.9,      # 180° 回弯头
    "elbow_180_open": 0.3,       # 180° 大半径回弯头
    # 三通
    "tee_run": 0.4,              # 三通（直通）
    "tee_branch": 1.0,           # 三通（支路）
    # 进出口
    "entrance_sharp": 0.5,       # 锐边入口
    "entrance_rounded": 0.04,    # 圆滑入口
    "entrance_reentrant": 0.8,   # 内伸入口
    "exit_sharp": 1.0,           # 锐边出口
    "exit_rounded": 0.8,         # 圆滑出口
    # 异径管
    "sudden_enlargement": 1.0,   # 突然扩大
    "sudden_contraction": 0.5,   # 突然缩小
    "gradual_enlargement": 0.3,  # 渐扩管
    "gradual_contraction": 0.1,  # 渐缩管
}


# ============================================================
# 摩擦因子计算
# ============================================================

def churchill_friction(Re: float, roughness: float = 0.0, diameter: float = 1.0) -> float:
    """Churchill 摩擦因子（全流区，层流到湍流）
    
    Churchill 方程覆盖层流、过渡区和湍流，无需判断流态。
    
    Args:
        Re: 雷诺数
        roughness: 绝对粗糙度 (m)
        diameter: 管道内径 (m)
        
    Returns:
        Darcy 摩擦因子 f_D
    """
    if Re <= 0:
        raise ValueError("雷诺数必须为正数")
    
    rel_rough = roughness / diameter if diameter > 0 else 0.0
    
    # 修正：采用安全下限防护，防止极端高 Re 下光滑管 rel_rough 导致 log 内部下溢零引发 ValueError 崩溃
    log_arg = max(1e-300, (7.0 / Re) ** 0.9 + 0.27 * rel_rough)
    A = (-2.457 * math.log(log_arg)) ** 16
    B = (37530.0 / Re) ** 16
    
    f = 8.0 * ((8.0 / Re) ** 12 + 1.0 / (A + B) ** 1.5) ** (1.0 / 12.0)
    
    return f


def darcy_weisbach_friction(Re: float, rel_rough: float = 0.0) -> float:
    """Darcy-Weisbach 摩擦因子（Swamee-Jain 近似，湍流区）
    
    适用于 Re > 4000 的湍流区。
    
    Args:
        Re: 雷诺数
        rel_rough: 相对粗糙度 ε/D
        
    Returns:
        Darcy 摩擦因子 f_D
    """
    if Re < 2300:
        # 层流区
        return 64.0 / Re
    else:
        # 修正：添加安全防护，防止 Re 数极高时 5.74/Re**0.9 发生下溢，在 rel_rough=0 时引发 log10 域计算崩溃
        log_arg = max(1e-300, rel_rough / 3.7 + 5.74 / Re ** 0.9)
        f = 0.25 / (math.log10(log_arg)) ** 2
        return f


# ============================================================
# 单相流体直管段压降
# ============================================================

def calculate_pipe_pressure_drop(
    Q: float,
    L: float,
    D: float = None,
    rho: float = None,
    mu: float = None,
    roughness: float = 4.5e-5,
    elevation_change: float = 0.0,
    method: str = "churchill",
    # ── 局部阻力参数 ──
    K_fittings: float = 0.0,          # 局部阻力系数之和 (无量纲)
    n_elbow_90: int = 0,              # 90°弯头个数 (K=0.3/个)
    n_elbow_45: int = 0,              # 45°弯头个数 (K=0.15/个)
    n_gate_valve: int = 0,            # 闸阀个数 (K=0.15/个)
    n_globe_valve: int = 0,           # 截止阀个数 (K=6.0/个)
    n_check_valve: int = 0,           # 止回阀个数 (K=2.0/个)
    # ── 可选：用户给定流体信息，自动查询物性 ──
    fluid: str = None,
    T: float = None,
    P: float = 101325.0,
    phase: str = "liquid",
    z: list = None,
    # ── 可选：自动选径参数 ──
    fluid_type: str = "liquid",       # 用于 select_pipe_diameter 的流体类型
    max_velocity: float = None,       # 覆盖默认最大流速
    min_velocity: float = None,       # 覆盖默认最小流速
) -> dict:
    """单相流体直管段压降计算（Darcy-Weisbach 方程）
    
    物性参数支持用户给定优先：
    - 若提供 rho/mu，直接使用
    - 若未提供但提供 fluid/T，自动从数据库查询
    
    管径 D 支持自动选型：
    - 若提供 D，直接使用该值
    - 若未提供 D，自动调用 select_pipe_diameter 选择合适公称直径
    
    局部阻力支持两种方式：
    - K_fittings: 直接给定局部阻力系数之和
    - n_elbow_90/n_gate_valve等: 按管件个数自动累加 K 值
    - 两者可叠加使用
    """
    from physics_engine.thermo_helper import get_fluid_density, get_fluid_viscosity
    from physics_engine.common import safe_float, safe_int
    
    # ── 防御性类型转换 ──
    Q = safe_float(Q, name='Q')
    L = safe_float(L, name='L')
    D = safe_float(D, name='D')
    rho = safe_float(rho, name='rho')
    mu = safe_float(mu, name='mu')
    roughness = safe_float(roughness, default=4.5e-5, name='roughness')
    elevation_change = safe_float(elevation_change, default=0.0, name='elevation_change')
    K_fittings = safe_float(K_fittings, default=0.0, name='K_fittings')
    n_elbow_90 = safe_int(n_elbow_90, default=0, name='n_elbow_90')
    n_elbow_45 = safe_int(n_elbow_45, default=0, name='n_elbow_45')
    n_gate_valve = safe_int(n_gate_valve, default=0, name='n_gate_valve')
    n_globe_valve = safe_int(n_globe_valve, default=0, name='n_globe_valve')
    n_check_valve = safe_int(n_check_valve, default=0, name='n_check_valve')
    T = safe_float(T, name='T')
    P = safe_float(P, default=101325.0, name='P')
    max_velocity = safe_float(max_velocity, name='max_velocity')
    min_velocity = safe_float(min_velocity, name='min_velocity')
    
    conversion_notes = []
    
    # ── 密度：用户给定值优先 ──
    if rho is not None:
        rho_val = float(rho)
        conversion_notes.append(f"[用户给定] 密度 rho = {rho_val} kg/m³")
    elif fluid is not None and T is not None:
        try:
            rho_val = get_fluid_density(fluid, T, P, phase=phase, z=z)
            conversion_notes.append(f"[自动查询] {fluid} 在 {T}K, {P}Pa 的密度 = {rho_val:.2f} kg/m³")
        except Exception as e:
            return {"error": f"密度查询失败：{e}。请手动提供 rho 参数。"}
    else:
        return {"error": "必须提供 rho 或 (fluid + T) 参数"}
    
    # ── 粘度：用户给定值优先 ──
    if mu is not None:
        mu_val = float(mu)
        conversion_notes.append(f"[用户给定] 粘度 mu = {mu_val} Pa·s")
    elif fluid is not None and T is not None:
        try:
            mu_val = get_fluid_viscosity(fluid, T, P, phase=phase, z=z)
            conversion_notes.append(f"[自动查询] {fluid} 在 {T}K, {P}Pa 的粘度 = {mu_val:.6f} Pa·s")
        except Exception as e:
            return {"error": f"粘度查询失败：{e}。请手动提供 mu 参数。"}
    else:
        return {"error": "必须提供 mu 或 (fluid + T) 参数"}
    
    # ── 管径：用户给定优先，未给定时自动选型 ──
    pipe_sizing_result = None
    if D is not None:
        D_val = float(D)
        conversion_notes.append(f"[用户给定] 管道内径 D = {D_val} m")
    else:
        # 自动选型
        pipe_sizing_result = select_pipe_diameter(
            Q, fluid_type=fluid_type,
            max_velocity=max_velocity, min_velocity=min_velocity,
        )
        D_val = pipe_sizing_result["ID_m"]
        conversion_notes.append(
            f"[自动选型] DN{pipe_sizing_result['DN']} "
            f"(ID={D_val*1000:.1f} mm, v={pipe_sizing_result['velocity_m_s']:.2f} m/s)"
        )
    
    g = 9.81
    
    # 流速
    A = math.pi * (D_val / 2) ** 2
    v = Q / A
    
    # 雷诺数
    Re = rho_val * v * D_val / mu_val
    
    # 摩擦因子
    if method == "churchill":
        f = churchill_friction(Re, roughness, D_val)
    else:
        f = darcy_weisbach_friction(Re, roughness / D_val)
    
    # 摩擦压降
    dP_friction = f * (L / D_val) * (rho_val * v ** 2 / 2)
    
    # ── 局部阻力计算 ──
    # 标准管件阻力系数 (工程常用值)
    K_elbow_90 = 0.3       # 90°标准弯头
    K_elbow_45 = 0.15      # 45°标准弯头
    K_gate_valve = 0.15    # 全开闸阀
    K_globe_valve = 6.0    # 截止阀
    K_check_valve = 2.0    # 止回阀
    
    K_local = float(K_fittings)
    K_local += n_elbow_90 * K_elbow_90
    K_local += n_elbow_45 * K_elbow_45
    K_local += n_gate_valve * K_gate_valve
    K_local += n_globe_valve * K_globe_valve
    K_local += n_check_valve * K_check_valve
    
    dP_local = K_local * (rho_val * v ** 2 / 2)
    
    if K_local > 0:
        fitting_parts = []
        if n_elbow_90: fitting_parts.append(f"90°弯头×{n_elbow_90}")
        if n_elbow_45: fitting_parts.append(f"45°弯头×{n_elbow_45}")
        if n_gate_valve: fitting_parts.append(f"闸阀×{n_gate_valve}")
        if n_globe_valve: fitting_parts.append(f"截止阀×{n_globe_valve}")
        if n_check_valve: fitting_parts.append(f"止回阀×{n_check_valve}")
        if K_fittings: fitting_parts.append(f"其他K={K_fittings}")
        conversion_notes.append(f"[局部阻力] {', '.join(fitting_parts)} → ΣK={K_local:.2f}")
    
    # 高程压降
    dP_elevation = rho_val * g * elevation_change
    
    # 总压降
    dP_total = dP_friction + dP_elevation + dP_local
    
    # 流态判断
    if Re < 2300:
        regime = "laminar"
    elif Re < 4000:
        regime = "transitional"
    else:
        regime = "turbulent"
    
    result = {
        "pressure_drop_Pa": dP_total,
        "friction_drop_Pa": dP_friction,
        "local_resistance_drop_Pa": dP_local,
        "elevation_drop_Pa": dP_elevation,
        "velocity_m_s": v,
        "Re": Re,
        "friction_factor": f,
        "flow_regime": regime,
        "K_fittings_total": K_local,
        "pressure_drop_per_m": dP_total / L if L > 0 else 0.0,
        "conversion_notes": conversion_notes,
    }
    
    # 如果是自动选型，附加选型信息
    if pipe_sizing_result is not None:
        result["pipe_sizing"] = pipe_sizing_result
        result["D_m"] = D_val
        result["DN"] = pipe_sizing_result["DN"]
    
    return result


# ============================================================
# 气液两相流压降（重构：Friedel 1979 严密经典关联式）
# ============================================================

def calculate_two_phase_pressure_drop(
    Q_liquid: float,
    Q_gas: float,
    D: float,
    L: float,
    rho_liquid: float,
    rho_gas: float,
    mu_liquid: float,
    mu_gas: float,
    sigma: float = 0.072,
    roughness: float = 4.5e-5,
) -> dict:
    """气液两相流直管段压降（基于 Friedel 经典学术关联式）
    
    标准 Friedel 方法：ΔP_tp = Φ_lo² · ΔP_lo
    Phi_lo_sq = E_term + 3.24 * F_term * H_term / (Fr**0.0454 * We**0.035)
    
    Args:
        Q_liquid: 液相体积流量 (m³/s)
        Q_gas: 气相体积流量 (m³/s)
        D: 管道内径 (m)
        L: 管道长度 (m)
        rho_liquid: 液相密度 (kg/m³)
        rho_gas: 气相密度 (kg/m³)
        mu_liquid: 液相黏度 (Pa·s)
        mu_gas: 气相黏度 (Pa·s)
        sigma: 表面张力 (N/m)
        roughness: 绝对粗糙度 (m)
        
    Returns:
        两相压降计算结果字典
    """
    g = 9.81
    
    # 质量流量
    m_liquid = rho_liquid * Q_liquid
    m_gas = rho_gas * Q_gas
    m_total = m_liquid + m_gas
    
    if m_total <= 0:
        raise ValueError("总质量流量必须为正数")
    
    # 质量含气率（干度）
    x = m_gas / m_total
    
    # 总质量流速 G, kg/(m²·s)
    A = math.pi * (D / 2) ** 2
    G = m_total / A  
    
    # 1. 均相密度 rho_H
    if x <= 0:
        rho_H = rho_liquid
    elif x >= 1.0:
        rho_H = rho_gas
    else:
        rho_H = 1.0 / (x / rho_gas + (1.0 - x) / rho_liquid)
    
    # 纯液相雷诺数和全流量液相压降 ΔP_lo
    Re_lo = G * D / mu_liquid
    f_lo = churchill_friction(Re_lo, roughness, D)
    dP_lo = f_lo * (L / D) * (G ** 2 / (2 * rho_liquid))
    
    # 纯气相雷诺数和全流量气相压降 ΔP_go
    Re_go = G * D / mu_gas
    f_go = churchill_friction(Re_go, roughness, D)
    dP_go = f_go * (L / D) * (G ** 2 / (2 * rho_gas))
    
    # 2. 计算标准 Friedel (1979) 关联式的三大核心校正项
    if 0.0 < x < 1.0:
        # E_term 项：代表密度比与摩擦因子的基准比对
        E_term = (1 - x) ** 2 + x ** 2 * (rho_liquid * f_go) / (rho_gas * f_lo) if f_lo > 0 else 1.0
        
        # F_term 项：代表相间两相干度的交互作用
        F_term = (x ** 0.78) * ((1 - x) ** 0.224)
        
        # H_term 项：粘度比与密度比交互修正
        mu_ratio = mu_gas / mu_liquid if mu_liquid > 0 else 1.0
        if mu_ratio < 1.0:
            H_term = ((rho_liquid / rho_gas) ** 0.91) * (mu_ratio ** 0.19) * ((1.0 - mu_ratio) ** 0.7)
        else:
            H_term = ((rho_liquid / rho_gas) ** 0.91) * (mu_ratio ** 0.19)
            
        # 弗劳德数 (Froude) 和 韦伯数 (Weber)
        Fr = G ** 2 / (g * D * rho_H ** 2) if (rho_H > 0 and D > 0) else 1e-5
        We = G ** 2 * D / (rho_H * sigma) if (rho_H > 0 and sigma > 0) else 1e-5
        
        # 3. 计算 Φ_lo² (Friedel 两相倍增因子)
        denom_tp = (Fr ** 0.0454) * (We ** 0.035)
        if denom_tp > 0:
            Phi_lo_sq = E_term + 3.24 * F_term * H_term / denom_tp
        else:
            Phi_lo_sq = E_term
            
        Phi_lo_sq = max(Phi_lo_sq, 1.0)
    elif x <= 0.0:
        Phi_lo_sq = 1.0
    else:
        # 纯气相极限状态
        Phi_lo_sq = (rho_liquid / rho_gas) * (f_go / f_lo) if f_lo > 0 else 1.0
    
    # 绝热两相阻力压降
    dP_tp = Phi_lo_sq * dP_lo
    
    # 流型简易估计 (基于干度分档)
    if x < 0.1:
        flow_regime = "bubble"
    elif x < 0.3:
        flow_regime = "slug"
    elif x < 0.7:
        flow_regime = "annular"
    else:
        flow_regime = "mist"
    
    return {
        "pressure_drop_Pa": dP_tp,
        "pressure_drop_per_m": dP_tp / L if L > 0 else 0.0,
        "two_phase_multiplier": Phi_lo_sq,
        "liquid_only_drop_Pa": dP_lo,
        "gas_only_drop_Pa": dP_go,
        "quality": x,
        "flow_regime_estimate": flow_regime,
        "mass_velocity": G,
        "Re_liquid": Re_lo,
        "Re_gas": Re_go,
    }


# ============================================================
# 管件局部阻力压降
# ============================================================

def calculate_pipe_fitting_pressure_drop(
    Q: float,
    D: float,
    rho: float,
    fittings: List[Dict[str, any]],
) -> dict:
    """管件局部阻力压降计算
    """
    # 流速
    A = math.pi * (D / 2) ** 2
    v = Q / A
    
    # 速度头
    velocity_head = rho * v ** 2 / 2
    
    # 空安全防护
    if not fittings:
        return {
            "pressure_drop_Pa": 0.0,
            "velocity_m_s": v,
            "velocity_head_Pa": velocity_head,
            "total_K": 0.0,
            "fittings_breakdown": [],
        }
        
    total_K = 0.0
    breakdown = []
    
    for fitting in fittings:
        ftype = fitting.get("type", "")
        count = fitting.get("count", 1)
        
        if ftype not in FITTING_K_VALUES:
            raise ValueError(f"未知管件类型: {ftype}。可用类型: {list(FITTING_K_VALUES.keys())}")
        
        K = FITTING_K_VALUES[ftype]
        K_total = K * count
        dP = K_total * velocity_head
        
        total_K += K_total
        breakdown.append({
            "type": ftype,
            "K_single": K,
            "count": count,
            "K_total": K_total,
            "pressure_drop_Pa": dP,
        })
    
    dP_total = total_K * velocity_head
    
    return {
        "pressure_drop_Pa": dP_total,
        "velocity_m_s": v,
        "velocity_head_Pa": velocity_head,
        "total_K": total_K,
        "fittings_breakdown": breakdown,
    }


# ============================================================
# 管径选择
# ============================================================

def select_pipe_diameter(
    Q: float,
    fluid_type: str = "liquid",
    max_velocity: Optional[float] = None,
    min_velocity: Optional[float] = None,
) -> dict:
    """根据流量选择标准管径（符合工艺标准流速限制的精确比选设计）
    
    Args:
        Q: 体积流量 (m³/s)
        fluid_type: 流体类型
        max_velocity: 最大允许流速 (m/s)
        min_velocity: 最小允许流速 (m/s)
        
    Returns:
        管道匹配设计结果字典
    """
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    Q = safe_float(Q, name='Q')
    max_velocity = safe_float(max_velocity, name='max_velocity')
    min_velocity = safe_float(min_velocity, name='min_velocity')
    
    # 获取推荐流速限制范围
    if fluid_type in _RECOMMENDED_VELOCITY:
        v_min, v_max = _RECOMMENDED_VELOCITY[fluid_type]
    else:
        v_min, v_max = _RECOMMENDED_VELOCITY["liquid"]
    
    if min_velocity is not None:
        v_min = min_velocity
    if max_velocity is not None:
        v_max = max_velocity
    
    # 设定目标设计流速
    v_target = (v_min + v_max) / 2
    
    # 修正：重构管径匹配策略。
    # 1. 优先在标准管道库中筛选实际运行流速严格落在推荐流速区间 [v_min, v_max] 内部的所有公称口径候选者。
    valid_candidates = []
    for dn, od, sched, id_mm in _STANDARD_PIPE_SIZES:
        id_m = id_mm / 1000.0
        A = math.pi * (id_m / 2) ** 2
        v_actual = Q / A
        if v_min <= v_actual <= v_max:
            valid_candidates.append(((dn, od, sched, id_mm), abs(v_actual - v_target)))
            
    if valid_candidates:
        # 2. 如果存在多个候选规格，选择实际流速最接近目标中值流速的那一个，保证最佳流体阻力与设备投资性价比
        best_match = min(valid_candidates, key=lambda x: x[1])[0]
    else:
        # 3. 如果因流量过大/过小，导致管道数据库中无一能落在理想流速范围内，退而求其次寻找实际流速最贴近目标值的管道
        best_match = min(_STANDARD_PIPE_SIZES, key=lambda s: abs((Q / (math.pi * ((s[3] / 1000.0) / 2) ** 2)) - v_target))
        
    dn, od, sched, id_mm = best_match
    id_m = id_mm / 1000.0
    
    # 计算最终实际流速
    A = math.pi * (id_m / 2) ** 2
    v_actual = Q / A
    within_range = v_min <= v_actual <= v_max
    
    return {
        "DN": dn,
        "OD_mm": od,
        "ID_mm": id_mm,
        "ID_m": id_m,
        "schedule": sched,
        "velocity_m_s": v_actual,
        "recommended_velocity_range": (v_min, v_max),
        "within_range": within_range,
    }


# ============================================================
# 管道总压降（直管 + 管件）
# ============================================================

def calculate_total_pipe_pressure_drop(
    Q: float,
    D: float,
    L: float,
    rho: float,
    mu: float,
    roughness: float = 4.5e-5,
    elevation_change: float = 0.0,
    fittings: Optional[List[Dict[str, any]]] = None,
) -> dict:
    """管道总压降计算（直管段 + 管件局部阻力）
    """
    # 直管段压降
    straight_result = calculate_pipe_pressure_drop(
        Q=Q, L=L, D=D, rho=rho, mu=mu,
        roughness=roughness, elevation_change=elevation_change,
    )
    
    dP_straight = straight_result["friction_drop_Pa"]
    dP_elevation = straight_result["elevation_drop_Pa"]
    
    # 管件压降
    if fittings:
        fitting_result = calculate_pipe_fitting_pressure_drop(Q, D, rho, fittings)
        dP_fitting = fitting_result["pressure_drop_Pa"]
    else:
        fitting_result = None
        dP_fitting = 0.0
    
    dP_total = dP_straight + dP_fitting + dP_elevation
    
    return {
        "total_pressure_drop_Pa": dP_total,
        "straight_pipe_drop_Pa": dP_straight,
        "fitting_drop_Pa": dP_fitting,
        "elevation_drop_Pa": dP_elevation,
        "velocity_m_s": straight_result["velocity_m_s"],
        "Re": straight_result["Re"],
        "friction_factor": straight_result["friction_factor"],
        "flow_regime": straight_result["flow_regime"],
        "fitting_details": fitting_result,
    }