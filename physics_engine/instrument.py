"""
仪表计算工具包 — 孔板流量计
包含：
- 孔板流量计计算（基于 ISO 5167-2 标准，严密处理雷诺数与流出系数迭代）
- 孔板孔径精确工艺设计（由目标流量反推孔径，基于连续寻根）
"""
from __future__ import annotations
import math
from typing import Optional, Dict


# ============================================================
# 孔板流量计（ISO 5167-2）
# ============================================================

def calculate_orifice_flow(
    dP: float,
    D_pipe: float,
    d_orifice: float = None,
    rho: float = None,
    mu: float = None,
    beta: Optional[float] = None,
    discharge_coefficient: Optional[float] = None,
    expansibility_factor: Optional[float] = None,
    tap_type: str = "corner",
    P1_absolute: Optional[float] = None,
    kappa: float = 1.4,
    # ── 可选：用户给定流体信息，自动查询物性 ──
    fluid: str = None,
    T: float = None,
    P: float = 101325.0,
    phase: str = "liquid",
    z: list = None,
) -> dict:
    """孔板流量计计算（ISO 5167 标准）
    
    由差压求流量（带雷诺数-流出系数自洽迭代）。
    物性参数支持用户给定优先，缺失时自动查询。
    """
    from physics_engine.thermo_helper import get_fluid_density, get_fluid_viscosity
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    dP = safe_float(dP, name='dP')
    D_pipe = safe_float(D_pipe, name='D_pipe')
    d_orifice = safe_float(d_orifice, name='d_orifice')
    rho = safe_float(rho, name='rho')
    mu = safe_float(mu, name='mu')
    beta = safe_float(beta, name='beta')
    discharge_coefficient = safe_float(discharge_coefficient, name='discharge_coefficient')
    expansibility_factor = safe_float(expansibility_factor, name='expansibility_factor')
    P1_absolute = safe_float(P1_absolute, name='P1_absolute')
    kappa = safe_float(kappa, default=1.4, name='kappa')
    T = safe_float(T, name='T')
    P = safe_float(P, default=101325.0, name='P')
    
    conversion_notes = []
    
    # ── 密度：用户给定值优先 ──
    if rho is not None:
        rho_val = float(rho)
        conversion_notes.append(f"[用户给定] 密度 rho = {rho_val} kg/m³")
    elif fluid is not None and T is not None:
        try:
            rho_val = get_fluid_density(fluid, T, P, phase=phase, z=z)
            conversion_notes.append(f"[自动查询] {fluid} 密度 = {rho_val:.2f} kg/m³")
        except Exception as e:
            return {"error": f"密度查询失败：{e}。请手动提供 rho。"}
    else:
        return {"error": "必须提供 rho 或 (fluid + T)"}
    
    # ── 粘度：用户给定值优先 ──
    if mu is not None:
        mu_val = float(mu)
        conversion_notes.append(f"[用户给定] 粘度 mu = {mu_val} Pa·s")
    elif fluid is not None and T is not None:
        try:
            mu_val = get_fluid_viscosity(fluid, T, P, phase=phase, z=z)
            conversion_notes.append(f"[自动查询] {fluid} 粘度 = {mu_val:.6f} Pa·s")
        except Exception as e:
            return {"error": f"粘度查询失败：{e}。请手动提供 mu。"}
    else:
        return {"error": "必须提供 mu 或 (fluid + T)"}
    
    # ── 孔径：给定 d_orifice 或用 beta 计算 ──
    if d_orifice is None and beta is not None:
        d_orifice = beta * D_pipe
    elif d_orifice is None:
        return {"error": "必须提供 d_orifice 或 beta"}
    if dP <= 0:
        raise ValueError("差压必须为正数")
    
    # 直径比
    if beta is None:
        beta = d_orifice / D_pipe
    
    if beta <= 0.1 or beta >= 0.75:
        # ISO 5167 标准推荐的 beta 范围通常在 0.1 ~ 0.75 之间
        pass
    if beta <= 0 or beta >= 1:
        raise ValueError("直径比 β 必须在 (0, 1) 范围内")
    
    # 渐进速度系数
    E = 1.0 / math.sqrt(1 - beta ** 4)  
    A_pipe = math.pi * (D_pipe / 2) ** 2
    
    # 1. 自动计算可膨胀性系数 ε
    if expansibility_factor is None:
        if P1_absolute is not None and P1_absolute > 0 and kappa > 0:
            # 气体/蒸汽可膨胀性系数 ISO 5167-2 标准公式
            epsilon = 1.0 - (0.41 + 0.35 * beta ** 4) * (dP / (kappa * P1_absolute))
            epsilon = max(0.66, min(1.0, epsilon))  # 标准物理下限安全保护
        else:
            # 默认液体
            epsilon = 1.0
    else:
        epsilon = expansibility_factor
        
    # 2. 迭代计算流出系数 C 与管道雷诺数 Re_D 的自洽解
    if discharge_coefficient is None:
        C = 0.6  # 工业初估值
        Re_D = 1e6  # 初始设定一个高雷诺数以启动迭代
        
        # 迭代收敛（通常 5 次内即可完全收敛）
        for _ in range(15):
            C_old = C
            C = _iso5167_discharge_coefficient(beta, D_pipe, tap_type, Re_D)
            
            # 计算当前的临时流量
            q_m_temp = C * E * epsilon * (math.pi / 4) * d_orifice ** 2 * math.sqrt(2 * rho_val * dP)
            v_temp = (q_m_temp / rho_val) / A_pipe if rho_val > 0 else 0.0
            
            # 更新雷诺数
            Re_D = rho_val * v_temp * D_pipe / mu_val if mu_val > 0 else 1e6
            Re_D = max(1e-5, Re_D)  # 防止除零或负数
            
            if abs(C - C_old) < 1e-6:
                break
    else:
        C = discharge_coefficient
        
    # 3. 计算最终物理流量
    q_m = C * E * epsilon * (math.pi / 4) * d_orifice ** 2 * math.sqrt(2 * rho_val * dP)
    q_v = q_m / rho_val
    v = q_v / A_pipe
    Re_D = rho_val * v * D_pipe / mu_val if mu_val > 0 else 0.0
    
    # 4. 精确计算永久压力损失 (ISO 5167-2 标准公式)
    sqrt_term = math.sqrt(1 - (beta ** 4) * (1 - C ** 2))
    denom_loss = sqrt_term + C * (beta ** 2)
    if denom_loss > 0:
        dP_loss = dP * (sqrt_term - C * (beta ** 2)) / denom_loss
    else:
        dP_loss = (1 - beta ** 2) * dP

    return {
        "mass_flow_kg_s": q_m,
        "volume_flow_m3_s": q_v,
        "beta": beta,
        "C": C,
        "epsilon": epsilon,
        "Re_D": Re_D,
        "velocity_m_s": v,
        "permanent_pressure_loss_Pa": dP_loss,
        "conversion_notes": conversion_notes,
    }


def _iso5167_discharge_coefficient(beta: float, D_pipe: float, tap_type: str = "corner", Re_D: float = 1e6) -> float:
    """ISO 5167-2 流出系数（Reader-Harris/Gallagher 公式，带雷诺数修正）
    """
    # 基础 Reader-Harris/Gallagher (R-H/G) 核心项
    C = 0.5961 + 0.0261 * beta ** 2 - 0.216 * beta ** 8
    
    # 考虑流道黏性边界层（雷诺数项）修正
    Re = max(100.0, Re_D)
    C += 0.000521 * ((1e6 * beta / Re) ** 0.7)
    
    # 针对不同取压方式的额外边界偏置修正
    if tap_type == "flange":
        # 法兰取压修正项 (简易工业近似)
        C += 0.0029 * beta ** 2.5 * (1e6 / Re) ** 0.3
    elif tap_type == "D_D/2":
        # D-D/2 径距取压修正
        C += 0.0015 * beta ** 2.5 * (10 ** 6 / Re) ** 0.3
    else:
        # 角接取压
        C += 0.0029 * beta ** 2.5 * (10 ** 6 / Re) ** 0.35 if Re > 0 else 0.0
        
    return C


def orifice_sizing(
    Q_target: float,
    D_pipe: float,
    dP_max: float,
    rho: float = None,
    mu: float = None,
    beta_min: float = 0.2,
    beta_max: float = 0.75,
    # ── 可选：用户给定流体信息，自动查询物性 ──
    fluid: str = None,
    T: float = None,
    P: float = 101325.0,
    phase: str = "liquid",
    z: list = None,
) -> dict:
    """孔板孔径设计（由目标流量反推孔径，基于连续二分法精确定位）
    
    物性参数支持用户给定优先，缺失时自动查询。
    
    Args:
        Q_target: 目标体积流量 (m³/s)
        D_pipe: 管道内径 (m)
        dP_max: 最大允许差压 (Pa)
        rho: 流体密度 (kg/m³)，给定后跳过查询
        mu: 流体动力黏度 (Pa·s)，给定后跳过查询
        beta_min: 最小直径比
        beta_max: 最大直径比
        fluid: 流体名称，用于自动查询
        T: 温度 (K)
        P: 压力 (Pa)，默认 101325
        phase: 相态: 'liquid'(默认)/'gas'
        z: 混合物摩尔分率列表
        
    Returns:
        孔径设计计算结果字典
    """
    from physics_engine.thermo_helper import get_fluid_density, get_fluid_viscosity
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    Q_target = safe_float(Q_target, name='Q_target')
    D_pipe = safe_float(D_pipe, name='D_pipe')
    dP_max = safe_float(dP_max, name='dP_max')
    rho = safe_float(rho, name='rho')
    mu = safe_float(mu, name='mu')
    beta_min = safe_float(beta_min, default=0.2, name='beta_min')
    beta_max = safe_float(beta_max, default=0.75, name='beta_max')
    T = safe_float(T, name='T')
    P = safe_float(P, default=101325.0, name='P')
    
    conversion_notes = []
    
    # ── 密度：用户给定值优先 ──
    if rho is not None:
        rho_val = float(rho)
        conversion_notes.append(f"[用户给定] 密度 rho = {rho_val} kg/m³")
    elif fluid is not None and T is not None:
        try:
            rho_val = get_fluid_density(fluid, T, P, phase=phase, z=z)
            conversion_notes.append(f"[自动查询] {fluid} 密度 = {rho_val:.2f} kg/m³")
        except Exception as e:
            return {"error": f"密度查询失败：{e}。请手动提供 rho。"}
    else:
        return {"error": "必须提供 rho 或 (fluid + T)"}
    
    # ── 粘度：用户给定值优先 ──
    if mu is not None:
        mu_val = float(mu)
        conversion_notes.append(f"[用户给定] 粘度 mu = {mu_val} Pa·s")
    elif fluid is not None and T is not None:
        try:
            mu_val = get_fluid_viscosity(fluid, T, P, phase=phase, z=z)
            conversion_notes.append(f"[自动查询] {fluid} 粘度 = {mu_val:.6f} Pa·s")
        except Exception as e:
            return {"error": f"粘度查询失败：{e}。请手动提供 mu。"}
    else:
        return {"error": "必须提供 mu 或 (fluid + T)"}
    
    q_m = rho_val * Q_target
    v = Q_target / (math.pi * (D_pipe / 2) ** 2)
    Re_D = rho_val * v * D_pipe / mu_val if mu_val > 0 else 1e6
    
    # 设定满量程设计时的目标设计差压等于系统允许的最大差压 dP_max
    dP_target = dP_max
    
    # 采用连续二分法替代步长为 1% 的粗糙离散循环，使 β 求解精度达到 1e-6 级别
    lo = beta_min
    hi = beta_max
    best_beta = 0.5 * (beta_min + beta_max)
    
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        C_temp = _iso5167_discharge_coefficient(mid, D_pipe, "corner", Re_D)
        E_temp = 1.0 / math.sqrt(1 - mid ** 4)
        d_temp = mid * D_pipe
        denom = C_temp * E_temp * (math.pi / 4) * d_temp ** 2
        
        if denom <= 0:
            lo = mid
            continue
            
        dP_calc = (q_m / denom) ** 2 / (2 * rho_val)
        
        if abs(dP_calc - dP_target) < 1e-2 or (hi - lo) < 1e-6:
            best_beta = mid
            break
            
        # β 增大 -> 孔径变大 -> 阻力与差压减小
        # 如果计算出的差压大于目标差压，说明阻力过大，孔径需要开得更大，应当提高 β 的下界
        if dP_calc > dP_target:
            lo = mid
        else:
            hi = mid
            best_beta = mid
            
    d_orifice = best_beta * D_pipe

    # 最终校核：在设计孔径下反算实际差压和流量
    C_final = _iso5167_discharge_coefficient(best_beta, D_pipe, "corner", Re_D)
    E_final = 1.0 / math.sqrt(1 - best_beta ** 4)
    A_orifice = math.pi / 4 * d_orifice ** 2
    dP_final = (q_m / (C_final * E_final * A_orifice)) ** 2 / (2 * rho_val)

    result = {
        "d_orifice_m": d_orifice,
        "d_orifice_mm": d_orifice * 1000,
        "beta": best_beta,
        "D_pipe_m": D_pipe,
        "Q_target_m3_s": Q_target,
        "rho_kg_m3": rho_val,
        "mu_Pa_s": mu_val,
        "Re_D": Re_D,
        "dP_max_Pa": dP_max,
        "dP_design_Pa": dP_final,
        "C_discharge": C_final,
        "E_velocity": E_final,
        "A_orifice_m2": A_orifice,
        "v_orifice_m_s": Q_target / A_orifice if A_orifice > 0 else 0,
        "valid": 0.2 <= best_beta <= 0.75,
        "conversion_notes": conversion_notes,
    }
    return result