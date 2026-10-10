"""
闪蒸罐工具包 — 闪蒸计算、气液分离、容积设计
包含：
- 闪蒸计算（Rachford-Rice 方程）
- 等温闪蒸、绝热闪蒸
- 气液平衡（VLE）计算
- 闪蒸罐容积设计（修复容积叠加逻辑）
- 分离效率校核
- 压力降计算
- 液位控制与停留时间
"""
from __future__ import annotations
import numpy as np
import math
from typing import Optional, List, Dict, Tuple

# 导入 thermo 辅助模块
from .thermo_helper import (
    USE_THERMO,
    get_chemical_properties,
    get_mixture_k_values,
)


# ============================================================
# 闪蒸计算核心 (Rachford-Rice)
# ============================================================

def rachford_rice_equation(psi: float, z: np.ndarray, K: np.ndarray) -> float:
    """Rachford-Rice 方程
    
    Σ (z_i * (K_i - 1) / (1 + ψ * (K_i - 1))) = 0
    
    Args:
        psi: 气相分率 (V/F)
        z: 进料组成 (摩尔分数)
        K: 平衡常数 K_i = y_i / x_i
        
    Returns:
        f(ψ) 值
    """
    if len(z) != len(K):
        raise ValueError("组成数组长度必须与K值数组长度相同")
    
    result = 0.0
    for i in range(len(z)):
        denom = 1 + psi * (K[i] - 1)
        if abs(denom) < 1e-12:
            return float('inf') if z[i] * (K[i] - 1) > 0 else float('-inf')
        result += z[i] * (K[i] - 1) / denom
    return result


def solve_rachford_rice(z: np.ndarray, K: np.ndarray, tol: float = 1e-8, max_iter: int = 100) -> dict:
    """求解 Rachford-Rice 方程得到气相分率
    
    Args:
        z: 进料组成 (摩尔分数)
        K: 平衡常数数组
        tol: 容差
        max_iter: 最大迭代次数
        
    Returns:
        {"psi": 气相分率, "iterations": 迭代次数, "converged": 是否收敛}
    """
    # 修正：进行数学区间的边界校验以处理单相区。
    # 1. 如果 f(0) <= 0，代表处于或低于泡点压力/温度（纯液相），物理根 psi = 0.0
    f_0 = sum(z[i] * (K[i] - 1.0) for i in range(len(z)))
    if f_0 <= 0:
        return {"psi": 0.0, "iterations": 0, "converged": True}
        
    # 2. 如果 f(1) >= 0，代表处于或高于露点压力/温度（纯气相），物理根 psi = 1.0
    try:
        f_1 = sum(z[i] * (K[i] - 1.0) / K[i] for i in range(len(z)) if K[i] > 1e-12)
    except ZeroDivisionError:
        f_1 = -1e9
        
    if f_1 >= 0:
        return {"psi": 1.0, "iterations": 0, "converged": True}

    # 3. 确定两相区的数学边界
    psi_low = 0.0
    psi_high = 1.0
    
    for i in range(len(K)):
        if K[i] > 1.0:
            psi_low = max(psi_low, -1.0 / (K[i] - 1.0) + 1e-12)
        else:
            if abs(1.0 - K[i]) > 1e-12:
                psi_high = min(psi_high, 1.0 / (1.0 - K[i]) - 1e-12)
    
    # 二分法求解物理根
    psi = 0.5 * (psi_low + psi_high)
    
    for iteration in range(max_iter):
        f_low = rachford_rice_equation(psi_low, z, K)
        f_psi = rachford_rice_equation(psi, z, K)
        
        if abs(f_psi) < tol:
            return {"psi": psi, "iterations": iteration + 1, "converged": True}
        
        if f_low * f_psi < 0:
            psi_high = psi
        else:
            psi_low = psi
        
        psi = 0.5 * (psi_low + psi_high)
    
    return {"psi": psi, "iterations": max_iter, "converged": False}


def isothermal_flash(feed_flow: float, z: np.ndarray, K: np.ndarray, 
                     tol: float = 1e-8) -> dict:
    """等温闪蒸计算
    
    Args:
        feed_flow: 进料摩尔流量 (mol/s)
        z: 进料组成 (摩尔分数)
        K: 平衡常数 (在闪蒸温度下)
        tol: 求解容差
        
    Returns:
        等温闪蒸结果字典
    """
    result = solve_rachford_rice(z, K, tol)
    psi = result["psi"]
    converged = result["converged"]
    
    V = feed_flow * psi
    L = feed_flow * (1.0 - psi)
    
    n = len(z)
    y = np.zeros(n)
    x = np.zeros(n)
    
    # 修正：单相状态下的组成极值保护，避免除零错误
    if psi <= 1e-12:
        x = np.array(z, dtype=float)
        y = np.array(z, dtype=float) * K
        y_sum = np.sum(y)
        if y_sum > 0:
            y = y / y_sum
        return {"psi": 0.0, "V": 0.0, "L": feed_flow, "y": y, "x": x, "converged": converged}
    elif psi >= 1.0 - 1e-12:
        y = np.array(z, dtype=float)
        for i in range(n):
            x[i] = z[i] / K[i] if K[i] > 1e-12 else 0.0
        x_sum = np.sum(x)
        if x_sum > 0:
            x = x / x_sum
        return {"psi": 1.0, "V": feed_flow, "L": 0.0, "y": y, "x": x, "converged": converged}
        
    for i in range(n):
        denom = 1 + psi * (K[i] - 1)
        x[i] = z[i] / denom
        y[i] = K[i] * x[i]
    
    x = x / np.sum(x)
    y = y / np.sum(y)
    
    return {
        "psi": psi,
        "V": V,
        "L": L,
        "y": y,
        "x": x,
        "converged": converged,
    }


def adiabatic_flash(feed_flow: float, z: np.ndarray, T_feed: float, P_feed: float,
                   Cp_vapor: float, Cp_liquid: float, delta_H_vap: float,
                   K_func, tol: float = 1e-6, max_iter: int = 50) -> dict:
    """绝热闪蒸计算（带阻尼及边界约束的割线迭代法）
    
    Args:
        feed_flow: 进料流量 (mol/s)
        z: 进料组成
        T_feed: 进料温度 (K)
        P_feed: 进料压力 (Pa)
        Cp_vapor: 气相摩尔热容 (J/(mol·K))
        Cp_liquid: 液相摩尔热容 (J/(mol·K))
        delta_H_vap: 蒸发焓 (J/mol)
        K_func: K值函数 K = f(T, P, z)
        tol: 容差
        max_iter: 最大迭代次数
        
    Returns:
        闪蒸结果字典
    """
    T = T_feed
    T_prev = T_feed
    residual_prev = None
    
    for iteration in range(max_iter):
        K = K_func(T, P_feed, z)
        flash_result = isothermal_flash(feed_flow, z, K, tol)
        
        V = flash_result["V"]
        L = flash_result["L"]
        
        # 绝热闪蒸热衡算：Q_sensible + Q_latent = 0
        Q_sensible = V * Cp_vapor * (T - T_feed) + L * Cp_liquid * (T - T_feed)
        Q_latent = V * delta_H_vap
        residual = Q_sensible + Q_latent
        
        if abs(residual) < tol * delta_H_vap * feed_flow:
            flash_result["T_flash"] = T
            flash_result["residual"] = residual
            flash_result["iterations"] = iteration + 1
            return flash_result
        
        # 修正：通过割线法计算温度导数，放弃粗糙且不稳定的固定比例导数
        if residual_prev is not None and abs(T - T_prev) > 1e-5:
            d_residual_dT = (residual - residual_prev) / (T - T_prev)
            if abs(d_residual_dT) > 1e-3:
                dT = -residual / d_residual_dT
            else:
                dT = -residual / (feed_flow * Cp_liquid + 1e-5)
        else:
            dT = -residual / (feed_flow * Cp_liquid + 1e-5)
            
        T_prev = T
        residual_prev = residual
        
        # 限制单步温变步长并添加阻尼因子以增强数值收敛韧性
        dT = max(-50.0, min(50.0, dT))
        T += dT * 0.8  
        
        # 绝对温度下限约束保护
        if T < 10.0:
            T = 10.0
            
    flash_result["T_flash"] = T
    flash_result["iterations"] = max_iter
    return flash_result


# ============================================================
# 平衡常数计算
# ============================================================

def henry_law_k_values(P_sat: np.ndarray, P: float) -> np.ndarray:
    """亨利定律计算K值：K_i = P_sat_i / P"""
    return P_sat / P


def antoine_k_values(T: float, P: float, A: np.ndarray = None, B: np.ndarray = None, C: np.ndarray = None,
                        use_thermo: bool = False, components: List[str] = None) -> np.ndarray:
    """用安托万方程或 thermo 计算K值
    
    log10(P_sat) = A - B/(T_C + C)  (注意：传统 Antoine 关联式温度 T_C 须为摄氏度，A, B, C需对应 mmHg 单位)
    """
    # 修正：直接安全提取导入的标准统一物性接口，防止 Chemical 依赖本地静态库出错
    if use_thermo and USE_THERMO and components is not None:
        try:
            return np.array(get_mixture_k_values(components, T, P))
        except Exception as e:
            pass
    
    if A is None or B is None or C is None:
        raise ValueError("必须提供 A, B, C 或在 thermo_helper 畅通时开启 use_thermo")
    
    n = len(A)
    K = np.zeros(n)
    T_C = T - 273.15  # 摄氏度转换
    for i in range(n):
        log10_P_sat = A[i] - B[i] / (T_C + C[i])
        P_sat = 10 ** log10_P_sat * 133.322  # mmHg to Pa
        K[i] = P_sat / P
    return K


# ============================================================
# 闪蒸罐容积设计
# ============================================================

def flash_drum_volume(vapor_flow: float, liquid_flow: float, 
                      vapor_density: float, liquid_density: float,
                      vapor_holdup_time: float = 10.0,
                      liquid_holdup_time: float = 300.0,
                      vapor_space_factor: float = 0.3,
                      liquid_level_factor: float = 0.6) -> dict:
    """闪蒸罐容积设计（符合容器共腔的标准设计法）
    
    Args:
        vapor_flow: 气相体积流量 (m³/s)
        liquid_flow: 液相体积流量 (m³/s)
        vapor_density: 气相密度 (kg/m³)
        liquid_density: 液相密度 (kg/m³)
        vapor_holdup_time: 气相停留时间 (s)
        liquid_holdup_time: 液相停留时间 (s)
        vapor_space_factor: 气相所占最小体积比例
        liquid_level_factor: 液相所占最高体积比例
        
    Returns:
        闪蒸罐主尺寸及容积字典
    """
    # 实际所需的最少气相和液相容积
    V_vapor_req = vapor_flow * vapor_holdup_time
    V_liquid_req = liquid_flow * liquid_holdup_time
    
    # 修正：根据标准工程设计方法，闪蒸罐中气液在同一壳体内共处，总容积基于液相容积系数与气相空间容积协同校核
    V_total = V_liquid_req / liquid_level_factor
    V_vapor = V_total - V_liquid_req
    
    # 气相空间比例校核，若实际气相空间不足则以气相设计基准反向扩大总容积
    if V_vapor / V_total < vapor_space_factor:
        V_total = V_vapor_req / vapor_space_factor
        V_vapor = V_vapor_req
        V_liquid_req = V_total - V_vapor
    
    # 假设推荐长径比为 3:1
    aspect_ratio = 3.0
    D = (4 * V_total / (aspect_ratio * math.pi)) ** (1/3)
    L_total = aspect_ratio * D
    
    # 液相真实高度与气相高度分配
    L_liquid = V_liquid_req / (math.pi * D**2 / 4)
    L_vapor = L_total - L_liquid
    
    # 额外附加 1.0m 的气液重力脱开分离高度（Disengagement Height）
    disengagement_height = 1.0
    L_total += disengagement_height
    L_vapor += disengagement_height
    
    return {
        "V_total": V_total,
        "V_vapor": V_vapor,
        "V_liquid": V_liquid_req,
        "D": D,
        "L_liquid": L_liquid,
        "L_vapor": L_vapor,
        "L_total": L_total,
        "aspect_ratio": L_total / D,
        "disengagement_height": disengagement_height,
        "tau_vapor": vapor_holdup_time,       # 修正：补充输出键以与校验模块完全对齐
        "tau_liquid": liquid_holdup_time,
    }


def flash_drum_diameter(vapor_flow: float, vapor_density: float, 
                       max_vapor_velocity: float = 3.0) -> dict:
    """基于气速计算闪蒸罐直径（避免夹带）"""
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    vapor_flow = safe_float(vapor_flow, name='vapor_flow')
    vapor_density = safe_float(vapor_density, name='vapor_density')
    max_vapor_velocity = safe_float(max_vapor_velocity, default=3.0, name='max_vapor_velocity')
    
    A_required = vapor_flow / max_vapor_velocity
    D = math.sqrt(4 * A_required / math.pi)
    u_actual = vapor_flow / (math.pi * D**2 / 4)
    
    return {
        "D": D,
        "u_actual": u_actual,
        "u_max": max_vapor_velocity,
        "A_cross": math.pi * D**2 / 4,
    }


# ============================================================
# 雾沫夹带计算
# ============================================================

def entrainment_fraction(vapor_velocity: float,
                        liquid_viscosity: Optional[float] = None,
                        surface_tension: Optional[float] = None,
                        liquid_density: Optional[float] = None,
                        gas_density: Optional[float] = None,
                        separator_diameter: Optional[float] = None,
                        strict: bool = False) -> dict:
    """估算雾沫夹带风险

    主要方法:
    1. Souders-Brown 气速校核法（标准工程实践）:
       u_max = K × sqrt((rho_L - rho_G) / rho_G)
       其中 K = 0.107（标准 Watkins 系数，适用于无除沫器的分离器）

    2. 当提供表面张力时，补充计算临界液滴直径:
       d_critical = sqrt(12 × sigma / (rho_L × u²))
       用于评估可能被夹带的液滴尺寸

    Args:
        vapor_velocity: 空塔气速 (m/s)
        liquid_viscosity: 液相黏度 (Pa.s)，可选（备用信息）
        surface_tension: 气液表面张力 (N/m)，可选
        liquid_density: 液相密度 (kg/m3)
        gas_density: 气相密度 (kg/m3)
        separator_diameter: 分离器内径 (m)，可选（备用信息）
        strict: 是否启用严格模式（预留参数，当前版本与默认行为一致）

    Returns:
        dict，包含:
        - method: 所用方法名称
        - velocity_ratio: u / u_max（Souders-Brown 速度比）
        - u_max_m_s: 最大允许气速 (m/s)
        - entrainment_risk: 夹带风险等级 ("low"/"medium"/"high"/"critical")
        - entrainment_fraction: 估算夹带分率（仅当可计算时）
        - critical_droplet_diameter_um: 临界液滴直径 (μm)，若提供表面张力
        - note: 文字说明
    """
    # ---- 主方法: Souders-Brown 气速校核法（标准工程实践）----
    if gas_density is not None and liquid_density is not None and gas_density > 0 and liquid_density > gas_density:
        K_sb = 0.107  # 标准 Watkins/Souders-Brown 系数（无除沫器）
        u_max = K_sb * math.sqrt((liquid_density - gas_density) / gas_density)
        u_max = max(u_max, 0.3)  # 工程下限
        u_max = min(u_max, 6.0)  # 工程上限
        velocity_ratio = vapor_velocity / u_max

        # 夹带风险评估
        if velocity_ratio < 0.6:
            risk_level = "low"
            ent_est = 0.0
            note = f"u/u_max={velocity_ratio:.4f}<0.6, 气速远低于夹带限, 无夹带风险"
        elif velocity_ratio < 0.8:
            risk_level = "medium"
            ent_est = 0.001
            note = f"u/u_max={velocity_ratio:.4f}, 气速偏低, 夹带风险低"
        elif velocity_ratio < 1.0:
            risk_level = "high"
            ent_est = 0.005
            note = f"u/u_max={velocity_ratio:.4f}, 接近夹带限, 需关注"
        else:
            risk_level = "critical"
            ent_est = 0.05
            note = f"u/u_max={velocity_ratio:.4f}>=1.0, 超允许气速, 存在严重夹带风险!"

        result = {
            "method": "Souders-Brown",
            "velocity_ratio": round(velocity_ratio, 4),
            "u_max_m_s": round(u_max, 4),
            "entrainment_risk": risk_level,
            "entrainment_fraction": ent_est,
            "note": note,
        }

        # ---- 补充: 临界液滴直径计算（若提供表面张力）----
        if surface_tension is not None and surface_tension > 0:
            # 临界液滴直径：在此尺寸以下的液滴可能被气流夹带
            # d_critical = sqrt(12 * sigma / (rho_L * u^2))
            # 注：这是基于液滴破碎的简化估算，非精确关联式
            if vapor_velocity > 0.01:  # 避免除零
                d_crit_m = math.sqrt(12 * surface_tension / (liquid_density * vapor_velocity**2))
                d_crit_um = d_crit_m * 1e6  # 转换为微米
                result["critical_droplet_diameter_um"] = round(d_crit_um, 1)
                result["note"] += f"\n临界液滴直径 d_crit ≈ {d_crit_um:.1f} μm（小于此尺寸的液滴可能被夹带）"

        return result

    # ---- 信息不足 ----
    kinetic_energy = 0.5 * (gas_density or 1.0) * vapor_velocity**2
    return {
        "method": "insufficient_data",
        "kinetic_energy_J_m3": round(kinetic_energy, 2),
        "velocity_ratio": None,
        "u_max_m_s": None,
        "entrainment_risk": "unknown",
        "entrainment_fraction": None,
        "note": "缺少密度数据，无法进行 Souders-Brown 校核。请提供 gas_density 和 liquid_density。",
    }


# ============================================================
# 压降计算
# ============================================================

def flash_drum_pressure_drop(inlet_vel: float, outlet_vel: float,
                             density: float, n_turns: float = 1.5) -> float:
    """闪蒸罐进出口压降（简化）"""
    delta_P = 0.5 * density * (inlet_vel**2 + outlet_vel**2) + \
             n_turns * 0.5 * density * inlet_vel**2 * 0.3  
    return delta_P


# ============================================================
# 设计校验
# ============================================================

def validate_flash_drum_design(design: dict) -> dict:
    """校验闪蒸罐设计"""
    checks = []
    warnings = []
    
    D = design.get("D", 0)
    L_total = design.get("L_total", 0)
    aspect_ratio = L_total / D if D > 0 else 0
    
    # 1. 长径比检查
    if 2 <= aspect_ratio <= 5:
        checks.append({"name": "aspect_ratio", "pass": True, "value": aspect_ratio, "limit": "2-5"})
    else:
        checks.append({"name": "aspect_ratio", "pass": False, "value": aspect_ratio, "limit": "2-5"})
        warnings.append(f"长径比 {aspect_ratio:.2f} 超出推荐范围 2-5")
    
    # 2. 气相速度检查
    u = design.get("u_actual", 0)
    u_max = design.get("u_max", 3.0)
    if u < u_max:
        checks.append({"name": "vapor_velocity", "pass": True, "value": u, "limit": f"<{u_max}"})
    else:
        checks.append({"name": "vapor_velocity", "pass": False, "value": u, "limit": f"<{u_max}"})
        warnings.append(f"气速 {u:.2f} m/s 过高，可能导致夹带")
    
    # 3. 停留时间检查（修正：支持两组常用键读取）
    tau_vapor = design.get("tau_vapor", design.get("vapor_holdup_time", 0))
    if tau_vapor >= 5:
        checks.append({"name": "vapor_residence_time", "pass": True, "value": tau_vapor, "limit": ">5s"})
    else:
        checks.append({"name": "vapor_residence_time", "pass": False, "value": tau_vapor, "limit": ">5s"})
    
    tau_liquid = design.get("tau_liquid", design.get("liquid_holdup_time", 0))
    if tau_liquid >= 180:  
        checks.append({"name": "liquid_residence_time", "pass": True, "value": tau_liquid, "limit": ">180s"})
    else:
        checks.append({"name": "liquid_residence_time", "pass": False, "value": tau_liquid, "limit": ">180s"})
    
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


# ============================================================
# 标准尺寸系列
# ============================================================

_STANDARD_DRUM_DIAMETERS = [0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0]

def select_standard_diameter(D_calculated: float) -> float:
    """选择标准罐径"""
    for D_std in _STANDARD_DRUM_DIAMETERS:
        if D_std >= D_calculated:
            return D_std
    return _STANDARD_DRUM_DIAMETERS[-1]  


# ============================================================
# 气液分离器尺寸设计（Watkins 法 + 停留时间法）
# ============================================================

def size_gas_liquid_separator(
    Q_gas: float,
    Q_liquid: float,
    rho_gas: float,
    rho_liquid: float,
    mu_gas: float = 2e-5,
    mu_liquid: Optional[float] = None,
    method: str = "watkins",
    residence_time_liquid: float = 300.0,
    residence_time_gas: float = 10.0,
    max_gas_velocity: Optional[float] = None,
    droplet_size: Optional[float] = None,
    D_existing: Optional[float] = None,
) -> dict:
    """气液分离器尺寸设计/校核
    
    droplet_size / mu_liquid 不提供时跳过 Stokes 沉降效率计算，
    仅输出基于 Watkins 气速校核的几何尺寸。
    
    D_existing 提供时，自动校核已有直径的夹带风险（速度比 + 风险等级）。
    """
    
    # ── 防御性类型转换：确保所有数值参数为 float（防止 LLM 输出字符串） ──
    Q_gas = float(Q_gas)
    Q_liquid = float(Q_liquid)
    rho_gas = float(rho_gas)
    rho_liquid = float(rho_liquid)
    mu_gas = float(mu_gas)
    residence_time_liquid = float(residence_time_liquid)
    residence_time_gas = float(residence_time_gas)
    if mu_liquid is not None:
        mu_liquid = float(mu_liquid)
    if max_gas_velocity is not None:
        max_gas_velocity = float(max_gas_velocity)
    if droplet_size is not None:
        droplet_size = float(droplet_size)
    if D_existing is not None:
        D_existing = float(D_existing)
    
    # ── method 别名兼容 ──
    if method in ("souders_brown", "souders-brown"):
        method = "watkins"
    
    # 提前定义好容积
    V_liquid = Q_liquid * residence_time_liquid
    V_gas = Q_gas * residence_time_gas
    
    if method == "watkins" or method == "both":
        if max_gas_velocity is None:
            K = 0.107  
            u_max = K * math.sqrt((rho_liquid - rho_gas) / rho_gas)
            u_max = max(u_max, 0.3)  
            u_max = min(u_max, 6.0)  
        else:
            u_max = max_gas_velocity
        
        A_gas = Q_gas / u_max
        D_watkins = math.sqrt(4 * A_gas / math.pi)
        A_total = math.pi * (D_watkins / 2) ** 2
        L_liquid = V_liquid / A_total if A_total > 0 else 0
        L_gas = residence_time_gas * u_max
        L_total = L_liquid + L_gas + 0.5  
        
        D_result = D_watkins
        L_result = L_total
    
    elif method == "residence_time":
        V_total = V_gas + V_liquid
        if Q_gas > Q_liquid * 10:
            aspect = 3.0
        else:
            aspect = 1.5
        
        D_rt = (4 * V_total / (aspect * math.pi)) ** (1/3)
        L_rt = aspect * D_rt
        
        D_result = D_rt
        L_result = L_rt
    
    else:
        return {"error": f"不支持的方法: {method}"}
    
    D_std = select_standard_diameter(D_result)
    A_std = math.pi * (D_std / 2) ** 2
    
    u_actual = Q_gas / A_std if A_std > 0 else 0
    liquid_level = V_liquid / A_std if A_std > 0 else 0
    
    # 修正：仅当 droplet_size 已显式提供时才执行 Stokes 沉降计算
    settling_result = {}
    if droplet_size is not None and droplet_size > 0:
        settling_result = calculate_stokes_settling_efficiency(
            droplet_diameter=droplet_size,
            rho_liquid=rho_liquid,
            rho_gas=rho_gas,
            mu_gas=mu_gas,
            separator_length=L_result,
            gas_velocity=u_actual
        )
    
    result = {
        "D": D_std,
        "L": L_result,
        "V_total": A_std * L_result,
        "method": method,
        "gas_velocity": u_actual,
        "liquid_level": liquid_level,
        "aspect_ratio": L_result / D_std,
        "separation_efficiency": settling_result.get("efficiency", None),
        "settling_detail": settling_result if settling_result else "droplet_size not provided, Stokes settling skipped",
        "V_liquid": V_liquid,
        "V_gas": V_gas,
    }
    
    # ── 已有直径校核（D_existing 提供时）──
    if D_existing is not None and D_existing > 0:
        A_existing = math.pi * (D_existing / 2) ** 2
        u_existing = Q_gas / A_existing if A_existing > 0 else 0
        
        # Souders-Brown 允许气速
        K = 0.107
        u_max_sb = K * math.sqrt((rho_liquid - rho_gas) / rho_gas)
        u_max_sb = max(u_max_sb, 0.3)
        u_max_sb = min(u_max_sb, 6.0)
        
        velocity_ratio = u_existing / u_max_sb if u_max_sb > 0 else 999
        
        if velocity_ratio < 0.6:
            risk = "low"
            risk_note = "远低于夹带风险阈值，分离器尺寸安全"
        elif velocity_ratio < 0.8:
            risk = "medium"
            risk_note = "接近风险区域，建议关注"
        elif velocity_ratio < 1.0:
            risk = "high"
            risk_note = "接近液泛，存在夹带风险"
        else:
            risk = "critical"
            risk_note = "已超过最大允许气速，必然夹带"
        
        result["verification"] = {
            "D_existing": D_existing,
            "A_cross_m2": round(A_existing, 4),
            "u_actual_m_s": round(u_existing, 4),
            "u_max_Souders_Brown_m_s": round(u_max_sb, 4),
            "velocity_ratio": round(velocity_ratio, 4),
            "velocity_ratio_percent": f"{velocity_ratio * 100:.2f}%",
            "risk_level": risk,
            "risk_note": risk_note,
            "pass": velocity_ratio < 0.8,
        }
    
    return result


# ============================================================
# 重力沉降分离效率
# ============================================================

def calculate_stokes_settling_efficiency(
    droplet_diameter: float,
    rho_liquid: float,
    rho_gas: float,
    mu_gas: float,
    separator_length: float,
    gas_velocity: float,
    separator_height: float = None,
) -> dict:
    """重力沉降分离效率计算"""
    g = 9.81
    
    if separator_height is None:
        separator_height = separator_length / 3
    
    u_terminal = g * droplet_diameter ** 2 * (rho_liquid - rho_gas) / (18 * mu_gas)
    Re_p = rho_gas * u_terminal * droplet_diameter / mu_gas
    
    if Re_p > 1.0:
        u_t = u_terminal
        for _ in range(10):
            Re_p = rho_gas * u_t * droplet_diameter / mu_gas
            if Re_p < 0.01:
                break
            C_d = 18.5 / Re_p ** 0.6 if Re_p < 500 else 0.44
            u_t_new = math.sqrt(4 * g * droplet_diameter * (rho_liquid - rho_gas) / (3 * C_d * rho_gas))
            u_t = 0.5 * (u_t + u_t_new)
        u_terminal = u_t
    
    t_settling = separator_height / u_terminal if u_terminal > 0 else float('inf')
    t_residence = separator_length / gas_velocity if gas_velocity > 0 else float('inf')
    
    if t_settling > 0 and t_settling < float('inf'):
        efficiency = 1 - math.exp(-t_residence / t_settling)
    else:
        efficiency = 0.0
    
    efficiency = max(0, min(efficiency, 1.0))
    
    u_t_critical = separator_height / t_residence if t_residence < float('inf') else 0
    d_critical = math.sqrt(18 * mu_gas * u_t_critical / (g * (rho_liquid - rho_gas))) if u_t_critical > 0 else 0
    
    if Re_p < 0.1:
        flow_regime = "stokes"
    elif Re_p < 500:
        flow_regime = "transition"
    else:
        flow_regime = "newton"
    
    return {
        "terminal_velocity": u_terminal,
        "settling_time": t_settling,
        "residence_time": t_residence,
        "efficiency": efficiency,
        "critical_droplet_size": d_critical,
        "flow_regime": flow_regime,
        "Re_particle": Re_p,
    }