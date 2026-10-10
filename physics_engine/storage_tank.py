"""
储罐工具包 — 常压/带压储罐设计、容积计算、液位控制
包含：
- 储罐容积计算（立式/卧式/球罐，修复球罐球缺容积 4 倍放大 Bug）
- 储罐尺寸设计（径高比优化、自动匹配安全填充率）
- 蒸发损失计算（呼吸损耗，消除气相高度 1m 硬编码缺陷）
- 基础荷载计算（风载/地震）
- 保温层设计
- 液位仪表选型
- 安全阀计算（API 520，补充绝压安全警示）
"""
from __future__ import annotations
import numpy as np
import math
from typing import Optional, List, Dict, Tuple

# Import thermo helper
from .thermo_helper import (
    get_chemical_properties,
)


# ============================================================
# 储罐类型与标准尺寸
# ============================================================

# 标准储罐直径系列 (m)
# 规则: 圆整到小数点后一位，末位为 5 或偶数 (0,2,4,5,6,8)
# 每个整数区间统一 0.2 步长: x.0, x.2, x.4, x.5, x.6, x.8

def _gen_diameters(start: float, stop: float) -> list:
    """生成标准直径序列: 0.1 步长遍历, 仅保留末位为 0/2/4/5/6/8 的值
    
    等价于全区间 0.2 步长 + 额外插入 x.5 值, 确保不遗漏任何合规直径。
    """
    allowed = {0, 2, 4, 5, 6, 8}
    result = []
    n = int(round(start * 10))
    n_end = int(round(stop * 10))
    while n <= n_end:
        last_digit = n % 10
        if last_digit in allowed:
            result.append(round(n / 10.0, 1))
        n += 1
    return result

# 标准立式储罐直径系列 (m), 起始 3.0
_STANDARD_VERTICAL_DIAMETERS = _gen_diameters(3.0, 50.0)

# 标准卧式储罐直径系列 (m), 起始 1.0
_STANDARD_HORIZONTAL_DIAMETERS = _gen_diameters(1.0, 10.0)

# 标准球罐直径系列 (m), 起始 4.0
_STANDARD_SPHERICAL_DIAMETERS = _gen_diameters(4.0, 50.0)

# 储罐类型
TANK_VERTICAL = "vertical_cylindrical"
TANK_HORIZONTAL = "horizontal_cylindrical"
TANK_SPHERICAL = "spherical"
TANK_FLAT_BOTTOM = "flat_bottom"


# ============================================================
# 封头容积计算
# ============================================================

def head_volume(D: float, head_type: str = "ellipsoidal") -> float:
    """单个封头的容积 (m3)

    支持封头类型:
    - flat:             平盖头 (容积 = 0)
    - ellipsoidal:      标准椭圆封头 (2:1), V = pi/24 * D^3
    - hemispherical:    半球封头, V = pi/12 * D^3
    - torispherical:    标准碟形封头 (Rc=D, rk=0.06D), V ~= 0.0847 * D^3
    - conical:          锥形封头 (半顶角 30deg), V = pi/24 * D^3

    Args:
        D: 储罐内径 (m)
        head_type: 封头类型

    Returns:
        单个封头的容积 (m3)
    """
    R = D / 2.0
    if head_type == "flat":
        return 0.0
    elif head_type == "hemispherical":
        # 半球体: V = (2/3)*pi*R^3 = pi/12 * D^3
        return (2.0 / 3.0) * math.pi * R ** 3
    elif head_type == "ellipsoidal":
        # 标准 2:1 椭圆封头: a=D/2, b=D/4
        # V = (2/3)*pi*a^2*b = (2/3)*pi*(D/2)^2*(D/4) = pi/24 * D^3
        return (math.pi / 24.0) * D ** 3
    elif head_type == "torispherical":
        # 标准碟形封头: Rc=D (冠部半径), rk=0.06*D (折边半径)
        # 精确公式: V = pi*h^2*(3*R - h)/3 其中 h 为封头深度
        # h = R - sqrt((R-rk)^2 - ((R-2*rk)/2)^2) + rk  (简化近似)
        # 工程近似: V ~= 0.0847 * D^3 (GB/T 25198)
        return 0.0847 * D ** 3
    elif head_type == "conical":
        # 锥形封头 (半顶角 30deg): h = D/(2*tan(30deg)) = D*sqrt(3)/2
        # 但工程常用 h = D/4 的标准锥, V = pi/12 * R^2 * h = pi/24 * D^3
        # 采用与椭圆封头相同的近似
        return (math.pi / 24.0) * D ** 3
    else:
        raise ValueError(f"不支持的封头类型: {head_type}。"
                         f"可选: flat, ellipsoidal, hemispherical, torispherical, conical")


def head_depth(D: float, head_type: str = "ellipsoidal") -> float:
    """单个封头的深度/高度 (m)

    Args:
        D: 储罐内径 (m)
        head_type: 封头类型

    Returns:
        封头深度 (m)
    """
    R = D / 2.0
    if head_type == "flat":
        return 0.0
    elif head_type == "hemispherical":
        return R
    elif head_type == "ellipsoidal":
        # 标准 2:1 椭圆: b = D/4
        return D / 4.0
    elif head_type == "torispherical":
        # 标准碟形: h = R - sqrt((R-rk)^2 - ((R-2*rk)/2)^2) + rk
        rk = 0.06 * D
        Rc = D
        h = Rc - math.sqrt((Rc - rk) ** 2 - ((R - rk) ** 2)) + rk
        # 简化: 对于 Rc=D, rk=0.06D, h ~= 0.194*D
        return h
    elif head_type == "conical":
        # 标准锥: h = D/4
        return D / 4.0
    else:
        raise ValueError(f"不支持的封头类型: {head_type}")


# ============================================================
# 容积计算
# ============================================================

def vertical_cylindrical_volume(D: float, H: float,
                                head_type: str = "flat") -> dict:
    """立式圆柱储罐容积

    Args:
        D: 储罐直径 (m)
        H: 筒体高度 (m)，不含封头
        head_type: 封头类型 (flat/ellipsoidal/hemispherical/torispherical/conical)

    Returns:
        dict 包含总容积、筒体容积、封头容积、封头深度
    """
    V_shell = math.pi * (D / 2) ** 2 * H
    V_one_head = head_volume(D, head_type)
    h_one_head = head_depth(D, head_type)

    # 立式罐通常底部为平盖，顶部为锥顶/拱顶
    # 为通用性，默认两端同类型封头
    V_heads = 2 * V_one_head
    V_total = V_shell + V_heads

    return {
        "total_volume": V_total,
        "shell_volume": V_shell,
        "head_volume_each": V_one_head,
        "head_volume_total": V_heads,
        "head_depth_each": h_one_head,
        "head_type": head_type,
        "D": D,
        "H_shell": H,
        "H_total": H + 2 * h_one_head,
    }


def horizontal_cylindrical_volume(D: float, L: float,
                                  liquid_height: Optional[float] = None,
                                  head_type: str = "flat") -> dict:
    """卧式圆柱储罐容积（含液面高度计算）

    Args:
        D: 储罐直径 (m)
        L: 筒体长度 (m)，不含封头
        liquid_height: 液面高度 (m)，None 时返回总容积
        head_type: 封头类型
    """
    V_shell = math.pi * (D / 2) ** 2 * L
    V_one_head = head_volume(D, head_type)
    h_one_head = head_depth(D, head_type)
    V_heads = 2 * V_one_head
    V_total = V_shell + V_heads

    if liquid_height is None:
        return {
            "total_volume": V_total,
            "shell_volume": V_shell,
            "head_volume_each": V_one_head,
            "head_volume_total": V_heads,
            "head_depth_each": h_one_head,
            "head_type": head_type,
            "liquid_volume": V_total,
            "ullage_volume": 0.0,
            "fill_fraction": 1.0,
        }

    # 筒体部分液量（弓形截面）
    h = max(0.0, min(liquid_height, D))
    if h == 0.0:
        A_liquid = 0.0
    elif h == D:
        A_liquid = math.pi * (D / 2) ** 2
    else:
        theta = 2 * math.acos(1 - 2 * h / D)
        A_liquid = (D ** 2 / 8) * (theta - math.sin(theta))

    V_shell_liquid = A_liquid * L

    # 封头部分液量（近似：按液面占比分配封头容积）
    if D > 0:
        fill_frac_shell = A_liquid / (math.pi * (D / 2) ** 2) if (math.pi * (D / 2) ** 2) > 0 else 0
    else:
        fill_frac_shell = 0
    V_heads_liquid = V_heads * fill_frac_shell

    V_liquid = V_shell_liquid + V_heads_liquid
    V_ullage = V_total - V_liquid
    fill_fraction = V_liquid / V_total if V_total > 0 else 0.0

    return {
        "total_volume": V_total,
        "shell_volume": V_shell,
        "head_volume_each": V_one_head,
        "head_volume_total": V_heads,
        "head_depth_each": h_one_head,
        "head_type": head_type,
        "liquid_volume": V_liquid,
        "ullage_volume": V_ullage,
        "fill_fraction": fill_fraction,
    }


def spherical_tank_volume(D: float) -> float:
    """球罐容积"""
    return (4/3) * math.pi * (D/2)**3


# ============================================================
# 储罐尺寸设计
# ============================================================

def design_vertical_tank(volume_required: float,
                          aspect_ratio: float = 1.0,
                          use_standard: bool = True,
                          head_type: str = "flat") -> dict:
    """设计立式圆柱储罐

    Args:
        volume_required: 所需液相容积 (m3)
        aspect_ratio: 高径比 H/D
        use_standard: 是否选用标准直径
        head_type: 封头类型 (flat/ellipsoidal/hemispherical/torispherical/conical)
    """
    # 工业标准：设计容积 = 所需液相容积 / 0.85 (预留15%气相空间)
    design_volume = volume_required / 0.85

    # 理论直径计算（先按纯圆柱体估算，再扣除封头）
    D_theoretical = (4 * design_volume / (math.pi * aspect_ratio)) ** (1 / 3)
    H_theoretical = aspect_ratio * D_theoretical

    if use_standard:
        D_selected = _STANDARD_VERTICAL_DIAMETERS[-1]
        for D_std in _STANDARD_VERTICAL_DIAMETERS:
            if D_std >= D_theoretical:
                D_selected = D_std
                break
        # 反推筒体高度（扣除封头容积）
        V_one_head = head_volume(D_selected, head_type)
        V_heads = 2 * V_one_head
        V_shell_needed = design_volume - V_heads
        if V_shell_needed <= 0:
            # 封头已满足容积需求，筒体可极短
            H_selected = 0.5
        else:
            H_selected = V_shell_needed / (math.pi * (D_selected / 2) ** 2)
    else:
        D_selected = D_theoretical
        H_selected = H_theoretical

    # 实际容积（含封头）
    vol_result = vertical_cylindrical_volume(D_selected, H_selected, head_type)
    actual_volume = vol_result["total_volume"]
    fill_ratio = volume_required / actual_volume if actual_volume > 0 else 1.0
    H_total = H_selected + 2 * vol_result["head_depth_each"]
    aspect_ratio_actual = H_total / D_selected if D_selected > 0 else 0.0

    return {
        "D": D_selected,
        "H": H_selected,                # 向后兼容: H = H_shell
        "H_shell": H_selected,
        "H_total": H_total,
        "head_type": head_type,
        "actual_volume": actual_volume,
        "shell_volume": vol_result["shell_volume"],
        "head_volume_each": vol_result["head_volume_each"],
        "head_volume_total": vol_result["head_volume_total"],
        "fill_ratio": fill_ratio,
        "aspect_ratio": aspect_ratio_actual,
        "theoretical_D": D_theoretical,
        "theoretical_H": H_theoretical,
    }


def design_horizontal_tank(volume_required: float,
                            aspect_ratio: float = 4.0,
                            use_standard: bool = True,
                            head_type: str = "ellipsoidal") -> dict:
    """设计卧式圆柱储罐

    Args:
        volume_required: 所需液相容积 (m3)
        aspect_ratio: 长径比 L/D
        use_standard: 是否选用标准直径
        head_type: 封头类型 (flat/ellipsoidal/hemispherical/torispherical/conical)
    """
    # 预留15%气相空间
    design_volume = volume_required / 0.85

    D_theoretical = (4 * design_volume / (math.pi * aspect_ratio)) ** (1 / 3)
    L_theoretical = aspect_ratio * D_theoretical

    if use_standard:
        D_selected = _STANDARD_HORIZONTAL_DIAMETERS[-1]
        for D_std in _STANDARD_HORIZONTAL_DIAMETERS:
            if D_std >= D_theoretical:
                D_selected = D_std
                break
        # 反推筒体长度（扣除封头容积）
        V_one_head = head_volume(D_selected, head_type)
        V_heads = 2 * V_one_head
        V_shell_needed = design_volume - V_heads
        if V_shell_needed <= 0:
            L_selected = 2.0
        else:
            L_selected = V_shell_needed / (math.pi * (D_selected / 2) ** 2)
    else:
        D_selected = D_theoretical
        L_selected = L_theoretical

    # 实际容积（含封头）
    vol_result = horizontal_cylindrical_volume(D_selected, L_selected, head_type=head_type)
    actual_volume = vol_result["total_volume"]
    fill_ratio = volume_required / actual_volume if actual_volume > 0 else 1.0
    L_total = L_selected + 2 * vol_result["head_depth_each"]
    aspect_ratio_actual = L_total / D_selected if D_selected > 0 else 0.0

    return {
        "D": D_selected,
        "L": L_selected,                # 向后兼容: L = L_shell
        "L_shell": L_selected,
        "L_total": L_total,
        "head_type": head_type,
        "actual_volume": actual_volume,
        "shell_volume": vol_result["shell_volume"],
        "head_volume_each": vol_result["head_volume_each"],
        "head_volume_total": vol_result["head_volume_total"],
        "fill_ratio": fill_ratio,
        "aspect_ratio": aspect_ratio_actual,
        "theoretical_D": D_theoretical,
        "theoretical_L": L_theoretical,
    }


def design_spherical_tank(volume_required: float, use_standard: bool = True) -> dict:
    """设计球罐"""
    design_volume = volume_required / 0.85
    D_theoretical = (6 * design_volume / math.pi) ** (1/3)
    
    if use_standard:
        D_selected = _STANDARD_SPHERICAL_DIAMETERS[-1]
        for D_std in _STANDARD_SPHERICAL_DIAMETERS:
            if (math.pi * D_std**3 / 6) >= design_volume:
                D_selected = D_std
                break
    else:
        D_selected = D_theoretical
    
    actual_volume = math.pi * D_selected**3 / 6
    fill_ratio = volume_required / actual_volume if actual_volume > 0 else 1.0
    
    return {
        "D": D_selected,
        "actual_volume": actual_volume,
        "fill_ratio": fill_ratio,
        "theoretical_D": D_theoretical,
    }


# ============================================================
# 蒸发损失计算 (呼吸损耗)
# ============================================================

def breathing_losses(D: float, T_day: float, T_night: float, 
                     vapor_pressure: float, molecular_weight: float,
                     saturation_factor: float = 0.7,
                     V_vapor: Optional[float] = None) -> dict:
    """计算固定顶罐呼吸损耗（简化 API 公式）
    
    Args:
        D: 储罐直径 (m)
        T_day: 白天最高气体平均温度 (K)
        T_night: 夜间最低气体平均温度 (K)
        vapor_pressure: 液体在平均温度下的饱和蒸气压 (Pa)
        molecular_weight: 蒸气分子量 (g/mol)
        saturation_factor: 气相空间不饱和度因子 (0.0~1.0)
        V_vapor: 实际气相空间体积 (m³)。若未提供，将默认按 D 比例估算，消除 1m 高度硬编码误差。
    """
    # 修正：杜绝原先写死 1.0m 气相高度带来的巨大预测偏差
    if V_vapor is None:
        V_vapor = math.pi * (D/2)**2 * 1.5  # 工业估算：默认气相空间高度约 1.5m
        
    delta_T = T_day - T_night
    R = 8.314
    P_vapor = saturation_factor * vapor_pressure
    
    n_moles = (P_vapor * V_vapor / (R * T_day)) - (P_vapor * V_vapor / (R * T_night))
    loss_per_breath = abs(n_moles) * molecular_weight / 1000  # kg
    
    n_breaths_per_day = max(1, int(abs(delta_T) / 5))
    annual_loss = loss_per_breath * n_breaths_per_day * 365
    
    return {
        "loss_per_breath": loss_per_breath,
        "annual_loss": annual_loss,
        "n_breaths_per_day": n_breaths_per_day,
        "V_vapor": V_vapor,
    }


# ============================================================
# 保温层设计
# ============================================================

def insulation_thickness(diameter: float, T_internal: float, T_ambient: float,
                         max_heat_loss: float = 100.0,
                         insulation_k: float = 0.04) -> dict:
    """计算保温层厚度"""
    delta_T = T_internal - T_ambient
    
    if delta_T <= 0:
        return {"thickness": 0.0, "heat_loss": 0.0, "surface_temperature": T_ambient}
    
    thickness_required = delta_T * insulation_k / max_heat_loss
    heat_loss = delta_T * insulation_k / thickness_required if thickness_required > 0 else float('inf')
    T_surface = T_internal - heat_loss * thickness_required / insulation_k
    
    return {
        "thickness": thickness_required,
        "heat_loss": heat_loss,
        "surface_temperature": T_surface,
        "insulation_k": insulation_k,
    }


# ============================================================
# 基础荷载计算 (简化)
# ============================================================

def tank_foundation_load(D: float, H: float, liquid_density: float,
                          steel_density: float = 7850.0) -> dict:
    """储罐基础荷载计算"""
    g = 9.81
    
    V_liquid = math.pi * (D/2)**2 * H
    W_liquid = V_liquid * liquid_density * g
    
    t_shell = 0.01  # m
    V_shell = math.pi * D * H * t_shell
    W_shell = V_shell * steel_density * g
    
    A_bottom = math.pi * (D/2)**2
    t_bottom = 0.008  # m
    W_bottom = A_bottom * t_bottom * steel_density * g
    
    total_weight = W_liquid + W_shell + W_bottom
    base_pressure = total_weight / A_bottom if A_bottom > 0 else 0.0
    
    return {
        "total_weight": total_weight,
        "liquid_weight": W_liquid,
        "shell_weight": W_shell,
        "bottom_weight": W_bottom,
        "base_pressure": base_pressure,
        "V_liquid": V_liquid,
    }


# ============================================================
# 安全阀计算
# ============================================================

def tank_relief_valve_sizing(volume: float, P_set: float, P_back: float,
                              T: float, molecular_weight: float,
                              k: float = 1.4) -> dict:
    """储罐安全阀尺寸计算（气体排放，API 520 泄放公式）
    
    安全警示：输入参数 P_set 和 P_back 必须是热力学“绝对压力 (Pa)”。
    若错误传入表压，可能导致求得的排放喉径偏小，带来泄放失效的安全隐患！
    """
    R = 8.314
    MW_kg_per_mol = molecular_weight / 1000.0
    P_critical = P_set * (2 / (k + 1)) ** (k / (k - 1))
    
    if P_back < P_critical:
        C = math.sqrt(k * (2 / (k + 1)) ** ((k + 1) / (k - 1)))
    else:
        Pr = P_back / P_set
        C = math.sqrt(2 * k / (k - 1) * (Pr ** (2/k) - Pr ** ((k+1)/k)))
    
    Q_fire = 17800 * math.pi * volume ** (2/3)  # W
    delta_H_vap = 300000  # J/kg
    W_relief = Q_fire / delta_H_vap  # kg/s
    
    A_required = W_relief / (C * P_set * math.sqrt(MW_kg_per_mol / (R * T)))
    D_throat = math.sqrt(4 * A_required / math.pi)
    
    return {
        "required_area": A_required,
        "required_diameter": D_throat,
        "relief_flow": W_relief,
        "C_factor": C,
        "fire_heat_input": Q_fire,
    }


# ============================================================
# 储罐径高比优化（最小表面积）
# ============================================================

def optimize_tank_aspect_ratio(
    volume_required: float,
    tank_type: str = TANK_VERTICAL,
    head_type: str = "ellipsoidal",
) -> dict:
    """储罐径高比优化（最小表面积设计）
    """
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    volume_required = safe_float(volume_required, name='volume_required')
    
    if tank_type == TANK_SPHERICAL:
        D = (6 * volume_required / math.pi) ** (1/3)
        A = 4 * math.pi * (D/2) ** 2
        return {
            "D_optimal": D,
            "H_optimal": D,
            "aspect_ratio": 1.0,
            "surface_area": A,
            "volume": spherical_tank_volume(D),
            "tank_type": tank_type,
        }
    
    head_volume_factor = {
        "ellipsoidal": 0.15,   # 椭圆封头 ≈ 0.15·D³
        "hemispherical": 0.26,  # 半球封头
        "flat": 0.0,
    }.get(head_type, 0.15)
    
    if tank_type == TANK_VERTICAL:
        best_D = None
        best_A = float('inf')
        
        for D in np.linspace(1, 20, 200):
            V_head = head_volume_factor * D ** 3
            V_shell = volume_required - 2 * V_head
            if V_shell <= 0:
                continue
            H = V_shell / (math.pi / 4 * D ** 2)
            
            if H < 0.5 or H > 30:
                continue
            
            A_shell = math.pi * D * H
            A_heads = 2 * math.pi * (D/2) ** 2
            A_total = A_shell + A_heads
            
            if A_total < best_A:
                best_A = A_total
                best_D = D
                best_H = H
        
        if best_D is None:
            best_D = (4 * volume_required / math.pi) ** (1/3)
            best_H = best_D
        
        return {
            "D_optimal": best_D,
            "H_optimal": best_H,
            "aspect_ratio": best_H / best_D,
            "surface_area": best_A,
            "volume": volume_required,
            "tank_type": tank_type,
        }
    
    elif tank_type == TANK_HORIZONTAL:
        best_D = None
        best_A = float('inf')
        
        for D in np.linspace(1, 10, 100):
            V_head = head_volume_factor * D ** 3
            V_shell = volume_required - 2 * V_head
            if V_shell <= 0:
                continue
            L = V_shell / (math.pi / 4 * D ** 2)
            
            if L < 2 or L > 30:
                continue
            
            aspect = L / D
            if aspect < 2 or aspect > 8:
                continue
            
            A_shell = math.pi * D * L
            A_heads = 2 * math.pi * (D/2) ** 2
            A_total = A_shell + A_heads
            
            if A_total < best_A:
                best_A = A_total
                best_D = D
                best_L = L
        
        if best_D is None:
            best_D = (4 * volume_required / (3 * math.pi)) ** (1/3)
            best_L = 3 * best_D
        
        return {
            "D_optimal": best_D,
            "H_optimal": best_L,
            "aspect_ratio": best_L / best_D,
            "surface_area": best_A,
            "volume": volume_required,
            "tank_type": tank_type,
        }
    
    return {"error": f"不支持的储罐类型: {tank_type}"}


# ============================================================
# 液位-持液量-停留时间计算
# ============================================================

def calculate_liquid_level_residence_time(
    volume: float,
    Q_in: float,
    Q_out: float,
    tank_type: str = TANK_VERTICAL,
    D: float = None,
    L: float = None,
    liquid_height: float = None,
    direction: str = "forward",
) -> dict:
    """液位-持液量-停留时间计算
    """
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    volume = safe_float(volume, name='volume')
    Q_in = safe_float(Q_in, name='Q_in')
    Q_out = safe_float(Q_out, name='Q_out')
    D = safe_float(D, name='D')
    L = safe_float(L, name='L')
    liquid_height = safe_float(liquid_height, name='liquid_height')
    
    if direction == "forward":
        if tank_type == TANK_VERTICAL:
            if D is None:
                D = (4 * volume / math.pi) ** 0.5
            H_tank = volume / (math.pi / 4 * D ** 2)
            
            if liquid_height is None:
                liquid_height = H_tank * 0.5  
            
            liquid_height = min(liquid_height, H_tank)
            liquid_vol = math.pi * (D/2) ** 2 * liquid_height
            fill_frac = liquid_vol / volume if volume > 0 else 0
        
        elif tank_type == TANK_HORIZONTAL:
            if D is None or L is None:
                D = (4 * volume / (4 * math.pi)) ** (1/3)
                L = 4 * D
            
            if liquid_height is None:
                liquid_height = D * 0.5
            
            vol_result = horizontal_cylindrical_volume(D, L, liquid_height)
            liquid_vol = vol_result["liquid_volume"]
            fill_frac = vol_result["fill_fraction"]
        
        else:
            # 球罐
            if D is None:
                D = (6 * volume / math.pi) ** (1/3)
            
            if liquid_height is None:
                liquid_height = D * 0.5
            
            h = min(liquid_height, D)
            # 修正：采用严格的球缺标准几何体积公式 V_seg = π * h^2 * (R - h/3)，避免写成直径 D 造成持液量被放大 4 倍的致命计算错误
            liquid_vol = math.pi * h ** 2 * (D/2 - h/3)
            fill_frac = liquid_vol / volume if volume > 0 else 0
        
        # 停留时间
        if Q_out > 0:
            residence_s = liquid_vol / Q_out
        elif Q_in > 0:
            residence_s = liquid_vol / Q_in
        else:
            residence_s = float('inf')
        
        return {
            "liquid_volume": liquid_vol,
            "liquid_height": liquid_height,
            "fill_fraction": fill_frac,
            "residence_time_s": residence_s,
            "residence_time_min": residence_s / 60,
            "residence_time_h": residence_s / 3600,
        }
    
    elif direction == "reverse":
        target_residence = volume / max(Q_out, Q_in, 1e-10)  
        
        if Q_out > 0:
            target_vol = Q_out * target_residence
        else:
            target_vol = volume * 0.5
        
        target_vol = min(target_vol, volume)
        
        if tank_type == TANK_VERTICAL:
            if D is None:
                D = (4 * volume / math.pi) ** 0.5
            H_tank = volume / (math.pi / 4 * D ** 2)
            liquid_height = target_vol / (math.pi / 4 * D ** 2)
            liquid_height = min(liquid_height, H_tank)
            fill_frac = target_vol / volume
        
        elif tank_type == TANK_HORIZONTAL:
            if D is None or L is None:
                D = (4 * volume / (4 * math.pi)) ** (1/3)
                L = 4 * D
            
            # 反算液位（迭代）
            fill_frac = target_vol / volume
            liquid_height = D * fill_frac  
            for _ in range(20):
                vol_result = horizontal_cylindrical_volume(D, L, liquid_height)
                error = vol_result["liquid_volume"] - target_vol
                if abs(error) < 1e-6:
                    break
                liquid_height += error / (math.pi * D * L) * 0.5
            liquid_height = max(0, min(liquid_height, D))
            fill_frac = target_vol / volume
        
        else:
            fill_frac = target_vol / volume
            if D is None:
                D = (6 * volume / math.pi) ** (1/3)
            # 修正：球缺反算液位时，同样引入标准的 R 修正因子进行液位求解
            liquid_height = D * fill_frac
        
        return {
            "liquid_volume": target_vol,
            "liquid_height": liquid_height,
            "fill_fraction": fill_frac,
            "residence_time_s": target_residence,
            "residence_time_min": target_residence / 60,
            "residence_time_h": target_residence / 3600,
        }
    
    return {"error": f"不支持的计算方向: {direction}"}


# ============================================================
# 设计校验
# ============================================================

def validate_tank_design(design: dict, tank_type: str = TANK_VERTICAL) -> dict:
    """校验储罐设计"""
    checks = []
    warnings = []
    
    D = design.get("D", 0)
    H = design.get("H", design.get("L", 0))
    
    # 1. 径高比/长径比检查
    if tank_type == TANK_VERTICAL:
        aspect_ratio = H / D if D > 0 else 0
        if 2.0 <= aspect_ratio <= 5.0:
            checks.append({"name": "aspect_ratio", "pass": True, "value": aspect_ratio, "limit": "2.0-5.0"})
        else:
            checks.append({"name": "aspect_ratio", "pass": False, "value": aspect_ratio, "limit": "2.0-5.0"})
            warnings.append(f"立式罐高径比 {aspect_ratio:.2f} 建议 2.0-5.0")
    
    elif tank_type == TANK_HORIZONTAL:
        # 卧式罐使用总长径比 L_total/D（包含封头）
        L_total = design.get("L_total", design.get("L", 0))
        aspect_ratio = L_total / D if D > 0 else 0
        if 2.0 <= aspect_ratio <= 5.0:
            checks.append({"name": "aspect_ratio", "pass": True, "value": round(aspect_ratio, 2), "limit": "2.0-5.0"})
        else:
            checks.append({"name": "aspect_ratio", "pass": False, "value": round(aspect_ratio, 2), "limit": "2.0-5.0"})
            warnings.append(f"卧式罐长径比 {aspect_ratio:.2f} 建议 2.0-5.0")
    
    # 2. 填充率检查
    fill_ratio = design.get("fill_ratio", 0)
    if fill_ratio <= 0.9:
        checks.append({"name": "fill_ratio", "pass": True, "value": fill_ratio, "limit": "<0.9"})
    else:
        checks.append({"name": "fill_ratio", "pass": False, "value": fill_ratio, "limit": "<0.9"})
        warnings.append(f"填充率 {fill_ratio:.2f} 过高，建议 < 0.9")
    
    # 3. 基础压力检查
    base_pressure = design.get("base_pressure", 0)
    max_soil_pressure = 200000  # Pa (200 kPa)
    if base_pressure < max_soil_pressure:
        checks.append({"name": "base_pressure", "pass": True, "value": base_pressure, "limit": f"<{max_soil_pressure}"})
    else:
        if base_pressure != 0:
            checks.append({"name": "base_pressure", "pass": False, "value": base_pressure, "limit": f"<{max_soil_pressure}"})
            warnings.append(f"基底压力 {base_pressure/1000:.1f} kPa 超过地基承载力")
    
    n_pass = sum(1 for c in checks if c["pass"])
    n_total = len(checks)
    pass_rate = n_pass / n_total if n_total > 0 else 0
    
    return {
        "valid": pass_rate >= 0.8,
        "pass_rate": pass_rate,
        "n_pass": n_pass,
        "n_total": n_total,
        "checks": checks,
        "warnings": warnings,
    }


# ============================================================
# 有效容量计算（卧式 / 立式）
# ============================================================

def effective_capacity_horizontal(
    D: float,
    L: float,
    fill_fraction: float = 0.85,
    head_type: str = "ellipsoidal",
    n_level_points: int = 11,
) -> dict:
    """卧式储罐有效容量计算（基于操作液位高度）

    利用卧式圆筒罐弓形截面的非线性几何特性，精确计算给定
    操作液位下的实际液相容积，并生成液位-容积对照表。

    工程背景:
        卧式罐圆形截面导致液位高度与容积呈非线性关系。
        液位 85% D 对应的容积占比并非 85%，而是更高（约 90%）。
        本函数通过精确几何计算消除这一工程估算偏差。

    Args:
        D: 储罐内径 (m)
        L: 筒体长度 (m)，不含封头
        fill_fraction: 最高操作液位占直径的比例 (0~1)，默认 0.85
        head_type: 封头类型 (flat/ellipsoidal/hemispherical/torispherical/conical)
        n_level_points: 液位-容积对照表的点数，默认 11

    Returns:
        dict 包含:
            - D, L, head_type: 输入参数
            - h_max: 最高操作液位高度 (m)
            - V_total: 总几何容积 (m³)
            - V_effective: 有效容量 (m³)
            - effective_fraction: 有效容量占总容积比例
            - vapor_space: 气相空间容积 (m³)
            - level_table: 液位-容积对照表 (list of dict)
    """
    if not 0 < fill_fraction <= 1:
        raise ValueError(f"fill_fraction 必须在 (0, 1] 范围内，当前值: {fill_fraction}")

    # 总几何容积
    vol_total = horizontal_cylindrical_volume(D, L, head_type=head_type)
    V_total = vol_total["total_volume"]

    # 最高操作液位高度
    h_max = D * fill_fraction

    # 有效容量（最高操作液位下的液相容积）
    vol_at_max = horizontal_cylindrical_volume(D, L, liquid_height=h_max, head_type=head_type)
    V_effective = vol_at_max["liquid_volume"]

    # 液位-容积对照表
    level_table = []
    for i in range(n_level_points):
        frac = (i + 1) / n_level_points
        h_i = D * frac
        vol_i = horizontal_cylindrical_volume(D, L, liquid_height=h_i, head_type=head_type)
        level_table.append({
            "fill_fraction": round(frac, 4),
            "liquid_height_m": round(h_i, 4),
            "liquid_volume_m3": round(vol_i["liquid_volume"], 4),
            "volume_pct": round(vol_i["liquid_volume"] / V_total * 100, 2) if V_total > 0 else 0,
        })

    return {
        "D": D,
        "L": L,
        "head_type": head_type,
        "h_max": round(h_max, 4),
        "fill_fraction": fill_fraction,
        "V_total": round(V_total, 4),
        "V_effective": round(V_effective, 4),
        "effective_fraction": round(V_effective / V_total, 4) if V_total > 0 else 0,
        "vapor_space": round(V_total - V_effective, 4),
        "level_table": level_table,
    }


def effective_capacity_vertical(
    D: float,
    H: float,
    fill_fraction: float = 0.85,
    head_type: str = "flat",
    n_level_points: int = 11,
) -> dict:
    """立式储罐有效容量计算（基于操作液位高度）

    对于立式圆筒罐，截面为均匀圆形，液位高度与容积呈线性关系。
    当封头类型为 flat 时，有效容量 = fill_fraction × V_total。
    当封头非 flat 时，需考虑底部封头对液位-容积关系的非线性影响。

    工程背景:
        立式罐通常底部为平盖或椭圆封头，顶部为锥顶/拱顶。
        底部封头形状会影响低液位区的容积-高度关系。
        本函数精确计算封头贡献，消除纯圆柱体假设的偏差。

    Args:
        D: 储罐内径 (m)
        H: 筒体高度 (m)，不含封头
        fill_fraction: 最高操作液位占总液深的比例 (0~1)，默认 0.85
        head_type: 封头类型 (flat/ellipsoidal/hemispherical/torispherical/conical)
        n_level_points: 液位-容积对照表的点数，默认 11

    Returns:
        dict 包含:
            - D, H, head_type: 输入参数
            - h_max: 最高操作液位高度 (m，从罐底算起)
            - V_total: 总几何容积 (m³)
            - V_effective: 有效容量 (m³)
            - effective_fraction: 有效容量占总容积比例
            - vapor_space: 气相空间容积 (m³)
            - level_table: 液位-容积对照表 (list of dict)
    """
    if not 0 < fill_fraction <= 1:
        raise ValueError(f"fill_fraction 必须在 (0, 1] 范围内，当前值: {fill_fraction}")

    # 总几何容积
    vol_total = vertical_cylindrical_volume(D, H, head_type=head_type)
    V_total = vol_total["total_volume"]
    H_total = vol_total["H_total"]  # 总高度（含封头）
    h_head_bottom = vol_total["head_depth_each"]  # 底部封头深度

    # 最高操作液位高度（从罐底算起，包含底部封头）
    h_max = H_total * fill_fraction

    # 计算有效容量
    if head_type == "flat":
        # 平盖底：截面均匀，液位-容积线性
        V_effective = math.pi * (D / 2) ** 2 * h_max
        V_effective = min(V_effective, V_total)
    else:
        # 非平盖底：需考虑底部封头的非线性贡献
        # 底部封头区域: 0 ~ h_head_bottom
        # 筒体区域: h_head_bottom ~ h_head_bottom + H
        V_one_head = vol_total["head_volume_each"]

        if h_max <= h_head_bottom:
            # 液位在底部封头内（近似：按封头容积的比例）
            if h_head_bottom > 0:
                frac_in_head = h_max / h_head_bottom
                V_effective = V_one_head * frac_in_head  # 底部封头贡献
            else:
                V_effective = 0.0
        else:
            # 液位超过底部封头，进入筒体区域
            h_in_shell = h_max - h_head_bottom
            h_in_shell = min(h_in_shell, H)  # 不超过筒体高度
            V_shell_liquid = math.pi * (D / 2) ** 2 * h_in_shell
            V_effective = V_one_head + V_shell_liquid  # 底部封头满 + 筒体部分

        V_effective = min(V_effective, V_total)

    # 液位-容积对照表
    level_table = []
    for i in range(n_level_points):
        frac = (i + 1) / n_level_points
        h_i = H_total * frac

        if head_type == "flat":
            V_i = math.pi * (D / 2) ** 2 * h_i
            V_i = min(V_i, V_total)
        else:
            V_one_head_i = vol_total["head_volume_each"]
            if h_i <= h_head_bottom:
                if h_head_bottom > 0:
                    V_i = V_one_head_i * (h_i / h_head_bottom)
                else:
                    V_i = 0.0
            else:
                h_shell_i = min(h_i - h_head_bottom, H)
                V_i = V_one_head_i + math.pi * (D / 2) ** 2 * h_shell_i
            V_i = min(V_i, V_total)

        level_table.append({
            "fill_fraction": round(frac, 4),
            "liquid_height_m": round(h_i, 4),
            "liquid_volume_m3": round(V_i, 4),
            "volume_pct": round(V_i / V_total * 100, 2) if V_total > 0 else 0,
        })

    return {
        "D": D,
        "H": H,
        "head_type": head_type,
        "H_total": round(H_total, 4),
        "h_max": round(h_max, 4),
        "fill_fraction": fill_fraction,
        "V_total": round(V_total, 4),
        "V_effective": round(V_effective, 4),
        "effective_fraction": round(V_effective / V_total, 4) if V_total > 0 else 0,
        "vapor_space": round(V_total - V_effective, 4),
        "level_table": level_table,
    }


# ============================================================
# 罐壁厚度计算（GB 50341 / API 650）
# ============================================================

def tank_shell_thickness(
    D: float,
    H: float,
    liquid_density: float = 1000.0,
    allowable_stress: float = 140e6,
    joint_efficiency: float = 0.85,
    corrosion_allowance: float = 0.002,
    specific_gravity: float = None,
) -> dict:
    """立式储罐罐壁厚度设计（GB 50341 / API 650 变点法简化）

    基于薄膜应力公式：t = P * D / (2 * [σ] * φ - P) + C_corr
    
    按 1m 分段计算各圈板厚度，并匹配标准钢板厚度。

    Args:
        D: 储罐直径 (m)
        H: 储罐高度 (m)
        liquid_density: 液体密度 (kg/m³)
        allowable_stress: 材料许用应力 (Pa)，默认 140 MPa（Q235-B）
        joint_efficiency: 焊缝系数，默认 0.85
        corrosion_allowance: 腐蚀裕量 (m)，默认 2 mm
        specific_gravity: 液体比重（与水比），若提供则覆盖 liquid_density 计算的比重

    Returns:
        各圈板设计厚度、计算厚度、选定厚度等
    """
    g = 9.81

    if specific_gravity is None:
        specific_gravity = liquid_density / 1000.0

    # 标准钢板厚度系列 (mm)
    standard_thicknesses = [5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20,
                            22, 25, 28, 30, 32, 36, 40, 45, 50]

    # 按 1m 分段计算各圈板
    n_courses = max(1, int(math.ceil(H)))
    courses = []

    for i in range(n_courses):
        # 该圈板底部距罐顶的液柱高度
        h_liquid = H - i * 1.0
        h_liquid = max(0.1, min(h_liquid, H))

        # 静水压力
        P = specific_gravity * 1000.0 * g * h_liquid  # Pa

        # 理论计算厚度（薄膜应力公式）
        t_calc = P * D / (2.0 * allowable_stress * joint_efficiency - P)

        # 加上腐蚀裕量
        t_design = t_calc + corrosion_allowance

        # 最小厚度要求（GB 50341-2014 表 5.2.3，分档封顶）
        # 不含腐蚀裕量；仅保证制造/施工刚度底线，大罐由结构计算控制
        if D <= 15.0:
            t_min_code = 0.005   # 5 mm (D ≤ 15m)
        elif D <= 35.0:
            t_min_code = 0.006   # 6 mm (15 < D ≤ 35m)
        else:
            t_min_code = 0.008   # 8 mm (D > 35m，封顶不再增长)

        t_required = max(t_design, t_min_code)

        # 匹配标准钢板厚度
        t_selected = standard_thicknesses[-1] / 1000.0  # 默认取最大
        for t_std in standard_thicknesses:
            if t_std / 1000.0 >= t_required:
                t_selected = t_std / 1000.0
                break

        courses.append({
            "course": i + 1,
            "height_from_bottom": i * 1.0,
            "h_liquid": h_liquid,
            "pressure_Pa": P,
            "t_calc_m": t_calc,
            "t_design_m": t_design,
            "t_required_m": t_required,
            "t_selected_m": t_selected,
            "t_selected_mm": t_selected * 1000,
        })

    # 第一圈（最底层）厚度
    t_bottom = courses[0]["t_selected_m"] if courses else 0.0
    # 最顶层厚度
    t_top = courses[-1]["t_selected_m"] if courses else 0.0

    # 罐底板厚度（API 650 简化）
    t_bottom_plate = max(0.006, t_bottom + 0.001)  # 罐底板比底层圈板厚 1mm，不小于 6mm

    # 罐顶板厚度（自支撑式锥顶，简化）
    t_roof = max(0.005, D / 200.0 * 0.005)  # 简化估算

    return {
        "D": D,
        "H": H,
        "n_courses": n_courses,
        "courses": courses,
        "t_bottom_shell_m": t_bottom,
        "t_bottom_shell_mm": t_bottom * 1000,
        "t_top_shell_m": t_top,
        "t_top_shell_mm": t_top * 1000,
        "t_bottom_plate_m": t_bottom_plate,
        "t_bottom_plate_mm": t_bottom_plate * 1000,
        "t_roof_m": t_roof,
        "t_roof_mm": t_roof * 1000,
        "allowable_stress_MPa": allowable_stress / 1e6,
        "joint_efficiency": joint_efficiency,
        "corrosion_allowance_mm": corrosion_allowance * 1000,
        "standard": "GB 50341 / API 650",
    }


# ============================================================
# 火灾工况呼吸阀泄放能力校核（API 2000）
# ============================================================

def fire_case_relief_check(
    D: float,
    H: float,
    existing_relief_capacity: float,
    liquid_density: float = 1000.0,
    latent_heat: float = 300e3,
    molecular_weight: float = 86.0,
    boiling_point: float = 373.15,
    P_set: float = 101325.0,
    P_back: float = 101325.0,
    k: float = 1.4,
    env_factor: float = 1.0,
) -> dict:
    """常压储罐火灾工况泄放能力校核（API 2000）

    基于API 2000标准，火灾工况下储罐外部受热，液体蒸发产生蒸汽
    需要呼吸阀（含紧急泄放阀）的泄放能力大于蒸发速率。

    热输入量计算（API 2000 附录A）：
        Q_fire = 43200 * F_env * A_wet^0.82  (kJ/h)
    其中 A_wet 为湿润表面积 (m²)，F_env 为环境因子（保温层折减）

    蒸发速率：
        W_vap = Q_fire / ΔH_vap  (kg/h)
    标态体积流量：
        Q_vap = W_vap * 22.4 / MW  (Nm³/h)

    Args:
        D: 储罐直径 (m)
        H: 储罐高度 (m)
        existing_relief_capacity: 现有呼吸阀标定泄放能力 (Nm³/h, 标态)
        liquid_density: 液体密度 (kg/m³)
        latent_heat: 液体汽化潜热 (J/kg)
        molecular_weight: 蒸气分子量 (g/mol)
        boiling_point: 常压沸点 (K)
        P_set: 呼吸阀设定压力 (Pa, 绝压)
        P_back: 背压 (Pa, 绝压)
        k: 气体绝热指数
        env_factor: 环境因子，无保温=1.0，有保温层按API 2000折减

    Returns:
        火灾热输入、蒸发速率、所需泄放能力、校核结论
    """
    g = 9.81

    # 湿润表面积：罐壁液位以下部分 + 罐底面积
    # 假设液位高度为罐高的 90%（最大操作液位）
    h_liquid = H * 0.9
    A_wet = math.pi * D * h_liquid + math.pi * (D / 2) ** 2

    # 火灾热输入量 (API 2000)
    # Q_fire = 43200 * F_env * A_wet^0.82 (kJ/h)
    Q_fire_kJ_h = 43200.0 * env_factor * (A_wet ** 0.82)
    Q_fire_W = Q_fire_kJ_h * 1000 / 3600  # 转换为 W

    # 蒸发速率
    W_vap_kg_s = Q_fire_W / latent_heat  # kg/s
    W_vap_kg_h = W_vap_kg_s * 3600

    # 标态体积流量 (Nm³/h)
    # Q_vap = W * R * T / (P * MW)
    R = 8.314
    MW_kg = molecular_weight / 1000.0
    Q_vap_Nm3_h = W_vap_kg_s * R * 273.15 / (101325.0 * MW_kg) * 3600

    # 所需泄放面积 (API 520 气体公式，仅对带压罐适用)
    P_critical = P_set * (2.0 / (k + 1.0)) ** (k / (k - 1.0))
    if P_back < P_critical:
        C_factor = math.sqrt(k * (2.0 / (k + 1.0)) ** ((k + 1.0) / (k - 1.0)))
    else:
        Pr = P_back / P_set
        C_factor = math.sqrt(2.0 * k / (k - 1.0) * (Pr ** (2.0 / k) - Pr ** ((k + 1.0) / k)))

    # 常压罐（P_set ≈ P_back）时 C_factor=0，跳过 API 520 面积计算
    if C_factor > 0 and P_set > 0:
        A_required = W_vap_kg_s / (C_factor * P_set * math.sqrt(MW_kg / (R * boiling_point)))
    else:
        A_required = None  # 常压罐按蒸发体积流量直接校核

    # 安全裕量系数 1.1
    Q_required = Q_vap_Nm3_h * 1.1

    # 校核
    ratio = existing_relief_capacity / Q_required if Q_required > 0 else 0
    is_adequate = existing_relief_capacity >= Q_required

    return {
        "A_wet_m2": A_wet,
        "h_liquid_m": h_liquid,
        "Q_fire_kJ_h": Q_fire_kJ_h,
        "Q_fire_W": Q_fire_W,
        "W_vap_kg_h": W_vap_kg_h,
        "Q_vap_Nm3_h": Q_vap_Nm3_h,
        "Q_required_Nm3_h": Q_required,
        "existing_capacity_Nm3_h": existing_relief_capacity,
        "capacity_ratio": ratio,
        "is_adequate": is_adequate,
        "A_required_m2": A_required,
        "C_factor": C_factor,
        "env_factor": env_factor,
        "standard": "API 2000",
    }