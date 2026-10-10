"""
填料床/固定床压降计算工具包
包含：
- Ergun 方程（经典形式）
- 13 种工业及学术界经典压降关联式可选（深度修正 Reichelt 壁效应与 Tallmadge 算法）
"""
from __future__ import annotations
import math
from typing import Optional, Dict


# ============================================================
# Ergun 方程（经典）
# ============================================================

def calculate_packed_bed_pressure_drop(
    Q: float,
    D_bed: float = None,
    L_bed: float = None,
    rho: float = None,
    mu: float = None,
    d_particle: float = None,
    epsilon: float = None,
    method: str = "ergun",
    # ── 散装填料专用参数 ──
    packing_Fp: float = None,      # 填料因子 (m⁻¹)，如 50mm PP阶梯环 Fp=98.4
    packing_a: float = None,       # 填料比表面积 (m²/m³)
    packing_name: str = None,      # 填料名称，如 'step_cascade_ring_50'，自动查表
    # ── 气液两相流参数 ──
    Q_liq: float = None,           # 液相体积流量 (m³/s)，提供时计算湿填料压降
    rho_liq: float = None,         # 液相密度 (kg/m³)
    mu_liq: float = None,          # 液相粘度 (Pa·s)
    sigma_liq: float = 0.0725,     # 气液表面张力 (N/m)，默认水
    flooding_factor: float = 0.7,  # 泛点因子（D_bed自动估算时）
    # ── 可选：用户给定流体信息，自动查询物性 ──
    fluid: str = None,
    T: float = None,
    P: float = 101325.0,
    phase: str = "gas",
    z: list = None,
) -> dict:
    """填料床/固定床压降计算（支持气液两相湿压降 + 自动塔径估算）
    
    Ergun 方程：
    ΔP/L = 150·(1-ε)²/ε³ · μ·u/d_p² + 1.75·(1-ε)/ε³ · ρ·u²/d_p
    
    物性参数支持用户给定优先：
    - 若提供 rho/mu，直接使用
    - 若未提供但提供 fluid/T，自动从数据库查询
    
    散装填料支持：
    - 提供 packing_Fp 时，自动换算当量直径 d_p = 6(1-ε)/Fp
    - 提供 packing_name 时，自动查表获取 Fp/ε/a
    
    塔径自动估算：
    - D_bed 未提供时，基于 Sherwood-Fair 泛点气速法自动估算
    - 需要提供 Q_liq 和 rho_liq 才能准确计算泛点
    
    气液两相湿压降：
    - 提供 Q_liq 时，采用 Leva 修正：ΔP_wet = ΔP_dry × (ε/ε_eff)³
    """
    from physics_engine.thermo_helper import get_fluid_density, get_fluid_viscosity
    from physics_engine.common import safe_float, safe_int
    
    # ── 防御性类型转换 ──
    Q = safe_float(Q, name='Q')
    D_bed = safe_float(D_bed, name='D_bed')
    L_bed = safe_float(L_bed, name='L_bed')
    rho = safe_float(rho, name='rho')
    mu = safe_float(mu, name='mu')
    d_particle = safe_float(d_particle, name='d_particle')
    epsilon = safe_float(epsilon, name='epsilon')
    packing_Fp = safe_float(packing_Fp, name='packing_Fp')
    packing_a = safe_float(packing_a, name='packing_a')
    Q_liq = safe_float(Q_liq, name='Q_liq')
    rho_liq = safe_float(rho_liq, name='rho_liq')
    mu_liq = safe_float(mu_liq, name='mu_liq')
    sigma_liq = safe_float(sigma_liq, default=0.0725, name='sigma_liq')
    flooding_factor = safe_float(flooding_factor, default=0.7, name='flooding_factor')
    T = safe_float(T, name='T')
    P = safe_float(P, default=101325.0, name='P')
    
    conversion_notes = []
    
    # ── 填料参数自动查表 ──
    if packing_name is not None:
        _packing = _PACKING_DATA.get(packing_name)
        if _packing is None:
            _key = packing_name.lower().replace(" ", "_").replace("-", "_")
            _packing = _PACKING_DATA.get(_key)
        if _packing is not None:
            if packing_Fp is None:
                packing_Fp = _packing["Fp"]
                conversion_notes.append(f"[填料查表] {packing_name}: Fp={packing_Fp} m⁻¹")
            if epsilon is None:
                epsilon = _packing["epsilon"]
                conversion_notes.append(f"[填料查表] {packing_name}: ε={epsilon}")
            if packing_a is None:
                packing_a = _packing["a"]
            if d_particle is None:
                d_particle = _packing["dp"]
                conversion_notes.append(f"[填料查表] {packing_name}: dp={d_particle} m")
    
    # epsilon 默认值
    if epsilon is None:
        epsilon = 0.90  # 散装填料典型值
        conversion_notes.append(f"[默认] 空隙率 ε={epsilon}")
    
    # ── 密度：用户给定值优先 ──
    if rho is not None:
        rho_val = float(rho)
        conversion_notes.append(f"[用户给定] 气相密度 rho = {rho_val} kg/m³")
    elif fluid is not None and T is not None:
        try:
            rho_val = get_fluid_density(fluid, T, P, phase=phase, z=z)
            conversion_notes.append(f"[自动查询] {fluid} 在 {T}K, {P}Pa 的气相密度 = {rho_val:.2f} kg/m³")
        except Exception as e:
            return {"error": f"密度查询失败：{e}。请手动提供 rho 参数。"}
    else:
        return {"error": "必须提供 rho 或 (fluid + T) 参数"}
    
    # ── 粘度：用户给定值优先 ──
    if mu is not None:
        mu_val = float(mu)
        conversion_notes.append(f"[用户给定] 气相粘度 mu = {mu_val} Pa·s")
    elif fluid is not None and T is not None:
        try:
            mu_val = get_fluid_viscosity(fluid, T, P, phase=phase, z=z)
            conversion_notes.append(f"[自动查询] {fluid} 在 {T}K, {P}Pa 的气相粘度 = {mu_val:.6f} Pa·s")
        except Exception as e:
            return {"error": f"粘度查询失败：{e}。请手动提供 mu 参数。"}
    else:
        return {"error": "必须提供 mu 或 (fluid + T) 参数"}
    
    # ── 颗粒直径：填料因子换算 ──
    # 关键：Ergun 方程必须使用水力当量直径，而非填料公称直径
    # 优先使用 a-based 公式（Billet-Schultes 标准）：d_p,equiv = 6(1-ε)/a
    # 备选使用 Fp-based 公式：d_p,equiv = 6(1-ε)/Fp
    dp_nominal = None  # 保留公称直径用于其他用途
    if packing_a is not None and packing_a > 0:
        # 优先：基于比表面积的当量直径（Billet-Schultes 标准）
        dp_nominal = float(d_particle) if d_particle is not None else None
        dp_val = 6 * (1 - epsilon) / packing_a
        conversion_notes.append(f"[当量直径] a={packing_a} m²/m³, ε={epsilon} → d_p,equiv=6(1-ε)/a={dp_val:.5f} m")
        if dp_nominal is not None:
            conversion_notes.append(f"[注意] 公称直径 dp={dp_nominal} m 不用于Ergun方程")
    elif packing_Fp is not None:
        # 备选：基于填料因子的当量直径
        dp_nominal = float(d_particle) if d_particle is not None else None
        dp_val = 6 * (1 - epsilon) / packing_Fp
        conversion_notes.append(f"[当量直径] Fp={packing_Fp} m⁻¹, ε={epsilon} → d_p,equiv=6(1-ε)/Fp={dp_val:.5f} m")
        if dp_nominal is not None:
            conversion_notes.append(f"[注意] 公称直径 dp={dp_nominal} m 不用于Ergun方程")
    elif d_particle is not None:
        dp_val = float(d_particle)
        conversion_notes.append(f"[颗粒直径] d_particle={dp_val} m")
    else:
        return {"error": "必须提供 d_particle、packing_Fp 或 packing_a 参数"}
    
    if packing_Fp is not None:
        conversion_notes.append(f"[填料参数] Fp={packing_Fp} m⁻¹, ε={epsilon}")
    
    # ── 塔径自动估算（基于泛点气速法）──
    tower_sizing_result = None
    if D_bed is not None:
        D_val = float(D_bed)
    else:
        # 使用 Sherwood-Fair 泛点气速法估算塔径
        # 需要液相信息
        if Q_liq is None or rho_liq is None:
            return {"error": "D_bed 未提供时，必须提供 Q_liq 和 rho_liq 以估算塔径"}
        
        # 液相粘度默认值
        mu_liq_val = float(mu_liq) if mu_liq is not None else 0.001
        
        # 调用泛点气速法估算塔径
        tower_sizing_result = select_packed_tower_diameter(
            Q_gas=float(Q),
            Q_liquid=float(Q_liq),
            rho_gas=rho_val,
            rho_liquid=float(rho_liq),
            mu_liquid=mu_liq_val,
            packing_name=packing_name if packing_name else "step_cascade_ring_50",
            flooding_factor=flooding_factor,
            F_packing=packing_Fp,
        )
        D_val = tower_sizing_result["D_selected_m"]
        conversion_notes.append(
            f"[自动选型] 泛点气速={tower_sizing_result['u_flood_m_s']:.2f} m/s, "
            f"操作气速={tower_sizing_result['u_actual_m_s']:.2f} m/s, "
            f"液泛率={tower_sizing_result['flooding_ratio']*100:.1f}%, "
            f"选定塔径 DN={D_val} m"
        )
    
    # L_bed 默认值
    if L_bed is None:
        L_bed_val = 3.0  # 典型填料层高度
        conversion_notes.append(f"[默认] 填料层高度 L_bed={L_bed_val} m")
    else:
        L_bed_val = float(L_bed)
    
    # 表观气速
    A = math.pi * (D_val / 2) ** 2
    u_G = float(Q) / A
    
    # 颗粒雷诺数
    Re_p = rho_val * u_G * dp_val / mu_val if mu_val > 0 else 1e6
    
    if epsilon <= 0 or epsilon >= 1:
        raise ValueError("空隙率必须在 (0, 1) 范围内")
    
    # ── 干床压降（Ergun 或其他关联式）──
    if method == "ergun":
        dP_L_dry = _ergun(u_G, rho_val, mu_val, dp_val, epsilon)
    elif method == "carman_kozeny":
        dP_L_dry = _carman_kozeny(u_G, rho_val, mu_val, dp_val, epsilon)
    elif method == "burke_plummer":
        dP_L_dry = _burke_plummer(u_G, rho_val, mu_val, dp_val, epsilon)
    elif method == "tallmadge":
        dP_L_dry = _tallmadge(u_G, rho_val, mu_val, dp_val, epsilon, Re_p)
    elif method == "macdonald":
        dP_L_dry = _macdonald(u_G, rho_val, mu_val, dp_val, epsilon, Re_p)
    elif method == "hicks":
        dP_L_dry = _hicks(u_G, rho_val, mu_val, dp_val, epsilon, Re_p)
    elif method == "eyeron_gunther":
        dP_L_dry = _eyeron_gunther(u_G, rho_val, mu_val, dp_val, epsilon)
    elif method == "morcom":
        dP_L_dry = _morcom(u_G, rho_val, mu_val, dp_val, epsilon, Re_p)
    elif method == "fahien":
        dP_L_dry = _fahien(u_G, rho_val, mu_val, dp_val, epsilon, Re_p)
    elif method == "reichelt":
        dP_L_dry = _reichelt(u_G, rho_val, mu_val, dp_val, epsilon, Re_p, D_val)
    elif method == "winterberg":
        dP_L_dry = _winterberg(u_G, rho_val, mu_val, dp_val, epsilon, Re_p)
    elif method == "nemec":
        dP_L_dry = _nemec(u_G, rho_val, mu_val, dp_val, epsilon, Re_p)
    elif method == "ouwer":
        dP_L_dry = _ouwer(u_G, rho_val, mu_val, dp_val, epsilon, Re_p)
    else:
        raise ValueError(f"未知关联式: {method}")
    
    # ── 湿床压降修正（Leva 修正，仅当提供 Q_liq 时）──
    dP_L_wet = dP_L_dry
    liquid_holdup = None
    if Q_liq is not None and float(Q_liq) > 0:
        u_L = float(Q_liq) / A  # 液相表观流速
        
        # Leva 持液量估算
        # 静态持液量：按填料类型
        if packing_name and ("raschig" in packing_name.lower()):
            h_static = 0.04
        elif packing_name and ("pall" in packing_name.lower() or "cmr" in packing_name.lower()):
            h_static = 0.025
        else:
            h_static = 0.03  # 阶梯环默认
        
        # 动态持液量
        h_dynamic = 0.12 * u_L ** 0.4 if u_L > 0 else 0
        
        # 总持液量，施加上限
        h_total = h_static + h_dynamic
        h_max = epsilon * 0.5
        h_total = min(h_total, h_max)
        
        # 有效空隙率
        eps_eff = epsilon - h_total
        
        # Leva 修正：ΔP_wet = ΔP_dry × (ε/ε_eff)³
        if eps_eff > 0 and eps_eff < epsilon:
            correction = (epsilon / eps_eff) ** 3
            dP_L_wet = dP_L_dry * correction
            liquid_holdup = {
                "h_static": h_static,
                "h_dynamic": h_dynamic,
                "h_total": h_total,
                "epsilon_eff": eps_eff,
                "correction_factor": correction,
                "u_L_m_s": u_L,
            }
            conversion_notes.append(
                f"[湿修正] 持液量 h={h_total:.4f}, ε_eff={eps_eff:.3f}, "
                f"修正因子={correction:.2f}"
            )
    
    dP_total = dP_L_wet * L_bed_val
    
    # Ergun 分解（黏性项 + 惯性项）
    viscous_term = 150 * (1 - epsilon) ** 2 / epsilon ** 3 * mu_val * u_G / dp_val ** 2
    inertial_term = 1.75 * (1 - epsilon) / epsilon ** 3 * rho_val * u_G ** 2 / dp_val
    
    result = {
        "total_pressure_drop_Pa": dP_total,
        "pressure_drop_per_m": dP_L_wet,
        "pressure_drop_dry_per_m": dP_L_dry,
        "superficial_gas_velocity": u_G,
        "Re_particle": Re_p,
        "d_particle_m": dp_val,
        "method": method,
        "viscous_term_Pa_m": viscous_term,
        "inertial_term_Pa_m": inertial_term,
        "conversion_notes": conversion_notes,
    }
    
    if tower_sizing_result is not None:
        result["tower_sizing"] = tower_sizing_result
    
    if liquid_holdup is not None:
        result["liquid_holdup"] = liquid_holdup
    
    return result


# --- 13 种关联式 ---

def _ergun(u, rho, mu, dp, eps):
    """经典 Ergun 方程"""
    viscous = 150 * (1 - eps) ** 2 / eps ** 3 * mu * u / dp ** 2
    inertial = 1.75 * (1 - eps) / eps ** 3 * rho * u ** 2 / dp
    return viscous + inertial


def _carman_kozeny(u, rho, mu, dp, eps):
    """Carman-Kozeny 方程（仅黏性项，适用于 Re_p < 10）"""
    return 180 * (1 - eps) ** 2 / eps ** 3 * mu * u / dp ** 2


def _burke_plummer(u, rho, mu, dp, eps):
    """Burke-Plummer 方程（仅惯性项，适用于 Re_p > 1000）"""
    return 1.75 * (1 - eps) / eps ** 3 * rho * u ** 2 / dp


def _tallmadge(u, rho, mu, dp, eps, Re_p):
    """Tallmadge 修正 Ergun（高 Re 区更准确）
    
    修正：还原 Tallmadge (1970) 原作中惯性阻力项随雷诺数 1/6 次幂衰减的非线性物理规律，
    并且将 (1-eps) 项的幂次精确修正为 7/6 (1.1667)。
    """
    viscous = 150 * (1 - eps) ** 2 / eps ** 3 * mu * u / dp ** 2
    if Re_p > 0:
        inertial = 4.2 * (1.0 - eps) ** 1.1667 / (eps ** 3 * Re_p ** (1.0 / 6.0)) * rho * u ** 2 / dp
    else:
        inertial = 4.2 * (1.0 - eps) ** 1.1667 / eps ** 3 * rho * u ** 2 / dp
    return viscous + inertial


def _macdonald(u, rho, mu, dp, eps, Re_p):
    """MacDonald 修正（根据 Re_p 调整系数）"""
    if Re_p < 10:
        C1, C2 = 180, 1.8
    else:
        C1, C2 = 150, 1.75
    return C1 * (1 - eps) ** 2 / eps ** 3 * mu * u / dp ** 2 + \
           C2 * (1 - eps) / eps ** 3 * rho * u ** 2 / dp


def _hicks(u, rho, mu, dp, eps, Re_p):
    """Hicks 关联式"""
    f = 3.5 / Re_p + 0.12 * Re_p ** (-0.1)
    return f * (1 - eps) / eps ** 3 * rho * u ** 2 / dp


def _eyeron_gunther(u, rho, mu, dp, eps):
    """Eyeron-Gunther 关联式"""
    return 150 * (1 - eps) ** 2 / eps ** 3 * mu * u / dp ** 2 + \
           1.75 * (1 - eps) / eps ** 3 * rho * u ** 2 / dp * (1 + 0.15 * (1 - eps))


def _morcom(u, rho, mu, dp, eps, Re_p):
    """Morcom 关联式"""
    f = 16 / Re_p + 3.5 / Re_p ** 0.5 + 0.12
    return f * (1 - eps) / eps ** 3 * rho * u ** 2 / dp


def _fahien(u, rho, mu, dp, eps, Re_p):
    """Fahien-Smith 关联式"""
    f = 150 / Re_p + 1.75
    return f * (1 - eps) ** 2 / eps ** 3 * rho * u ** 2 / dp


def _reichelt(u, rho, mu, dp, eps, Re_p, D_bed):
    """Reichelt 关联式（严格考虑容器壁摩擦对剪切效应影响的物理模型）
    
    修正：接收 D_bed 替换掉有 bug 的 locals/dir() 判断，
    并精确复原 Reichelt 1972 经典的 Aw^2 和 Bw 雷诺数修正因子。
    """
    if D_bed > 0 and dp > 0:
        ratio = D_bed / dp
        aw = 1.0 + 2.0 / (3.0 * ratio * (1.0 - eps))
        bw = (1.5 / ratio + 0.88) ** 2
        f = 150 * aw ** 2 / Re_p + aw * bw
    else:
        f = 150 / Re_p + 1.75
        
    return f * (1 - eps) ** 2 / eps ** 3 * rho * u ** 2 / dp


def _winterberg(u, rho, mu, dp, eps, Re_p):
    """Winterberg 关联式"""
    f = 150 / Re_p + 2.0
    return f * (1 - eps) ** 2 / eps ** 3 * rho * u ** 2 / dp


def _nemec(u, rho, mu, dp, eps, Re_p):
    """Nemec-Levec 关联式"""
    f = 150 / Re_p + 1.75 * (1 + 0.5 / (1 - eps))
    return f * (1 - eps) ** 2 / eps ** 3 * rho * u ** 2 / dp


def _ouwer(u, rho, mu, dp, eps, Re_p):
    """Ouwerkerk 关联式"""
    f = 150 / Re_p + 1.75 * Re_p ** 0.1
    return f * (1 - eps) ** 2 / eps ** 3 * rho * u ** 2 / dp


# ============================================================
# Billet-Schultes 严格模型 (1993, 1999)
# 工业标准散装填料压降计算方法
# 子模型: 干床压降关联 + 持液量(静态+动态) + 湿床修正
# ============================================================

def calculate_packed_bed_pressure_drop_bs(
    u_G: float,
    u_L: float,
    rho_G: float,
    rho_L: float,
    mu_G: float,
    mu_L: float,
    sigma: float,
    epsilon: float,
    a: float,
    Fp: float,
    dp_nominal: float = None,
) -> dict:
    """Billet-Schultes 严格模型计算散装填料层压降

    参考文献:
      Billet R, Schultes M. "Prediction of mass transfer columns with dumped
      and arranged packings - updated summary of the calculation method
      Billet/Schultes". Chem Eng Res Des, 1999, 77(A6): 498-510.

    模型框架:
      1. 干床压降: ΔP_dry/Z = Ψ₀ · (a/ε³) · (ρ_G·u_G²/2)
         其中 Ψ₀ = B₀/Re_G + C₀/√Re_G (层流/湍流分区)
      2. 持液量: h_total = h_static + h_dynamic
         h_static = C_st / Fr_L^(1/3)  (表面张力/重力/几何)
         h_dynamic = C_dyn · (Re_L)^(1/3)  (液相雷诺数关联)
      3. 湿床压降: ΔP_wet = ΔP_dry · (ε/(ε-h_total))^(25/12)

    Args:
        u_G: 气相表观流速 (m/s)
        u_L: 液相表观流速 (m/s)
        rho_G: 气相密度 (kg/m³)
        rho_L: 液相密度 (kg/m³)
        mu_G: 气相动力黏度 (Pa·s)
        mu_L: 液相动力黏度 (Pa·s)
        sigma: 气液表面张力 (N/m)
        epsilon: 填料空隙率 (-)
        a: 填料比表面积 (m²/m³)
        Fp: 填料因子 (m⁻¹)
        dp_nominal: 填料公称直径 (m), 可选

    Returns:
        压降计算结果字典，包含干/湿压降、持液量分项、雷诺数等
    """
    g = 9.81

    # 当量直径
    dp_eq = 6 * (1 - epsilon) / a

    # ---- 气相修正雷诺数 ----
    # Re_G,mod = u_G · ρ_G / (a · μ_G)  (Billet-Schultes 定义)
    Re_G_mod = u_G * rho_G / (a * mu_G) if mu_G > 0 else 1e6

    # ---- 1. 干床压降 ----
    # 阻力系数 Ψ₀ = B₀/Re + C₀/√Re
    # 层流区 (Re < 2000): B₀=160, C₀=3.5  (对应 Ergun 150→160, 1.75→3.5)
    # 湍流区 (Re >= 2000): B₀=8.0,  C₀=1.5  (高Re下阻力系数更低)
    if Re_G_mod < 2000:
        B0, C0 = 160.0, 3.5
        flow_regime = "laminar"
    else:
        B0, C0 = 8.0, 1.5
        flow_regime = "turbulent"

    Psi0 = B0 / Re_G_mod + C0 / math.sqrt(Re_G_mod)
    dP_dry = Psi0 * (a / epsilon ** 3) * (rho_G * u_G ** 2 / 2)

    # ---- 2. 持液量模型 ----
    # 液相修正雷诺数
    Re_L_mod = u_L * rho_L / (a * mu_L) if mu_L > 0 else 1e-10

    # Froude 数 (液相)
    Fr_L = a * u_L ** 2 / g if g > 0 else 0

    # 静态持液量 h_static
    # h_static = C_st · Fr_L^(1/3) · (σ_crit/σ)^(1/3)
    # 其中 C_st ≈ 0.055 (50mm 散装填料典型值)
    # σ_crit 为参考表面张力, 取水 0.0725 N/m
    # 当 σ 较大时, 静态持液量减小 (液体更容易流下)
    C_st = 0.055
    sigma_ref = 0.0725  # 参考表面张力 (水, 20°C), N/m
    sigma_ratio = min(sigma_ref / max(sigma, 0.01), 2.0)  # 限制极端值
    h_static = C_st * Fr_L ** (1.0 / 3.0) * sigma_ratio ** (1.0 / 3.0)

    # 动态持液量 h_dynamic
    # h_dynamic = C_dyn · Re_L^(1/3)
    # C_dyn 需按填料尺寸修正: 大填料理应持液更少
    # 基准: C_dyn,ref = 0.035 (对应 25mm 填料)
    # 修正: C_dyn = C_dyn,ref × (25mm/dp_nominal)^0.5
    # 若无公称直径, 从 a 和 ε 反推: dp ≈ 6(1-ε)/a
    C_dyn_ref = 0.035  # 25mm 填料基准
    if dp_nominal and dp_nominal > 0:
        dp_ref_for_correction = dp_nominal
    else:
        dp_ref_for_correction = dp_eq
    size_correction = (0.025 / dp_ref_for_correction) ** 0.5
    size_correction = max(0.3, min(size_correction, 2.0))  # 限制极端值
    C_dyn = C_dyn_ref * size_correction
    h_dynamic = C_dyn * Re_L_mod ** (1.0 / 3.0) if Re_L_mod > 0 else 0

    # 总持液量, 施加上限约束
    h_total = h_static + h_dynamic
    h_max = epsilon * 0.5  # 安全上限: 不超过空隙率的一半
    h_total = min(h_total, h_max)

    # 有效空隙率
    eps_eff = epsilon - h_total

    # ---- 3. 湿床压降 ----
    # Billet-Schultes 修正: ΔP_wet = ΔP_dry · (ε/ε_eff)^(25/12)
    # 指数 25/12 ≈ 2.083, 介于 Ergun 理论值 3 与工程简化值 2 之间
    if eps_eff > 0 and eps_eff < epsilon:
        correction = (epsilon / eps_eff) ** (25.0 / 12.0)
        dP_wet = dP_dry * correction
    else:
        correction = 1.0
        dP_wet = dP_dry

    # 安全警告
    if h_total >= h_max:
        warning = f"持液量已达上限 ({h_max:.4f}), 计算结果可能偏保守"
    elif correction > 5.0:
        warning = f"湿修正因子 {correction:.1f} 过大, 可能接近液泛"
    else:
        warning = None

    return {
        "pressure_drop_dry_Pa_m": dP_dry,
        "pressure_drop_wet_Pa_m": dP_wet,
        "correction_factor": correction,
        "Psi0": Psi0,
        "Re_G_modified": Re_G_mod,
        "Re_L_modified": Re_L_mod,
        "Fr_L": Fr_L,
        "flow_regime": flow_regime,
        "B0": B0,
        "C0": C0,
        "h_static": h_static,
        "h_dynamic": h_dynamic,
        "h_total": h_total,
        "epsilon_eff": eps_eff,
        "dp_equivalent_m": dp_eq,
        "warning": warning,
    }


# ============================================================
# 填料塔直径选择
# ============================================================

# 常见散装填料参数
# Fp 单位: m⁻¹ (SI)，原始数据通常为 ft⁻¹，需乘以 3.281 换算
# 注释中标注原始 ft⁻¹ 值供工程对照
# 数据来源: Koch-Glitsch, NJC, 化工填料手册等厂家数据
_PACKING_DATA = {
    #                            Fp(m⁻¹)  = Fp(ft⁻¹) × 3.281
    # 50mm 鲍尔环 (PP): Fp≈66-80 m⁻¹, ε≈0.94-0.95
    "pall_ring_50":          {"dp": 0.050, "epsilon": 0.95, "Fp":  72.2, "a":  95.0},  # 22 ft⁻¹
    "pall_ring_38":          {"dp": 0.038, "epsilon": 0.94, "Fp": 114.8, "a": 130.0},  # 35 ft⁻¹
    "pall_ring_25":          {"dp": 0.025, "epsilon": 0.93, "Fp": 196.9, "a": 200.0},  # 60 ft⁻¹
    # 50mm 拉西环 (PP): Fp≈500-600 m⁻¹, ε≈0.85-0.90
    "raschig_ring_50":       {"dp": 0.050, "epsilon": 0.85, "Fp": 557.8, "a":  95.0},  # 170 ft⁻¹
    "raschig_ring_25":       {"dp": 0.025, "epsilon": 0.78, "Fp":1640.5, "a": 190.0},  # 500 ft⁻¹
    # 50mm 矩鞍环 (PP): Fp≈180-210 m⁻¹, ε≈0.90-0.92
    "intalox_saddle_50":     {"dp": 0.050, "epsilon": 0.91, "Fp": 196.9, "a": 125.0},  # 60 ft⁻¹
    "intalox_saddle_25":     {"dp": 0.025, "epsilon": 0.90, "Fp": 984.3, "a": 250.0},  # 300 ft⁻¹
    # 50mm 阶梯环/CMR (PP): Fp≈90-110 m⁻¹, ε≈0.90-0.92
    "cascade_mini_ring_50":  {"dp": 0.050, "epsilon": 0.92, "Fp":  95.1, "a": 115.0},  # 29 ft⁻¹
    "cmr_50":                {"dp": 0.050, "epsilon": 0.92, "Fp":  95.1, "a": 115.0},  # 29 ft⁻¹
    "step_cascade_ring_50":  {"dp": 0.050, "epsilon": 0.90, "Fp":  98.4, "a": 120.0},  # 30 ft⁻¹
    # 50mm Hilflow/NMR环 (PP): Fp≈80-100 m⁻¹, ε≈0.92-0.94
    "hilflow_ring_50":       {"dp": 0.050, "epsilon": 0.93, "Fp":  85.3, "a":  90.0},  # 26 ft⁻¹
    "nmr_50":                {"dp": 0.050, "epsilon": 0.93, "Fp":  85.3, "a":  90.0},  # 26 ft⁻¹
}

# 标准塔径系列 (m)
_STANDARD_TOWER_DIAMETERS = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4,
                             1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0]


def select_packed_tower_diameter(
    Q_gas: float,
    Q_liquid: float,
    rho_gas: float,
    rho_liquid: float,
    mu_liquid: float = 0.001,
    packing_name: str = "pall_ring_50",
    flooding_factor: float = 0.7,
    F_packing: float = None,
) -> dict:
    """填料塔直径选择（基于GPDC泛点气速法）

    使用 Sherwood-Fair 关联式估算泛点气速，按设计因子（通常 0.6~0.8 倍泛点气速）
    确定操作气速，进而计算所需塔径并匹配标准系列。

    Sherwood-Fair 泛点关联式：
        Y_f = exp[A - B*(log X)^2]
    其中：
        X = (L/G) * sqrt(ρ_g / ρ_l)    (流动参数)
        Y = u_f² * Fp / g * (ρ_g / ρ_l) * (μ_l / μ_w)^0.2  (容量参数)
    
    泛点气速：
        u_f = sqrt(Y_f * g * ρ_l / (Fp * ρ_g * (μ_l/μ_w)^0.2))

    Args:
        Q_gas: 气相体积流量 (m³/s)
        Q_liquid: 液相体积流量 (m³/s)
        rho_gas: 气相密度 (kg/m³)
        rho_liquid: 液相密度 (kg/m³)
        mu_liquid: 液相黏度 (Pa·s)，默认 0.001
        packing_name: 填料名称，参考 _PACKING_DATA
        flooding_factor: 泛点因子，设计操作气速/泛点气速，默认 0.7
        F_packing: 填料因子 Fp，若提供则覆盖 packing_name 中的值

    Returns:
        塔径设计结果，包括泛点气速、操作气速、计算塔径、标准塔径等
    """
    from physics_engine.common import safe_float
    
    # ── 防御性类型转换 ──
    Q_gas = safe_float(Q_gas, name='Q_gas')
    Q_liquid = safe_float(Q_liquid, name='Q_liquid')
    rho_gas = safe_float(rho_gas, name='rho_gas')
    rho_liquid = safe_float(rho_liquid, name='rho_liquid')
    mu_liquid = safe_float(mu_liquid, default=0.001, name='mu_liquid')
    flooding_factor = safe_float(flooding_factor, default=0.7, name='flooding_factor')
    F_packing = safe_float(F_packing, name='F_packing')
    
    g = 9.81
    mu_w = 0.001  # 水在 20°C 下的黏度参考值

    # 获取填料参数
    if F_packing is not None:
        Fp = F_packing
        epsilon = 0.95
        dp = 0.050
    else:
        packing = _PACKING_DATA.get(packing_name)
        if packing is None:
            # 尝试模糊匹配
            key = packing_name.lower().replace(" ", "_").replace("-", "_")
            packing = _PACKING_DATA.get(key)
        if packing is None:
            raise ValueError(
                f"未知填料类型: {packing_name}。"
                f"可用类型: {list(_PACKING_DATA.keys())}"
            )
        Fp = packing["Fp"]
        epsilon = packing["epsilon"]
        dp = packing["dp"]

    # 质量流量
    m_gas = rho_gas * Q_gas
    m_liquid = rho_liquid * Q_liquid

    # 流动参数 X = (L/G) * sqrt(ρ_g / ρ_l)
    if m_gas > 0:
        X = (m_liquid / m_gas) * math.sqrt(rho_gas / rho_liquid)
    else:
        X = 0.01

    X = max(X, 1e-4)

    # Sherwood-Fair 泛点关联式
    # Eckert GPDC 泛点曲线（Kister/Stichlmair 经典拟合）:
    # log10(Y_f) = -1.668 - 1.085*log10(X) - 0.297*(log10(X))^2
    log_X = math.log10(X)
    log_Y_flood = -1.668 - 1.085 * log_X - 0.297 * log_X ** 2
    Y_flood = 10 ** log_Y_flood

    # 泛点气速
    viscosity_correction = (mu_liquid / mu_w) ** 0.2
    # Y = u_f² * Fp / g * (ρ_g / ρ_l) * viscosity_correction
    # u_f = sqrt(Y_f * g * ρ_l / (Fp * ρ_g * viscosity_correction))
    u_flood = math.sqrt(
        Y_flood * g * rho_liquid / (Fp * rho_gas * viscosity_correction)
    )

    # 操作气速
    u_operating = u_flood * flooding_factor

    # 泛点气速合理性校验（基于填料尺寸的经验范围）
    # 散装填料泛点气速经验上限：
    #   dp >= 50mm: 3-5 m/s
    #   dp 25-50mm: 2-4 m/s  
    #   dp < 25mm: 1.5-3 m/s
    if dp >= 0.050:
        u_flood_max_expected = 5.0
    elif dp >= 0.025:
        u_flood_max_expected = 4.0
    else:
        u_flood_max_expected = 3.0
    
    flood_velocity_warning = None
    if u_flood > u_flood_max_expected:
        flood_velocity_warning = (
            f"泛点气速 {u_flood:.2f} m/s 超出 {dp*1000:.0f}mm 填料的经验上限 "
            f"({u_flood_max_expected} m/s)，请检查输入参数或填料类型"
        )

    # 计算塔径
    if u_operating > 0:
        A_required = Q_gas / u_operating
        D_calculated = math.sqrt(4 * A_required / math.pi)
    else:
        D_calculated = 1.0

    # 匹配标准塔径
    D_selected = _STANDARD_TOWER_DIAMETERS[-1]
    for D_std in _STANDARD_TOWER_DIAMETERS:
        if D_std >= D_calculated:
            D_selected = D_std
            break

    # 实际操作参数
    A_actual = math.pi * (D_selected / 2) ** 2
    u_actual = Q_gas / A_actual if A_actual > 0 else 0
    flooding_ratio = u_actual / u_flood if u_flood > 0 else 0

    # 压降估算（Robbins 或 Leva 修正）
    # ΔP/Z ≈ C1 * Fp^0.7 * (G²/(ρ_g*ρ_l)) * μ_l^0.1
    # 简化估算：操作点压降约 200-400 Pa/m（60-80% 液泛）
    if flooding_ratio < 0.5:
        dp_per_m_estimate = 150.0
    elif flooding_ratio < 0.7:
        dp_per_m_estimate = 300.0
    elif flooding_ratio < 0.9:
        dp_per_m_estimate = 600.0
    else:
        dp_per_m_estimate = 1000.0

    return {
        "D_calculated_m": D_calculated,
        "D_selected_m": D_selected,
        "u_flood_m_s": u_flood,
        "u_operating_m_s": u_operating,
        "u_actual_m_s": u_actual,
        "flooding_ratio": flooding_ratio,
        "flooding_factor_design": flooding_factor,
        "flow_parameter_X": X,
        "capacity_parameter_Y": Y_flood,
        "Fp": Fp,
        "packing_name": packing_name,
        "packing_dp": dp,
        "packing_epsilon": epsilon,
        "A_cross_section_m2": A_actual,
        "dp_per_m_estimate_Pa": dp_per_m_estimate,
        "is_safe": flooding_ratio < 0.95,
        "flood_velocity_warning": flood_velocity_warning,
        "standard_diameters": _STANDARD_TOWER_DIAMETERS,
    }


# ============================================================
# 填料塔液泛率计算（校验已有塔）
# ============================================================

def packed_bed_flooding_check(
    D_tower: float,
    G_mass: float,
    L_mass: float,
    rho_gas: float,
    rho_liquid: float,
    mu_liquid: float = 0.001,
    packing_name: str = "pall_ring_50",
    F_packing: float = None,
) -> dict:
    """填料塔液泛率校核（基于GPDC泛点关联式）

    对已有填料塔，计算当前操作点相对于液泛点的液泛率，评估安全裕度。

    液泛率 = u_actual / u_flood * 100%

    安全操作区间：液泛率 < 80%（一般设计取 60-80%）

    Args:
        D_tower: 塔径 (m)
        G_mass: 气相质量流速 (kg/(m²·s))，或总气相质量流量 (kg/s)
        L_mass: 液相质量流速 (kg/(m²·s))，或总液相质量流量 (kg/s)
        rho_gas: 气相密度 (kg/m³)
        rho_liquid: 液相密度 (kg/m³)
        mu_liquid: 液相黏度 (Pa·s)
        packing_name: 填料名称
        F_packing: 填料因子，若提供则覆盖

    Returns:
        液泛率、泛点气速、操作气速、安全评估等
    """
    g = 9.81
    mu_w = 0.001

    # 获取填料参数
    if F_packing is not None:
        Fp = F_packing
    else:
        packing = _PACKING_DATA.get(packing_name)
        if packing is None:
            key = packing_name.lower().replace(" ", "_").replace("-", "_")
            packing = _PACKING_DATA.get(key)
        if packing is None:
            raise ValueError(
                f"未知填料类型: {packing_name}。"
                f"可用类型: {list(_PACKING_DATA.keys())}"
            )
        Fp = packing["Fp"]

    # 判断输入是质量流速还是质量流量
    # 如果 G_mass 或 L_mass 数值较大，可能是质量流量，需转换为流速
    A_tower = math.pi * (D_tower / 2) ** 2

    # 统一处理：如果 G_mass / A_tower 与 G_mass 差别大，说明是流量
    # 简化处理：直接按质量流速处理（更符合工程惯例）
    G_flux = G_mass  # kg/(m²·s)
    L_flux = L_mass  # kg/(m²·s)

    # 流动参数 X = (L/G) * sqrt(ρ_g / ρ_l)
    if G_flux > 0:
        X = (L_flux / G_flux) * math.sqrt(rho_gas / rho_liquid)
    else:
        X = 0.01
    X = max(X, 1e-4)

    # 泛点容量参数 Y_f（Eckert GPDC）
    # log10(Y_f) = -1.668 - 1.085*log10(X) - 0.297*(log10(X))^2
    log_X = math.log10(X)
    log_Y_flood = -1.668 - 1.085 * log_X - 0.297 * log_X ** 2
    Y_flood = 10 ** log_Y_flood

    # 泛点气速
    viscosity_correction = (mu_liquid / mu_w) ** 0.2
    u_flood = math.sqrt(
        Y_flood * g * rho_liquid / (Fp * rho_gas * viscosity_correction)
    )

    # 实际气速
    u_actual = G_flux / rho_gas if rho_gas > 0 else 0

    # 液泛率
    flooding_ratio = u_actual / u_flood if u_flood > 0 else float('inf')

    # 安全评估
    if flooding_ratio < 0.5:
        safety_level = "very_safe"
        recommendation = "操作点远离液泛点，可考虑提升处理量"
    elif flooding_ratio < 0.7:
        safety_level = "safe"
        recommendation = "操作点在安全范围内"
    elif flooding_ratio < 0.8:
        safety_level = "acceptable"
        recommendation = "接近液泛点上限，注意监测"
    elif flooding_ratio < 0.95:
        safety_level = "risky"
        recommendation = "液泛率偏高，建议降低气相负荷"
    else:
        safety_level = "dangerous"
        recommendation = "已超过液泛点，必须立即降负荷"

    # 提升 20% 处理量后的液泛率
    flooding_ratio_120 = flooding_ratio * 1.2
    if flooding_ratio_120 < 0.8:
        can_increase_20 = True
    else:
        can_increase_20 = False

    # 操作点压降估算
    if flooding_ratio < 0.5:
        dp_per_m = 150.0
    elif flooding_ratio < 0.7:
        dp_per_m = 300.0
    elif flooding_ratio < 0.9:
        dp_per_m = 600.0
    else:
        dp_per_m = 1000.0

    return {
        "D_tower_m": D_tower,
        "u_actual_m_s": u_actual,
        "u_flood_m_s": u_flood,
        "flooding_ratio": flooding_ratio,
        "flooding_ratio_percent": flooding_ratio * 100,
        "flooding_ratio_120_percent": flooding_ratio_120 * 100,
        "can_increase_20": can_increase_20,
        "safety_level": safety_level,
        "recommendation": recommendation,
        "flow_parameter_X": X,
        "capacity_parameter_Y_flood": Y_flood,
        "Fp": Fp,
        "packing_name": packing_name,
        "dp_per_m_estimate_Pa": dp_per_m,
        "G_flux": G_flux,
        "L_flux": L_flux,
    }