"""
反应器工具包 — 固定床、流化床、CSTR 的核心计算
对应研究内容一：固定床反应器设计、流化床反应器设计

包含：
- 反应速率计算、转化率求解
- 固定床 (PFR/PBR) 设计（单筒与列管式）
- 流化床 (FBR) 设计（鼓泡/湍流流化）
- CSTR 设计（单个与串联）
- 传热、传质、压降计算
- 约束校验与自动优化建议
"""
from __future__ import annotations
import numpy as np
import math
from typing import Optional, List, Dict, Tuple

from .common import _REQUIRED, check_required_params


# ============================================================
# 反应动力学与热力学
# ============================================================

def reaction_rate(
    C_A: float,
    k: float,
    n: float = 1.0,
    reversible: bool = False,
    K_eq: Optional[float] = None,
    C_B: float = 0.0,
) -> float:
    """反应速率计算 (r = k * C_A^n)

    Args:
        C_A: 反应物A浓度 (mol/m³)
        k: 速率常数
        n: 反应级数
        reversible: 是否可逆反应
        K_eq: 平衡常数
        C_B: 产物B浓度

    Returns:
        反应速率 r (mol/(m³·s))
    """
    # 保护 C_A 免于因数值微调产生微小负数导致幂运算报错
    C_A_safe = max(0.0, C_A)
    if reversible and K_eq is not None:
        k_r = k / K_eq
        return k * (C_A_safe ** n) - k_r * C_B
    else:
        return k * (C_A_safe ** n)


def reaction_enthalpy(
    delta_H_rxn: float,
    conversion: float,
    molar_flow: float,
) -> float:
    """反应放热量计算（单一反应，简化接口）
    """
    if delta_H_rxn is None:
        return 0.0
    return -delta_H_rxn * molar_flow * conversion


def reaction_enthalpy_multi(
    species: List[Dict],
    T: float = 298.15,
    basis_flow: Optional[float] = None,
    basis_species: Optional[str] = None,
    conversion: float = 1.0,
) -> dict:
    """多组分反应热计算（支持多反应物/生成物）
    """
    T_ref = 298.15
    R_gas = 8.314

    if not species:
        raise ValueError("species 列表不能为空")

    # ── 1. 标准反应焓 (298.15 K) ──────────────────────────────
    delta_H_std = sum(sp["nu"] * sp["dHf"] for sp in species)

    # ── 2. Kirchhoff 温度修正 ──────────────────────────────────
    delta_H_T_corr = 0.0
    if abs(T - T_ref) > 0.01:
        for sp in species:
            coeffs = sp.get("Cp_coeffs", None)
            if coeffs is None or len(coeffs) < 1:
                continue  
            a = coeffs[0] if len(coeffs) > 0 else 0.0
            b = coeffs[1] if len(coeffs) > 1 else 0.0
            c = coeffs[2] if len(coeffs) > 2 else 0.0
            d = coeffs[3] if len(coeffs) > 3 else 0.0
            dH_i = (
                a * (T - T_ref)
                + b / 2 * (T**2 - T_ref**2)
                + c / 3 * (T**3 - T_ref**3)
                + d / 4 * (T**4 - T_ref**4)
            )
            delta_H_T_corr += sp["nu"] * dH_i

    delta_H_rxn = delta_H_std + delta_H_T_corr

    # ── 3. 各组分贡献明细 ─────────────────────────────────────
    contributions = []
    for sp in species:
        role = "product" if sp["nu"] > 0 else "reactant"
        contributions.append({
            "name": sp["name"],
            "nu": sp["nu"],
            "dHf": sp["dHf"],
            "contribution_J_mol": sp["nu"] * sp["dHf"],
            "role": role,
        })

    # ── 4. 基准组分确定 & 放热功率 Q ─────────────────────────
    Q = None
    if basis_flow is not None and basis_flow > 0:
        if basis_species is not None:
            basis = next((s for s in species if s["name"] == basis_species), None)
            if basis is None:
                raise ValueError(f"basis_species '{basis_species}' 不在 species 列表中")
            nu_basis = abs(basis["nu"])
        else:
            reactants = [s for s in species if s["nu"] < 0]
            if not reactants:
                nu_basis = 1.0
            else:
                nu_basis = max(abs(s["nu"]) for s in reactants)

        Q = -delta_H_rxn * (basis_flow * conversion) / nu_basis

    if delta_H_rxn is None:
        reaction_type = "unknown"
    elif delta_H_rxn < -100:
        reaction_type = "exothermic"
    elif delta_H_rxn > 100:
        reaction_type = "endothermic"
    else:
        reaction_type = "thermoneutral"

    return {
        "delta_H_rxn": delta_H_rxn,
        "delta_H_rxn_kJ": delta_H_rxn / 1000.0 if delta_H_rxn is not None else None,
        "Q": Q,
        "T_ref": T_ref,
        "T": T,
        "delta_H_std": delta_H_std,
        "delta_H_T_correction": delta_H_T_corr,
        "species_contributions": contributions,
        "reaction_type": reaction_type,
    }


def adiabatic_temperature_rise(
    delta_H_rxn: float,
    conversion: float,
    F_A0: float,
    Cp_mix: float,
) -> float:
    """绝热温升计算"""
    if delta_H_rxn is None or Cp_mix is None or Cp_mix <= 0:
        return 0.0
    dT = (-delta_H_rxn * conversion) / Cp_mix
    return dT


def estimate_reactor_U(reactor_type: str, phase: str = "liquid",
                       k_fluid: float = None, rho_gas: float = None,
                       Cp_gas: float = None, mu: float = None,
                       u_s: float = None, d_particle: float = None) -> dict:
    """估算反应器总传热系数 U (W/m²·K)
    """
    rtype = reactor_type.lower().replace(" ", "_")

    if rtype == "fluidized_bed":
        U = 200.0
        dominant = "床层-壁面颗粒对流"
        method = "流化床典型范围 100-400 W/m²·K，取中值"
    elif rtype == "cstr":
        if phase in ("gas", "vapor"):
            U = 80.0
            dominant = "气侧对流"
            method = "CSTR气相典型范围 30-150 W/m²·K"
        else:
            U = 350.0
            dominant = "液体侧对流"
            method = "CSTR液相搅拌典型范围 200-500 W/m²·K"
    else:  
        if all(v is not None for v in [k_fluid, rho_gas, mu, u_s, d_particle, Cp_gas]):
            try:
                Re = rho_gas * u_s * d_particle / mu
                Pr = Cp_gas * mu / k_fluid
                Nu_w = 0.17 * Re**0.79 * Pr**0.33
                h_w = Nu_w * k_fluid / d_particle
                h_cool = 500.0
                k_wall, t_wall = 50.0, 0.005
                U = 1.0 / (1.0/h_w + t_wall/k_wall + 1.0/h_cool)
                dominant = f"床层侧对流 h_w={h_w:.1f} W/m²·K (Re={Re:.0f}, Pr={Pr:.2f})"
                method = "Wakao壁面关联式 Nu_w=0.17·Re^0.79·Pr^0.33 + 串联热阻"
            except Exception:
                U = 100.0
                dominant = "默认经验值（物性计算失败）"
                method = "固定床/PFR 典型范围 50-200 W/m²·K"
        else:
            if phase in ("gas", "vapor"):
                U = 100.0
                dominant = "气侧对流（未提供流速/物性）"
                method = "固定床/PFR 气相典型范围 50-200 W/m²·K"
            else:
                U = 300.0
                dominant = "液侧对流（未提供流速/物性）"
                method = "固定床/PFR 液相典型范围 100-500 W/m²·K"

    return {
        "U": U,
        "dominant_resistance": dominant,
        "method": method,
        "note": f"提供 k_fluid/rho_gas/mu/u_s/d_particle/Cp_gas 可获取更精确的 U 估算",
    }


def selectivity_analysis(
    C_A0: float,
    C_D: float,
    C_U: float,
) -> dict:
    """选择性分析"""
    conv = (C_A0 - (C_D + C_U)) / C_A0 if C_A0 > 0 else 0.0
    select = C_D / (C_D + C_U) if (C_D + C_U) > 0 else 0.0
    yield_d = C_D / C_A0 if C_A0 > 0 else 0.0
    return {
        "selectivity": select,
        "yield": yield_d,
        "conversion": conv,
    }


def calc_concentration(
    mass_flow: Optional[float] = None,
    molar_flow: Optional[float] = None,
    volumetric_flow: Optional[float] = None,
    MW: Optional[float] = None,
    density: Optional[float] = None,
    P: Optional[float] = None,
    T: Optional[float] = None,
    phase: str = "liquid",
) -> dict:
    """通过流量、分子量、密度计算浓度
    """
    R_gas = 8.314  

    MW_kg = None
    if MW is not None:
        MW_kg = MW / 1000.0 if MW > 1.0 else MW  

    phase = phase.lower()

    if molar_flow is not None and volumetric_flow is not None and volumetric_flow > 0:
        C = molar_flow / volumetric_flow
        method = "C = ṅ / Q_v"
        if MW_kg is not None:
            mass_flow = mass_flow or (molar_flow * MW_kg)
            density = density or (mass_flow / volumetric_flow)

    elif mass_flow is not None and MW_kg is not None and volumetric_flow is not None and volumetric_flow > 0:
        molar_flow = mass_flow / MW_kg
        C = molar_flow / volumetric_flow
        method = "C = (ṁ / MW) / Q_v"
        density = density or (mass_flow / volumetric_flow)

    elif mass_flow is not None and MW_kg is not None and density is not None and density > 0:
        volumetric_flow = mass_flow / density
        molar_flow = mass_flow / MW_kg
        C = molar_flow / volumetric_flow  
        method = "C = ρ / MW  (Q_v = ṁ / ρ)"

    elif molar_flow is not None and MW_kg is not None and density is not None and density > 0:
        mass_flow = mass_flow or (molar_flow * MW_kg)
        volumetric_flow = mass_flow / density
        C = molar_flow / volumetric_flow
        method = "C = ṅ / Q_v  (Q_v = ṁ / ρ,  ṁ = ṅ × MW)"

    elif phase == "gas" and P is not None and T is not None and T > 0:
        C = P / (R_gas * T)
        method = "C = P / (RT)  [ideal gas]"
        if molar_flow is not None:
            volumetric_flow = molar_flow / C if C > 0 else None

    else:
        raise ValueError(
            "参数不足，无法计算浓度。\n"
            "液相至少需要以下组合之一：\n"
            "  (a) molar_flow + volumetric_flow\n"
            "  (b) mass_flow + MW + volumetric_flow\n"
            "  (c) mass_flow + MW + density\n"
            "  (d) molar_flow + MW + density\n"
            "气相需要：phase='gas', P, T（以及可选 molar_flow）"
        )

    return {
        "C": C,
        "C_kmol_m3": C / 1000.0,
        "C_mol_L": C / 1000.0,
        "molar_flow": molar_flow,
        "volumetric_flow": volumetric_flow,
        "mass_flow": mass_flow,
        "density": density,
        "MW_kg_mol": MW_kg,
        "method": method,
    }


def calc_concentration_multi(
    components: List[Dict],
    total_volumetric_flow: Optional[float] = None,
    total_molar_flow: Optional[float] = None,
    total_mass_flow: Optional[float] = None,
    P: Optional[float] = None,
    T: Optional[float] = None,
    phase: str = "liquid",
) -> dict:
    """多组分混合物各组分浓度计算
    """
    R_gas = 8.314
    phase = phase.lower()

    comps = []
    for raw in components:
        c = dict(raw)
        mw = c.get("MW")
        if mw is not None:
            c["MW"] = mw / 1000.0 if mw > 1.0 else mw  
        comps.append(c)

    if phase == "gas" and P is not None and T is not None and T > 0:
        C_total = P / (R_gas * T)
        total_molar_flow_calc = total_molar_flow  

        results = []
        for c in comps:
            name = c.get("name", "?")
            y_i = c.get("mole_fraction")
            n_i = c.get("molar_flow")
            MW_i = c.get("MW")

            if n_i is None and y_i is not None and total_molar_flow is not None:
                n_i = y_i * total_molar_flow

            if y_i is not None:
                C_i = y_i * C_total
            elif n_i is not None and total_volumetric_flow is not None and total_volumetric_flow > 0:
                C_i = n_i / total_volumetric_flow
                y_i = C_i / C_total if C_total > 0 else None
            else:
                C_i = None

            results.append({
                "name": name,
                "C": C_i,
                "C_mol_L": C_i / 1000.0 if C_i is not None else None,
                "mole_fraction": y_i,
                "molar_flow": n_i,
                "MW_kg_mol": MW_i,
            })

        if total_volumetric_flow is None and total_molar_flow is not None and C_total > 0:
            total_volumetric_flow = total_molar_flow / C_total

        return {
            "components": results,
            "C_total": C_total,
            "total_molar_flow": total_molar_flow_calc,
            "total_volumetric_flow": total_volumetric_flow,
            "phase": phase,
        }

    for c in comps:
        mw = c.get("MW")
        if c.get("molar_flow") is not None:
            continue  
        if c.get("mole_fraction") is not None and total_molar_flow is not None:
            c["molar_flow"] = c["mole_fraction"] * total_molar_flow
        elif c.get("mass_flow") is not None and mw is not None and mw > 0:
            c["molar_flow"] = c["mass_flow"] / mw
        elif (c.get("mass_fraction") is not None and total_mass_flow is not None
              and mw is not None and mw > 0):
            c["molar_flow"] = (c["mass_fraction"] * total_mass_flow) / mw

    total_molar_flow_calc = total_molar_flow  
    if total_molar_flow_calc is None:
        _flows = [c["molar_flow"] for c in comps if c.get("molar_flow") is not None]
        if _flows:
            total_molar_flow_calc = sum(_flows)

    if total_volumetric_flow is None and total_molar_flow_calc is not None:
        _mass_flows = []
        for c in comps:
            n_i = c.get("molar_flow")
            mw_i = c.get("MW")
            if n_i is not None and mw_i is not None:
                _mass_flows.append(n_i * mw_i)
        if len(_mass_flows) == len(comps) and total_mass_flow is None:
            total_mass_flow_est = sum(_mass_flows)
            rho_est = 900.0
            total_volumetric_flow = total_mass_flow_est / rho_est
            _vflow_estimated = True
        else:
            _vflow_estimated = False
    else:
        _vflow_estimated = False

    results = []
    for c in comps:
        name = c.get("name", "?")
        n_i = c.get("molar_flow")
        MW_i = c.get("MW")
        y_i = c.get("mole_fraction")

        if n_i is not None and total_molar_flow_calc and total_molar_flow_calc > 0:
            y_i = n_i / total_molar_flow_calc

        if n_i is not None and total_volumetric_flow is not None and total_volumetric_flow > 0:
            C_i = n_i / total_volumetric_flow
        elif y_i is not None and total_molar_flow_calc and total_volumetric_flow is not None and total_volumetric_flow > 0:
            n_i = y_i * total_molar_flow_calc
            C_i = n_i / total_volumetric_flow
        else:
            C_i = None

        results.append({
            "name": name,
            "C": C_i,
            "C_mol_L": C_i / 1000.0 if C_i is not None else None,
            "mole_fraction": y_i,
            "molar_flow": n_i,
            "MW_kg_mol": MW_i,
        })

    C_total = (
        total_molar_flow_calc / total_volumetric_flow
        if total_volumetric_flow and total_molar_flow_calc
        else None
    )

    ret = {
        "components": results,
        "C_total": C_total,
        "total_molar_flow": total_molar_flow_calc,
        "total_volumetric_flow": total_volumetric_flow,
        "phase": phase,
    }
    if _vflow_estimated:
        ret["volumetric_flow_note"] = (
            "total_volumetric_flow 未提供，已用各组分 MW 汇总质量流量并假设混合液密度 "
            "900 kg/m³ 进行估算，浓度结果仅供参考，建议传入实测体积流量或密度。"
        )
    return ret


# ============================================================
# 固定床反应器 (PFR/PBR) — 核心计算
# ============================================================

def pfr_design(
    F_A0: float,
    X_target: float,
    k: float,
    C_A0: Optional[float] = None,
    n: float = 1.0,
    epsilon: float = 0.0,
    C_B0: Optional[float] = None,  # 双分子反应支持
) -> dict:
    """PFR/PBR 设计方程积分
    
    支持两种模式:
    1. 单组分: r = k * C_A^n (默认, C_B0=None)
    2. 双分子: r = k * C_A * C_B (C_B0 提供时)
    """
    if X_target <= 0 or F_A0 <= 0 or k <= 0:
        return {"V": 0.0, "X": X_target, "tau": 0.0}

    X_target = min(X_target, 0.999)
    if C_A0 is None or C_A0 <= 0:
        C_A0 = 1.0

    # ── 双分子反应模式: r = k * C_A * C_B ──
    if C_B0 is not None and C_B0 > 0:
        M = C_B0 / C_A0  # 浓度比
        
        # 转化率上限受限于 B 的化学计量
        X_max = min(M, 1.0) if M < 1 else 1.0
        if X_target >= X_max:
            X_target = X_max * 0.999
        
        if abs(epsilon) < 1e-10:  # 液相反应，解析解
            if abs(M - 1.0) < 1e-6:
                # M = 1 时: V = F_A0 / (k * C_A0²) * X / (1-X)
                V = (F_A0 / (k * C_A0**2)) * (X_target / (1 - X_target))
            else:
                # M ≠ 1 时: V = F_A0 / (k * C_A0² * (M-1)) * ln[(M-X)/(M*(1-X))]
                V = (F_A0 / (k * C_A0**2 * (M - 1))) * math.log((M - X_target) / (M * (1 - X_target)))
            tau = V * C_A0 / F_A0 if F_A0 > 0 else 0.0
        else:
            # 气相反应，数值积分
            # r = k * C_A0² * (1-X)(M-X) / (1+εX)²
            N = 1000
            X_arr = np.linspace(0, X_target, N + 1)
            dX = X_target / N
            integrand = []
            for x in X_arr:
                if x >= 0.9999:
                    integrand.append(0.0)
                else:
                    CA = C_A0 * (1 - x) / (1 + epsilon * x) if (1 + epsilon * x) > 0 else 0
                    CB = (C_B0 - C_A0 * x) / (1 + epsilon * x) if (1 + epsilon * x) > 0 else 0
                    if CA > 0 and CB > 0:
                        r = k * CA * CB
                        integrand.append(1.0 / r if r > 0 else 0.0)
                    else:
                        integrand.append(0.0)
            integral = dX * (0.5 * integrand[0] + sum(integrand[1:-1]) + 0.5 * integrand[-1])
            V = F_A0 * integral
            tau = V * C_A0 / F_A0 if F_A0 > 0 else 0.0
        
        return {"V": V, "X": X_target, "tau": tau, "bimolecular": True, "M_ratio": M}

    # ── 单组分模式: r = k * C_A^n ──
    # 梯形法则数值积分
    N = 1000
    X_arr = np.linspace(0, X_target, N + 1)
    dX = X_target / N

    if abs(epsilon) < 1e-10:  # 液相反应，解析解
        if abs(n - 1) < 1e-6:
            V = (F_A0 / (k * C_A0)) * (-np.log(1 - X_target))
            tau = V * C_A0 / F_A0 if F_A0 > 0 else 0.0
        else:
            factor = 1.0 / ((n - 1) * (C_A0 ** n))
            term = 1.0 / ((1 - X_target) ** (n - 1)) - 1.0
            V = (F_A0 / k) * factor * term
            tau = V * C_A0 / F_A0 if F_A0 > 0 else 0.0
    else:
        # 气相反应，数值积分
        integrand_full = [((1 + epsilon * x) / (1 - x)) ** n if x < 0.9999 else 0.0 for x in X_arr]
        integral = dX * (0.5 * integrand_full[0] + sum(integrand_full[1:-1]) + 0.5 * integrand_full[-1])
        V = (F_A0 / (k * (C_A0 ** n))) * integral
        tau = V * C_A0 / F_A0 if F_A0 > 0 else 0.0

    return {"V": V, "X": X_target, "tau": tau}


def pfr_heat_transfer(
    delta_H_rxn: float,
    F_A0: float,
    X: float,
    U: float,
    A: float,
    T_in: float,
    T_wall: float,
    Cp_mix: float,
    T_coolant: Optional[float] = None,
) -> dict:
    """PFR 传热计算
    """
    if T_coolant is None:
        T_coolant = T_wall

    Q_gen = abs(delta_H_rxn) * F_A0 * X
    T_adia = T_in + adiabatic_temperature_rise(delta_H_rxn, X, F_A0, Cp_mix)

    if U > 0 and A > 0:
        dT_hot = T_adia - T_coolant
        dT_cold = T_in - T_coolant
        if dT_hot > 0.1 and dT_cold > 0.1 and abs(dT_hot - dT_cold) > 0.01:
            dT_lm = (dT_hot - dT_cold) / math.log(dT_hot / dT_cold)
        elif dT_hot > 0:
            dT_lm = dT_hot
        else:
            dT_lm = 0.0

        Q_removed = U * A * dT_lm
        Q_net = Q_gen - Q_removed  

        T_out = T_in + Q_net / (F_A0 * Cp_mix) if (F_A0 * Cp_mix) > 0 else T_in
        A_required = Q_gen / (U * dT_lm) if (U * dT_lm) > 0 else float('inf')
        NTU = U * A / (F_A0 * Cp_mix) if (F_A0 * Cp_mix) > 0 else 0.0
    else:
        Q_removed = 0.0
        T_out = T_adia
        A_required = float('inf')
        NTU = 0.0

    return {
        "Q_gen": Q_gen,
        "Q_removed": Q_removed,
        "T_out": T_out,
        "T_adia": T_adia,
        "A_heat_transfer": A,
        "A_required_adiabatic": A_required,
        "NTU": NTU,
        "delta_T_wall": T_wall - T_adia if T_wall > T_adia else 0,
        "cooling_ratio": Q_removed / Q_gen if Q_gen > 0 else 0.0,
    }


def multitubular_design(
    F_A0: float,
    X_target: float,
    k: float,
    C_A0: float,
    n: float,
    d_tube_inner: float,
    L_tube: float,
    U: float,
    T_coolant: float,
    rho_cat: float,
    W_cat_total: Optional[float] = None,
) -> dict:
    """列管式固定床设计"""
    V_single_tube = math.pi * (d_tube_inner ** 2) / 4 * L_tube

    N_tubes_low = 1
    N_tubes_high = 100000
    N_tubes = 1

    for _ in range(50):
        N_tubes = (N_tubes_low + N_tubes_high) // 2
        F_per_tube = F_A0 / N_tubes
        result = pfr_design(F_per_tube, X_target, k, C_A0, n, epsilon=0.0)

        if result["V"] <= V_single_tube:
            N_tubes_high = N_tubes
        else:
            N_tubes_low = N_tubes

        if N_tubes_high - N_tubes_low <= 1:
            break

    N_tubes = N_tubes_high
    F_per_tube = F_A0 / N_tubes
    result = pfr_design(F_per_tube, X_target, k, C_A0, n, epsilon=0.0)

    V_total = N_tubes * V_single_tube
    A_total = N_tubes * math.pi * d_tube_inner * L_tube

    if W_cat_total is None:
        W_cat_total = V_total * rho_cat

    return {
        "N_tubes": N_tubes,
        "L_tube": L_tube,
        "d_tube_inner": d_tube_inner,
        "V_total": V_total,
        "A_total": A_total,
        "W_cat_total": W_cat_total,
        "Q_total": V_total * abs(k * C_A0),
        "X_estimated": X_target,
        "per_tube_volume": V_single_tube,
        "tau_per_tube": result.get("tau", 0),
    }


def fixed_bed_outlet_conditions(
    F_A0: float,
    X: float,
    T_in: float,
    P_in: float,
    delta_H_rxn: float,
    Cp_mix: float,
    epsilon: float = 0.0,
) -> dict:
    """床层出口条件计算"""
    T_out = T_in + (-delta_H_rxn * X) / Cp_mix if Cp_mix > 0 else T_in
    P_out = P_in  
    F_out = F_A0 * (1 + epsilon * X)

    return {
        "T_out": T_out,
        "P_out": P_out,
        "F_out": F_out,
        "X_out": X,
    }


def fixed_bed_catalyst_inventory(
    V_bed: float,
    rho_cat: float,
    particle_diameter: float,
    voidage: float = 0.4,
) -> dict:
    """催化剂库存计算"""
    W_cat = V_bed * rho_cat * (1 - voidage)
    N_particles = (V_bed * (1 - voidage)) / (math.pi * particle_diameter**3 / 6)
    S_BET = 6 / particle_diameter * (1 - voidage) / voidage

    return {
        "W_cat": W_cat,
        "N_particles": N_particles,
        "S_BET": S_BET,
        "V_particles": V_bed * (1 - voidage),
        "V_void": V_bed * voidage,
    }


def estimate_bed_void_fraction(
    D_particle: float,
    D_bed: float = None,
    packing_type: str = "random_spheres",
    sphericity: float = 1.0,
) -> float:
    """床层空隙率估算

    根据颗粒直径、床层直径和填充方式估算固定床的空隙率。

    Args:
        D_particle: 催化剂颗粒直径 (m)
        D_bed: 床层直径 (m)，用于考虑壁面效应
        packing_type: 填充类型
            - "random_spheres": 随机填充球形颗粒 (ε ≈ 0.35-0.45)
            - "ordered_spheres": 规则填充球形颗粒 (ε ≈ 0.26-0.30)
            - "cylinders": 圆柱形颗粒 (ε ≈ 0.35-0.40)
            - "rings": 环形颗粒 (ε ≈ 0.40-0.50)
        sphericity: 球形度 (0-1)，1.0 为完美球形

    Returns:
        epsilon: 床层空隙率 (0-1)

    References:
        - Dixon, A. G. (1988). Correlations for wall and particle shape effects on fixed bed bulk voidage.
        - Benenati, R. F., & Brosilow, C. B. (1962). Void fraction distribution in packed beds.
    """
    if D_particle <= 0:
        raise ValueError("颗粒直径必须大于 0")

    # 基础空隙率（无穷大床层，无壁面效应）
    if packing_type == "random_spheres":
        epsilon_inf = 0.40  # 随机填充球形颗粒
    elif packing_type == "ordered_spheres":
        epsilon_inf = 0.26  # 规则填充（FCC/HCP）
    elif packing_type == "cylinders":
        epsilon_inf = 0.37  # 圆柱形颗粒
    elif packing_type == "rings":
        epsilon_inf = 0.45  # 环形颗粒
    else:
        epsilon_inf = 0.40  # 默认值

    # 球形度修正：非球形颗粒空隙率更高
    if sphericity < 1.0:
        epsilon_inf += 0.1 * (1.0 - sphericity)

    # 壁面效应修正（Dixon 关联式）
    # 当 D_bed/D_particle < 10 时，壁面效应显著
    if D_bed is not None and D_bed > 0:
        ratio = D_bed / D_particle
        if ratio < 50:
            # Dixon (1988): ε = ε_inf + (1 - ε_inf) × 0.25 × (d_p/D)
            # 壁面附近空隙率增大
            epsilon = epsilon_inf + (1.0 - epsilon_inf) * 0.25 / ratio
        else:
            epsilon = epsilon_inf
    else:
        # 无床层直径信息，使用基础值
        epsilon = epsilon_inf

    # 限制在合理范围内
    epsilon = max(0.25, min(0.60, epsilon))

    return epsilon


def fixed_bed_dimensions(V_bed: float, D_bed: float) -> dict:
    """床层几何尺寸"""
    A_bed = math.pi * D_bed**2 / 4
    L_bed = V_bed / A_bed
    H_D_ratio = L_bed / D_bed

    return {
        "D_bed": D_bed,
        "L_bed": L_bed,
        "H_D_ratio": H_D_ratio,
        "A_cross": A_bed,
        "V_bed": V_bed,
    }


def fixed_bed_pressure_drop_ergun(
    L_bed: float,
    D_particle: float,
    epsilon: float,
    rho_gas: float,
    mu: float,
    u_s: float,
    sphericity: float = 0.8,
) -> dict:
    """Ergun方程计算固定床压降"""
    # 修正：根据 Ergun 经典阻力方程，针对非球形颗粒引入有效当量直径 dp_eff 参与粘性项和惯性项压降阻力计算，防球形度参数失效
    dp_eff = sphericity * D_particle
    a_s = 6 / dp_eff
    Re_p = rho_gas * u_s * dp_eff / (mu * (1 - epsilon)) if mu > 0 else 1e6

    term1 = 150 * mu * (1 - epsilon)**2 / (epsilon**3 * dp_eff**2) * u_s
    term2 = 1.75 * rho_gas * (1 - epsilon) / (epsilon**3 * dp_eff) * u_s**2
    dP_dL = term1 + term2

    delta_P = dP_dL * L_bed

    return {
        "delta_P": delta_P,
        "dP_dL": dP_dL,
        "term_viscous": term1 * L_bed,
        "term_kinetic": term2 * L_bed,
        "Re_particle": Re_p,
        "u_s": u_s,
        "epsilon": epsilon,
    }


def fixed_bed_residence_time(
    V_bed: float,
    F_total: float,
    P: float,
    T: float,
    epsilon: float = 0.4,
    R: float = 8.314,
) -> dict:
    """固定床停留时间计算"""
    v_volumetric = F_total * R * T / P
    V_void = V_bed * epsilon
    tau = V_void / v_volumetric if v_volumetric > 0 else 0.0

    return {
        "tau": tau,
        "v_volumetric": v_volumetric,
        "V_void": V_void,
    }


# ============================================================
# 流化床反应器 (FBR) — 核心计算
# ============================================================

def minimum_fluidization_velocity(
    d_particle: float,
    rho_particle: float,
    rho_gas: float,
    mu: float,
    epsilon_mf: float = 0.45,
    sphericity: float = 0.8,
    g: float = 9.81,
) -> dict:
    """最小流化速度计算 (Wen-Yu 方程)"""
    Ar = (d_particle**3 * rho_gas * (rho_particle - rho_gas) * g) / (mu**2) if mu > 0 else 1e-12
    Re_mf = np.sqrt(33.7**2 + 0.0408 * Ar) - 33.7
    u_mf = Re_mf * mu / (rho_gas * d_particle) if (rho_gas * d_particle) > 0 else 0.0

    return {
        "u_mf": u_mf,
        "Re_mf": Re_mf,
        "Ar": Ar,
        "epsilon_mf": epsilon_mf,
    }


def terminal_velocity(
    d_particle: float,
    rho_particle: float,
    rho_gas: float,
    mu: float,
    sphericity: float = 0.8,
    g: float = 9.81,
) -> dict:
    """终端速度计算"""
    u_t = 0.1 

    for _ in range(50):
        Re_t = rho_gas * u_t * d_particle / mu if mu > 0 else 1e6
        # 保护：防止 Re_t 趋于 0 时，Re_t 倒数除零溢出崩溃
        Re_t = max(1e-10, Re_t)
        
        if Re_t <= 1000:
            C_d = 24 / Re_t * (1 + 0.15 * Re_t**0.687)
        else:
            C_d = 0.44
        C_d = C_d / (sphericity**(2/3)) if sphericity > 0 else C_d

        u_t_new = np.sqrt(4 * d_particle * g * (rho_particle - rho_gas) / (3 * rho_gas * C_d))

        if abs(u_t_new - u_t) < 1e-6:
            u_t = u_t_new
            break
        u_t = u_t_new

    Re_t = rho_gas * u_t * d_particle / mu if mu > 0 else 1e6

    return {
        "u_terminal": u_t,
        "Re_terminal": Re_t,
        "C_d": C_d,
    }


def check_fluidization_regime(
    u_actual: float,
    u_mf: float,
    u_t: float,
) -> dict:
    """判断流化状态"""
    u_u_mf = u_actual / u_mf if u_mf > 0 else float("inf")

    if u_actual < u_mf:
        regime = "fixed_bed"
    elif u_u_mf <= 5:
        regime = "bubbling_fluidized"
    elif u_u_mf <= 15:
        regime = "turbulent_fluidized"
    else:
        regime = "fast_fluidized"

    is_bubbling = 1.0 < u_u_mf <= 5.0
    is_turbulent = 5.0 < u_u_mf <= 15.0

    carryover = False
    if u_actual > 0.8 * u_t:
        carryover = True
        regime += " (carryover_risk)"

    return {
        "regime": regime,
        "u_u_mf": u_u_mf,
        "is_bubbling": is_bubbling,
        "is_turbulent": is_turbulent,
        "carryover_risk": carryover,
    }


def fluidization_velocity_limits(u_mf: float, u_t: float) -> dict:
    """流化速度范围"""
    u_bubbling_max = min(5 * u_mf, 0.8 * u_t)
    u_turbulent_max = min(15 * u_mf, 0.8 * u_t)

    return {
        "bubbling_range": (u_mf, u_bubbling_max),
        "turbulent_range": (5 * u_mf, u_turbulent_max),
        "safe_max": 0.8 * u_t,
        "u_mf": u_mf,
        "u_t": u_t,
    }


def fluidized_bed_outlet_conditions(
    F_in: float,
    T_in: float,
    P_in: float,
    conversion: float,
    delta_H_rxn: float,
    Cp_mix: float,
    heat_transfer: Dict,
    epsilon: float = 0.0,
) -> dict:
    """流化床出口条件"""
    T_out = T_in + (-delta_H_rxn * conversion) / Cp_mix if Cp_mix > 0 else T_in

    # 修正：同时读取 Q_removed 和 Q 键名，解决两端接口传参错配导致冷却温降永远不生效的 Bug
    Q_removed = heat_transfer.get("Q_removed", heat_transfer.get("Q", 0))
    if Q_removed > 0 and F_in * Cp_mix > 0:
        T_out -= Q_removed / (F_in * Cp_mix)

    F_out = F_in * (1 + epsilon * conversion) 
    P_out = P_in  

    return {
        "T_out": T_out,
        "P_out": P_out,
        "F_out": F_out,
        "X_out": conversion,
    }


def fluidized_bed_catalyst_inventory(
    D_bed: float,
    H_dense: float,
    rho_cat: float,
    epsilon_mf: float = 0.45,
) -> dict:
    """流化床催化剂质量计算"""
    V_dense = math.pi * D_bed**2 / 4 * H_dense
    W_cat = V_dense * rho_cat * (1 - epsilon_mf)

    return {
        "W_cat": W_cat,
        "V_dense": V_dense,
        "V_expanded": None, 
        "H_dense": H_dense,
    }


def fluidized_bed_residence_time(
    H_dense: float,
    u_mf: float,
    u_actual: float,
) -> dict:
    """流化床停留时间"""
    tau = H_dense / u_mf if u_mf > 0 else 0.0

    return {
        "tau_gas": H_dense / u_actual if u_actual > 0 else float("inf"),
        "tau_particles": tau,
    }


def fluidized_bed_heat_transfer_compat(
    D_bed: float,
    H_dense: float,
    T_bed: float,
    T_wall: float,
    k_gas: float,
    mu: float,
    rho_gas: float,
    u_mf: float,
    d_particle: float,
    Cp_gas: float,
) -> dict:
    """流化床传热系数估算 (兼容包装)"""
    # 修正：自适应摩尔比热容与质量比热容不匹配造成的 Prandtl 数数量级错误，
    # 常见气相介质的 Prandtl 数高度聚集在 0.7 左右，若失真则自适应回退至 0.7 物理常数。
    Pr = Cp_gas * mu / k_gas if k_gas > 0 else 0.7
    if Pr < 0.1 or Pr > 2.0:
        Pr = 0.7
        
    Re_mf = rho_gas * u_mf * d_particle / mu if mu > 0 else 0.0

    Nu = 2 + 0.589 * Re_mf**0.5 * Pr**(1/3) / (1 + 0.0001 * Re_mf**1.25)**0.6
    h_fw = Nu * k_gas / d_particle if d_particle > 0 else 0.0
    A = math.pi * D_bed * H_dense

    Q_max = h_fw * A * (T_bed - T_wall)

    return {
        "h_fw": h_fw,
        "Nu": Nu,
        "Re_mf": Re_mf,
        "Pr": Pr,
        "A": A,
        "Q": Q_max,
    }


# ============================================================
# CSTR (全混釜) — 核心计算
# ============================================================

def cstr_design(
    F_A0: float,
    X_target: float,
    k: float,
    C_A0: float,
    n: float = 1.0,
    num_tanks: int = 1,
    C_B0: Optional[float] = None,
) -> dict:
    """CSTR设计（单个或多个等体积串联）

    支持两种动力学模型:
      1. 单组分近似: r = k * C_A^n  (默认, C_B0=None)
      2. 双分子反应: r = k * C_A * C_B  (传入 C_B0 后自动启用)

    Args:
        F_A0: 关键组分A的摩尔进料 (mol/s)
        X_target: 目标转化率 (0~1)
        k: 反应速率常数
        C_A0: 入口浓度 C_A0 (mol/m3)
        n: 反应级数 (仅单组分模式使用)
        num_tanks: 串联釜数
        C_B0: 反应物B的入口浓度 (mol/m3), 传入后启用双分子动力学 r=k*CA*CB

    Returns:
        设计结果字典
    """
    if X_target <= 0 or X_target >= 1 or F_A0 <= 0 or k <= 0 or C_A0 <= 0:
        return {"V_total": 0.0, "V_per_tank": 0.0, "tau": 0.0, "num_tanks": num_tanks}

    # ── 双分子反应模式: r = k * C_A * C_B ──
    if C_B0 is not None and C_B0 > 0:
        return _cstr_design_bimolecular(F_A0, X_target, k, C_A0, C_B0, num_tanks)

    # ── 单组分模式: r = k * C_A^n ──
    if num_tanks == 1:
        if n == 1:
            tau = X_target / (k * (1 - X_target))
        else:
            tau = (X_target / (k * (C_A0 ** (n - 1)))) / ((1 - X_target) ** n)
            
        r_A = k * (C_A0 * (1 - X_target)) ** n
        V_single = F_A0 * X_target / r_A if r_A > 0 else 0.0

        V_total = V_single
        tau_space = V_total / (F_A0 / C_A0) if (F_A0 / C_A0) > 0 else 0.0

    else:
        if n == 1:
            if X_target >= 1.0:
                X_target = 0.999
            tau = (1.0 / k) * ((1.0 / (1.0 - X_target)) ** (1.0 / num_tanks) - 1.0)
            volumetric_flow = F_A0 / C_A0 
            V_single = tau * volumetric_flow
            V_total = V_single * num_tanks
            tau_space = V_total / volumetric_flow
        else:
            r_A_avg = k * (C_A0 * (1 - X_target)) ** n
            V_total = F_A0 * X_target / r_A_avg
            V_single = V_total / num_tanks
            volumetric_flow = F_A0 / C_A0
            tau_space = V_total / volumetric_flow

    return {
        "V_total": V_total,
        "V_per_tank": V_single,
        "tau": tau_space,
        "num_tanks": num_tanks,
        "C_A_out": C_A0 * (1 - X_target),
    }


def _cstr_design_bimolecular(
    F_A0: float,
    X_target: float,
    k: float,
    C_A0: float,
    C_B0: float,
    num_tanks: int = 1,
) -> dict:
    """双分子二级反应 CSTR 精确设计 (A + B -> P, r = k*CA*CB)

    单釜: 直接代入 CSTR 设计方程
      V = F_A0 * X / (k * C_A_out * C_B_out)
      其中 C_A_out = C_A0*(1-X), C_B_out = C_B0 - C_A0*X

    多釜串联: 逐釜数值求解 (二分法)
      每釜满足: V = F_A0*(X_i - X_{i-1}) / (k * C_A_i * C_B_i)
      等体积约束下通过二分法求公共 V
    """
    v0 = F_A0 / C_A0  # 体积流量 (m3/s)

    # 转化率上限: 受限于 B 的化学计量
    X_max = min(C_B0 / C_A0, 1.0) if C_A0 > 0 else 1.0
    if X_target >= X_max:
        X_target = X_max * 0.999

    # ── 单釜 ──
    if num_tanks == 1:
        CA_out = C_A0 * (1 - X_target)
        CB_out = C_B0 - C_A0 * X_target

        if CA_out <= 0 or CB_out <= 0:
            return {"V_total": 0.0, "V_per_tank": 0.0, "tau": 0.0,
                    "num_tanks": num_tanks, "C_A_out": CA_out,
                    "error": "出口浓度为零或负, 转化率超出极限"}

        r_out = k * CA_out * CB_out
        V = F_A0 * X_target / r_out if r_out > 0 else float('inf')
        tau = V / v0 if v0 > 0 else 0.0

        return {
            "V_total": V,
            "V_per_tank": V,
            "tau": tau,
            "num_tanks": 1,
            "C_A_out": CA_out,
            "C_B_out": CB_out,
            "kinetics": "bimolecular",
        }

    # ── 多釜串联: 等体积逐釜求解 ──
    # 策略: 二分搜索公共体积 V, 使得 N 釜后 X_N >= X_target
    # 对给定 V, 每釜通过二分法求解出口转化率

    def _solve_tank(X_in, V_tank):
        """给定入口转化率和釜体积, 求出口转化率"""
        # 二分搜索 X_out in [X_in, X_max)
        lo, hi = X_in, X_max * 0.9999
        for _ in range(200):
            mid = (lo + hi) / 2
            CA = C_A0 * (1 - mid)
            CB = C_B0 - C_A0 * mid
            if CA <= 0 or CB <= 0:
                hi = mid
                continue
            r = k * CA * CB
            X_required = r * V_tank / F_A0 + X_in
            if X_required > mid:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-12:
                break
        return (lo + hi) / 2

    # 二分搜索公共体积 V
    V_lo, V_hi = 1e-10, F_A0 * 100 / (k * 1e-6)  # 宽泛上下界
    for _ in range(300):
        V_mid = (V_lo + V_hi) / 2
        X = 0.0
        for _t in range(num_tanks):
            X = _solve_tank(X, V_mid)
        if X >= X_target:
            V_hi = V_mid
        else:
            V_lo = V_mid
        if V_hi - V_lo < 1e-12:
            break

    V_per = (V_lo + V_hi) / 2
    V_total = V_per * num_tanks
    tau = V_total / v0 if v0 > 0 else 0.0

    # 计算各釜中间转化率
    X_intermediates = []
    X = 0.0
    for _t in range(num_tanks):
        X = _solve_tank(X, V_per)
        X_intermediates.append(X)

    CA_final = C_A0 * (1 - X)
    CB_final = C_B0 - C_A0 * X

    return {
        "V_total": V_total,
        "V_per_tank": V_per,
        "tau": tau,
        "num_tanks": num_tanks,
        "C_A_out": CA_final,
        "C_B_out": CB_final,
        "X_intermediates": X_intermediates,
        "kinetics": "bimolecular",
    }


# ============================================================
# 催化剂物性估算
# ============================================================

def estimate_catalyst_properties(cat_type: str = _REQUIRED) -> dict:
    """根据催化剂类型返回参考物性数据"""
    missing = check_required_params(locals(), {
        "cat_type": "催化剂类型 (zeolite/alumina/metal)，必须显式指定",
    })
    if missing:
        return missing

    reference_data = {
        "zeolite": {
            "particle_diameter": 0.003,
            "particle_density": 1500.0,
            "bulk_density": 900.0,
            "bed_voidage": 0.40,
            "sphericity": 0.85,
        },
        "alumina": {
            "particle_diameter": 0.005,
            "particle_density": 3600.0,
            "bulk_density": 1800.0,
            "bed_voidage": 0.38,
            "sphericity": 0.95,
        },
        "metal": {
            "particle_diameter": 0.002,
            "particle_density": 8000.0,
            "bulk_density": 4800.0,
            "bed_voidage": 0.35,
            "sphericity": 0.90,
        },
    }

    result = reference_data.get(cat_type)
    if result is None:
        return {
            "error": "unknown_catalyst_type",
            "message": f"未知催化剂类型 '{cat_type}'。支持的类型: {list(reference_data.keys())}",
            "note": "以上为参考值，实际设计请使用测量数据或用户指定的准确值",
        }
    result["_note"] = "以上为文献参考值（数量级参考），实际设计请使用测量值或用户指定值"
    return result


# ============================================================
# 设计验证
# ============================================================

# ============================================================
# 目标转化率可达性约束检查
# ============================================================

def check_conversion_feasibility(
    reactor_type: str,
    X_target: float,
    k: float,
    C_A0: float,
    n: float,
    F_A0: float,
    V: float,
    tau: float,
    C_B0: Optional[float] = None,
    V_max_practical: float = 10000.0,
) -> dict:
    """目标转化率可达性约束检查

    在反应器设计完成后，综合判断目标转化率 X_target 在当前动力学条件下
    是否能够工程实现。检查维度：

    1. 化学计量限制（双分子反应）：X_target 不能超过 B 耗尽对应的极限转化率
    2. 体积爆炸检测：所需反应器体积是否超出工程合理范围
    3. 出口速率衰减：r_exit / r_inlet 过低意味着末段转化需要不成比例的巨大体积
    4. 停留时间合理性：tau 过大说明动力学上难以达到目标转化率

    Args:
        reactor_type: 反应器类型 (pfr/cstr/fixed_bed/...)
        X_target: 目标转化率 (0~1)
        k: 反应速率常数
        C_A0: 进料浓度 (mol/m³)
        n: 反应级数
        F_A0: A 的摩尔进料 (mol/s)
        V: 计算所得反应器体积 (m³)
        tau: 计算所得停留时间 (s)
        C_B0: 反应物B的进料浓度 (mol/m³)，双分子反应时提供
        V_max_practical: 工程上可接受的最大体积 (m³)，默认 10000

    Returns:
        dict: 可达性检查结果
    """
    checks = []
    is_feasible = True
    warnings = []
    x_max_feasible = None  # 当前条件下可达的最大转化率估计

    # ── 检查 1: 双分子反应化学计量限制 ──
    if C_B0 is not None and C_B0 > 0 and C_A0 > 0:
        M = C_B0 / C_A0  # 浓度比
        X_stoich_limit = min(M, 1.0)  # B 耗尽时 A 的最大转化率
        if X_target >= X_stoich_limit:
            is_feasible = False
            x_max_feasible = X_stoich_limit * 0.95  # 建议留 5% 余量
            checks.append({
                "name": "stoichiometric_limit",
                "status": "fail",
                "message": (
                    f"目标转化率 X={X_target:.4f} 超出化学计量极限 "
                    f"X_max=M(C_B0/C_A0)={M:.4f}。"
                    f"B 组分不足以支撑该转化率，请降低 X_target 或增大 C_B0。"
                ),
            })
        elif X_target >= 0.95 * X_stoich_limit:
            is_feasible = False
            checks.append({
                "name": "stoichiometric_limit",
                "status": "warn",
                "message": (
                    f"目标转化率 X={X_target:.4f} 接近化学计量极限 "
                    f"X_max={X_stoich_limit:.4f}（已达 {X_target/X_stoich_limit*100:.1f}%），"
                    f"所需体积将急剧增大，建议留有余量。"
                ),
            })
            x_max_feasible = X_stoich_limit * 0.95
        else:
            checks.append({
                "name": "stoichiometric_limit",
                "status": "pass",
                "message": f"化学计量余量充足：X_target={X_target:.4f} < X_max={X_stoich_limit:.4f}",
            })

    # ── 检查 2: 体积爆炸检测 ──
    if V <= 0:
        is_feasible = False
        checks.append({
            "name": "volume_feasibility",
            "status": "fail",
            "message": f"计算体积 V={V:.2f} m³ 无效（≤0），转化率不可达。",
        })
    elif V > V_max_practical:
        is_feasible = False
        checks.append({
            "name": "volume_feasibility",
            "status": "fail",
            "message": (
                f"所需体积 V={V:.1f} m³ 超出工程上限 {V_max_practical:.0f} m³。"
                f"目标转化率 X={X_target:.4f} 在当前动力学条件下不可行，"
                f"建议：提高反应温度(增大k)、降低 X_target、或采用多级反应器。"
            ),
        })
    elif V > V_max_practical * 0.5:
        checks.append({
            "name": "volume_feasibility",
            "status": "warn",
            "message": (
                f"所需体积 V={V:.1f} m³ 偏大（接近工程上限 {V_max_practical:.0f} m³），"
                f"转化率可达但经济性较差。"
            ),
        })
    else:
        checks.append({
            "name": "volume_feasibility",
            "status": "pass",
            "message": f"反应器体积 V={V:.2f} m³ 在合理范围内",
        })

    # ── 检查 3: 出口速率衰减比 ──
    r_inlet = k * (C_A0 ** n) if C_A0 > 0 and k > 0 else 0.0
    C_A_exit = C_A0 * (1 - X_target)
    r_exit = k * (max(C_A_exit, 0) ** n) if C_A_exit > 0 and k > 0 else 0.0

    if r_inlet > 0:
        rate_ratio = r_exit / r_inlet
        # rate_ratio = (1 - X_target)^n
        if rate_ratio < 1e-4:
            is_feasible = False
            checks.append({
                "name": "rate_decay",
                "status": "fail",
                "message": (
                    f"出口反应速率仅为入口的 {rate_ratio:.2e} 倍（r_exit/r_inlet），"
                    f"高转化率下反应几乎停滞。"
                    f"当前 X_target={X_target:.4f} 在动力学上极难达到，"
                    f"建议降低转化率或提高反应温度。"
                ),
            })
            # 估算更合理的最大转化率：使 r_exit/r_inlet >= 0.01
            # (1-X)^n >= 0.01 => X <= 1 - 0.01^(1/n)
            if n > 0:
                x_max_rate = 1.0 - 0.01 ** (1.0 / n)
                if x_max_feasible is None or x_max_rate < x_max_feasible:
                    x_max_feasible = x_max_rate
        elif rate_ratio < 0.01:
            checks.append({
                "name": "rate_decay",
                "status": "warn",
                "message": (
                    f"出口反应速率降至入口的 {rate_ratio:.4f} 倍，"
                    f"末段反应推动力极低，体积效率差。"
                    f"建议考虑适当降低 X_target 或采用分段进料。"
                ),
            })
        else:
            checks.append({
                "name": "rate_decay",
                "status": "pass",
                "message": f"出口速率保持良好：r_exit/r_inlet = {rate_ratio:.4f}",
            })

    # ── 检查 4: 停留时间合理性 ──
    # 工业典型范围：液相 10s~10000s，气相 0.1s~100s
    # 这里用宽松阈值：tau > 100000s (~27h) 视为不合理
    tau_max = 100000.0
    if tau > tau_max:
        is_feasible = False
        checks.append({
            "name": "residence_time",
            "status": "fail",
            "message": (
                f"所需停留时间 tau={tau:.0f} s（{tau/3600:.1f} h）过长，"
                f"超出工程合理范围。转化率 X={X_target:.4f} 在动力学上不可行。"
            ),
        })
    elif tau > tau_max * 0.1:
        checks.append({
            "name": "residence_time",
            "status": "warn",
            "message": (
                f"停留时间 tau={tau:.0f} s（{tau/3600:.1f} h）偏长，"
                f"请确认是否可接受。"
            ),
        })
    else:
        checks.append({
            "name": "residence_time",
            "status": "pass",
            "message": f"停留时间 tau={tau:.2f} s 在合理范围内",
        })

    # ── 汇总 ──
    n_fail = sum(1 for c in checks if c["status"] == "fail")
    n_warn = sum(1 for c in checks if c["status"] == "warn")

    if x_max_feasible is not None:
        x_max_feasible = round(min(x_max_feasible, 0.999), 4)

    suggestion = None
    if not is_feasible:
        suggestion = (
            f"目标转化率 X_target={X_target:.4f} 在当前条件下不可达。"
        )
        if x_max_feasible is not None:
            suggestion += f"建议将 X_target 降至 ≤ {x_max_feasible:.4f}，"
        suggestion += "或尝试：(1) 提高反应温度以增大 k；(2) 增大进料浓度 C_A0；(3) 采用多级反应器串联。"

    return {
        "is_feasible": is_feasible,
        "checks": checks,
        "n_fail": n_fail,
        "n_warn": n_warn,
        "X_target": X_target,
        "X_max_feasible": x_max_feasible,
        "r_exit_r_inlet_ratio": round(r_exit / r_inlet, 6) if r_inlet > 0 else None,
        "V_m3": V,
        "tau_s": tau,
        "suggestion": suggestion,
    }


def validate_reactor_design(design_params: dict, results: dict) -> dict:
    """反应器设计约束校验"""
    checks = []
    n_total = 0
    n_pass = 0

    reactor_type = design_params.get("reactor_type", "fixed_bed")
    
    V = results.get("V", results.get("V_total", 0))
    n_total += 1
    if V > 0:
        n_pass += 1
        checks.append({"name": "volume_positive", "status": "pass", "message": f"体积 V={V:.2f} m³ > 0"})
    else:
        checks.append({"name": "volume_positive", "status": "fail", "message": f"体积 V={V:.2f} m³ 无效"})

    delta_P = results.get("delta_P", results.get("pressure_drop", 0))
    n_total += 1
    max_dP = design_params.get("max_pressure_drop", 150000) 
    if delta_P <= max_dP:
        n_pass += 1
        checks.append({"name": "pressure_drop", "status": "pass", "message": f"压降 {delta_P/1000:.1f} kPa <= {max_dP/1000:.0f} kPa"})
    else:
        checks.append({"name": "pressure_drop", "status": "fail", "message": f"压降 {delta_P/1000:.1f} kPa 过高"})

    T_out = results.get("T_out", results.get("T_exit", 0))
    n_total += 1
    if 273 < T_out < 1273: 
        n_pass += 1
        checks.append({"name": "temperature_range", "status": "pass", "message": f"出口温度 {T_out:.1f} K 在合理范围"})
    else:
        checks.append({"name": "temperature_range", "status": "warn" if T_out != 0 else "fail", "message": f"出口温度 {T_out:.1f} K 可能不合理"})
        if T_out != 0:
            n_pass += 0.5 

    if reactor_type == "fluidized_bed":
        u_ratio = results.get("u_u_mf", 0)
        n_total += 1
        if 2.0 <= u_ratio <= 15.0:
            n_pass += 1
            checks.append({"name": "fluidization_regime", "status": "pass", "message": f"u/u_mf={u_ratio:.1f} 在推荐范围内"})
        else:
            checks.append({"name": "fluidization_regime", "status": "warn", "message": f"u/u_mf={u_ratio:.1f} 可能不在最优流化范围"})
            n_pass += 0.5

    L_D = results.get("L_D", results.get("H_D_ratio", 0))
    if L_D > 0:
        n_total += 1
        if 1.0 <= L_D <= 20.0:
            n_pass += 1
            checks.append({"name": "L_D_ratio", "status": "pass", "message": f"长径比 L/D={L_D:.1f} 合理"})
        else:
            checks.append({"name": "L_D_ratio", "status": "warn", "message": f"长径比 L/D={L_D:.1f} 过大或过小"})
            n_pass += 0.5

    pass_rate = n_pass / n_total if n_total > 0 else 0
    critical_failures = [c["name"] for c in checks if c["status"] == "fail" and c["name"] in ["volume_positive", "pressure_drop"]]

    return {
        "checks": checks,
        "n_pass": n_pass,
        "n_total": n_total,
        "pass_rate": pass_rate,
        "critical_failures": critical_failures,
        "is_valid": pass_rate >= 0.8 and len(critical_failures) == 0,
    }


# ============================================================
# 主设计函数及高阶物化推演
# ============================================================

def catalyst_deactivation(
    time_on_stream: float,
    deactivation_model: str = "exponential",
    k0: float = 1.0,
    kd: float = 0.01,
    order: float = 1.0,
) -> dict:
    """Catalyst deactivation models"""
    import math
    if deactivation_model == "exponential":
        kt = k0 * math.exp(-kd * time_on_stream)
    elif deactivation_model == "power":
        kt = k0 / ((1 + kd * time_on_stream) ** order)
    elif deactivation_model == "linear":
        kt = max(0, k0 * (1 - kd * time_on_stream))
    else:
        kt = k0
    
    return {
        "k_current": float(kt),
        "k0": float(k0),
        "relative_activity": float(kt / k0 if k0 > 0 else 1.0),
        "time_on_stream": time_on_stream,
        "model": deactivation_model,
    }


def calculate_thiele_modulus(
    k: float,
    D_eff: float,
    L_c: float,
    n: float = 1.0,
) -> float:
    """Calculate Thiele modulus (Default to first order equivalence)"""
    import math
    phi = L_c * math.sqrt(k / D_eff)
    return float(phi)


def effectiveness_factor(
    Thiele_modulus: float,
    pellet_shape: str = "sphere",
) -> float:
    """Calculate effectiveness factor eta"""
    import math
    phi = Thiele_modulus
    if phi < 1e-6:
        return 1.0
    
    if pellet_shape == "sphere":
        if phi < 0.5:
            eta = 1.0 / (1.0 + phi**2 / 9.0)
        else:
            eta = (3.0 / phi) * (1.0 / math.tanh(phi) - 1.0 / phi)
    elif pellet_shape == "cylinder":
        if phi < 0.5:
            eta = 1.0 / (1.0 + phi**2 / 16.0)
        else:
            eta = (2.0 / phi) * math.exp(phi - 1.0)
    elif pellet_shape == "slab":
        eta = 1.0 / math.cosh(phi) if phi < 100 else 0.0
    else:
        eta = 1.0 / (1.0 + phi**2 / 9.0)
    
    return float(eta)


def pore_diffusion_correction(
    D_AB: float,
    T: float,
    P: float,
    pore_diameter: float = 1e-7,
    tortuosity: float = 3.0,
    porosity: float = 0.4,
) -> float:
    """Pore diffusion correction (Bosanquet equation)"""
    import math
    M_avg = 0.029
    D_knudsen = (pore_diameter / 3) * math.sqrt(8 * 8.314 * T / (math.pi * M_avg))
    D_eff_bos = (D_AB * D_knudsen) / (D_AB + D_knudsen)
    D_eff = (porosity / tortuosity) * D_eff_bos
    
    return float(D_eff)


def nonisothermal_pfr(
    F_A0: float,
    X_target: float,
    k_ref: float,
    E_a: float,
    delta_H_rxn: float,
    T_in: float,
    C_A0: float,
    U: float,
    A: float,  
    T_coolant: float,
    rho_cp: float,
    n: float = 1.0,
    num_segments: int = 100,
) -> dict:
    """Non-isothermal PFR design (segmented integration)"""
    import numpy as np
    import math
    
    R_gas = 8.314
    X_arr = np.linspace(0, X_target, num_segments + 1)
    dX = X_target / num_segments
    
    T_arr = [T_in]
    V_arr = [0.0]
    
    v0 = F_A0 / C_A0 if C_A0 > 0 else 1.0
    
    for i in range(num_segments):
        X_i = X_arr[i]
        X_next = X_arr[i + 1]
        X_avg = (X_i + X_next) / 2
        
        T_i = T_arr[-1]
        k = k_ref * math.exp(-E_a / R_gas * (1 / T_i - 1 / 300))
        C_A = C_A0 * (1 - X_avg)
        r_A = k * (C_A ** n)
        
        dV = F_A0 * dX / r_A if r_A > 0 else 0
        
        dT_adiabatic = (-delta_H_rxn * C_A0 / rho_cp) * dX if rho_cp > 0 else 0
        dT_cooling = U * A * (T_coolant - T_i) * dV / (v0 * rho_cp) if rho_cp > 0 else 0
        
        dT = dT_adiabatic + dT_cooling
        T_next = T_i + dT
        
        T_arr.append(T_next)
        V_arr.append(V_arr[-1] + dV)
    
    return {
        "X_profile": X_arr.tolist(),
        "T_profile": T_arr,
        "V_profile": V_arr,
        "V_total": V_arr[-1],
        "T_out": T_arr[-1],
        "T_max": max(T_arr),
        "T_min": min(T_arr),
    }


def design_reactor(
    reactor_type: str,
    F_A0: float,
    X_target: float,
    delta_H_rxn: float,
    k: float,
    C_A0: float,
    n: float = _REQUIRED,
    rho_cat: float = _REQUIRED,
    particle_diameter: float = _REQUIRED,
    epsilon: float = _REQUIRED,
    T_in: float = _REQUIRED,
    P_in: float = _REQUIRED,
    D_bed: Optional[float] = None,
    L_bed: Optional[float] = None,
    d_tube_inner: float = _REQUIRED,
    L_tube: float = _REQUIRED,
    U: float = _REQUIRED,
    T_coolant: float = _REQUIRED,
    num_tanks: int = _REQUIRED,
    rho_gas: float = _REQUIRED,
    mu: float = _REQUIRED,
    Cp_mix: float = _REQUIRED,
    C_B0: Optional[float] = None,
    **kwargs,
) -> dict:
    """通用反应器设计主函数

    Args:
        ... (原有参数)
        C_B0: 反应物B的入口浓度 (mol/m3), 传入后 CSTR 启用双分子动力学 r=k*CA*CB
    """
    missing = check_required_params(locals(), {
        "n": "反应级数（无量纲），如一级反应 n=1，由设备层传入",
        "rho_cat": "催化剂密度 (kg/m³)，由设备层从物性引擎或用户输入获取",
        "particle_diameter": "催化剂粒径 (m)，由设备层从物性引擎或用户输入获取",
        "epsilon": "床层空隙率（无量纲），由设备层从物性引擎或用户输入获取",
        "T_in": "反应器入口温度 (K)，由设备层计算提供",
        "P_in": "反应器入口压力 (Pa)，由设备层计算提供",
        "d_tube_inner": "列管内径 (m)，由设备层传入",
        "L_tube": "列管长度 (m)，由设备层传入",
        "U": "总传热系数 (W/m²·K)，由设备层计算提供",
        "T_coolant": "冷却介质温度 (K)，由设备层计算提供",
        "num_tanks": "CSTR 串联釜数，由设备层传入",
        "rho_gas": "气相密度 (kg/m³)，由设备层从物性引擎获取",
        "mu": "流体动力粘度 (Pa·s)，由设备层从物性引擎获取",
        "Cp_mix": "混合物热容 (J/(mol_A·K))，由设备层从物性引擎获取",
    })
    if missing:
        return missing

    reactor_type = reactor_type.lower().replace(" ", "_")
    _type_aliases = {
        "fbr": "fluidized_bed", "fluidized": "fluidized_bed",
        "fixed": "fixed_bed", "packed_bed": "fixed_bed",
        "plug_flow": "pfr", "cstr": "cstr", "cstr_tank": "cstr",
    }
    reactor_type = _type_aliases.get(reactor_type, reactor_type)

    results = {
        "reactor_type": reactor_type,
        "F_A0": F_A0,
        "X_target": X_target,
        "delta_H_rxn": delta_H_rxn,
    }

    # ==========================
    # 1. 固定床 (PFR/单段床)
    # ==========================
    if reactor_type in ["fixed_bed", "pfr", "plug_flow"]:
        kin_result = pfr_design(F_A0, X_target, k, C_A0, n, epsilon=0.0, C_B0=C_B0)
        results.update(kin_result)
        V = kin_result["V"]

        if V > 0:
            if D_bed is None and L_bed is None:
                L_D_target = 4.0
                D_bed = (4 * V / (math.pi * L_D_target)) ** (1/3)
                L_bed = L_D_target * D_bed
            elif D_bed is None and L_bed is not None:
                D_bed = math.sqrt(4 * V / (math.pi * L_bed))
            elif L_bed is None and D_bed is not None:
                L_bed = V / (math.pi * D_bed**2 / 4)
        else:
            if D_bed is None: D_bed = 0.3
            if L_bed is None: L_bed = 1.0

        geo = fixed_bed_dimensions(V, D_bed)
        results.update(geo)

        A_cross = geo["A_cross"]
        v_volumetric = F_A0 / C_A0 if C_A0 > 0 else 0.0
        u_actual = v_volumetric / A_cross if A_cross > 0 else 0.1

        press = fixed_bed_pressure_drop_ergun(
            L_bed, particle_diameter, epsilon, rho_gas, mu,
            u_s=u_actual, sphericity=0.8,
        )
        results.update(press)

        cat = fixed_bed_catalyst_inventory(V, rho_cat, particle_diameter, epsilon)
        results.update(cat)

        res = fixed_bed_residence_time(V, F_A0, P_in, T_in, epsilon)
        results.update(res)

        A_ht = math.pi * D_bed * L_bed if D_bed and L_bed else geo["A_cross"]

        heat = pfr_heat_transfer(
            delta_H_rxn, F_A0, X_target, U, A_ht, T_in, T_coolant, Cp_mix, T_coolant
        )
        results.update(heat)

        outlet = fixed_bed_outlet_conditions(F_A0, X_target, T_in, P_in, delta_H_rxn, Cp_mix, epsilon)
        results.update(outlet)
        # 修正：防止单纯进行绝热计算的 outlet 数据直接覆盖掉已考虑换热冷却的正确 T_out
        results["T_out"] = heat.get("T_out", outlet["T_out"])

        results["sub_type"] = "single_bed_pfr"

    # ==========================
    # 2. 多管固定床 (列管式)
    # ==========================
    elif reactor_type == "multitubular":
        mt = multitubular_design(F_A0, X_target, k, C_A0, n, d_tube_inner, L_tube, U, T_coolant, rho_cat)
        results.update(mt)

        V = mt["V_total"]
        geo = fixed_bed_dimensions(V, d_tube_inner * (mt["N_tubes"]**0.5) * 1.1)

        v_volumetric = F_A0 / C_A0 if C_A0 > 0 else 0.0
        A_tube_total = mt["N_tubes"] * math.pi * d_tube_inner**2 / 4
        u_actual = v_volumetric / A_tube_total if A_tube_total > 0 else 0.1

        press = fixed_bed_pressure_drop_ergun(L_tube, particle_diameter, epsilon, rho_gas, mu, u_s=u_actual)

        results.update({
            "V_bed": V,
            "D_bed": geo["D_bed"],
            "L_bed": L_tube,
            "delta_P": press["delta_P"],
            "sub_type": "multitubular",
        })

    # ==========================
    # 3. 流化床 (FBR)
    # ==========================
    elif reactor_type == "fluidized_bed":
        u_mf_res = minimum_fluidization_velocity(particle_diameter, rho_cat, rho_gas, mu, epsilon_mf=epsilon)
        u_t_res = terminal_velocity(particle_diameter, rho_cat, rho_gas, mu)

        u_actual = u_mf_res["u_mf"] * 3.0  
        regime = check_fluidization_regime(u_actual, u_mf_res["u_mf"], u_t_res["u_terminal"])

        v_volumetric = F_A0 / C_A0 if C_A0 > 0 else 1.0

        if D_bed is None:
            A_req = v_volumetric / u_actual if u_actual > 0 else 0.2
            D_bed = math.sqrt(4 * A_req / math.pi)
        if L_bed is None:
            L_bed = 2.0 * D_bed

        cat = fluidized_bed_catalyst_inventory(D_bed, L_bed, rho_cat, epsilon)
        results.update(cat)

        V_bed = math.pi * D_bed**2 / 4 * L_bed
        res = fluidized_bed_residence_time(L_bed, u_mf_res["u_mf"], u_actual)
        
        ht = fluidized_bed_heat_transfer_compat(
            D_bed, L_bed, T_in, T_coolant, kwargs.get("k_gas", 0.025), mu, rho_gas,
            u_mf_res["u_mf"], particle_diameter, Cp_mix
        )

        Q_gen = -delta_H_rxn * F_A0 * X_target
        Q_removed = Q_gen  

        outlet = fluidized_bed_outlet_conditions(F_A0, T_in, P_in, X_target, delta_H_rxn, Cp_mix, ht)

        results.update({
            "u_mf": u_mf_res["u_mf"],
            "u_terminal": u_t_res["u_terminal"],
            "u_actual": u_actual,
            "u_u_mf": regime["u_u_mf"],
            "regime": regime["regime"],
            "V_bed": V_bed,
            "H_dense": L_bed,
            "D_bed": D_bed,
            "delta_P": (rho_cat - rho_gas) * (1 - epsilon) * 9.81 * L_bed,  
            "Q_gen": Q_gen,
            "Q_removed": Q_removed,
            "tau_gas": res["tau_gas"],
            "T_out": outlet["T_out"],
            "sub_type": "bubbling_fluidized",
        })

    # ==========================
    # 4. CSTR
    # ==========================
    elif reactor_type == "cstr":
        cstr_res = cstr_design(F_A0, X_target, k, C_A0, n, num_tanks, C_B0=C_B0)
        results.update(cstr_res)

        V = cstr_res["V_total"]

        if num_tanks == 1:
            if D_bed is None:
                D_bed = (4 * V / (math.pi * 2)) ** (1/3)  
            H_design = 2 * D_bed
        else:
            V_single = cstr_res["V_per_tank"]
            if D_bed is None:
                D_bed = (4 * V_single / (math.pi * 2)) ** (1/3)
            H_design = 2 * D_bed

        results.update({
            "D_tank": D_bed,
            "H_tank": H_design,
            "V_single": V,
            "delta_P": 0.0, 
            "sub_type": f"{num_tanks}_tank_cstr",
        })

        outlet = {
            "T_out": T_in + adiabatic_temperature_rise(delta_H_rxn, X_target, F_A0, Cp_mix),
            "P_out": P_in,
            "F_out": F_A0 * (1 + X_target), 
            "X_out": X_target,
        }
        results.update(outlet)

    else:
        raise ValueError(f"Unknown reactor type: {reactor_type}")

    validation = validate_reactor_design(
        {
            "reactor_type": reactor_type,
            "max_pressure_drop": kwargs.get("max_pressure_drop", 300000),
        },
        results,
    )
    results["validation"] = validation

    return results


# ============================================================
# 反应器移热能力校核
# ============================================================

def reactor_heat_removal_check(
    Q_reaction: float,
    U: float,
    A_heat: float,
    delta_T: float,
    safety_factor: float = 1.2,
) -> dict:
    """反应器移热能力校核

    对比反应放热速率与夹套/盘管换热能力，判断反应器是否热失控。

    Q_removed = U * A * ΔT

    当 Q_removed >= safety_factor * Q_reaction 时，移热能力充足。

    Args:
        Q_reaction: 反应放热速率 (W)，正值表示放热
        U: 总传热系数 (W/(m²·K))
        A_heat: 换热面积 (m²)，夹套或盘管
        delta_T: 反应物料与冷却介质间的平均温差 (K)
        safety_factor: 安全系数，默认 1.2（移热能力需为放热的 1.2 倍）

    Returns:
        移热能力校核结果，包含是否充足、裕量、建议等
    """
    # 换热能力
    Q_removed = U * A_heat * delta_T

    # 所需移热速率（含安全裕量）
    Q_required = Q_reaction * safety_factor

    # 移热裕量
    Q_margin = Q_removed - Q_reaction
    ratio = Q_removed / Q_reaction if Q_reaction > 0 else float('inf')

    # 判断
    is_adequate = Q_removed >= Q_required

    # 所需换热面积
    A_required = Q_required / (U * delta_T) if (U * delta_T) > 0 else float('inf')

    # 所需温差
    delta_T_required = Q_required / (U * A_heat) if (U * A_heat) > 0 else float('inf')

    # 建议措施
    suggestions = []
    if not is_adequate:
        if A_required > A_heat * 1.5:
            suggestions.append("换热面积严重不足，建议增设盘管或扩大夹套")
        elif A_required > A_heat:
            suggestions.append(f"换热面积不足，需增加至 {A_required:.2f} m²")
        if delta_T_required > delta_T:
            suggestions.append(f"可降低冷却介质温度至增大温差（需 ΔT ≥ {delta_T_required:.1f} K）")
        if U < 200:
            suggestions.append("传热系数偏低，建议改善搅拌或清理换热面")
        if not suggestions:
            suggestions.append("移热能力略不足，建议适当降低进料量或反应温度")

    return {
        "Q_reaction_W": Q_reaction,
        "Q_removed_W": Q_removed,
        "Q_required_W": Q_required,
        "Q_margin_W": Q_margin,
        "removal_ratio": ratio,
        "safety_factor": safety_factor,
        "is_adequate": is_adequate,
        "U": U,
        "A_heat": A_heat,
        "delta_T": delta_T,
        "A_required_m2": A_required,
        "delta_T_required_K": delta_T_required,
        "suggestions": suggestions if not is_adequate else ["移热能力充足"],
    }