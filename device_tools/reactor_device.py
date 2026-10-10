"""反应器设备级 MCP 工具 — 支持副反应、多级反应、零/一/二级动力学的完整设计流水线

改进点：
1. 修复了 n, T_in, P_in 在参数列表丢失导致的 UnboundLocalError 漏洞。
2. 支持组分/混合比例（components & z & F_total）传入，自动反推精确的反应物 C_A0 与进料 F_A0。
3. 切断 **kwargs 向物理引擎底层的传递，实施绝对安全防护。
4. 补全各反应器类型的完整参数（几何、传热、压降、停留时间、催化剂库存）。
5. 完善单位换算（温度/压力/流量）与参数校验。
"""
from __future__ import annotations
import math
from typing import Optional, List, Dict, Any

from physics_engine.reactor import (
    design_reactor,
    validate_reactor_design,
    check_conversion_feasibility,
    reaction_rate,
    reaction_enthalpy,
    adiabatic_temperature_rise,
    selectivity_analysis,
    pfr_design,
    cstr_design,
    minimum_fluidization_velocity,
    terminal_velocity,
    check_fluidization_regime,
    fixed_bed_pressure_drop_ergun,
    fixed_bed_dimensions,
    fixed_bed_catalyst_inventory,
    fixed_bed_residence_time,
    catalyst_deactivation,
    calculate_thiele_modulus,
    effectiveness_factor,
    estimate_reactor_U,
)
from physics_engine.common import _REQUIRED, check_required_params
from tools.base import ToolRegistry
from tools.registry import UniversalToolRegistry, ToolCapability


# ============================================================
# 单位换算工具
# ============================================================

def _normalize_temperature(T: float) -> float:
    """≤200 视为 °C 自动转 K"""
    return T + 273.15 if T <= 200 else T


def _normalize_pressure(P: float) -> float:
    """≤500 视为 bar → Pa；≤5000 视为 kPa → Pa"""
    if P <= 0:
        return 101325.0
    if P <= 500:
        return P * 1e5
    if P <= 5000:
        return P * 1e3
    return P


def _normalize_flow(F: float) -> float:
    """"1000 视为 mol/h → mol/s"""
    return F / 3600.0 if F > 1000 else F


def _normalize_concentration(C: float) -> float:
    """≤0 时给默认值 1000 mol/m³"""
    return C if C > 0 else 1000.0


# ============================================================
# 副反应与多反应体系
# ============================================================

def _solve_parallel_side_reaction(
    F_A0: float,
    X_target: float,
    k_main: float,
    k_side: float,
    n_main: float,
    n_side: float,
    C_A0: float,
    reactor_type: str,
) -> dict:
    """平行副反应：A→D（主）+ A→U（副）"""
    import numpy as np
    N = 500
    X_arr = np.linspace(0, X_target, N + 1)
    dX = X_target / N

    S_D_vals = []
    for X in X_arr:
        C_A = C_A0 * (1 - X)
        if C_A <= 0:
            S_D_vals.append(0.0)
            continue
        r_main = k_main * (C_A ** n_main)
        r_side = k_side * (C_A ** n_side)
        r_tot = r_main + r_side
        S_D_vals.append(r_main / r_tot if r_tot > 0 else 1.0)

    S_D_avg = (dX * (0.5 * S_D_vals[0] + sum(S_D_vals[1:-1]) + 0.5 * S_D_vals[-1])) / X_target

    F_D = F_A0 * X_target * S_D_avg        
    F_U = F_A0 * X_target * (1 - S_D_avg)  
    yield_D = X_target * S_D_avg            

    return {
        "S_D_avg": S_D_avg,
        "yield_D": yield_D,
        "yield_U": X_target * (1 - S_D_avg),
        "F_D": F_D,
        "F_U": F_U,
        "k_ratio": k_main / k_side if k_side > 0 else float("inf"),
        "selectivity_note": (
            "n_main > n_side：高浓度有利于主反应，建议 PFR/低转化率"
            if n_main > n_side else
            "n_main < n_side：低浓度有利于主反应，建议 CSTR/高转化率"
            if n_main < n_side else
            "反应级数相同，选择性仅由速率常数之比决定"
        ),
    }


def _solve_series_side_reaction(
    F_A0: float,
    X_target: float,
    k1: float,
    k2: float,
    C_A0: float,
    n1: float = 1.0,
    n2: float = 1.0,
) -> dict:
    """串联副反应：A→D→U"""
    import numpy as np

    if abs(n1 - 1.0) < 1e-6 and abs(n2 - 1.0) < 1e-6 and abs(k1 - k2) > 1e-10:
        tau_opt = math.log(k2 / k1) / (k2 - k1) if k2 > k1 else math.log(k1 / k2) / (k1 - k2)
        C_A_opt = C_A0 * math.exp(-k1 * tau_opt)
        C_D_max = C_A0 * k1 / (k2 - k1) * (math.exp(-k1 * tau_opt) - math.exp(-k2 * tau_opt))
        X_D_opt = 1 - C_A_opt / C_A0

        tau = -math.log(1 - X_target) / k1 if X_target < 1 else 1e6
        C_D_out = C_A0 * k1 / (k2 - k1) * (math.exp(-k1 * tau) - math.exp(-k2 * tau))
        C_U_out = C_A0 - C_A0 * math.exp(-k1 * tau) - C_D_out
    else:
        N = 1000
        dt = 10.0 / max(k1, k2) / N
        C_A, C_D, C_U = C_A0, 0.0, 0.0
        C_D_max, X_D_opt = 0.0, 0.0
        for i in range(N * 100):
            r1 = k1 * max(C_A, 0) ** n1
            r2 = k2 * max(C_D, 0) ** n2
            dC_A = -r1 * dt
            dC_D = (r1 - r2) * dt
            dC_U = r2 * dt
            C_A += dC_A
            C_D += dC_D
            C_U += dC_U
            if C_D > C_D_max:
                C_D_max = C_D
                X_D_opt = 1 - C_A / C_A0
            X_cur = 1 - C_A / C_A0
            if X_cur >= X_target:
                break
        C_D_out = C_D
        C_U_out = C_U
        tau_opt = X_D_opt / k1 if k1 > 0 else 0.0

    yield_D = C_D_out / C_A0 if C_A0 > 0 else 0.0

    return {
        "C_D_max": C_D_max,
        "X_D_opt": X_D_opt,
        "C_D_out": C_D_out,
        "C_U_out": C_U_out,
        "yield_D": yield_D,
        "series_note": (
            f"D 浓度在 X={X_D_opt:.3f} 处达到最大值 {C_D_max:.2f} mol/m³；"
            f"当前 X={X_target:.3f} 下 D 收率 {yield_D:.3f}。"
            + ("建议降低转化率以提升 D 收率。" if X_target > X_D_opt + 0.05 else "")
        ),
    }


def _compute_nth_order_volume(
    reactor_type: str,
    F_A0: float,
    X_target: float,
    k: float,
    C_A0: float,
    n: float,
    C_B0: float = None,  # 双分子反应支持
) -> dict:
    """零/一/二/任意级反应体积解析解（液相）
    
    支持两种模式:
    1. 单组分: r = k * C_A^n (默认, C_B0=None)
    2. 双分子: r = k * C_A * C_B (C_B0 提供时)
    """
    X = min(X_target, 0.9999)
    v0 = F_A0 / C_A0  

    # ── 双分子反应模式: r = k * C_A * C_B ──
    if C_B0 is not None and C_B0 > 0:
        M = C_B0 / C_A0  # 浓度比
        
        # 转化率上限受限于 B 的化学计量
        X_max = min(M, 1.0) if M < 1 else 1.0
        if X >= X_max:
            X = X_max * 0.999
        
        if reactor_type.lower() in ("pfr", "fixed_bed", "plug_flow"):
            # PFR 双分子反应解析解
            # V = F_A0 / (k * C_A0²) * ∫dX / [(1-X)(M-X)]
            if abs(M - 1.0) < 1e-6:
                # M = 1 时: V = F_A0 / (k * C_A0²) * X / (1-X)
                V = (F_A0 / (k * C_A0**2)) * (X / (1 - X))
            else:
                # M ≠ 1 时: V = F_A0 / (k * C_A0² * (M-1)) * ln[(M-X)/(M*(1-X))]
                V = (F_A0 / (k * C_A0**2 * (M - 1))) * math.log((M - X) / (M * (1 - X)))
            tau = V * C_A0 / F_A0 if F_A0 > 0 else 0.0
        else:
            # CSTR 双分子反应
            # V = F_A0 * X / (k * C_A0² * (1-X) * (M-X))
            CA_out = C_A0 * (1 - X)
            CB_out = C_B0 - C_A0 * X
            if CA_out <= 0 or CB_out <= 0:
                V, tau = 0.0, 0.0
            else:
                r_out = k * CA_out * CB_out
                V = F_A0 * X / r_out if r_out > 0 else 0.0
                tau = V / v0
        
        return {
            "V": V,
            "tau": tau,
            "n_order": 2,
            "order_label": "二级(双分子)",
            "bimolecular": True,
            "M_ratio": M,
        }

    # ── 单组分模式: r = k * C_A^n ──
    if reactor_type.lower() in ("pfr", "fixed_bed", "plug_flow"):
        if abs(n) < 1e-6:          
            V = F_A0 * X / (k * C_A0)
            tau = X / k
        elif abs(n - 1.0) < 1e-6: 
            V = (v0 / k) * (-math.log(1 - X))
            tau = -math.log(1 - X) / k
        elif abs(n - 2.0) < 1e-6: 
            V = (v0 / (k * C_A0)) * (X / (1 - X))
            tau = X / (k * C_A0 * (1 - X))
        else:                       
            import numpy as np
            N = 2000
            X_arr = np.linspace(0, X, N + 1)
            integrand = 1.0 / (k * (C_A0 * (1 - X_arr)) ** n + 1e-300)
            V = float(F_A0 * np.trapz(integrand, X_arr))
            tau = V * C_A0 / F_A0

    else:  # CSTR
        C_A_out = C_A0 * (1 - X)
        r_out = k * max(C_A_out, 0) ** n
        if r_out <= 0:
            V, tau = 0.0, 0.0
        else:
            V = F_A0 * X / r_out
            tau = V / v0

    return {
        "V": V,
        "tau": tau,
        "n_order": n,
        "order_label": {0: "零级", 1: "一级", 2: "二级"}.get(round(n), f"{n:.1f}级"),
    }


# ============================================================
# 主设计函数 (已完美补齐参数声明)
# ============================================================

def _full_reactor_design(
    reactor_type: str = _REQUIRED,
    X_target: float = _REQUIRED,
    k: float = _REQUIRED,
    delta_H_rxn: float = _REQUIRED,
    T_in: float = _REQUIRED,
    F_A0: float = _REQUIRED,
    C_A0: float = _REQUIRED,
    n: float = _REQUIRED,
    P_in: float = _REQUIRED,
    C_B0: Optional[float] = None,  # 双分子反应：反应物B的入口浓度 (mol/m³)
    side_reaction_type: Optional[str] = None,
    k_side: float = 0.0,
    n_side: float = 1.0,
    k2: float = 0.0,
    n2: float = 1.0,
    rho_cat: float = None,
    particle_diameter: float = None,
    epsilon: float = None,
    cat_type: str = None,
    D_bed: Optional[float] = None,
    L_bed: Optional[float] = None,
    d_tube_inner: float = None,
    L_tube: float = None,
    U: float = None,
    T_coolant: float = None,
    rho_gas: float = _REQUIRED,
    mu: float = _REQUIRED,
    Cp_mix: float = _REQUIRED,
    num_tanks: Optional[float] = 1,
    D_eff: Optional[float] = None,
    components: Optional[list[str]] = None,
    z: Optional[list[float]] = None,
    F_total: Optional[float] = None,
    stoichiometry: Optional[Dict[str, float]] = None,
    **kwargs,
) -> dict:
    warnings = []
    conversion_log = []
    notes_missing = []

    # ── 阶段0：归一化 + 反应器类型识别 ──
    is_cstr = reactor_type and reactor_type.lower() in ("cstr", "cstr_tank")
    T_in = _normalize_temperature(T_in)
    P_in = _normalize_pressure(P_in)

    # ── 阶段1：混合物物性推算（在校验之前执行！） ──
    if components is not None and z is not None and len(components) == len(z):
        from physics_engine.thermo_helper import (
            get_fluid_density, get_fluid_MW, get_mixture_properties, get_formation_enthalpy,
            calc_liquid_density, calc_gas_density,
        )
        z = [float(zi) for zi in z]

        # 1a. 流量推算
        if F_total is not None and F_A0 is _REQUIRED:
            F_A0 = _normalize_flow(F_total) * z[0]
            conversion_log.append(f"[流量推算] F_A0={F_A0:.5f} mol/s")

        # 1b. 密度→C_A0 推算（逐组分液相+气相EOS+理想气体回退）
        if C_A0 is _REQUIRED:
            R_gas = 8.314  # J/(mol·K)
            V_m_sum = 0.0
            density_sources = []  # 记录每个组分的密度来源

            for i, name in enumerate(components):
                try:
                    mw_i = get_fluid_MW(name)
                except Exception:
                    mw_i = None

                rho_liq_T = None   # 液相密度 @ T_in
                rho_liq_298 = None  # 液相密度 @ 298.15 K
                rho_gas_i = None

                # 获取液相密度 @ T_in
                try:
                    val = calc_liquid_density(name, T_in, P_in)
                    if val and val > 100:
                        rho_liq_T = val
                except Exception:
                    pass

                # 获取液相密度 @ 298.15 K（回退）
                if rho_liq_T is None:
                    try:
                        val = calc_liquid_density(name, 298.15, P_in)
                        if val and val > 100:
                            rho_liq_298 = val
                    except Exception:
                        pass

                # 获取气相密度（EOS @ T_in, P_in）
                try:
                    val = calc_gas_density(name, T_in, P_in)
                    if val and val > 0:
                        rho_gas_i = val
                except Exception:
                    pass

                # 理想气体回退
                if rho_gas_i is None and mw_i is not None and T_in > 0:
                    rho_gas_i = P_in * (mw_i / 1000.0) / (R_gas * T_in)

                # 相态判定：
                # - 液相@T_in 存在 → 该温度下确为液相，使用它
                # - 仅有液相@298 + 气相@T → 实际操作温度下为气相，用气相密度
                # - 仅有气相 → 使用气相
                rho_i = None
                source = ""
                if rho_liq_T is not None:
                    rho_i, source = rho_liq_T, "液相密度@T_in"
                elif rho_gas_i is not None and rho_gas_i > 0:
                    # T_in 下有气相密度：该组分在操作条件下为气相
                    rho_i, source = rho_gas_i, "气相EOS@T_in"
                elif rho_liq_298 is not None:
                    # 仅 298K 液相密度可用，T_in 下无液/气数据，回退到 298K 液相
                    rho_i, source = rho_liq_298, "液相密度@298K(回退)"

                if rho_i is not None and rho_i > 0 and mw_i is not None:
                    V_m_sum += z[i] * (mw_i / 1000.0) / rho_i
                    density_sources.append(f"{name}: ρ={rho_i:.2f}({source})")
                else:
                    # 最终回退：MW 缺失时假设 80 g/mol + 理想气体
                    mw_fb = mw_i if mw_i is not None else 80.0
                    rho_fb = P_in * (mw_fb / 1000.0) / (R_gas * T_in) if T_in > 0 else 1.2
                    if rho_fb > 0:
                        V_m_sum += z[i] * (mw_fb / 1000.0) / rho_fb
                    density_sources.append(f"{name}: ρ≈{rho_fb:.2f}(回退EOS,MW={mw_fb:.0f})")

            if V_m_sum > 0:
                C_A0 = z[0] / V_m_sum
                notes_missing.append(f"C_A0 由逐组分密度加权推算：{C_A0:.2f} mol/m³")
                conversion_log.append("[C_A0推算] " + "; ".join(density_sources))
            else:
                C_A0 = 1000.0
                warnings.append("C_A0 无法从密度推算，使用缺省值 1000 mol/m³")

        # 1c. 混合物物性推算
        try:
            mp = get_mixture_properties(components, z, T_in, P_in)
            if rho_gas is _REQUIRED and mp.get("rho_L"):
                rho_gas = mp["rho_L"]; notes_missing.append(f"rho_gas 推算：{rho_gas:.1f} kg/m³")
            elif rho_gas is _REQUIRED and mp.get("rho_V"):
                rho_gas = mp["rho_V"]; notes_missing.append(f"rho_gas 推算：{rho_gas:.3f} kg/m³(气)")
            if mu is _REQUIRED and mp.get("mu_L"):
                mu = mp["mu_L"]; notes_missing.append(f"mu 推算：{mu:.5f} Pa·s")
            if Cp_mix is _REQUIRED and mp.get("Cp_L"):
                Cp_molar_mix = mp["Cp_L"]  # J/(mol·K) 摩尔加权混合热容
                # 换算为“每摩尔 A”基准: Cp_eff = Cp_mix × F_total / F_A0
                if F_A0 is not _REQUIRED and F_A0 > 0 and F_total is not None:
                    F_total_norm = _normalize_flow(F_total)
                    Cp_mix = Cp_molar_mix * F_total_norm / F_A0
                elif z[0] > 0:
                    Cp_mix = Cp_molar_mix / z[0]
                else:
                    Cp_mix = Cp_molar_mix
                notes_missing.append(f"Cp_mix 推算：{Cp_mix:.1f} J/(mol_A·K)（混合摩尔热容 {Cp_molar_mix:.1f} J/(mol·K)）")
        except Exception:
            pass

    # ── 阶段1d：反应焓推算（赫斯定律）—— 独立于 components，只要有 stoichiometry 就执行 ──
    if delta_H_rxn is _REQUIRED and stoichiometry:
        try:
            from physics_engine.thermo_helper import get_formation_enthalpy
            dH_rxn = 0.0
            dH_log = []
            all_ok = True
            for comp_name, nu_i in stoichiometry.items():
                hf = get_formation_enthalpy(comp_name, T_in, P_in, phase="gas")
                if hf is None:
                    hf = get_formation_enthalpy(comp_name, T_in, P_in, phase="liquid")
                if hf is None:
                    warnings.append(f"无法获取 '{comp_name}' 的生成焓")
                    all_ok = False; break
                dH_rxn += nu_i * hf
                dH_log.append(f"  {comp_name}(ν={nu_i}): ΔHf°={hf/1000:.2f} kJ/mol")
            if all_ok:
                delta_H_rxn = dH_rxn
                rtype = "放热" if dH_rxn < 0 else "吸热"
                notes_missing.append(f"[赫斯定律] ΔH_rxn = {dH_rxn/1000:.2f} kJ/mol ({rtype})\n" + "\n".join(dH_log))
        except Exception as e:
            warnings.append(f"反应焓推算失败: {e}")

    # ── 阶段2：必填参数校验（在推算之后！按反应器类型动态调整） ──
    required_params = {
        "reactor_type": "反应器类型", "X_target": "目标转化率", "k": "速率常数",
        "T_in": "入口温度", "n": "反应级数", "P_in": "入口压力",
    }
    # delta_H_rxn 不列入必填参数——缺失时走兜底逻辑(设为0.0，热平衡不可用)
    # rho_gas, mu, Cp_mix 也不列入必填——仅用于压降/传热计算，核心体积计算不需要
    if C_A0 is _REQUIRED:
        required_params["C_A0"] = "进料浓度；或提供 components+z 自动推算"
    if F_A0 is _REQUIRED:
        required_params["F_A0"] = "进料流量；或提供 F_total+z 自动推算"
    if num_tanks is _REQUIRED and is_cstr:
        required_params["num_tanks"] = "CSTR串联个数"

    missing = check_required_params(locals(), required_params)
    if missing:
        return missing

    # ── 阶段3：兜底填充 ──
    if C_A0 is _REQUIRED or C_A0 is None:
        C_A0 = 1000.0; notes_missing.append("C_A0 使用缺省值 1000 mol/m³")
    else:
        C_A0 = _normalize_concentration(C_A0)
    if F_A0 is _REQUIRED or F_A0 is None:
        F_A0 = 0.01; notes_missing.append("F_A0 使用缺省值 0.01 mol/s")
    else:
        F_A0 = _normalize_flow(F_A0)
    if delta_H_rxn is _REQUIRED:
        delta_H_rxn = 0.0; notes_missing.append("delta_H_rxn 未提供，热平衡不可用")
    if num_tanks is _REQUIRED or num_tanks is None:
        num_tanks = 1; notes_missing.append("num_tanks 默认为 1")
    # 物性参数兜底（仅用于压降/传热，核心体积计算不需要）
    if rho_gas is _REQUIRED or rho_gas is None:
        rho_gas = 800.0; notes_missing.append("rho_gas 使用缺省值 800 kg/m³（液相典型值）")
    if mu is _REQUIRED or mu is None:
        mu = 0.001; notes_missing.append("mu 使用缺省值 0.001 Pa·s（液相典型值）")
    if Cp_mix is _REQUIRED or Cp_mix is None:
        Cp_mix = 75.0; notes_missing.append("Cp_mix 使用缺省值 75 J/(mol·K)")

    if F_A0 <= 0:
        warnings.append("F_A0 必须 > 0，已重置为 0.01 mol/s。")
        F_A0 = 0.01

    if not (0 < X_target < 1):
        warnings.append(f"X_target={X_target} 超出 (0,1)，已截断。")
        X_target = min(max(X_target, 0.001), 0.999)

    if k <= 0:
        warnings.append("速率常数 k 必须 > 0，已重置为 0.1。")
        k = 0.1

    n = float(n)
    if n < 0:
        warnings.append(f"反应级数 n={n} 为负数，请确认。")

    # ── 阶段2b: 默认值参数追踪（在 _clean 之前记录哪些参数用户未提供） ──
    _USER_DEFAULTED_PARAMS = {}
    _DEFAULT_VALUES = {
        "rho_cat": 2000.0, "particle_diameter": 0.005, "epsilon": 0.4,
        "cat_type": "general", "d_tube_inner": 0.038, "L_tube": 6.0,
        "U": 100.0, "T_coolant": T_in,
    }
    for _pname, _pdefault in _DEFAULT_VALUES.items():
        _pval = locals().get(_pname)
        if _pval is None or _pval is _REQUIRED:
            _USER_DEFAULTED_PARAMS[_pname] = _pdefault

    # ── 阶段3b: _REQUIRED 哨兵清理 + defaulted_params 追踪 ──
    defaulted_params = {}
    def _clean(v, default=0.0, name=""):
        nonlocal defaulted_params
        if v is _REQUIRED:
            defaulted_params[name] = default
            return default
        if v is None:
            defaulted_params[name] = default
            return default
        return v
    rho_cat       = _clean(rho_cat,       2000.0,  "rho_cat")
    particle_diameter = _clean(particle_diameter, 0.005, "particle_diameter")
    epsilon       = _clean(epsilon,       0.4,    "epsilon")
    cat_type      = _clean(cat_type,      "general", "cat_type")
    d_tube_inner  = _clean(d_tube_inner,  0.038,  "d_tube_inner")
    L_tube        = _clean(L_tube,        6.0,    "L_tube")
    U             = _clean(U,             100.0,  "U")
    T_coolant     = _clean(T_coolant,     T_in,   "T_coolant")
    rho_gas       = _clean(rho_gas,       1.2,    "rho_gas")
    mu            = _clean(mu,            2e-5,   "mu")
    Cp_mix        = _clean(Cp_mix,        200.0,  "Cp_mix")  # J/(mol_A·K) 每摩尔A基准
    delta_H_rxn   = _clean(delta_H_rxn,   0.0,    "delta_H_rxn")
    num_tanks     = int(_clean(num_tanks,  1,      "num_tanks"))

    # ── 构建结果字典（在 _clean 之后，确保 input_parameters 是实际值） ──
    result: Dict[str, Any] = {
        "reactor_type": reactor_type,
        "input_parameters": {
            "F_A0_mol_s": F_A0, "X_target": X_target, "T_in_K": T_in, "P_in_Pa": P_in,
            "C_A0_mol_m3": C_A0, "n_order": n, "k": k,
            "delta_H_rxn_J_mol": delta_H_rxn,
            "rho_cat_kg_m3": rho_cat, "particle_diameter_m": particle_diameter,
            "epsilon": epsilon, "cat_type": cat_type,
            "d_tube_inner_m": d_tube_inner, "L_tube_m": L_tube,
            "U_W_m2K": U, "T_coolant_K": T_coolant,
            "rho_gas_kg_m3": rho_gas, "mu_Pas": mu, "Cp_mix_J_molA_K": Cp_mix,
            "num_tanks": num_tanks,
        },
        "warnings": warnings,
    }

    # 4. 催化剂失活修正
    k_effective = k
    time_on_stream = kwargs.get("time_on_stream")
    if time_on_stream is not None and time_on_stream > 0:
        deact = catalyst_deactivation(
            time_on_stream=time_on_stream,
            deactivation_model=kwargs.get("deactivation_model", "exponential"),
            k0=k,
            kd=kwargs.get("kd", 0.01),
            order=kwargs.get("deact_order", 1.0),
        )
        k_effective = deact["k_current"]
        result["catalyst_deactivation"] = deact
        conversion_log.append(f"催化剂失活修正：k {k:.4f} → {k_effective:.4f} （相对活性 {deact['relative_activity']:.3f}）")

    # 5. 催化剂效率因子（Thiele 模量）修正
    eta = 1.0
    if D_eff is not None and D_eff > 0:
        L_c = particle_diameter / 6.0  
        phi = calculate_thiele_modulus(k_effective, D_eff, L_c, n)
        eta = effectiveness_factor(phi, pellet_shape=kwargs.get("pellet_shape", "sphere"))
        k_effective = k_effective * eta
        result["thiele_analysis"] = {
            "L_c": L_c,
            "Thiele_modulus": phi,
            "effectiveness_factor": eta,
            "k_apparent": k_effective,
        }
        conversion_log.append(f"Thiele φ={phi:.3f}，效率因子 η={eta:.3f}，有效速率常数 k*={k_effective:.4f}")

    # 6. 零/一/二级反应解析体积
    nth_order_result = _compute_nth_order_volume(
        reactor_type=reactor_type,
        F_A0=F_A0, X_target=X_target,
        k=k_effective, C_A0=C_A0, n=n,
        C_B0=C_B0,  # 双分子反应支持
    )
    result["nth_order_analysis"] = nth_order_result
    if nth_order_result.get("bimolecular"):
        conversion_log.append(f"双分子二级反应解析解(A+B→P): V={nth_order_result['V']:.4f} m³, M=C_B0/C_A0={nth_order_result['M_ratio']:.3f}")
    else:
        conversion_log.append(f"{nth_order_result['order_label']}反应解析解：V={nth_order_result['V']:.4f} m³")

    # 7. 副反应计算
    if side_reaction_type == "parallel" and k_side > 0:
        side = _solve_parallel_side_reaction(
            F_A0=F_A0, X_target=X_target,
            k_main=k_effective, k_side=k_side,
            n_main=n, n_side=n_side,
            C_A0=C_A0, reactor_type=reactor_type,
        )
        result["parallel_side_reaction"] = side
        result["selectivity_D"] = side["S_D_avg"]
        result["yield_D"] = side["yield_D"]
        conversion_log.append(f"平行副反应：平均选择性 S_D={side['S_D_avg']:.3f}，主产物收率 Y_D={side['yield_D']:.3f}")

    elif side_reaction_type == "series" and (k2 > 0 or k_side > 0):
        _k2 = k2 if k2 > 0 else k_side
        series = _solve_series_side_reaction(
            F_A0=F_A0, X_target=X_target,
            k1=k_effective, k2=_k2,
            C_A0=C_A0, n1=n, n2=n2,
        )
        result["series_side_reaction"] = series
        result["yield_D"] = series["yield_D"]
        conversion_log.append(f"串联副反应：最优转化率 X_opt={series['X_D_opt']:.3f}，当前 D 收率={series['yield_D']:.3f}")

    # 8. 自动估算总传热系数 U（若用户未提供）
    if "U" in defaulted_params:
        # 估算空塔气速用于 Wakao 关联式
        V_est = nth_order_result.get("V", 0)
        if V_est > 0 and C_A0 > 0:
            A_cross_est = math.pi * (4 * V_est / math.pi / 4) ** (2/3) / 4  # 粗略估算
            u_s_est = (F_A0 / C_A0) / A_cross_est if A_cross_est > 0 else 0.1
        else:
            u_s_est = 0.1

        # 尝试获取流体导热系数
        k_fluid = None
        try:
            from physics_engine.thermo_helper import calc_thermal_conductivity_liquid, calc_thermal_conductivity_gas
            k_fluid = calc_thermal_conductivity_liquid(components[0], T_in) if components else None
        except Exception:
            pass

        # Cp 从 J/(mol_A·K) 转为 J/(kg·K) 用于 Wakao 关联式中的 Pr 计算
        Cp_mass = Cp_mix / (rho_gas * V_est / (F_A0 / C_A0)) if (rho_gas > 0 and V_est > 0 and F_A0 > 0 and C_A0 > 0) else 1000.0

        u_est = estimate_reactor_U(
            reactor_type=reactor_type,
            phase="gas" if rho_gas < 50 else "liquid",
            k_fluid=k_fluid,
            rho_gas=rho_gas,
            Cp_gas=Cp_mass,
            mu=mu,
            u_s=u_s_est,
            d_particle=particle_diameter,
        )
        U = u_est["U"]
        result["U_estimation"] = u_est
        notes_missing.append(f"U 自动估算：{U:.0f} W/m²·K（{u_est['method']}）")
        del defaulted_params["U"]

    # 9. 调用底层物理引擎 (不传 kwargs 防止底层崩溃)
    try:
        core = design_reactor(
            reactor_type=reactor_type,
            F_A0=F_A0, X_target=X_target,
            delta_H_rxn=delta_H_rxn,
            k=k_effective,       
            C_A0=C_A0, n=n,
            rho_cat=rho_cat,
            particle_diameter=particle_diameter,
            epsilon=epsilon,
            T_in=T_in, P_in=P_in,
            D_bed=D_bed, L_bed=L_bed,
            d_tube_inner=d_tube_inner,
            L_tube=L_tube,
            U=U, T_coolant=T_coolant,
            num_tanks=num_tanks,
            rho_gas=rho_gas,
            mu=mu, Cp_mix=Cp_mix,
            C_B0=C_B0,  # 双分子反应支持
        )
        result.update(core)
    except Exception as e:
        result["core_design_error"] = str(e)
        warnings.append(f"design_reactor 物理引擎计算失败：{e}。将以解析解作为回退方案。")
        V_fallback = nth_order_result["V"]
        if V_fallback > 0:
            L_D = 4.0
            D_calc = (4 * V_fallback / (math.pi * L_D)) ** (1 / 3)
            L_calc = L_D * D_calc
            result.update({
                "V": V_fallback,
                "D_bed": D_calc,
                "L_bed": L_calc,
                "H_D_ratio": L_D,
                "tau": nth_order_result["tau"],
            })

    # 10. 统一物理量纲与补全
    V = result.get("V", result.get("V_total", result.get("V_bed", 0.0)))
    if "H_D_ratio" not in result:
        D_r = result.get("D_bed", result.get("D_tank", 0))
        L_r = result.get("L_bed", result.get("H_tank", 0))
        if D_r > 0 and L_r > 0: result["H_D_ratio"] = L_r / D_r

    if "T_adia" not in result:
        dT_adia = adiabatic_temperature_rise(delta_H_rxn, X_target, F_A0, Cp_mix)
        result["T_adia"] = T_in + dT_adia
        result["delta_T_adiabatic"] = dT_adia

    if "Q_gen" not in result:
        result["Q_gen"] = abs(delta_H_rxn) * F_A0 * X_target

    W_cat = result.get("W_cat", 0)
    if W_cat > 0 and F_A0 > 0:
        result["WHSV"] = F_A0 / W_cat   

    # 10b. 流动与传质参数补全
    v0 = F_A0 / C_A0 if C_A0 > 0 else 0.0  # 体积流量 m³/s
    D_bed_val = result.get("D_bed", result.get("D_tank", 0))
    A_cross_val = math.pi * D_bed_val ** 2 / 4 if D_bed_val > 0 else 0.0
    u_s = v0 / A_cross_val if A_cross_val > 0 else 0.0  # 空塔气速
    u_interstitial = u_s / epsilon if epsilon > 0 else 0.0  # 颗粒间实际流速
    G_mass = rho_gas * u_s if rho_gas > 0 else 0.0  # 质量流速 kg/(m²·s)
    Re_p = rho_gas * u_s * particle_diameter / mu if (rho_gas > 0 and u_s > 0 and mu > 0) else 0.0

    result["flow_and_transport"] = {
        "v0_m3_s": round(v0, 6),
        "u_superficial_m_s": round(u_s, 4),
        "u_interstitial_m_s": round(u_interstitial, 4),
        "G_mass_kg_m2s": round(G_mass, 3),
        "Re_particle": round(Re_p, 2),
        "Re_note": "层流(Re<10)" if Re_p < 10 else "过渡区(10≤Re<1000)" if Re_p < 1000 else "湍流(Re≥1000)",
        "A_cross_m2": round(A_cross_val, 4),
    }

    # 多管反应器补全管数
    if result.get("sub_type") == "multitubular":
        result["flow_and_transport"]["N_tubes"] = result.get("N_tubes")
        result["flow_and_transport"]["d_tube_inner_m"] = d_tube_inner
        result["flow_and_transport"]["L_tube_m"] = L_tube
        A_tube_single = math.pi * d_tube_inner ** 2 / 4
        N_t = result.get("N_tubes", 0)
        if N_t and N_t > 0:
            A_tube_total = N_t * A_tube_single
            u_tube = v0 / A_tube_total if A_tube_total > 0 else 0.0
            result["flow_and_transport"]["u_tube_m_s"] = round(u_tube, 4)
            result["flow_and_transport"]["A_tube_total_m2"] = round(A_tube_total, 4)

    # 10c. 目标转化率可达性约束检查
    tau_for_check = nth_order_result.get("tau", 0.0)
    feasibility = check_conversion_feasibility(
        reactor_type=reactor_type,
        X_target=X_target,
        k=k_effective,
        C_A0=C_A0,
        n=n,
        F_A0=F_A0,
        V=V,
        tau=tau_for_check,
        C_B0=C_B0,
    )
    result["conversion_feasibility"] = feasibility
    if not feasibility["is_feasible"]:
        warnings.append(
            f"⚠ 目标转化率 X_target={X_target:.4f} 在当前条件下不可达！"
            + (f" 建议降至 ≤ {feasibility['X_max_feasible']:.4f}。" if feasibility["X_max_feasible"] else "")
        )
        conversion_log.append(
            f"[可达性检查] 不可行：{feasibility['n_fail']} 项约束未通过"
        )
    else:
        conversion_log.append(
            f"[可达性检查] 通过（{feasibility['n_warn']} 项警告）"
        )

    # 11. 约束校验
    validation_input = {
        "V": V,
        "delta_P": result.get("delta_P", 0),
        "T_out": result.get("T_out", T_in),
        "H_D_ratio": result.get("H_D_ratio", 0),
        "u_u_mf": result.get("u_u_mf", 0),
    }
    validation = validate_reactor_design(
        {"reactor_type": reactor_type, "max_pressure_drop": kwargs.get("max_pressure_drop", 300000)},
        validation_input,
    )
    result["validation"] = validation

    # 12. 工程优化建议
    suggestions = []
    # 转化率可达性建议
    if not feasibility["is_feasible"] and feasibility.get("suggestion"):
        suggestions.append(feasibility["suggestion"])
    dP = result.get("delta_P", 0)
    if dP > 200000: suggestions.append(f"压降 {dP/1000:.1f} kPa 偏高，建议增大颗粒粒径或减小床长。")
    LD = result.get("H_D_ratio", 0)
    if LD > 20: suggestions.append(f"长径比 L/D={LD:.1f} 过大，建议采用多管设计。")
    if 0 < LD < 1: suggestions.append(f"长径比 L/D={LD:.2f} 过小，建议增加高度。")
    if eta < 0.5:
        suggestions.append(f"内扩散严重限制（η={eta:.2f}），建议降低颗粒直径 dp 约 {particle_diameter*500:.1f} mm。")

    result["suggestions"] = suggestions
    result["conversion_log"] = conversion_log
    result["success"] = True

    # ── 阶13: 默认值参数提醒 + 缺失参数分类提示 ──
    if defaulted_params:
        result["defaulted_params"] = defaulted_params
    
    if notes_missing:
        result["notes_missing_params"] = notes_missing
    
    # ════════════════════════════════════════════════════════
    # 用户默认值参数提醒报告（8 个可动态调整的参数）
    # ════════════════════════════════════════════════════════
    if _USER_DEFAULTED_PARAMS:
        _UNIT_MAP = {
            "rho_cat": "kg/m³", "particle_diameter": "m", "epsilon": "",
            "cat_type": "", "d_tube_inner": "m", "L_tube": "m",
            "U": "W/m²·K", "T_coolant": "K",
        }
        _DESC_MAP = {
            "rho_cat": "催化剂堆密度（典型氧化铝/硅胶载体 ~1200-2500 kg/m³）",
            "particle_diameter": "催化剂粒径（常见球形颗粒 3-6 mm）",
            "epsilon": "床层空隙率（随机填充 ~0.35-0.45，规则填充 ~0.40-0.50）",
            "cat_type": "催化剂类型",
            "d_tube_inner": "反应管内径（工业常用 25-50 mm）",
            "L_tube": "反应管长（工业常用 3-12 m）",
            "U": "总传热系数（气相反应 ~50-200，液相 ~200-800 W/m²·K）",
            "T_coolant": "冷却/加热介质温度",
        }
        report_lines = [
            "═" * 60,
            "ℹ️ 以下参数使用了工程默认值（用户未显式提供）：",
            "═" * 60,
        ]
        for pname, pval in _USER_DEFAULTED_PARAMS.items():
            unit = _UNIT_MAP.get(pname, "")
            desc = _DESC_MAP.get(pname, "")
            val_str = f"{pval}" if isinstance(pval, str) else f"{pval:g}"
            unit_str = f" {unit}" if unit else ""
            report_lines.append(f"  ▸ {pname} = {val_str}{unit_str}")
            report_lines.append(f"    {desc}")
            report_lines.append(f"    💡 建议根据实际工况提供此参数以提高计算精度")
        report_lines.append("═" * 60)
        report_lines.append("💡 以上参数均为可动态调整的设计参数，可根据实际工况修改后重新计算")
        result["user_defaulted_params_report"] = "\n".join(report_lines)
        result["user_defaulted_params"] = _USER_DEFAULTED_PARAMS
    
    # 分类提示：哪些参数使用了默认值，提供后可做哪些计算
    param_hints = []

    _CATALYST_HINTS = {
        "rho_cat": ("催化剂颗粒密度 (kg/m³)", "计算床层压降、催化剂用量、WHSV"),
        "particle_diameter": ("催化剂颗粒直径 (m)", "计算 Ergun 压降、Thiele 模量、效率因子"),
        "epsilon": ("床层空隙率 (0~1)", "计算床层压降、最小流化速度"),
        "cat_type": ("催化剂类型名称", "标记催化剂类型，便于选型参考"),
    }
    _GEOMETRY_HINTS = {
        "d_tube_inner": ("反应管内径 (m)", "计算管内流速、传热面积、压降"),
        "L_tube": ("反应管长度 (m)", "计算传热面积、压降、停留时间"),
    }
    _HEAT_HINTS = {
        "U": ("总传热系数 (W/m²·K)", "计算非等温操作的热交换面积"),
        "T_coolant": ("冷却/加热介质温度 (K)", "计算传热温差和热负荷"),
    }
    _FLUID_HINTS = {
        "rho_gas": ("流体密度 (kg/m³)", "计算压降、Re、质量流速；提供 components+z 可自动推算"),
        "mu": ("流体动力粘度 (Pa·s)", "计算 Re、压降；提供 components+z 可自动推算"),
        "Cp_mix": ("混合物热容 (J/(mol_A·K))", "计算绝热温升、热平衡；提供 components+z 可自动推算"),
    }

    for hints_dict, cat_label, purpose_label in [
        (_CATALYST_HINTS, "催化剂参数", "催化剂设计"),
        (_GEOMETRY_HINTS, "几何参数", "反应器几何设计"),
        (_HEAT_HINTS, "传热参数", "热管理设计"),
        (_FLUID_HINTS, "流体物性", "压降与热平衡计算"),
    ]:
        missing_in_cat = {k: v for k, v in hints_dict.items()
                          if k in defaulted_params or k in _USER_DEFAULTED_PARAMS}
        if missing_in_cat:
            items = "\n".join(f"  - {k}: {desc}（当前默认值 {defaulted_params.get(k, _USER_DEFAULTED_PARAMS.get(k, '?'))}，提供后可{purpose}）"
                              for k, (desc, purpose) in missing_in_cat.items())
            param_hints.append(f"【{purpose_label}】\n{items}")

    if param_hints:
        result["suggestions"] = suggestions + param_hints
    if delta_H_rxn == 0.0 and "delta_H_rxn" in defaulted_params:
        result["suggestions"].append(
            "⚠ 反应焓 delta_H_rxn 使用默认值 0.0，热平衡和绝热温升计算不可用。"
            "提供 stoichiometry 参数可自动推算反应焓。"
        )

    return result


# ============================================================
# Executor 适配器
# ============================================================

class _FullDesignHandler:
    def execute(self, params: dict) -> dict:
        try:
            return _full_reactor_design(**params)
        except TypeError as e:
            return {
                "success": False,
                "error": f"参数错误：{e}。请检查参数名是否与工具定义一致。",
            }
        except Exception as e:
            return {"success": False, "error": str(e)}

_REACTOR_HANDLER = _FullDesignHandler()


# ============================================================
# 注册中心
# ============================================================

_PARAM_DESCRIPTIONS = {
    "F_A0":      "(必填) 反应物 A 的摩尔流量 mol/s，无默认值",
    "X_target":  "(必填) 目标转化率 (0~1)，无默认值",
    "reactor_type": "(必填) 反应器类型：PFR / fixed_bed / multitubular / FBR / fluidized_bed / CSTR，无默认值",
    "delta_H_rxn":  "(必填) 反应焓 J/mol，放热为负，无默认值",
    "k":            "(必填) 主反应速率常数，无默认值。二级反应常数单位为 m³/(mol·s)，若是 L/(mol·s) 请乘以 0.001 换算",
    "C_A0":         "(必填) 反应物A的进料浓度 mol/m³，无默认值。双分子反应中A为限量组分或用户指定的基准组分",
    "C_B0":         "(选填) 反应物B的进料浓度 mol/m³。双分子二级反应(r=k*CA*CB)必须提供！启用双分子动力学模式",
    "n":            "(必填) 主反应级数（支持 0/1/2/任意实数），无默认值",
    "T_in":         "(必填) 入口温度 K（≤200 自动按 °C 换算），无默认值",
    "P_in":         "(必填) 入口压力 Pa（≤500 按 bar，≤5000 按 kPa 换算），无默认值",
    "components":   "(选填) 进料组分名称列表（英文名），如 ['acetic_acid', 'ethanol']。提供后系统将自动高精度反推 C_A0",
    "z":            "(选填) 进料摩尔分数比例，必须与 components 对应，如 [0.2, 0.8]",
    "F_total":      "(选填) 进料总摩尔流量 mol/s，提供后结合 z 自动算出 A 的进料量 F_A0",
    # 副反应
    "side_reaction_type": "副反应类型：'parallel'（平行）/ 'series'（串联）/ None，默认 None",
    "k_side":      "副反应速率常数，默认 0",
    "n_side":      "平行副反应级数，默认 1.0",
    "k2":          "串联副反应 D→U 速率常数，默认 0",
    "n2":          "串联副反应 D→U 级数，默认 1.0",
    # 催化剂
    "rho_cat":          "(必填) 催化剂堆密度 kg/m³，无默认值",
    "particle_diameter":"(必填) 催化剂颗粒直径 m，无默认值",
    "epsilon":          "(必填) 床层空隙率，无默认值",
    "cat_type":         "(必填) 催化剂类型名称，无默认值",
    "D_eff":            "有效扩散系数 m²/s；提供后自动计算 Thiele 模量和效率因子",
    # 几何
    "d_tube_inner": "(必填) 反应管/单管内径 m，无默认值",
    "L_tube":       "(必填) 反应管长度 m，无默认值",
    # 传热
    "U":          "(必填) 总传热系数 W/m²·K，无默认值",
    "T_coolant":  "(必填) 冷却介质温度 K，无默认值",
    # 流体物性
    "rho_gas":  "(必填) 气体密度 kg/m³，无默认值",
    "mu":       "(必填) 流体动力粘度 Pa·s，无默认值",
    "Cp_mix":   "(必填) 混合物热容 J/(mol_A·K)（每摩尔A基准），无默认值",
    # CSTR
    "num_tanks": "(必填) CSTR 串联个数，无默认值",
}


def register_reactor_device(registry: ToolRegistry) -> ToolRegistry:
    registry.register_function(
        name="reactor_design",
        description=(
            "反应器完整设计。支持 PFR/固定床/列管/流化床/CSTR，"
            "支持零/一/二/任意级反应。支持自动基于组分与摩尔比例反推真实的反应浓度(C_A0)和流量(F_A0)。"
        ),
        func=_REACTOR_HANDLER,
        param_descriptions=_PARAM_DESCRIPTIONS,
        category="device",
        tags=["reactor", "PFR", "CSTR", "FBR", "固定床", "反应器"],
    )
    return registry


def register_to_universal(registry: UniversalToolRegistry) -> None:
    capability = ToolCapability(
        name="reactor_design",
        description=(
            "反应器完整设计。支持 PFR/固定床/列管式/流化床/CSTR，"
            "支持零/一/二/任意级反应。**支持双分子二级反应(r=k*CA*CB)**：提供 C_B0 启用双分子动力学模式。"
            "提供 components+z+F_total 可自动推算 C_A0/F_A0/流体物性；"
            "提供 stoichiometry 可通过赫斯定律自动推算反应焓 delta_H_rxn。"
        ),
        category="device_design",
        parameters={
            "reactor_type": {
                "type": "string", "required": True,
                "description": "(必填) 反应器类型：PFR / fixed_bed / multitubular / FBR / fluidized_bed / CSTR",
            },
            "X_target": {
                "type": "number", "required": True,
                "description": "(必填) 目标转化率 0~1",
            },
            "k": {
                "type": "number", "required": True,
                "description": "(必填) 主反应速率常数。注意二级反应常数单位为 m³/(mol·s)，若用户提供的是 L/(mol·s)，请乘以 0.001 换算！",
            },
            "T_in": {
                "type": "number", "required": True,
                "description": "(必填) 入口温度 K（≤200 自动按 °C 换算）",
            },
            "n": {
                "type": "number", "required": True,
                "description": "(必填) 主反应级数（0=零级，1=一级，2=二级）",
            },
            "P_in": {
                "type": "number", "required": True,
                "description": "(必填) 入口压力 Pa（≤500 按 bar，≤5000 按 kPa 换算）",
            },
            "delta_H_rxn": {
                "type": "number", "required": False,
                "description": "(选填) 反应焓 J/mol（放热为负）。若未提供但有 stoichiometry，系统自动通过赫斯定律推算",
            },
            "F_A0": {
                "type": "number", "required": False,
                "description": "(选填) A 的进料摩尔流量 mol/s。若未提供但有 F_total+z，系统自动推算",
            },
            "C_A0": {
                "type": "number", "required": False,
                "description": "(选填) 反应物A的进料浓度 mol/m³。若未提供但有 components+z，系统自动从混合物密度推算。双分子反应中A为基准组分",
            },
            "C_B0": {
                "type": "number", "required": False,
                "description": "(选填) 反应物B的进料浓度 mol/m³。**双分子二级反应(r=k*CA*CB)必须提供**，启用双分子动力学模式。如甲醇+乙酸反应，A=甲醇，B=乙酸",
            },
            "stoichiometry": {
                "type": "object", "required": False,
                "description": "(选填) 化学计量系数字典，负数=反应物，正数=生成物。如 {'methanol': -1, 'acetic acid': -1, 'methyl acetate': 1, 'water': 1}。提供后系统通过赫斯定律自动推算 delta_H_rxn",
            },
            "components": {
                "type": "list", "required": False,
                "description": "(选填) 组分英文名列表（空格分隔），如 ['acetic acid', 'ethanol']。提供后系统结合密度自动推算 C_A0、F_A0、流体物性",
            },
            "z": {
                "type": "list", "required": False,
                "description": "(选填) 进料摩尔分数比例，必须与 components 对应，如 [0.4, 0.6]",
            },
            "F_total": {
                "type": "number", "required": False,
                "description": "(选填) 总进料摩尔流量 mol/s，提供后结合 z[0] 自动计算 F_A0",
            },
            "num_tanks": {
                "type": "integer", "required": False,
                "description": "(选填) CSTR 串联级数，默认 1",
            },
            "rho_cat": {
                "type": "number", "required": False,
                "description": "(选填) 催化剂堆密度 kg/m³，默认 2000",
            },
            "particle_diameter": {
                "type": "number", "required": False,
                "description": "(选填) 催化剂颗粒直径 m，默认 0.005",
            },
            "epsilon": {
                "type": "number", "required": False,
                "description": "(选填) 床层空隙率，默认 0.4",
            },
            "cat_type": {
                "type": "string", "required": False,
                "description": "(选填) 催化剂类型名称，默认 'general'",
            },
            "d_tube_inner": {
                "type": "number", "required": False,
                "description": "(选填) 反应管/单管内径 m，默认 0.038",
            },
            "L_tube": {
                "type": "number", "required": False,
                "description": "(选填) 反应管长度 m，默认 6.0",
            },
            "U": {
                "type": "number", "required": False,
                "description": "(选填) 总传热系数 W/m²·K，默认 100",
            },
            "T_coolant": {
                "type": "number", "required": False,
                "description": "(选填) 冷却介质温度 K，默认等于 T_in",
            },
            "rho_gas": {
                "type": "number", "required": False,
                "description": "(选填) 气体密度 kg/m³。提供 components+z 可自动推算",
            },
            "mu": {
                "type": "number", "required": False,
                "description": "(选填) 流体动力粘度 Pa·s。提供 components+z 可自动推算",
            },
            "Cp_mix": {
                "type": "number", "required": False,
                "description": "(选填) 混合物热容 J/(mol_A·K)（每摩尔A基准）。提供 components+z 可自动推算",
            },
            "side_reaction_type": {
                "type": "string", "required": False,
                "description": "(选填) 副反应类型：parallel / series",
            },
            "k_side": {
                "type": "number", "required": False,
                "description": "(选填) 平行副反应速率常数，默认 0",
            },
            "k2": {
                "type": "number", "required": False,
                "description": "(选填) 串联副反应 D→U 速率常数，默认 0",
            },
            "D_eff": {
                "type": "number", "required": False,
                "description": "(选填) 有效导热/扩散系数 m²/s",
            },
        },
        tags=["reactor", "反应器", "PFR", "CSTR", "FBR"],
    )
    registry.register(capability, _full_reactor_design)