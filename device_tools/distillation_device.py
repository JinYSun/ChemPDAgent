"""精馏塔设备级 MCP 工具 — 完整 FUG 设计流水线"""
import math
from physics_engine.distillation import (
    design_distillation_column, calculate_column_diameter, get_component_properties, get_saturation_pressure,
    predict_mixture_viscosity, predict_mixture_molecular_weight,
    calculate_tray_hydraulics, bubble_point_temperature, dew_point_temperature,
    calculate_q_factor, underwood_minimum_reflux,
    calculate_condenser_duty, calculate_reboiler_duty, calculate_column_heat_duty,
)
from physics_engine.thermo_helper import get_fluid_density, get_fluid_Cp, calc_enthalpy_vaporization
from physics_engine.common import _REQUIRED, check_required_params
from tools.base import ToolRegistry
from tools.registry import UniversalToolRegistry, ToolCapability

# ============================================================
# 精馏塔标准参数集（工程常用值）
# ============================================================
# 回流比倍数 R/Rmin — 常用范围
STANDARD_RR_FACTOR = {"default": 2.0, "range": [1.2, 1.5, 2.0, 2.5, 3.0], "note": "常用 1.5-2.5，推荐 2.0"}
# 塔板类型
STANDARD_TRAY_TYPE = {"default": "sieve", "options": ["sieve", "valve"], "note": "sieve=筛板, valve=浮阀"}
# 液泛分率 — 常用范围
STANDARD_FLOODING_FRACTION = {"default": 0.80, "range": [0.70, 0.75, 0.80, 0.82, 0.85], "note": "常用 0.75-0.82，推荐 0.80"}
# 塔板间距 (m) — 常用规格
STANDARD_TRAY_SPACING = {"default": 0.60, "range": [0.45, 0.50, 0.60, 0.75, 0.90], "note": "常用 0.5-0.6m，推荐 0.60"}


def _full_design(components: list, z: list,
                 F: float = _REQUIRED,
                 x_D_spec: float = _REQUIRED, x_B_spec: float = _REQUIRED,
                 P: float = _REQUIRED,
                 q: float = None,
                 T_feed: float = None,
                 RR_factor: float = 2.0, tray_type: str = 'sieve',
                 tray_efficiency: float = None,
                 flooding_fraction: float = 0.80,
                 tray_spacing: float = 0.60,
                 flood_method: str = "simple",
                 # ── 可选：用户给定物性/设计参数（给定后跳过内部计算/查询） ──
                 alpha: float = None,        # 相对挥发度
                 rho_V: float = None,        # 气相密度 kg/m³
                 rho_L: float = None,        # 液相密度 kg/m³
                 mu_L: float = None,         # 液相粘度 cP
                 Cp_feed: float = None,      # 进料比热容 J/(mol·K)
                 **kwargs) -> dict:
    # 检查必填参数：所有核心设计参数和操作条件必须由用户提供
    missing = check_required_params(locals(), {
        "F": "进料流量 (kmol/h)",
        "P": "操作压力 (Pa)，如常压为 101325",
        "x_D_spec": "塔顶轻关键组分摩尔分数要求 (0~1)，如 0.99",
        "x_B_spec": "塔釜轻关键组分摩尔分数要求 (0~1)，如 0.02",
    })
    if missing:
        return missing
    q_auto_notes = []
    n = len(components)
    if n < 2: return {"error": "需要至少2个组分", "design_summary": {}}
    comp_props = [get_component_properties(c) for c in components]
    MWs = [p.get("MW", 50) for p in comp_props]
    Tbs = [p.get("Tb", 373) for p in comp_props]

    # === 温度计算 ===
    # 进料泡点/露点（用于 q 值计算，基于进料组成 z）
    T_feed_bub = bubble_point_temperature(P, components, z, T_guess=max(Tbs))
    T_feed_dew = dew_point_temperature(P, components, z, T_guess=min(Tbs))

    # 塔顶温度 = 塔顶液相组成 x_D 的泡点（全凝器）
    # 塔底温度 = 塔底液相组成 x_B 的泡点
    x_D_comp = [x_D_spec] + [0.0] * (n - 2) + [1.0 - x_D_spec] if n > 1 else [x_D_spec]
    x_B_comp = [x_B_spec] + [0.0] * (n - 2) + [1.0 - x_B_spec] if n > 1 else [x_B_spec]
    T_top = bubble_point_temperature(P, components, x_D_comp, T_guess=min(Tbs))
    T_bot = bubble_point_temperature(P, components, x_B_comp, T_guess=max(Tbs))
    T_avg = (T_bot.get("Tb", 390) + T_top.get("Tb", 360)) / 2

    # 自动计算 q 值（当用户未显式提供时）—— 使用进料组成 z 的泡点/露点
    if q is None:
        T_bub = T_feed_bub.get("Tb", 0)
        T_dew = T_feed_dew.get("Td", 0)
        if T_feed is not None and T_feed > 0 and T_bub > 0 and T_dew > 0:
            # 获取进料混合物的 Cp 和汽化焓
            # ── 进料比热容：用户给定值优先 ──
            if Cp_feed is not None:
                Cp_feed_val = float(Cp_feed)
                q_auto_notes.append(f"[用户给定] 进料比热容 Cp_feed = {Cp_feed_val:.2f} J/(mol·K)")
            else:
                try:
                    Cp_feed_val = get_fluid_Cp(components, T_feed, P, phase="liquid", z=z)
                except Exception:
                    Cp_feed_val = 0.0
            try:
                Hvap_vals = [calc_enthalpy_vaporization(c, T_feed) for c in components]
                Hvap_vals = [v for v in Hvap_vals if v is not None and v > 0]
                latent_heat_feed = sum(Hvap_vals[i] * z[i] for i in range(min(len(Hvap_vals), len(z)))) if Hvap_vals else 0.0
            except Exception:
                latent_heat_feed = 0.0
            q_result = calculate_q_factor(
                feed_temperature=T_feed, bubble_point=T_bub, dew_point=T_dew,
                Cp_feed=Cp_feed_val, latent_heat_feed=latent_heat_feed,
            )
            q = q_result["q"]
            q_auto_notes.append(f"[q自动计算] T_feed={T_feed:.1f}K, T_bub={T_bub:.1f}K, T_dew={T_dew:.1f}K → q={q:.3f} ({q_result.get('feed_state', '')})")
            if q_result.get("warnings"):
                q_auto_notes.extend(q_result["warnings"])
        else:
            q = 1.0
            q_auto_notes.append("[q默认] 未提供进料温度或泡露点数据，默认 q=1.0（饱和液体进料）")
    Psat_lk = get_saturation_pressure(components[0], T_avg) * 1000
    Psat_hk = get_saturation_pressure(components[-1], T_avg) * 1000
    # ── 相对挥发度：用户给定值优先 ──
    if alpha is not None:
        alpha_LK_HK = float(alpha)
        q_auto_notes.append(f"[用户给定] 相对挥发度 alpha = {alpha_LK_HK:.3f}")
    else:
        alpha_LK_HK = Psat_lk / Psat_hk if Psat_hk > 0 else 2.0
    x_LK_dist, x_HK_dist = x_D_spec, 1.0 - x_D_spec
    x_LK_btm, x_HK_btm = x_B_spec, 1.0 - x_B_spec
    # ── 分区物性计算 ────────────────────────────────────
    # 精馏段代表组成 ≈ x_D，提馏段代表组成 ≈ x_B，进料组成 z 保留作基准
    def _rho_V_sec(zs, T, P_sec=None):
        _P = P_sec if P_sec is not None else P
        return sum(get_fluid_density(components[i], T, _P, phase="gas") * zs[i] for i in range(n))

    def _rho_L_sec(zs, T):
        val = 0.0
        for i in range(n):
            try:
                val += get_fluid_density(components[i], T, P, phase="liquid") * zs[i]
            except Exception:
                val += 0.0
        return val

    MW_RS = sum(MWs[i] * x_D_comp[i] for i in range(n))
    MW_SS = sum(MWs[i] * x_B_comp[i] for i in range(n))
    MW_avg = sum(MWs[i] * z[i] for i in range(n))

    if rho_V is not None:
        rho_V_val = float(rho_V)
        rho_V_RS = rho_V_SS = rho_V_val
        q_auto_notes.append(f"[用户给定] 气相密度 rho_V = {rho_V_val:.2f} kg/m³")
    else:
        rho_V_RS = _rho_V_sec(x_D_comp, T_top.get("Tb", T_avg))
        rho_V_SS = _rho_V_sec(x_B_comp, T_bot.get("Tb", T_avg))
        rho_V_val = _rho_V_sec(z, T_avg)  # 进料基准（兼容设计塔径）

    if rho_L is not None:
        rho_L_val = float(rho_L)
        rho_L_RS = rho_L_SS = rho_L_val
        q_auto_notes.append(f"[用户给定] 液相密度 rho_L = {rho_L_val:.2f} kg/m³")
    else:
        rho_L_RS = _rho_L_sec(x_D_comp, T_top.get("Tb", T_avg))
        rho_L_SS = _rho_L_sec(x_B_comp, T_bot.get("Tb", T_avg))
        rho_L_val = _rho_L_sec(z, T_avg)  # 进料基准
    # 物料衡算 — 同时满足 x_D_spec 和 x_B_spec 约束
    F_kmolh = F
    # 联立求解: F = D + B, F*z[0] = D*x_D + B*x_B
    # => D = F * (z[0] - x_B) / (x_D - x_B)
    if x_D_spec > x_B_spec:
        D_kmolh = F_kmolh * (z[0] - x_B_spec) / (x_D_spec - x_B_spec)
    else:
        D_kmolh = F_kmolh * z[0] / x_D_spec if x_D_spec > 0 else F_kmolh * 0.5
    D_mol = D_kmolh / 3600.0  # 塔顶产品流量 mol/s

    # 相对挥发度（用于 Underwood 和 FUG 计算）
    Psat_hk = get_saturation_pressure(components[-1], T_avg) * 1000
    alpha_i = []
    for i in range(n):
        Psat_i = get_saturation_pressure(components[i], T_avg) * 1000
        alpha_i.append(Psat_i / Psat_hk if Psat_hk > 0 else 1.0)
    x_di = [0.0] * n
    x_di[0] = x_D_spec
    if n > 1: x_di[-1] = 1.0 - x_D_spec

    # Underwood 最小回流比 → 实际操作回流比
    underwood_pre = underwood_minimum_reflux(
        alpha_i=alpha_i, z_i=z, x_di=x_di, q=q,
    )
    R_min = underwood_pre.get("R_min", 1.0)
    R_actual = R_min * RR_factor
    # 塔内气液流量（分区计算，各段用各自 MW）
    # 精馏段: V_RS = D*(R+1), L_RS = D*R
    V_mol = D_mol * (R_actual + 1.0)
    V_mass = V_mol * MW_RS   # kg/s (D_mol in kmol/s × MW g/mol = kg/s)
    L_mass = D_mol * R_actual * MW_RS
    # 提馏段: L_SS = L_RS + q*F, V_SS = V_RS - (1-q)*F  （全部 kmol/s）
    F_kmol_s = F_kmolh / 3600.0  # kmol/s（与 D_mol 同单位）
    L_mol_SS = D_mol * R_actual + q * F_kmol_s
    V_mol_SS = V_mol - (1.0 - q) * F_kmol_s
    V_mass_SS = V_mol_SS * MW_SS   # kg/s
    L_mass_SS = L_mol_SS * MW_SS
    # 混合物粘度：用户给定值优先
    if mu_L is not None:
        mu_L_cP = float(mu_L)
        q_auto_notes.append(f"[用户给定] 液相粘度 mu_L = {mu_L_cP:.3f} cP")
    else:
        mu_result = predict_mixture_viscosity(components, z, T_avg)
        mu_L_cP = mu_result.get("mu_l", 0.3)
    # 塔板效率：若外部指定则直接使用，否则自动按 O'Connell 法估算
    if tray_efficiency is not None and 0 < tray_efficiency <= 1:
        eff_used = tray_efficiency
        eff_note = "（外部指定）"
    else:
        from physics_engine.distillation import estimate_tray_efficiency
        # mu_L_cP 单位是 cP，estimate_tray_efficiency 期望 Pa·s，需要转换
        eff_used = estimate_tray_efficiency(alpha_LK_HK, mu_L_cP / 1000.0)
        if isinstance(eff_used, dict):
            eff_used = eff_used.get("E0", eff_used.get("efficiency", 0.6))
        eff_note = "（O'Connell 自动估算）"
    # 重关键组分物料平衡：B = F - D，x_B_LK 由约束确定
    B_kmolh = F_kmolh - D_kmolh
    if B_kmolh <= 0:
        B_kmolh = F_kmolh * (1 - z[0] / x_D_spec) if x_D_spec > 0 else F_kmolh * 0.5
        D_kmolh = F_kmolh - B_kmolh
    # 各组分在产品中的流量
    comp_dist = []
    comp_bottom = []
    for i in range(n):
        fi_feed = F_kmolh * z[i]
        if i == 0:  # 轻关键
            fi_dist = D_kmolh * x_D_spec
            fi_bot = fi_feed - fi_dist
        elif i == n - 1:  # 重关键
            fi_bot = B_kmolh * (1 - x_B_spec)
            fi_dist = fi_feed - fi_bot
        else:
            fi_dist = fi_feed * 0.5
            fi_bot = fi_feed - fi_dist
        comp_dist.append({"component": components[i], "feed_kmolh": round(fi_feed, 3),
                          "distillate_kmolh": round(max(0, fi_dist), 3),
                          "bottoms_kmolh": round(max(0, fi_bot), 3)})
    # 实际产品组成（从组分流量计算）
    LK_dist = D_kmolh * x_D_spec
    LK_bot = F_kmolh * z[0] - LK_dist
    x_D_actual = round(LK_dist / D_kmolh, 4) if D_kmolh > 0 else x_D_spec
    x_B_actual = round(LK_bot / B_kmolh, 4) if B_kmolh > 0 else x_B_spec
    
    # === 热负荷计算 ===
    # 塔顶产品流量 mol/s
    D_mol_s = D_kmolh * 1000.0 / 3600.0
    B_mol_s = B_kmolh * 1000.0 / 3600.0
    
    # 塔顶组成列表
    x_D_list = [x_D_spec] + [0.0] * (n - 2) + [1.0 - x_D_spec] if n > 1 else [x_D_spec]
    # 塔底组成列表
    x_B_list = [x_B_spec] + [0.0] * (n - 2) + [1.0 - x_B_spec] if n > 1 else [x_B_spec]
    
    # ── 塔径设计基准（分段一致物性）──────────────────────────
    # 各段用各自的 V_mass/ρV/ρL 独立核算所需塔径（Souders-Brown，净面积基准），
    # 取所需塔径较大者为控制段，以该段的流量与物性定径 —— 避免混用进料组成物性
    _dia_RS = calculate_column_diameter(
        V=V_mass, rho_vapor=rho_V_RS, rho_liquid=rho_L_RS,
        tray_spacing=tray_spacing, flooding_fraction=flooding_fraction,
        L_mass=L_mass, flood_method=flood_method,
    )
    _dia_SS = calculate_column_diameter(
        V=V_mass_SS, rho_vapor=rho_V_SS, rho_liquid=rho_L_SS,
        tray_spacing=tray_spacing, flooding_fraction=flooding_fraction,
        L_mass=L_mass_SS, flood_method=flood_method,
    )
    if (_dia_SS.get("diameter_calculated") or 0) > (_dia_RS.get("diameter_calculated") or 0):
        _gov_section = "stripping"
        V_design, L_design = V_mass_SS, L_mass_SS
        rho_V_design, rho_L_design = rho_V_SS, rho_L_SS
    else:
        _gov_section = "rectifying"
        V_design, L_design = V_mass, L_mass
        rho_V_design, rho_L_design = rho_V_RS, rho_L_RS
    q_auto_notes.append(
        f"[塔径基准] 分段独立核算: D_calc(精馏段)={_dia_RS.get('diameter_calculated', 0):.3f} m, "
        f"D_calc(提馏段)={_dia_SS.get('diameter_calculated', 0):.3f} m → 控制段={_gov_section}，以该段流量与物性定径"
    )

    # 热负荷移至塔底压力修正后计算（见下方 T_bot 修正块）
    result = design_distillation_column(
        alpha_LK_HK=alpha_LK_HK, x_LK_dist=x_LK_dist, x_HK_dist=x_HK_dist,
        x_LK_btm=x_LK_btm, x_HK_btm=x_HK_btm, z_i=z, alpha_i=alpha_i, x_di=x_di,
        q=q, R_ratio=RR_factor, V_mass=V_design, L_mass=L_design,
        rho_vapor=rho_V_design, rho_liquid=rho_L_design, tray_efficiency=eff_used,
        flooding_fraction=flooding_fraction, tray_spacing=tray_spacing,
        P_top=P, flood_method=flood_method,
    )
    # 提取物理引擎子结果用于暴露中间参数
    diameter_result = result.get("diameter_result", {})
    hydraulics_result = result.get("hydraulics_result", {})
    fenske_detail = result.get("fenske", {})
    underwood_detail = result.get("underwood", {})
    gilliland_detail = result.get("gilliland", {})

    # ── 分区水力学（用各段自己的 ρ_V, ρ_L, V_mass, L_mass 独立计算） ──
    D_col_val = result.get("D_column", 1.0)
    hydraulics_RS = calculate_tray_hydraulics(
        V_mass=V_mass, L_mass=L_mass,
        rho_vapor=rho_V_RS, rho_liquid=rho_L_RS,
        sigma=None, mu_liquid=None, D_column=D_col_val,
        tray_spacing=tray_spacing, flood_method=flood_method,
    )
    hydraulics_SS = calculate_tray_hydraulics(
        V_mass=V_mass_SS, L_mass=L_mass_SS,
        rho_vapor=rho_V_SS, rho_liquid=rho_L_SS,
        sigma=None, mu_liquid=None, D_column=D_col_val,
        tray_spacing=tray_spacing, flood_method=flood_method,
    )

    # ── 塔底压力修正：T_bot 应在实际底部压力下计算 ──
    # 修正后的塔内物理塔板数（已扣除再沸器）与各段实际板数（均不含进料板）
    _N_act = result.get("N_actual_trays", 0)
    _N_RS = result.get("N_rectifying_actual") or 0
    _N_SS = result.get("N_stripping_actual") or 0
    _dP_RS = hydraulics_RS.get("pressure_drop_per_tray", 0) or 0
    _dP_SS = hydraulics_SS.get("pressure_drop_per_tray", 0) or 0
    # 总压降按修正后的物理塔板数累加：ΔP_total = (N_RS+1)×dP_RS + N_SS×dP_SS（进料板计入精馏段侧）
    if _N_act and _N_RS + _N_SS + 1 == _N_act:
        _total_dP = (_N_RS + 1) * _dP_RS + _N_SS * _dP_SS
    elif _N_act:
        # 分段板数缺失时按统一单板压降累加
        _total_dP = (result.get("pressure_drop_per_tray", 0) or 0) * _N_act
    else:
        _total_dP = 0
    _dP_avg = _total_dP / _N_act if _N_act else 0
    _P_bot = P + _total_dP
    T_bot_corrected = bubble_point_temperature(_P_bot, components, x_B_comp, T_guess=max(Tbs))
    T_top_K = T_top.get("Tb", 360)
    T_bot_K = T_bot_corrected.get("Tb", T_bot.get("Tb", 390))
    T_avg_corrected = (T_bot_K + T_top_K) / 2.0
    # 用修正后 T_bot 和 P_bot 更新提馏段物性
    rho_V_SS = _rho_V_sec(x_B_comp, T_bot_K, P_sec=_P_bot)
    rho_L_SS = _rho_L_sec(x_B_comp, T_bot_K)
    # 用修正后物性重算提馏段水力学
    hydraulics_SS = calculate_tray_hydraulics(
        V_mass=V_mass_SS, L_mass=L_mass_SS,
        rho_vapor=rho_V_SS, rho_liquid=rho_L_SS,
        sigma=None, mu_liquid=None, D_column=D_col_val,
        tray_spacing=tray_spacing, flood_method=flood_method,
    )
    _dP_SS = hydraulics_SS.get("pressure_drop_per_tray", 0) or 0
    # 用修正后 dP_SS 更新总压降和 P_bot（板数口径与首次计算一致）
    if _N_act and _N_RS + _N_SS + 1 == _N_act:
        _total_dP = (_N_RS + 1) * _dP_RS + _N_SS * _dP_SS
    elif _N_act:
        _total_dP = (result.get("pressure_drop_per_tray", 0) or 0) * _N_act
    if _N_act:
        _dP_avg = _total_dP / _N_act
        _P_bot = P + _total_dP
    q_auto_notes.append(
        f"[P_bot修正] T_bubble(x_B) 在 P_bot={_P_bot:.0f} Pa 下重新计算: "
        f"T_bot={T_bot_K:.2f} K (原 P_top={P:.0f} Pa 下 T_bot={T_bot.get('Tb',0):.2f} K)"
    )

    # === 热负荷（使用修正后的 T_bot） ===
    heat_duty_result = calculate_column_heat_duty(
        D_mol_s=D_mol_s, B_mol_s=B_mol_s, R_actual=R_actual,
        components=components, x_D=x_D_list, x_B=x_B_list,
        T_top=T_top_K, T_bottom=T_bot_K,
        q=q, P=P,
    )
    condenser_detail = heat_duty_result.get("condenser_result", {})
    reboiler_detail = heat_duty_result.get("reboiler_result", {})

    # 气液摩尔流量
    V_molar = V_mol  # mol/s  (已在前计算)
    L_molar = D_mol * R_actual  # mol/s
    F_mol = F_kmolh * 1000.0 / 3600.0  # mol/s

    # 塔截面面积
    D_col = result.get("D_column", 0)
    A_column = math.pi * D_col ** 2 / 4 if D_col > 0 else 0

    # 全塔总压降（加权：N_RS×dP_RS + N_SS×dP_SS）
    dP_per_tray_avg = _dP_avg
    N_act = result.get("N_actual_trays", 0)
    total_pressure_drop = _total_dP if _total_dP > 0 else None

    return {
        "material_balance": {
            "F_kmolh": round(F_kmolh, 3), "D_kmolh": round(D_kmolh, 3), "B_kmolh": round(B_kmolh, 3),
            "F_mol_s": round(F_mol, 4), "D_mol_s": round(D_kmolh * 1000 / 3600, 4), "B_mol_s": round(B_kmolh * 1000 / 3600, 4),
            "x_D_LK": x_D_actual, "x_B_LK": x_B_actual,
            "component_distribution": comp_dist,
        },
        "design_summary": {
            "N_min": result.get("N_min"), "R_min": round(result.get("R_min", 0), 4),
            "R": round(result.get("R_actual", 0), 4), "R_Rmin_ratio": round(RR_factor, 2),
            "N_theoretical": result.get("N_theoretical"),
            "N_actual": result.get("N_actual_trays"), "feed_stage": result.get("feed_stage"),
            "N_rectifying": result.get("N_rectifying"), "N_stripping": result.get("N_stripping"),
            "feed_stage_actual": result.get("feed_stage_actual"),
            "N_rectifying_actual": result.get("N_rectifying_actual"),
            "N_stripping_actual": result.get("N_stripping_actual"),
            "tray_count_note": "N_actual 为塔内物理塔板数（已扣除再沸器1级）；feed_stage/N_rectifying/N_stripping 为理论值（含再沸器级），*_actual 为实际值（均不含进料板，N_RS+1+N_SS=N_actual）",
            "D_column": result.get("D_column"),
            "D_column_calculated": diameter_result.get("diameter_calculated"),
            "H_column": result.get("H_column"),
            "H_D_ratio": result.get("H_D_ratio"),
            "tray_efficiency": eff_used, "efficiency_note": eff_note,
            "alpha_LK_HK": round(alpha_LK_HK, 2),
            "q_state": result.get("q_state", ""),
        },
        "flow_rates": {
            "rectifying": {
                "V_molar_mol_s": round(V_molar, 4),
                "L_molar_mol_s": round(L_molar, 4),
                "V_mass_kg_s": round(V_mass, 4),
                "L_mass_kg_s": round(L_mass, 4),
                "MW_g_mol": round(MW_RS, 1),
            },
            "stripping": {
                "V_molar_mol_s": round(V_mol_SS, 4),
                "L_molar_mol_s": round(L_mol_SS, 4),
                "V_mass_kg_s": round(V_mass_SS, 4),
                "L_mass_kg_s": round(L_mass_SS, 4),
                "MW_g_mol": round(MW_SS, 1),
            },
            "A_column_m2": round(A_column, 4),
            "D_column_m": round(D_col, 3),
        },
        "temperatures": {
            "top_K": round(T_top_K, 2),
            "top_C": round(T_top_K - 273.15, 1),
            "bottom_K": round(T_bot_K, 2),
            "bottom_C": round(T_bot_K - 273.15, 1),
            "average_K": round(T_avg_corrected, 2),
            "average_C": round(T_avg_corrected - 273.15, 1),
            "feed_bubble_K": round(T_feed_bub.get("Tb", 0), 2),
            "feed_dew_K": round(T_feed_dew.get("Td", 0), 2),
        },
        "properties": {
            "rectifying": {
                "rho_V_kg_m3": round(rho_V_RS, 3),
                "rho_L_kg_m3": round(rho_L_RS, 1),
                "MW_g_mol": round(MW_RS, 1),
                "composition_basis": "x_D",
            },
            "stripping": {
                "rho_V_kg_m3": round(rho_V_SS, 3),
                "rho_L_kg_m3": round(rho_L_SS, 1),
                "MW_g_mol": round(MW_SS, 1),
                "composition_basis": "x_B",
            },
            "feed": {
                "rho_V_kg_m3": round(rho_V_val, 3),
                "rho_L_kg_m3": round(rho_L_val, 1),
                "MW_g_mol": round(MW_avg, 1),
                "composition_basis": "z_feed",
            },
            "mu_L_cP": round(mu_L_cP, 3),
        },
        "fenske_underwood_gilliland": {
            "fenske": {
                "N_min": fenske_detail.get("N_min"),
                "log_alpha": round(fenske_detail.get("log_alpha_term", 0), 4),
            },
            "underwood": {
                "R_min": round(underwood_detail.get("R_min", 0), 4),
                "theta": round(underwood_detail.get("theta", 0), 4),
                "q": round(underwood_detail.get("q", q), 4),
            },
            "gilliland": {
                "N_theoretical": gilliland_detail.get("N_actual"),
                "R_ratio": round(gilliland_detail.get("R_ratio", RR_factor), 3),
            },
        },
        "hydraulics": {
            "rectifying": {
                "V_mass_kg_s": round(V_mass, 4),
                "rho_V_kg_m3": round(rho_V_RS, 4),
                "rho_L_kg_m3": round(rho_L_RS, 2),
                "A_col_m2": round(A_column, 4),
                "A_net_m2": round(A_column * 0.90, 4),
                "A_hole_m2": round(A_column * 0.10, 4),
                "u_operating_m_s": hydraulics_RS.get("vapor_velocity"),
                "u_flooding_m_s": hydraulics_RS.get("u_flood"),
                "u_hole_m_s": hydraulics_RS.get("hole_velocity"),
                "flooding_percent": hydraulics_RS.get("flooding_percent"),
                "F_factor": hydraulics_RS.get("F_factor"),
                "F_hole": hydraulics_RS.get("F_hole"),
                "pressure_drop_per_tray_Pa": hydraulics_RS.get("pressure_drop_per_tray"),
                "composition_basis": "x_D",
            },
            "stripping": {
                "V_mass_kg_s": round(V_mass_SS, 4),
                "rho_V_kg_m3": round(rho_V_SS, 4),
                "rho_L_kg_m3": round(rho_L_SS, 2),
                "A_col_m2": round(A_column, 4),
                "A_net_m2": round(A_column * 0.90, 4),
                "A_hole_m2": round(A_column * 0.10, 4),
                "u_operating_m_s": hydraulics_SS.get("vapor_velocity"),
                "u_flooding_m_s": hydraulics_SS.get("u_flood"),
                "u_hole_m_s": hydraulics_SS.get("hole_velocity"),
                "flooding_percent": hydraulics_SS.get("flooding_percent"),
                "F_factor": hydraulics_SS.get("F_factor"),
                "F_hole": hydraulics_SS.get("F_hole"),
                "pressure_drop_per_tray_Pa": hydraulics_SS.get("pressure_drop_per_tray"),
                "composition_basis": "x_B (P_bot 修正后)",
            },
            "C_sb": diameter_result.get("C_factor"),
            "C_sb_note": "base coefficient (no surface tension correction); same C_sb used for all sections",
            "area_basis": "u_op/u_flood/FF 基于净面积 A_net = 0.90×A_col（扣除降液管）；u_hole 基于开孔面积 A_hole = 0.10×A_col",
            "net_area_fraction": 0.90,
            "total_pressure_drop_Pa": round(total_pressure_drop, 0) if total_pressure_drop else None,
            "dP_breakdown": f"({round(_N_RS)}+1)×{round(_dP_RS,0)} + {round(_N_SS)}×{round(_dP_SS,0)} = {round(_total_dP,0)} Pa (进料板计入精馏段侧; hydraulics only, no condenser/reboiler losses)",
            "P_bot_Pa": round(_P_bot, 0),
            "weeping_rectifying": hydraulics_RS.get("weeping_check", {}).get("weeping"),
            "weeping_stripping": hydraulics_SS.get("weeping_check", {}).get("weeping"),
        },
        "column_sizing": {
            "governing_section": _gov_section,
            "D_calc_rectifying_m": round(_dia_RS.get("diameter_calculated", 0) or 0, 4),
            "D_calc_stripping_m": round(_dia_SS.get("diameter_calculated", 0) or 0, 4),
            "u_max_governing_m_s": diameter_result.get("u_max"),
            "u_design_m_s": diameter_result.get("u_design"),
            "u_design_ratio": flooding_fraction,
            "basis": f"控制段={_gov_section}: 以该段 V_mass/ρV/ρL 定径 (Souders-Brown, 净面积 A_net=0.90×A_col 基准)",
            "rho_V_governing_kg_m3": round(rho_V_design, 4),
            "rho_L_governing_kg_m3": round(rho_L_design, 2),
            "note": "各段用各自 V_mass/ρV/ρL 独立核算所需塔径，取较大者（控制段）定径；u_design = u_max × flooding_fraction",
        },
        "heat_duty": {
            "Q_condenser_W": round(heat_duty_result.get("Q_condenser_W", 0), 0),
            "Q_condenser_kW": round(heat_duty_result.get("Q_condenser_kW", 0), 2),
            "Q_reboiler_W": round(heat_duty_result.get("Q_reboiler_W", 0), 0),
            "Q_reboiler_kW": round(heat_duty_result.get("Q_reboiler_kW", 0), 2),
            "Q_total_W": round(heat_duty_result.get("Q_total_W", 0), 0),
            "Q_total_kW": round(heat_duty_result.get("Q_total_kW", 0), 2),
            "condenser": {
                "V_top_mol_s": round(condenser_detail.get("V_top_mol_s", 0), 4),
                "L_top_mol_s": round(condenser_detail.get("L_top_mol_s", 0), 4),
                "Hvap_mix_J_mol": round(condenser_detail.get("Hvap_mix_J_mol", 0), 1),
            },
            "reboiler": {
                "V_strip_mol_s": round(reboiler_detail.get("V_strip_mol_s", 0), 4),
                "L_strip_mol_s": round(reboiler_detail.get("L_strip_mol_s", 0), 4),
                "Hvap_mix_J_mol": round(reboiler_detail.get("Hvap_mix_J_mol", 0), 1),
            },
        },
        "validation": {
            "checks": {
                "N_min_positive": {"value": result.get("N_min"), "pass": (result.get("N_min", 0) or 0) > 0, "criterion": "N_min > 0"},
                "R_ratio": {"value": round(RR_factor, 2), "pass": 1.1 <= RR_factor <= 3.0, "criterion": "1.1 ≤ R/R_min ≤ 3.0"},
                "efficiency": {"value": round(eff_used * 100, 1), "pass": 30 <= eff_used * 100 <= 90, "criterion": "E_o = 30-90%"},
                "flooding": {
                    "value": max(hydraulics_RS.get("flooding_percent", 0), hydraulics_SS.get("flooding_percent", 0)),
                    "section": "rectifying" if (hydraulics_RS.get("flooding_percent", 0) or 0) >= (hydraulics_SS.get("flooding_percent", 0) or 0) else "stripping",
                    "pass": max(hydraulics_RS.get("flooding_percent", 0) or 0, hydraulics_SS.get("flooding_percent", 0) or 0) <= 85,
                    "criterion": "Flooding fraction ≤ 85%",
                },
                "pressure_drop": {
                    "value": max(hydraulics_RS.get("pressure_drop_per_tray", 0) or 0, hydraulics_SS.get("pressure_drop_per_tray", 0) or 0),
                    "section": "rectifying" if (hydraulics_RS.get("pressure_drop_per_tray", 0) or 0) >= (hydraulics_SS.get("pressure_drop_per_tray", 0) or 0) else "stripping",
                    "pass": max(hydraulics_RS.get("pressure_drop_per_tray", 0) or 0, hydraulics_SS.get("pressure_drop_per_tray", 0) or 0) <= 1500,
                    "criterion": "Tray ΔP ≤ 1500 Pa",
                },
                "weeping": {
                    "value": min(hydraulics_RS.get("F_hole", 999) or 999, hydraulics_SS.get("F_hole", 999) or 999),
                    "section": "rectifying" if (hydraulics_RS.get("F_hole", 999) or 999) <= (hydraulics_SS.get("F_hole", 999) or 999) else "stripping",
                    "pass": min(hydraulics_RS.get("F_hole", 999) or 999, hydraulics_SS.get("F_hole", 999) or 999) >= 5,
                    "criterion": "F_hole ≥ 5 (no weeping)",
                },
                "H_D_ratio": {"value": round(result.get("H_D_ratio", 0), 2), "pass": 2 <= (result.get("H_D_ratio", 0) or 0) <= 30, "criterion": "H/D = 2-30"},
                "column_diameter": {"value": result.get("D_column"), "pass": (result.get("D_column", 0) or 0) >= 0.3, "criterion": "D ≥ 0.3 m"},
            },
            "n_pass": sum(1 for c in [
                (result.get("N_min", 0) or 0) > 0,
                1.1 <= RR_factor <= 3.0,
                30 <= eff_used * 100 <= 90,
                max(hydraulics_RS.get("flooding_percent", 0) or 0, hydraulics_SS.get("flooding_percent", 0) or 0) <= 85,
                max(hydraulics_RS.get("pressure_drop_per_tray", 0) or 0, hydraulics_SS.get("pressure_drop_per_tray", 0) or 0) <= 1500,
                min(hydraulics_RS.get("F_hole", 999) or 999, hydraulics_SS.get("F_hole", 999) or 999) >= 5,
                2 <= (result.get("H_D_ratio", 0) or 0) <= 30,
                (result.get("D_column", 0) or 0) >= 0.3,
            ] if c),
            "n_total": 8,
        },
        "warnings": result.get("warnings", []) + q_auto_notes + heat_duty_result.get("warnings", []) + (
            [f"[T_bot修正] 泡点在 P_bot={_P_bot:.0f} Pa 下重新计算: T_bot={T_bot_K:.2f} K"]
        ),
        "converged": result.get("valid", True),
        "q_info": {"q": round(q, 4), "source": "user" if not q_auto_notes else "auto"},
    }


def register_distillation_device(registry: ToolRegistry) -> ToolRegistry:
    registry.register_function(name="distillation_column_design",
        description="精馏塔完整设计 (FUG法)。若遇到液泛问题需要调整塔径，请调小 flooding_fraction，严禁传入 D_column。",
        func=_full_design,
        param_descriptions={
            "components": "(必填) 组分名称列表，如 ['benzene', 'toluene']",
            "z": "(必填) 进料摩尔分数列表，必须与 components 对应",
            "F": "(必填) 进料流量 (kmol/h)，无默认值",
            "P": "(必填) 操作压力 (Pa)，无默认值（常压为 101325）",
            "x_D_spec": "(必填) 塔顶轻关键组分摩尔分数要求 (0~1)，无默认值",
            "x_B_spec": "(必填) 塔釜轻关键组分摩尔分数要求 (0~1)，无默认值",
            "q": "(选填) 进料热状态 (1=饱和液体, 0=饱和蒸气)。不提供时自动根据 T_feed 计算",
            "T_feed": "(选填) 进料温度 (K)。当 q 未提供时，用于自动计算 q 值",
            "RR_factor": "(选填) 回流比倍数 (R/Rmin)，默认 2.0，可选: 1.2/1.5/2.0/2.5/3.0",
            "tray_type": "(选填) 塔板类型，默认 'sieve'(筛板)，可选: 'sieve'/'valve'(浮阀)",
            "tray_efficiency": "(选填) 塔板效率 (0-1)，不填则自动按 O'Connell 法估算",
            "flooding_fraction": "(选填) 目标液泛分率，默认 0.80，可选: 0.70/0.75/0.80/0.82/0.85",
            "tray_spacing": "(选填) 塔板间距(m)，默认 0.60，可选: 0.45/0.50/0.60/0.75/0.90",
            "flood_method": "(选填) 液泛关联式: 'simple'(简化式，默认) 或 'fair'(Fair 1961 通用关联式，含 F_LV 修正)"
        },
        category="device", tags=["distillation", "FUG", "精馏塔"])
    return registry


def register_to_universal(registry: UniversalToolRegistry) -> None:
    """注册到通用注册中心 — 自动提取能力描述"""
    capability = ToolCapability(
        name="distillation_column_design",
        description="精馏塔完整设计 (FUG法)。可计算理论板数、回流比、进料位置、塔径等。若遇液泛需调大塔径，请调低 flooding_fraction，严禁传入 D_column。",
        category="device_design",
        parameters={
            "components": {"type": "list", "required": True, "description": "(必填) 组分英文名列表，如 ['benzene', 'toluene']"},
            "z": {"type": "list", "required": True, "description": "(必填) 进料摩尔分数列表，如 [0.4, 0.6]"},
            "F": {"type": "number", "required": True, "description": "(必填) 进料流量 kmol/h，无默认值"},
            "P": {"type": "number", "required": True, "description": "(必填) 操作压力 Pa，无默认值（常压=101325）"},
            "x_D_spec": {"type": "number", "required": True, "description": "(必填) 塔顶轻关键组分要求 0~1，无默认值"},
            "x_B_spec": {"type": "number", "required": True, "description": "(必填) 塔釜轻关键组分要求 0~1，无默认值"},
            "q": {"type": "number", "required": False, "description": "(选填) 进料热状态: 1=饱和液体, 0=饱和蒸气。不提供时自动根据T_feed计算"},
            "T_feed": {"type": "number", "required": False, "description": "(选填) 进料温度 K。当q未提供时用于自动计算q值"},
            "RR_factor": {"type": "number", "required": False, "description": "(选填) 回流比倍数 R/Rmin，默认 2.0，可选: 1.2/1.5/2.0/2.5/3.0"},
            "tray_type": {"type": "string", "required": False, "description": "(选填) 塔板类型: 'sieve'(筛板)/'valve'(浮阀)，默认 sieve"},
            "tray_efficiency": {"type": "number", "required": False, "description": "(选填) 塔板效率 (0-1)，不填则按 O'Connell 自动估算"},
            "flooding_fraction": {"type": "number", "required": False, "description": "(选填) 设计液泛分率，默认 0.80，可选: 0.70/0.75/0.80/0.82/0.85"},
            "tray_spacing": {"type": "number", "required": False, "description": "(选填) 塔板间距 m，默认 0.60，可选: 0.45/0.50/0.60/0.75/0.90"},
            "flood_method": {"type": "string", "required": False, "description": "(选填) 液泛关联式: 'simple'(简化式，默认) 或 'fair'(Fair 1961 通用关联式，含 F_LV 修正)"},
            # ── 可选物性参数 ──
            "alpha": {"type": "number", "required": False, "description": "(选填) 相对挥发度，给定后跳过内部计算"},
            "rho_V": {"type": "number", "required": False, "description": "(选填) 气相密度 kg/m³，给定后跳过物性查询"},
            "rho_L": {"type": "number", "required": False, "description": "(选填) 液相密度 kg/m³，给定后跳过物性查询"},
            "mu_L": {"type": "number", "required": False, "description": "(选填) 液相粘度 cP，给定后跳过粘度估算"},
            "Cp_feed": {"type": "number", "required": False, "description": "(选填) 进料比热容 J/(mol·K)，给定后跳过物性查询"},
        },
        examples=[
            {
                "description": "苯-甲苯精馏塔设计",
                "input": {
                    "components": ["benzene", "toluene"],
                    "z": [0.4, 0.6],
                    "F": 100,
                    "P": 101325,
                    "T_feed": 373.15,
                },
            },
        ],
        tags=["distillation", "FUG", "精馏塔", "分离"],
    )
    registry.register(capability, _full_design)