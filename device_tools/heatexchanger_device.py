"""换热器设备级 MCP 工具 — 完整管壳式换热器设计流水线"""
from physics_engine.heatexchanger import (
    design_heat_exchanger, lmtd_with_correction, heat_duty_from_flow,
    estimate_U, suggest_configuration, fouling_resistance,
    get_default_coolant,
)
from physics_engine.thermo_helper import (
    get_fluid_density, get_fluid_viscosity, get_fluid_Cp, get_fluid_thermal_conductivity,
    get_molecular_weight,
)
from physics_engine.common import _REQUIRED, check_required_params
from tools.base import ToolRegistry
from tools.registry import UniversalToolRegistry, ToolCapability

# ============================================================
# 换热管标准尺寸集（工程常用值）
# ============================================================
# 换热管外径 (m) — 国标/ASME 常用规格
STANDARD_TUBE_OD = [0.015, 0.019, 0.025, 0.032, 0.038]  # 15/19/25/32/38 mm
DEFAULT_TUBE_OD = 0.019  # 19 mm (¾" OD), 工业最常用

# 换热管长度 (m) — 常用标准长度
STANDARD_TUBE_LENGTH = [2.0, 3.0, 4.0, 6.0]  # m
DEFAULT_TUBE_LENGTH = 3.0  # 3 m, 工业最常用

_MOLAR_MASS = {
    "water": 0.018015, "methanol": 0.03204, "ethanol": 0.04607,
    "benzene": 0.07811, "toluene": 0.09214, "acetone": 0.05808,
    "ethylene": 0.02805, "propane": 0.04410, "butane": 0.05812,
    "ammonia": 0.01703, "co2": 0.04401, "nitrogen": 0.02802, "air": 0.02897,
}

# ============================================================
# 行业标准流速约束 (GB/T 151, SH/T 3121, TEMA)
# ============================================================
# 管程流速约束 (m/s)
TUBE_VELOCITY_CONSTRAINTS = {
    "liquid_min": 0.5,       # 液体最低流速 (防沉积)
    "liquid_optimal": 1.0,   # 液体最优流速
    "liquid_max": 3.0,       # 液体最高流速 (防冲蚀)
    "gas_min": 5.0,          # 气体最低流速
    "gas_optimal": 10.0,     # 气体最优流速
    "gas_max": 30.0,         # 气体最高流速
}

# 壳程流速约束 (m/s)
SHELL_VELOCITY_CONSTRAINTS = {
    "liquid_min": 0.3,       # 液体最低流速
    "liquid_optimal": 0.8,   # 液体最优流速
    "liquid_max": 2.0,       # 液体最高流速
    "gas_min": 3.0,          # 气体最低流速
    "gas_optimal": 8.0,      # 气体最优流速
    "gas_max": 20.0,         # 气体最高流速
}

# ============================================================
# 折流板约束 (GB/T 151, TEMA)
# ============================================================
BAFFLE_CONSTRAINTS = {
    "B_D_min": 0.2,          # 折流板间距/壳径 最小值
    "B_D_max": 1.0,          # 折流板间距/壳径 最大值
    "B_D_optimal_min": 0.3,  # 最优范围下限
    "B_D_optimal_max": 0.6,  # 最优范围上限
}

# ============================================================
# 热膨胀约束
# ============================================================
THERMAL_EXPANSION_CONSTRAINTS = {
    "alpha_steel": 1.2e-5,   # 碳钢线膨胀系数 (1/K)
    "alpha_ss": 1.6e-5,      # 不锈钢线膨胀系数 (1/K)
    "alpha_cu": 1.7e-5,      # 铜线膨胀系数 (1/K)
    "expansion_joint_max_mm": 10.0,  # 膨胀节最大补偿量 (mm)
}

_FLUID_ALIAS = {
    "甲醇": "methanol", "水": "water", "乙醇": "ethanol", "酒精": "ethanol",
    "苯": "benzene", "甲苯": "toluene", "丙酮": "acetone", "乙烯": "ethylene",
    "丙烷": "propane", "丁烷": "butane", "氨": "ammonia", "二氧化碳": "co2",
    "氮": "nitrogen", "氮气": "nitrogen", "空气": "air",
    "冷却水": "water", "清水": "water", "纯水": "water",
}

def _normalize_fluid_name(fluid: str) -> str:
    if not isinstance(fluid, str) or not fluid.strip():
        return ""  # 空字符串将在后续物性查询中报错，不静默回退
    fluid = fluid.strip()
    return _FLUID_ALIAS.get(fluid, fluid.lower())

def _celsius_to_kelvin(t_c: float) -> float:
    if t_c > 200: return float(t_c)
    return float(t_c) + 273.15

def _ensure_kelvin(value, param_name: str, notes: list) -> float:
    value = float(value)
    if value <= 200:
        converted = value + 273.15
        notes.append(f"[单位换算] {param_name}: {value} °C → {converted:.2f} K")
        return converted
    return value

def _ensure_pascal(value, param_name: str, notes: list) -> float:
    value = float(value)
    if value <= 0:
        notes.append(f"[参数错误] {param_name}: 压力值 {value} 无效，压力必须为正数")
        return 0.0  # 将在后续热物性查询时触发异常，不静默回退
    if value <= 500:
        converted = value * 100000.0
        notes.append(f"[单位换算] {param_name}: {value} bar → {converted:.0f} Pa")
        return converted
    if value <= 2000:
        converted = value * 1000.0
        notes.append(f"[单位换算] {param_name}: {value} kPa → {converted:.0f} Pa")
        return converted
    return value

def _ensure_mass_flow(value, fluid: str, param_name: str, notes: list) -> float:
    value = float(value)
    if value <= 0:
        notes.append(f"[参数错误] {param_name}: 流量值 {value} 无效，流量必须为正数")
        return 0.0  # 将在后续热物性查询时触发异常，不静默回退
    if value > 1000:
        converted = value / 3600.0
        notes.append(f"[单位换算] {param_name}: {value} kg/h → {converted:.4f} kg/s（推测单位）")
        return converted
    return value

def _validate_temperature_cross(T_hot_in: float, T_hot_out: float,
                                  T_cold_in: float, T_cold_out: float,
                                  notes: list) -> dict:
    warnings = []
    if T_hot_in <= T_cold_in:
        warnings.append(
            f"热侧入口 {T_hot_in:.1f} K ≤ 冷侧入口 {T_cold_in:.1f} K，存在温度交叉！请检查冷/热侧参数是否颠倒。"
        )
    if T_hot_out <= T_cold_in:
        warnings.append(
            f"热侧出口 {T_hot_out:.1f} K ≤ 冷侧入口 {T_cold_in:.1f} K，换热驱动力不足，建议降低冷侧出口温度或提高热侧入口温度。"
        )
    return {"valid": len(warnings) == 0, "warnings": warnings}


def _full_he_design(
    T_hot_in: float = _REQUIRED,
    T_hot_out: float = _REQUIRED,
    T_cold_in: float = _REQUIRED,
    m_hot: float = _REQUIRED,
    T_cold_out: float = None,
    m_cold: float = None,
    hot_fluid: str = _REQUIRED,
    cold_fluid: str = _REQUIRED,
    P_hot: float = _REQUIRED,
    P_cold: float = _REQUIRED,
    d_o: float = None,
    L: float = None,
    tube_layout: str = "auto",
    # ── 可选：用户给定物性/设计参数（给定后跳过内部计算/查询） ──
    U: float = None,
    K: float = None,            # U 的别名（中国工程习惯用 K）
    Cp_hot: float = None,       # 热侧比热容 J/(kg·K)
    Cp_cold: float = None,      # 冷侧比热容 J/(kg·K)
    rho_hot: float = None,      # 热侧密度 kg/m³
    rho_cold: float = None,     # 冷侧密度 kg/m³
    mu_hot: float = None,       # 热侧粘度 Pa·s
    mu_cold: float = None,      # 冷侧粘度 Pa·s
    k_hot: float = None,        # 热侧导热系数 W/(m·K)
    k_cold: float = None,       # 冷侧导热系数 W/(m·K)
    **kwargs,
) -> dict:
    # 检查必填参数（核心工艺参数必须由用户提供，d_o/L 有工程默认值）
    missing = check_required_params(locals(), {
        "T_hot_in": "热流体入口温度 (K 或 °C，≤200 自动按 °C 换算)",
        "T_hot_out": "热流体出口温度 (K 或 °C)",
        "T_cold_in": "冷流体入口温度 (K 或 °C)",
        "m_hot": "热流体质量流量 (kg/s 或 kg/h，>1000 自动按 kg/h 换算)",
        "hot_fluid": "热流体名称，如 'water', 'methanol', 'benzene'",
        "cold_fluid": "冷流体名称，如 'water'(冷却水), 'methanol'",
        "P_hot": "热侧操作压力 (Pa，或 bar/kPa 自动换算)",
        "P_cold": "冷侧操作压力 (Pa，或 bar/kPa 自动换算)",
    })
    if missing:
        return missing
    conversion_notes: list[str] = []
    auto_coolant_info: dict = {}
    _user_defaulted_params = {}

    # 检测无效的嵌套字典参数（planner 可能错误传递）
    _invalid_nested_params = [k for k in kwargs if isinstance(kwargs.get(k), dict)]
    if _invalid_nested_params:
        conversion_notes.append(
            f"[参数警告] 检测到无效嵌套字典参数: {_invalid_nested_params}。"
            f"工具不接受嵌套字典，请使用扁平参数如 Cp_hot/rho_hot/mu_hot/k_hot。"
        )

    hot_fluid  = _normalize_fluid_name(hot_fluid)
    cold_fluid_raw = cold_fluid
    cold_fluid = _normalize_fluid_name(cold_fluid)

    # 记录用户是否显式提供了冷侧参数
    _user_provided_T_cold_out = T_cold_out is not None
    _user_provided_m_cold = m_cold is not None

    T_hot_in  = _ensure_kelvin(T_hot_in,  "T_hot_in",  conversion_notes)
    T_hot_out = _ensure_kelvin(T_hot_out, "T_hot_out", conversion_notes)
    T_cold_in = _ensure_kelvin(T_cold_in, "T_cold_in", conversion_notes)
    if _user_provided_T_cold_out:
        T_cold_out = _ensure_kelvin(T_cold_out, "T_cold_out", conversion_notes)

    P_hot  = _ensure_pascal(P_hot,  "P_hot",  conversion_notes)
    P_cold = _ensure_pascal(P_cold, "P_cold", conversion_notes)

    m_hot  = _ensure_mass_flow(m_hot,  hot_fluid,  "m_hot",  conversion_notes)
    if _user_provided_m_cold:
        m_cold = _ensure_mass_flow(m_cold, cold_fluid, "m_cold", conversion_notes)

    # d_o / L 默认值处理（工程标准值集，用户未提供时使用默认并追踪）
    if d_o is None or d_o is _REQUIRED:
        d_o = DEFAULT_TUBE_OD
        _user_defaulted_params["d_o"] = f"{d_o} m ({int(d_o*1000)} mm，标准可选值: {[f'{v*1000:.0f}mm' for v in STANDARD_TUBE_OD]})"
        conversion_notes.append(f"[默认值] d_o 使用工程标准值 {d_o} m ({int(d_o*1000)} mm)，标准规格: {STANDARD_TUBE_OD}")
    d_o = float(d_o)

    if L is None or L is _REQUIRED:
        L = DEFAULT_TUBE_LENGTH
        _user_defaulted_params["L"] = f"{L} m（标准可选值: {STANDARD_TUBE_LENGTH} m）"
        conversion_notes.append(f"[默认值] L 使用工程标准值 {L} m，标准长度: {STANDARD_TUBE_LENGTH}")
    L = float(L)

    # ─── 冷却介质自动选择：仅当用户完全未指定冷侧参数时触发 ───
    _cold_unspecified = (
        not _user_provided_T_cold_out   # 用户未指定 T_cold_out
        and not _user_provided_m_cold   # 用户未指定 m_cold
        and abs(T_cold_in - 300.0) < 5.0
        and cold_fluid_raw.strip().lower() in ("water", "冷却水", "循环水", "")
    )
    if _cold_unspecified:
        coolant = get_default_coolant(T_hot_out_K=T_hot_out)
        T_cold_in  = coolant["T_cold_in_K"]
        T_cold_out = coolant["T_cold_out_K"]
        cold_fluid = coolant["fluid_name"]
        P_cold     = coolant["P_cold_Pa"]
        auto_coolant_info = coolant
        conversion_notes.append(
            f"[冷却介质自动选择] 未指定冷侧参数，根据热侧出口温度 {T_hot_out:.1f} K"
            f" 自动选择：{coolant['coolant_type']}。"
            f" 冷侧：{T_cold_in:.1f} K → {T_cold_out:.1f} K"
            f" 说明：{coolant['note']}"
        )

    cross_check = _validate_temperature_cross(
        T_hot_in, T_hot_out, T_cold_in, T_cold_out, conversion_notes
    )
    all_warnings: list[str] = list(cross_check["warnings"])

    # T_cold_out 未提供时估算（仅在自动冷却介质未触发时生效）
    if T_cold_out is None or abs(T_cold_out - T_cold_in) < 0.1:
        delta_hot = abs(T_hot_in - T_hot_out)
        T_cold_out = T_cold_in + delta_hot * 0.8
        conversion_notes.append(
            f"[参数估算] T_cold_out 未提供，估算为 {T_cold_out:.2f} K"
        )
        all_warnings.append(
            f"冷侧出口温度未提供，自动估算为 {T_cold_out:.2f} K，请确认"
        )

    T_h_avg = (T_hot_in  + T_hot_out)  / 2
    T_c_avg = (T_cold_in + T_cold_out) / 2

    # ── 物性参数：用户给定值优先，缺失时才查询/计算 ──
    _user_provided_props = []  # 记录用户显式提供的物性参数
    MW_h = None  # 分子量（仅当需要 Cp 单位换算时才查询）
    MW_c = None

    # 比热容 Cp (J/(kg·K))
    if Cp_hot is not None:
        # 自动单位换算：如果值 < 100，推测单位为 kJ/(kg·K)，转为 J/(kg·K)
        if Cp_hot < 100:
            Cp_h = Cp_hot * 1000.0
            conversion_notes.append(f"[用户给定] Cp_hot = {Cp_hot} kJ/(kg·K) → {Cp_h} J/(kg·K)")
        else:
            Cp_h = float(Cp_hot)
            conversion_notes.append(f"[用户给定] Cp_hot = {Cp_h} J/(kg·K)")
        _user_provided_props.append("Cp_hot")
    if Cp_cold is not None:
        if Cp_cold < 100:
            Cp_c = Cp_cold * 1000.0
            conversion_notes.append(f"[用户给定] Cp_cold = {Cp_cold} kJ/(kg·K) → {Cp_c} J/(kg·K)")
        else:
            Cp_c = float(Cp_cold)
            conversion_notes.append(f"[用户给定] Cp_cold = {Cp_c} J/(kg·K)")
        _user_provided_props.append("Cp_cold")

    # 密度 rho (kg/m³)
    if rho_hot is not None:
        rho_h = float(rho_hot)
        conversion_notes.append(f"[用户给定] rho_hot = {rho_h} kg/m³")
        _user_provided_props.append("rho_hot")
    if rho_cold is not None:
        rho_c = float(rho_cold)
        conversion_notes.append(f"[用户给定] rho_cold = {rho_c} kg/m³")
        _user_provided_props.append("rho_cold")

    # 粘度 mu (Pa·s)
    if mu_hot is not None:
        mu_h = float(mu_hot)
        conversion_notes.append(f"[用户给定] mu_hot = {mu_h} Pa·s")
        _user_provided_props.append("mu_hot")
    if mu_cold is not None:
        mu_c = float(mu_cold)
        conversion_notes.append(f"[用户给定] mu_cold = {mu_c} Pa·s")
        _user_provided_props.append("mu_cold")

    # 导热系数 k (W/(m·K))
    if k_hot is not None:
        k_h = float(k_hot)
        conversion_notes.append(f"[用户给定] k_hot = {k_h} W/(m·K)")
        _user_provided_props.append("k_hot")
    if k_cold is not None:
        k_c = float(k_cold)
        conversion_notes.append(f"[用户给定] k_cold = {k_c} W/(m·K)")
        _user_provided_props.append("k_cold")

    # 查询未提供的物性参数（每个属性独立查询+独立容错+工程兜底默认值）
    _any_query_needed = (
        'Cp_hot' not in _user_provided_props or
        'Cp_cold' not in _user_provided_props or
        'rho_hot' not in _user_provided_props or
        'rho_cold' not in _user_provided_props or
        'mu_hot' not in _user_provided_props or
        'mu_cold' not in _user_provided_props or
        'k_hot' not in _user_provided_props or
        'k_cold' not in _user_provided_props
    )

    # ── 水的工程默认物性（查询失败时的兜底值） ──
    _WATER_DEFAULTS = {
        "Cp": 4180.0,       # J/(kg·K) — 水在 300-373K 的典型值
        "rho": 980.0,       # kg/m³ — 水在 300-373K 的典型值
        "mu": 0.0005,       # Pa·s — 水在 340K 附近的典型值
        "k": 0.62,          # W/(m·K) — 水在 340K 附近的典型值
    }
    _ORGANIC_LIQUID_DEFAULTS = {
        "Cp": 2200.0,       # J/(kg·K) — 有机液体典型值
        "rho": 800.0,       # kg/m³
        "mu": 0.001,        # Pa·s
        "k": 0.15,          # W/(m·K)
    }

    def _get_defaults(fluid_name: str) -> dict:
        """根据流体名称返回工程默认物性"""
        fn = fluid_name.lower().strip()
        if fn in ("water", "冷却水", "循环水", "纯水", "清水"):
            return _WATER_DEFAULTS
        return _ORGANIC_LIQUID_DEFAULTS

    if _any_query_needed:
        # 每个物性独立查询，互不影响
        if 'Cp_hot' not in _user_provided_props:
            try:
                Cp_h_molar = get_fluid_Cp(hot_fluid, T_h_avg, P_hot, phase="liquid")
            except Exception as e:
                defaults = _get_defaults(hot_fluid)
                Cp_h_molar = defaults["Cp"] * 0.018015  # J/(kg·K) → J/(mol·K) 近似
                conversion_notes.append(f"[兜底默认值] Cp_hot 查询失败({e})，使用 {defaults['Cp']} J/(kg·K)")

        if 'Cp_cold' not in _user_provided_props:
            try:
                Cp_c_molar = get_fluid_Cp(cold_fluid, T_c_avg, P_cold, phase="liquid")
            except Exception as e:
                defaults = _get_defaults(cold_fluid)
                Cp_c_molar = defaults["Cp"] * 0.018015
                conversion_notes.append(f"[兜底默认值] Cp_cold 查询失败({e})，使用 {defaults['Cp']} J/(kg·K)")

        if 'rho_hot' not in _user_provided_props:
            try:
                rho_h = get_fluid_density(hot_fluid, T_h_avg, P_hot, phase="liquid")
            except Exception as e:
                defaults = _get_defaults(hot_fluid)
                rho_h = defaults["rho"]
                conversion_notes.append(f"[兜底默认值] rho_hot 查询失败({e})，使用 {rho_h} kg/m³")

        if 'rho_cold' not in _user_provided_props:
            try:
                rho_c = get_fluid_density(cold_fluid, T_c_avg, P_cold, phase="liquid")
            except Exception as e:
                defaults = _get_defaults(cold_fluid)
                rho_c = defaults["rho"]
                conversion_notes.append(f"[兜底默认值] rho_cold 查询失败({e})，使用 {rho_c} kg/m³")

        if 'mu_hot' not in _user_provided_props:
            try:
                mu_h = get_fluid_viscosity(hot_fluid, T_h_avg, P_hot, phase="liquid")
            except Exception as e:
                defaults = _get_defaults(hot_fluid)
                mu_h = defaults["mu"]
                conversion_notes.append(f"[兜底默认值] mu_hot 查询失败({e})，使用 {mu_h} Pa·s")

        if 'mu_cold' not in _user_provided_props:
            try:
                mu_c = get_fluid_viscosity(cold_fluid, T_c_avg, P_cold, phase="liquid")
            except Exception as e:
                defaults = _get_defaults(cold_fluid)
                mu_c = defaults["mu"]
                conversion_notes.append(f"[兜底默认值] mu_cold 查询失败({e})，使用 {mu_c} Pa·s")

        if 'k_hot' not in _user_provided_props:
            try:
                k_h = get_fluid_thermal_conductivity(hot_fluid, T_h_avg, P_hot, phase="liquid")
            except Exception as e:
                defaults = _get_defaults(hot_fluid)
                k_h = defaults["k"]
                conversion_notes.append(f"[兜底默认值] k_hot 查询失败({e})，使用 {k_h} W/(m·K)")

        if 'k_cold' not in _user_provided_props:
            try:
                k_c = get_fluid_thermal_conductivity(cold_fluid, T_c_avg, P_cold, phase="liquid")
            except Exception as e:
                defaults = _get_defaults(cold_fluid)
                k_c = defaults["k"]
                conversion_notes.append(f"[兜底默认值] k_cold 查询失败({e})，使用 {k_c} W/(m·K)")

    # Cp 单位转换: J/(mol·K) → J/(kg·K)（仅当用户未直接给定质量比热时）
    if 'Cp_hot' not in _user_provided_props or 'Cp_cold' not in _user_provided_props:
        # 仅查询需要转换的流体的分子量
        if 'Cp_hot' not in _user_provided_props:
            try:
                MW_h = get_molecular_weight(hot_fluid)
            except Exception:
                MW_h = _MOLAR_MASS.get(hot_fluid, 18.015)  # 水为默认
                conversion_notes.append(f"[兜底默认值] MW_hot 查询失败，使用 {MW_h} g/mol")
        if 'Cp_cold' not in _user_provided_props:
            try:
                MW_c = get_molecular_weight(cold_fluid)
            except Exception:
                MW_c = _MOLAR_MASS.get(cold_fluid, 18.015)
                conversion_notes.append(f"[兜底默认值] MW_cold 查询失败，使用 {MW_c} g/mol")
        if 'Cp_hot' not in _user_provided_props and MW_h:
            Cp_h = Cp_h_molar / (MW_h / 1000.0)  # J/(kg·K)
        if 'Cp_cold' not in _user_provided_props and MW_c:
            Cp_c = Cp_c_molar / (MW_c / 1000.0)  # J/(kg·K)

    # m_cold 未提供时由热平衡估算: m_cold = Q_hot / (Cp_c * ΔT_cold)
    if m_cold is None:
        Q_hot_est = abs(float(Cp_h) * m_hot * (T_hot_in - T_hot_out))
        delta_T_cold = abs(T_cold_out - T_cold_in)
        if delta_T_cold < 0.1 or Cp_c < 1.0:
            return {
                "error": f"无法估算冷侧流量：ΔT_cold={delta_T_cold:.2f}K 或 Cp_c={Cp_c:.1f} J/(kg·K) 异常",
                "conversion_notes": conversion_notes,
                "converged": False,
            }
        m_cold = Q_hot_est / (float(Cp_c) * delta_T_cold)
        conversion_notes.append(
            f"[流量估算] m_cold 未提供，由热平衡估算为 {m_cold:.4f} kg/s"
            f" (Q={Q_hot_est:.0f} W, ΔT_cold={delta_T_cold:.1f} K)"
        )

    Q_h = abs(float(Cp_h) * m_hot  * (T_hot_in  - T_hot_out))
    Q_c = abs(float(Cp_c) * m_cold * (T_cold_out - T_cold_in))

    if Q_h > 1e-6 and Q_c > 1e-6:
        imbalance = abs(Q_h - Q_c) / max(Q_h, Q_c) * 100
        if imbalance > 20:
            all_warnings.append(f"热平衡偏差 {imbalance:.1f}%，可能存在参数错误，建议核实流量。")

    try:
        lmtd_result = lmtd_with_correction(
            T_hot_in=T_hot_in, T_hot_out=T_hot_out,
            T_cold_in=T_cold_in, T_cold_out=T_cold_out,
        )
        delta_T_m = lmtd_result.get("LMTD") or lmtd_result.get("LMTD_corrected") or 10.0
    except Exception:
        delta_T_m = abs(T_hot_in - T_cold_out + T_hot_out - T_cold_in) / 2

    config_result = suggest_configuration(T_hot_in=T_hot_in, T_hot_out=T_hot_out, T_cold_in=T_cold_in, T_cold_out=T_cold_out)

    # ── 总传热系数 U：用户给定值优先 ──
    # K 是 U 的别名（中国工程习惯用 K 表示总传热系数）
    _user_U = U if U is not None else K
    if _user_U is not None:
        U_est = float(_user_U)
        conversion_notes.append(f"[用户给定] 总传热系数 U = {U_est} W/(m²·K)")
        _user_provided_props.append("U")
    else:
        u_est_result = estimate_U(f"{hot_fluid}/{cold_fluid}")
        U_est = u_est_result.get("U_typical") if isinstance(u_est_result, dict) else None
        if U_est is None:
            conversion_notes.append("[参数缺失] 无法从数据库获取 U 估算值，使用简化传热模型")

    Q = max(Q_h, Q_c)
    d_i_calc = max(d_o * (1 - 2 * 0.083), d_o * 0.75)

    # 管排列方式自动选择
    if tube_layout == "auto":
        # 默认三角排列（换热效率高）；若冷侧为易结垢流体，选正方排列（便于机械清洗）
        _fouling_fluids = {"slurry", "泥浆", "crude", "原油", "wastewater", "废水", "cooling_water"}
        if cold_fluid.strip().lower() in _fouling_fluids or hot_fluid.strip().lower() in _fouling_fluids:
            tube_layout = "square"
            conversion_notes.append("[管排列] 检测到易结垢流体，自动选择正方排列 (square) 以便机械清洗")
        else:
            tube_layout = "triangular"
            conversion_notes.append("[管排列] 自动选择三角排列 (triangular)，换热效率更高")

    # ============================================================
    # 约束驱动迭代设计 (Constraint-Driven Iterative Design)
    # ============================================================
    # 调整策略:
    #   1. 管程流速偏低 → 增加管程数 (2→4→6→8)
    #   2. 壳程流速偏低 → 减小折流板间距
    #   3. 长径比 L/D 过大 → 增加管程数 (减少管数→增大壳径)
    #   4. 折流板间距比超限 → 调整折流板间距
    # ============================================================
    MAX_ITERATIONS = 10
    tube_passes_list = [2, 4, 6, 8]
    baffle_spacing_options = [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60]
    
    best_result = None
    best_pass_count = 0
    adjust_log = []
    
    # 初始参数
    current_tube_passes = 2
    current_baffle_spacing = None  # None 表示使用默认值
    current_tube_length = L  # 当前管长
    current_force_passes = False
    
    # 标准管长选项 (m)
    tube_length_options = [0.5, 0.8, 1.0, 1.2, 1.8, 1.5, 2.2, 2.5, 2.8, 3.0, 3.2, 3.5, 3.8, 4.0, 4.2, 4.5, 4.8, 5.0, 5.5, 6.0, 6.8, 7.0, 8.0, 9.0]
    
    for iteration in range(MAX_ITERATIONS):
        # 调用设计函数（使用当前参数）
        try:
            if current_baffle_spacing is not None:
                result = design_heat_exchanger(
                    T_hot_in=T_hot_in, T_hot_out=T_hot_out, T_cold_in=T_cold_in, T_cold_out=T_cold_out,
                    m_hot=m_hot, m_cold=m_cold, Cp_hot=float(Cp_h), Cp_cold=float(Cp_c),
                    rho_hot=float(rho_h), rho_cold=float(rho_c), mu_hot=float(mu_h), mu_cold=float(mu_c),
                    k_hot=float(k_h), k_cold=float(k_c), d_o=d_o, d_i=d_i_calc, tube_length=current_tube_length,
                    tube_layout=tube_layout, tube_passes=current_tube_passes,
                    baffle_spacing=current_baffle_spacing,
                    force_tube_passes=current_force_passes,
                    force_baffle_spacing=current_baffle_spacing,
                )
            else:
                result = design_heat_exchanger(
                    T_hot_in=T_hot_in, T_hot_out=T_hot_out, T_cold_in=T_cold_in, T_cold_out=T_cold_out,
                    m_hot=m_hot, m_cold=m_cold, Cp_hot=float(Cp_h), Cp_cold=float(Cp_c),
                    rho_hot=float(rho_h), rho_cold=float(rho_c), mu_hot=float(mu_h), mu_cold=float(mu_c),
                    k_hot=float(k_h), k_cold=float(k_c), d_o=d_o, d_i=d_i_calc, tube_length=current_tube_length,
                    tube_layout=tube_layout, tube_passes=current_tube_passes,
                    force_tube_passes=current_force_passes,
                )
        except Exception as e:
            A_req = Q / (U_est * max(float(delta_T_m), 0.1))
            result = {
                "Q": Q, "LMTD": float(delta_T_m), "U": U_est, "A_required": A_req, "A_design": A_req * 1.15,
                "N_tubes": None, "D_shell": None, "valid": False, "warnings": [f"计算失败({e})"],
            }
            all_warnings.append(f"设计计算异常：{e}，已回退至简化估算")
            break
        
        # --- 检查所有约束 ---
        h_tube_detail_iter = result.get("h_tube_detail", {})
        h_shell_detail_iter = result.get("h_shell_detail", {})
        u_tube_iter = h_tube_detail_iter.get("velocity", 0)
        u_shell_iter = h_shell_detail_iter.get("velocity", 0)
        D_shell_iter = result.get("D_shell", 0)
        # 使用实际请求的折流板间距（而非引擎返回值）进行迭代比较
        baffle_sp_iter = current_baffle_spacing if current_baffle_spacing is not None else result.get("baffle_spacing_m", 0)
        L_tube_iter = result.get("L_tube", current_tube_length)
        L_D_iter = L_tube_iter / D_shell_iter if D_shell_iter > 0 else 99
        
        pass_count = 0
        total_checks = 5
        fail_reasons = []
        
        # 1. 管程流速
        if TUBE_VELOCITY_CONSTRAINTS["liquid_min"] <= u_tube_iter <= TUBE_VELOCITY_CONSTRAINTS["liquid_max"]:
            pass_count += 1
        else:
            fail_reasons.append(f"管程流速{u_tube_iter:.3f}m/s")
        
        # 2. 壳程流速
        if SHELL_VELOCITY_CONSTRAINTS["liquid_min"] <= u_shell_iter <= SHELL_VELOCITY_CONSTRAINTS["liquid_max"]:
            pass_count += 1
        else:
            fail_reasons.append(f"壳程流速{u_shell_iter:.3f}m/s")
        
        # 3. 长径比
        if 3 <= L_D_iter <= 15:
            pass_count += 1
        else:
            fail_reasons.append(f"L/D={L_D_iter:.1f}")
        
        # 4. 折流板间距比
        if D_shell_iter > 0 and baffle_sp_iter > 0:
            B_D = baffle_sp_iter / D_shell_iter
            if BAFFLE_CONSTRAINTS["B_D_min"] <= B_D <= BAFFLE_CONSTRAINTS["B_D_max"]:
                pass_count += 1
            else:
                fail_reasons.append(f"B/D={B_D:.3f}")
        else:
            pass_count += 1  # 无法计算时跳过
        
        # 5. 管束振动 (壳程流速 < 5 m/s)
        if u_shell_iter < 5.0:
            pass_count += 1
        else:
            fail_reasons.append(f"振动风险(u_shell={u_shell_iter:.2f})")
        
        # 记录最优结果
        if pass_count > best_pass_count:
            best_pass_count = pass_count
            best_result = result
        
        # 全部通过 → 退出
        if pass_count == total_checks:
            adjust_log.append(f"[迭代{iteration+1}] 所有约束通过 ✓ ({pass_count}/{total_checks})")
            break
        
        # --- 参数调整策略（允许多参数同时调整）---
        adjusted = False
        
        # 策略1: 管程流速偏低 → 增加管程数
        if u_tube_iter < TUBE_VELOCITY_CONSTRAINTS["liquid_min"]:
            for new_passes in tube_passes_list:
                if new_passes > current_tube_passes:
                    adjust_log.append(f"[迭代{iteration+1}] 管程流速{u_tube_iter:.3f}<0.5, 管程数 {current_tube_passes}→{new_passes}")
                    current_tube_passes = new_passes
                    current_force_passes = True
                    adjusted = True
                    break
        
        # 策略2: 壳程流速偏低 → 减小折流板间距
        if u_shell_iter < SHELL_VELOCITY_CONSTRAINTS["liquid_min"]:
            found_smaller = False
            for new_baffle in baffle_spacing_options:
                if new_baffle < baffle_sp_iter and new_baffle >= 0.10:
                    adjust_log.append(f"[迭代{iteration+1}] 壳程流速{u_shell_iter:.3f}<0.3, 折流板间距 {baffle_sp_iter:.2f}→{new_baffle:.2f}")
                    current_baffle_spacing = new_baffle
                    adjusted = True
                    found_smaller = True
                    break
            if not found_smaller:
                adjust_log.append(f"[迭代{iteration+1}] 壳程流速{u_shell_iter:.3f}<0.3, 折流板间距已为最小值{baffle_sp_iter:.2f}m，无法继续调整")
        
        # 策略3: 长径比过大 → 增加管程数 或 减小管长
        if L_D_iter > 15:
            # 先尝试增加管程数
            for new_passes in tube_passes_list:
                if new_passes > current_tube_passes:
                    adjust_log.append(f"[迭代{iteration+1}] L/D={L_D_iter:.1f}>15, 管程数 {current_tube_passes}→{new_passes}")
                    current_tube_passes = new_passes
                    current_force_passes = True
                    adjusted = True
                    break
            else:
                # 如果管程数已最大，尝试减小管长
                for new_length in tube_length_options:
                    if new_length < current_tube_length:
                        adjust_log.append(f"[迭代{iteration+1}] L/D={L_D_iter:.1f}>15, 管长 {current_tube_length:.1f}→{new_length:.1f}m")
                        current_tube_length = new_length
                        adjusted = True
                        break
        
        # 策略4: 折流板间距比超限 → 调整折流板间距
        if D_shell_iter > 0 and baffle_sp_iter > 0:
            B_D = baffle_sp_iter / D_shell_iter
            if B_D < BAFFLE_CONSTRAINTS["B_D_min"]:
                # 间距过小 → 增大
                for new_baffle in baffle_spacing_options:
                    if new_baffle / D_shell_iter >= BAFFLE_CONSTRAINTS["B_D_min"]:
                        adjust_log.append(f"[迭代{iteration+1}] B/D={B_D:.3f}<0.2, 折流板间距→{new_baffle:.2f}")
                        current_baffle_spacing = new_baffle
                        adjusted = True
                        break
            elif B_D > BAFFLE_CONSTRAINTS["B_D_max"]:
                # 间距过大 → 减小
                for new_baffle in reversed(baffle_spacing_options):
                    if new_baffle / D_shell_iter <= BAFFLE_CONSTRAINTS["B_D_max"]:
                        adjust_log.append(f"[迭代{iteration+1}] B/D={B_D:.3f}>1.0, 折流板间距→{new_baffle:.2f}")
                        current_baffle_spacing = new_baffle
                        adjusted = True
                        break
        
        if not adjusted:
            adjust_log.append(f"[迭代{iteration+1}] 无法进一步优化，当前 {pass_count}/{total_checks}")
            break
    else:
        adjust_log.append(f"[迭代{MAX_ITERATIONS}] 达到最大迭代次数，使用最优结果 ({best_pass_count}/{total_checks})")
        if best_result and best_pass_count < total_checks:
            result = best_result
    
    if adjust_log:
        all_warnings.extend(adjust_log)

    all_warnings.extend(result.get("warnings", []))

    # 提取物理引擎子结果用于暴露中间参数
    h_tube_detail = result.get("h_tube_detail", {})
    h_shell_detail = result.get("h_shell_detail", {})
    lmtd_detail = result.get("lmtd_result", {})
    u_detail = result.get("u_result", {})

    # ── ε-NTU 校核 ──
    from physics_engine.heatexchanger import effectiveness_ntu
    C_hot = m_hot * Cp_h   # W/K
    C_cold = m_cold * Cp_c  # W/K
    U_val = float(result.get("U", U_est)) if result.get("U") else (U_est or 0)
    A_val = float(result.get("A_design", 0) or 0)
    UA = U_val * A_val
    ntu_result = effectiveness_ntu(C_hot=C_hot, C_cold=C_cold, UA=UA, T_hot_in=T_hot_in, T_cold_in=T_cold_in)

    # 管内侧直径
    d_i_val = result.get("d_i", d_i_calc)
    # 管程数
    tube_passes = result.get("tube_passes", 2)
    # 管间距
    pitch = 1.25 * d_o

    ret = {
        "design_summary": {
            "Q_kW":             round(result.get("Q", Q) / 1000, 1),
            "LMTD_K":           round(float(result.get("LMTD", delta_T_m)), 1),
            "LMTD_raw_K":       round(float(lmtd_detail.get("LMTD_raw", delta_T_m)), 1),
            "F_factor":         result.get("F"),
            "F_quality":        result.get("F_quality"),
            # 用户给定 U 时优先使用用户值，否则使用计算值
            "U_W_m2K":          round(float(U_est if "U" in _user_provided_props else result.get("U", U_est)), 0),
            "A_required_m2":    round(float(Q / (U_est * max(float(result.get("LMTD", delta_T_m)), 0.1))), 2) if "U" in _user_provided_props else round(float(result.get("A_required", Q / (U_est * max(float(delta_T_m), 0.1)))), 2),
            "A_design_m2":      round(float(result.get("A_design", 0) or 0), 2) if "U" not in _user_provided_props else round(float(Q / (U_est * max(float(result.get("LMTD", delta_T_m)), 0.1))) * 1.15, 2),
            "N_tubes":          result.get("N_tubes"),
            "N_tube_passes":    tube_passes,
            "D_shell_m":        result.get("D_shell"),
            "L_tube_m":         L,
            "d_o_m":            result.get("d_o", d_o),
            "d_i_m":            round(d_i_val, 5),
            "pitch_m":          round(pitch, 5),
            "overspec_percent": result.get("overspec_percent"),
            "L_D_ratio":        result.get("length_diameter_ratio"),
            "approach_temp_min_K": lmtd_detail.get("approach_temp_min"),
            "recommended_config": config_result.get("recommended_config"),
            "tube_layout":      result.get("tube_layout", tube_layout),
            "N_baffles":        result.get("N_baffles"),
            "baffle_spacing_m": result.get("baffle_spacing_m"),
            "baffle_cut":       result.get("baffle_cut"),
            "fouling_fraction": result.get("fouling_fraction"),
        },
        "heat_balance": {
            "Q_hot_kW":  round(Q_h / 1000, 1),
            "Q_cold_kW": round(Q_c / 1000, 1),
            "imbalance_pct": round(abs(Q_h - Q_c) / max(Q_h, Q_c) * 100 if max(Q_h, Q_c) > 0 else 0, 1),
            "C_hot_W_K": round(C_hot, 1),
            "C_cold_W_K": round(C_cold, 1),
            "C_min_W_K": round(ntu_result.get("C_min", min(C_hot, C_cold)), 1),
        },
        "epsilon_ntu": {
            "epsilon": round(ntu_result.get("epsilon", 0), 4),
            "NTU": round(ntu_result.get("NTU", 0), 3),
            "C_r": round(ntu_result.get("C_r", 0), 4),
            "Q_max_kW": round(ntu_result.get("Q_max", 0) / 1000, 1),
            "Q_actual_kW": round(ntu_result.get("Q_actual", 0) / 1000, 1),
        },
        "U_decomposition": {
            "U_overall": round(u_detail.get("U", U_val), 1),
            "R_total": round(u_detail.get("R_total", 0), 6),
            "R_conv_tube": round(u_detail.get("R_tube", 0), 6),
            "R_conv_shell": round(u_detail.get("R_shell", 0), 6),
            "R_wall": round(u_detail.get("R_wall", 0), 6),
            "R_fouling_tube": round(u_detail.get("R_fouling_tube", 0), 6),
            "R_fouling_shell": round(u_detail.get("R_fouling_shell", 0), 6),
            "fouling_fraction_pct": round(u_detail.get("fouling_fraction", 0), 1),
        },
        "pressure_drop": {
            "tube_Pa":  result.get("delta_P_tube"),
            "shell_Pa": result.get("delta_P_shell"),
        },
        "heat_transfer": {
            "h_tube_W_m2K":  result.get("h_tube"),
            "h_shell_W_m2K": result.get("h_shell"),
            "u_tube_m_s":    h_tube_detail.get("velocity"),
            "u_shell_m_s":   h_shell_detail.get("velocity"),
            "Re_tube":       h_tube_detail.get("Re"),
            "Re_shell":      h_shell_detail.get("Re"),
            "Pr_tube":       h_tube_detail.get("Pr"),
            "Pr_shell":      h_shell_detail.get("Pr"),
            "Nu_tube":       h_tube_detail.get("Nu"),
            "Nu_shell":      h_shell_detail.get("Nu"),
            "d_e_shell_m":   result.get("d_e_shell"),
            "tube_flow_regime":  h_tube_detail.get("flow_regime"),
            "shell_flow_regime": h_shell_detail.get("flow_regime"),
        },
        "properties": {
            "Cp_hot_J_kgK": round(float(Cp_h), 0), "Cp_cold_J_kgK": round(float(Cp_c), 0),
            "rho_hot_kg_m3": round(float(rho_h), 1), "rho_cold_kg_m3": round(float(rho_c), 1),
            "mu_hot_Pas": round(float(mu_h), 6), "mu_cold_Pas": round(float(mu_c), 6),
            "k_hot_W_mK": round(float(k_h), 4), "k_cold_W_mK": round(float(k_c), 4),
            "MW_hot": round(MW_h, 2) if MW_h else None, "MW_cold": round(MW_c, 2) if MW_c else None,
        },
        "input_after_conversion": {
            "T_hot_in_K": T_hot_in, "T_hot_out_K": T_hot_out, "T_cold_in_K": T_cold_in, "T_cold_out_K": T_cold_out,
            "P_hot_Pa": P_hot, "P_cold_Pa": P_cold, "m_hot_kgs": m_hot, "m_cold_kgs": m_cold,
            "hot_fluid": hot_fluid, "cold_fluid": cold_fluid,
        },
        "auto_coolant": auto_coolant_info if auto_coolant_info else None,
        "conversion_notes": conversion_notes,
        "validation":  result.get("constraint_checks", {}),
        "warnings":    all_warnings,
        "converged":   result.get("valid", True),
    }

    # ============================================================
    # 行业标准约束校核 (GB/T 151, SH/T 3121, TEMA)
    # ============================================================
    industry_checks = {}
    
    # --- 1. 管程流速校核 ---
    u_tube = h_tube_detail.get("velocity", 0)
    tube_constraints = TUBE_VELOCITY_CONSTRAINTS
    tube_vel_check = {
        "value_m_s": round(u_tube, 3),
        "min_m_s": tube_constraints["liquid_min"],
        "max_m_s": tube_constraints["liquid_max"],
        "optimal_m_s": tube_constraints["liquid_optimal"],
        "pass": tube_constraints["liquid_min"] <= u_tube <= tube_constraints["liquid_max"],
        "status": "ok" if tube_constraints["liquid_min"] <= u_tube <= tube_constraints["liquid_max"] else "fail",
    }
    if u_tube < tube_constraints["liquid_min"]:
        tube_vel_check["warning"] = f"管程流速 {u_tube:.3f} m/s < {tube_constraints['liquid_min']} m/s，存在沉积风险"
    elif u_tube > tube_constraints["liquid_max"]:
        tube_vel_check["warning"] = f"管程流速 {u_tube:.3f} m/s > {tube_constraints['liquid_max']} m/s，存在冲蚀风险"
    industry_checks["tube_velocity"] = tube_vel_check
    
    # --- 2. 壳程流速校核 ---
    u_shell = h_shell_detail.get("velocity", 0)
    shell_constraints = SHELL_VELOCITY_CONSTRAINTS
    shell_vel_check = {
        "value_m_s": round(u_shell, 3),
        "min_m_s": shell_constraints["liquid_min"],
        "max_m_s": shell_constraints["liquid_max"],
        "optimal_m_s": shell_constraints["liquid_optimal"],
        "pass": shell_constraints["liquid_min"] <= u_shell <= shell_constraints["liquid_max"],
        "status": "ok" if shell_constraints["liquid_min"] <= u_shell <= shell_constraints["liquid_max"] else "fail",
    }
    if u_shell < shell_constraints["liquid_min"]:
        shell_vel_check["warning"] = f"壳程流速 {u_shell:.3f} m/s < {shell_constraints['liquid_min']} m/s，传热效率低"
    elif u_shell > shell_constraints["liquid_max"]:
        shell_vel_check["warning"] = f"壳程流速 {u_shell:.3f} m/s > {shell_constraints['liquid_max']} m/s，压降过大"
    industry_checks["shell_velocity"] = shell_vel_check
    
    # --- 3. 折流板间距比校核 (B/D_shell) ---
    baffle_spacing = result.get("baffle_spacing", 0)
    D_shell = result.get("D_shell", 0)
    if D_shell > 0 and baffle_spacing > 0:
        B_D_ratio = baffle_spacing / D_shell
        baffle_constraints = BAFFLE_CONSTRAINTS
        baffle_check = {
            "B_m": round(baffle_spacing, 4),
            "D_shell_m": round(D_shell, 4),
            "B_D_ratio": round(B_D_ratio, 3),
            "min": baffle_constraints["B_D_min"],
            "max": baffle_constraints["B_D_max"],
            "optimal_range": f"{baffle_constraints['B_D_optimal_min']}-{baffle_constraints['B_D_optimal_max']}",
            "pass": baffle_constraints["B_D_min"] <= B_D_ratio <= baffle_constraints["B_D_max"],
            "status": "ok" if baffle_constraints["B_D_min"] <= B_D_ratio <= baffle_constraints["B_D_max"] else "fail",
        }
        if B_D_ratio < baffle_constraints["B_D_min"]:
            baffle_check["warning"] = f"B/D={B_D_ratio:.3f} < {baffle_constraints['B_D_min']}，间距过小，压降异常"
        elif B_D_ratio > baffle_constraints["B_D_max"]:
            baffle_check["warning"] = f"B/D={B_D_ratio:.3f} > {baffle_constraints['B_D_max']}，间距过大，存在流体死区"
        industry_checks["baffle_spacing_ratio"] = baffle_check
    
    # --- 4. 管束振动风险校核 ---
    # 简化判据：壳程流速与管束固有频率的关系
    # 当壳程流速过高时，可能诱发管束振动
    vibration_check = {
        "u_shell_m_s": round(u_shell, 3),
        "risk_level": "low",
        "pass": True,
    }
    if u_shell > 5.0:
        vibration_check["risk_level"] = "medium"
        vibration_check["warning"] = f"壳程流速 {u_shell:.2f} m/s > 5 m/s，存在中等振动风险"
    if u_shell > 10.0:
        vibration_check["risk_level"] = "high"
        vibration_check["pass"] = False
        vibration_check["warning"] = f"壳程流速 {u_shell:.2f} m/s > 10 m/s，振动风险高，建议增加折流板或降低流速"
    industry_checks["tube_bundle_vibration"] = vibration_check
    
    # --- 5. 热膨胀差校核 ---
    # ΔL = α × L × ΔT
    # ΔT = T_hot_avg - T_cold_avg (管壳温差)
    T_hot_avg = (T_hot_in + T_hot_out) / 2
    T_cold_avg = (T_cold_in + T_cold_out) / 2
    delta_T = abs(T_hot_avg - T_cold_avg)
    L_tube = result.get("L_tube", L or DEFAULT_TUBE_LENGTH)
    
    # 默认碳钢
    alpha = THERMAL_EXPANSION_CONSTRAINTS["alpha_steel"]
    delta_L = alpha * L_tube * delta_T * 1000  # mm
    
    expansion_joint_max = THERMAL_EXPANSION_CONSTRAINTS["expansion_joint_max_mm"]
    thermal_check = {
        "T_hot_avg_K": round(T_hot_avg, 1),
        "T_cold_avg_K": round(T_cold_avg, 1),
        "delta_T_K": round(delta_T, 1),
        "alpha_1_K": alpha,
        "L_tube_m": round(L_tube, 2),
        "delta_L_mm": round(delta_L, 2),
        "expansion_joint_max_mm": expansion_joint_max,
        "pass": delta_L <= expansion_joint_max,
        "status": "ok" if delta_L <= expansion_joint_max else "warning",
    }
    if delta_L > expansion_joint_max:
        thermal_check["warning"] = f"热膨胀差 ΔL={delta_L:.2f} mm > {expansion_joint_max} mm，需要膨胀节"
    industry_checks["thermal_expansion"] = thermal_check
    
    # 添加到输出
    ret["industry_standard_checks"] = industry_checks
    
    # 添加警告
    for check_name, check_data in industry_checks.items():
        if not check_data.get("pass", True):
            warning_msg = check_data.get("warning", f"{check_name} 校核失败")
            all_warnings.append(f"[行业标准] {warning_msg}")

    # ── 默认值参数提醒报告 ──
    if _user_defaulted_params:
        report_lines = [
            "═" * 60,
            "ℹ️ 以下换热器几何参数使用了工程标准值（用户未显式提供）：",
            "═" * 60,
        ]
        for pname, pval in _user_defaulted_params.items():
            report_lines.append(f"  ▸ {pname} = {pval}")
            if pname == "d_o":
                report_lines.append(f"    换热管外径（工业常用 15-38 mm，¾\" OD 最常用）")
            elif pname == "L":
                report_lines.append(f"    换热管长（工业常用 2-6 m，3m 最常用）")
            report_lines.append(f"    💡 建议根据实际工况提供此参数以提高计算精度")
        report_lines.append("═" * 60)
        report_lines.append("💡 以上参数均为可动态调整的标准规格，可根据实际工况修改后重新计算")
        ret["user_defaulted_params_report"] = "\n".join(report_lines)
        ret["user_defaulted_params"] = _user_defaulted_params

    return ret

class _HexDesignHandler:
    def execute(self, params: dict) -> dict:
        try:
            return _full_he_design(**params)
        except TypeError as e:
            return {"success": False, "error": f"参数错误：{e}。请检查参数名是否与工具定义一致。"}
        except Exception as e:
            return {"success": False, "error": str(e)}

_HEX_HANDLER = _HexDesignHandler()

def register_heatexchanger_device(registry: ToolRegistry) -> ToolRegistry:
    registry.register_function(
        name="shell_and_tube_hex_design",
        description="管壳式换热器完整设计。自动处理单位换算（°C↔K、bar↔Pa、kg/h↔kg/s）。",
        func=_HEX_HANDLER,
        param_descriptions={
            "T_hot_in":  "(必填) 热流体入口温度（K 或 °C），无默认值",
            "T_hot_out": "(必填) 热流体出口温度（K 或 °C），无默认值",
            "T_cold_in": "(必填) 冷流体入口温度（K 或 °C），无默认值",
            "m_hot":     "(必填) 热流体质量流量（kg/s 或 kg/h），无默认值",
            "T_cold_out":"(选填) 冷流体出口温度，不填自动从热平衡估算",
            "m_cold":    "(选填) 冷流体质量流量，不填自动从热平衡计算",
            "hot_fluid": "(必填) 热流体名称，如 'water','methanol'，无默认值",
            "cold_fluid":"(必填) 冷流体名称，如 'water'(冷却水)，无默认值",
            "P_hot":     "(必填) 热侧压力 (Pa/bar/kPa)，无默认值",
            "P_cold":    "(必填) 冷侧压力 (Pa/bar/kPa)，无默认值",
            "d_o":       "(选填) 换热管外径(m)，标准可选值: 0.015/0.019(默认)/0.025/0.032/0.038 m，不填自动用默认 19mm",
            "L":         "(选填) 换热管长(m)，标准可选值: 2.0/3.0(默认)/4.0/6.0 m，不填自动用默认 3.0m",
            "tube_layout": "(选填) 管排列方式: 'triangular'(三角)/'square'(正方)/'auto'(自动，默认三角)，默认 'auto'",
        },
        category="device", tags=["heat_exchanger", "shell_and_tube", "LMTD"]
    )
    return registry

def register_to_universal(registry: UniversalToolRegistry) -> None:
    capability = ToolCapability(
        name="shell_and_tube_hex_design",
        description="管壳式换热器完整设计。自动处理单位换算，计算热负荷、LMTD、传热面积等。支持用户给定物性参数（U/Cp/rho/mu/k），给定后跳过内部查询。",
        category="device_design",
        parameters={
            "T_hot_in": {"type": "number", "required": True, "description": "(必填) 热流体入口温度 K/°C，无默认值"},
            "T_hot_out": {"type": "number", "required": True, "description": "(必填) 热流体出口温度 K/°C，无默认值"},
            "T_cold_in": {"type": "number", "required": True, "description": "(必填) 冷流体入口温度 K/°C，无默认值"},
            "m_hot": {"type": "number", "required": True, "description": "(必填) 热流体质量流量 kg/s 或 kg/h，无默认值"},
            "T_cold_out": {"type": "number", "required": False, "description": "(选填) 冷流体出口温度，不填自动从热平衡估算"},
            "m_cold": {"type": "number", "required": False, "description": "(选填) 冷流体质量流量，不填自动从热平衡计算"},
            "hot_fluid": {"type": "string", "required": True, "description": "(必填) 热流体名称，如 'water','methanol'，无默认值"},
            "cold_fluid": {"type": "string", "required": True, "description": "(必填) 冷流体名称，如 'water'(冷却水)，无默认值"},
            "P_hot": {"type": "number", "required": True, "description": "(必填) 热侧压力 Pa/bar/kPa，无默认值"},
            "P_cold": {"type": "number", "required": True, "description": "(必填) 冷侧压力 Pa/bar/kPa，无默认值"},
            "d_o": {"type": "number", "required": False, "description": "(选填) 管外径 m，标准可选值: 0.015/0.019(默认)/0.025/0.032/0.038，不填自动用 19mm"},
            "L": {"type": "number", "required": False, "description": "(选填) 管长 m，标准可选值: 2.0/3.0(默认)/4.0/6.0，不填自动用 3.0m"},
            "tube_layout": {"type": "string", "required": False, "description": "(选填) 管排列: 'triangular'/'square'/'auto'，默认 'auto'"},
            # ── 可选物性参数（用户给定后跳过内部查询） ──
            "U": {"type": "number", "required": False, "description": "(选填) 总传热系数 W/(m²·K)，给定后跳过内部估算。K 的别名"},
            "K": {"type": "number", "required": False, "description": "(选填) 总传热系数 W/(m²·K)，中国工程习惯用 K 表示。U 的别名"},
            "Cp_hot": {"type": "number", "required": False, "description": "(选填) 热侧比热容 J/(kg·K) 或 kJ/(kg·K)（<100 自动按 kJ 换算），给定后跳过物性查询"},
            "Cp_cold": {"type": "number", "required": False, "description": "(选填) 冷侧比热容 J/(kg·K) 或 kJ/(kg·K)，给定后跳过物性查询"},
            "rho_hot": {"type": "number", "required": False, "description": "(选填) 热侧密度 kg/m³，给定后跳过物性查询"},
            "rho_cold": {"type": "number", "required": False, "description": "(选填) 冷侧密度 kg/m³，给定后跳过物性查询"},
            "mu_hot": {"type": "number", "required": False, "description": "(选填) 热侧粘度 Pa·s，给定后跳过物性查询"},
            "mu_cold": {"type": "number", "required": False, "description": "(选填) 冷侧粘度 Pa·s，给定后跳过物性查询"},
            "k_hot": {"type": "number", "required": False, "description": "(选填) 热侧导热系数 W/(m·K)，给定后跳过物性查询"},
            "k_cold": {"type": "number", "required": False, "description": "(选填) 冷侧导热系数 W/(m·K)，给定后跳过物性查询"},
        },
        examples=[
            {
                "description": "甲醇冷却器",
                "input": {"T_hot_in": 60, "T_hot_out": 35, "T_cold_in": 25, "m_hot": 1.0, "hot_fluid": "methanol"}
            },
            {
                "description": "用户给定总传热系数和比热容",
                "input": {"T_hot_in": 150, "T_hot_out": 90, "T_cold_in": 30, "T_cold_out": 80, "m_hot": 5000, "K": 500, "Cp_hot": 2.5, "hot_fluid": "oil", "cold_fluid": "water", "P_hot": 101325, "P_cold": 101325}
            }
        ],
        tags=["heat_exchanger", "换热器", "单位换算", "用户给定物性"],
    )
    registry.register(capability, _full_he_design)