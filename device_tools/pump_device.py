"""泵选型完整设计 MCP 工具"""
import math
from physics_engine.pump import pump_selection_full, validate_pump_selection, suggest_pump_type
from physics_engine.common import _REQUIRED, check_required_params
from tools.base import ToolRegistry
from tools.registry import UniversalToolRegistry, ToolCapability

# ============================================================
# 泵标准参数集（工程常用值）
# ============================================================
# 电机转速 (rpm) — 常用规格
STANDARD_RPM = [2900, 1450, 960, 730]
DEFAULT_RPM = 2900  # 4极电机，最常用
_RPM_NOTE = {
    2900: "4极电机(50Hz)，最常用",
    1450: "8极电机(50Hz)，中大流量",
    960: "12极电机(50Hz)，大流量低转速",
    730: "16极电机(50Hz)，超大流量",
}

# ============================================================
# 比转速约束（硬编码）
# ============================================================
# 中国工程单位制: ns_CN = n(rpm) × √Q(m³/h) / H(m)^0.75
# 约束范围: ns_CN ∈ [30, 300]，最优范围 [80, 200]
# 
# 换算关系:
#   ns_CN = 3600 × ns_SI
#   ns_US = 2730 × ns_SI
#   ns_CN ≈ 1.32 × ns_US
#
# SI无量纲比转速约束 (rps, m³/s, m)
# 公式: ns = n(rps) × √Q(m³/s) / H(m)^0.75
PUMP_CONSTRAINTS = {
    "ns_min": 0.0083,      # 比转速下限 (SI) ≈ 中国标准 30
    "ns_max": 0.083,       # 比转速上限 (SI) ≈ 中国标准 300
    "ns_optimal_min": 0.022,  # 最优范围下限 (SI) ≈ 中国标准 80
    "ns_optimal_max": 0.056,  # 最优范围上限 (SI) ≈ 中国标准 200
}

# 中国工程单位制约束 (rpm, m³/h, m) - 仅供工程参考
PUMP_CONSTRAINTS_CN = {
    "ns_CN_min": 30,       # 离心泵最低比转速
    "ns_CN_max": 300,      # 离心泵最高比转速
    "ns_CN_optimal_min": 80,   # 最优范围下限
    "ns_CN_optimal_max": 200,  # 最优范围上限
}

# ============================================================
# 比转速单位转换关系
# ============================================================
# SI (无量纲): ns = n(rps) × √Q(m³/s) / H(m)^0.75
# 中国:        ns_CN = n(rpm) × √Q(m³/h) / H(m)^0.75
# US (传统):   Ns = N(rpm) × √Q(gpm) / H(ft)^0.75
# 
# 换算关系:
#   ns_CN = 3600 × ns_SI  (精确: 60 × √3600 = 3600)
#   ns_US = 2730 × ns_SI  (精确: 60 × √15850.3 / 3.28084^0.75 ≈ 2730)
#   ns_CN = 1.32 × ns_US
#
# 常用参考值 (三种单位制对照):
#   ns_SI=0.0083 → ns_CN=30   → ns_US=23   (离心泵下限)
#   ns_SI=0.022  → ns_CN=80   → ns_US=60   (最优范围下限)
#   ns_SI=0.056  → ns_CN=200  → ns_US=150  (最优范围上限)
#   ns_SI=0.083  → ns_CN=300  → ns_US=230  (离心泵上限)
NS_TO_NS_CN = 3600.0  # SI → 中国 转换系数
NS_TO_NS_US = 2730.0  # SI → US 转换系数

def _full_design(Q: float = _REQUIRED, H: float = _REQUIRED,
                 fluid: str = _REQUIRED, T: float = _REQUIRED,
                 NPSHa: float = _REQUIRED, pump_type: str = None,
                 P_suction: float = None, suction_head: float = 0.0,
                 n: float = None, d_pipe: float = None,
                 # ── 可选：用户给定物性/设计参数（给定后跳过内部查询） ──
                 rho: float = None,        # 流体密度 kg/m³
                 mu: float = None,         # 流体动力粘度 Pa·s
                 Pv: float = None,         # 饱和蒸气压 Pa
                 efficiency: float = None, # 泵效率 (0-1)
                 **kwargs) -> dict:
    # 检查必填参数（核心工艺参数必须由用户提供，n 有工程默认值）
    missing = check_required_params(locals(), {
        "Q": "设计流量 (m³/h)",
        "H": "设计扬程 (m)",
        "fluid": "流体名称，如 'water', 'methanol'",
        "T": "流体温度 (K)",
        "NPSHa": "有效汽蚀余量 (m)",
    })
    if missing:
        return missing

    # n 默认值处理
    _user_defaulted_params = {}
    conversion_notes = []
    if n is None or n is _REQUIRED:
        n = DEFAULT_RPM
        _user_defaulted_params["n"] = f"{n} rpm ({_RPM_NOTE.get(n, '')})，标准可选: {STANDARD_RPM}"
        conversion_notes.append(f"[默认值] 电机转速 n 使用 {n} rpm ({_RPM_NOTE.get(n, '')})")
    n = int(n)
    Q_m3s = Q / 3600.0  # m³/h → m³/s

    # 调用物理引擎进行完整选型
    pump_kwargs = dict(Q_required=Q_m3s, H_required=H, fluid=fluid, T=T, pump_type=pump_type)
    if P_suction is not None:
        pump_kwargs["P_suction"] = P_suction
    if suction_head != 0.0:
        pump_kwargs["suction_head"] = suction_head
    if n != 2900:
        pump_kwargs["n"] = n
    result = pump_selection_full(**pump_kwargs)
    if not isinstance(result, dict):
        return {"error": "泵选型计算返回非字典结果", "converged": False}

    # ── 补全所有输出参数 ──
    result['Q_m3h'] = Q
    result['Q_m3s'] = round(Q_m3s, 6)
    result['H_m'] = H
    result['NPSHa_user'] = NPSHa
    result['rpm'] = n

    # ============================================================
    # 比转速约束校验与自动重新选型
    # ============================================================
    adjust_notes = []
    n_rps = n / 60.0
    ns = 0.0
    if H > 0:
        ns = n_rps * math.sqrt(Q_m3s) / (H ** 0.75)
    
    # 检查比转速是否在合理范围内
    ns_min = PUMP_CONSTRAINTS["ns_min"]
    ns_max = PUMP_CONSTRAINTS["ns_max"]
    ns_optimal_min = PUMP_CONSTRAINTS["ns_optimal_min"]
    ns_optimal_max = PUMP_CONSTRAINTS["ns_optimal_max"]
    
    if ns < ns_min or ns > ns_max:
        # 比转速超出范围，尝试不同转速重新选型
        adjust_notes.append(f"[调整] 初始比转速 ns={ns:.3f} 超出范围({ns_min}-{ns_max})，尝试调整转速")
        
        # 策略：对于低比转速（高扬程小流量），选择最低标准转速以获得更高效率
        # 对于高比转速（大流量低扬程），选择最高标准转速
        if ns < ns_min:
            # 低比转速：选择最低转速 (730 rpm) 以改善效率
            best_rpm = min(STANDARD_RPM)
            adjust_notes.append(f"[调整] 低比转速工况，选择最低标准转速 {best_rpm} rpm 以改善效率")
        else:
            # 高比转速：选择最高转速 (2900 rpm)
            best_rpm = max(STANDARD_RPM)
            adjust_notes.append(f"[调整] 高比转速工况，选择最高标准转速 {best_rpm} rpm")
        
        if best_rpm != n:
            # 重新选型
            pump_kwargs_new = dict(Q_required=Q_m3s, H_required=H, fluid=fluid, T=T, pump_type=pump_type)
            if P_suction is not None:
                pump_kwargs_new["P_suction"] = P_suction
            if suction_head != 0.0:
                pump_kwargs_new["suction_head"] = suction_head
            pump_kwargs_new["n"] = best_rpm
            
            result_new = pump_selection_full(**pump_kwargs_new)
            if isinstance(result_new, dict):
                result = result_new
                n = best_rpm
                result['rpm'] = n
                
                # 重新计算比转速
                n_rps = n / 60.0
                ns = n_rps * math.sqrt(Q_m3s) / (H ** 0.75) if H > 0 else 0
                adjust_notes.append(f"[调整] 重新选型后比转速 ns={ns:.3f}")
                
                # 检查是否满足约束
                if ns < ns_min or ns > ns_max:
                    adjust_notes.append(f"[警告] 调整后比转速仍超出范围，该工况需特殊选型")
        else:
            adjust_notes.append(f"[信息] 当前转速 {n} rpm 已是边界值，无法进一步优化")

    # 比转速结果记录
    result['specific_speed'] = round(ns, 4)  # SI (无量纲)
    result['specific_speed_cn'] = round(ns * NS_TO_NS_CN, 1)  # 中国单位制
    result['specific_speed_us'] = round(ns * NS_TO_NS_US, 1)  # US单位制
    result['specific_speed_note'] = (
        "低比转速(ns_CN<30)" if ns * NS_TO_NS_CN < 30 else
        "中低比转速(30≤ns_CN<80)" if ns * NS_TO_NS_CN < 80 else
        "中比转速(80≤ns_CN<200)" if ns * NS_TO_NS_CN < 200 else
        "中高比转速(200≤ns_CN<300)" if ns * NS_TO_NS_CN < 300 else
        "高比转速(ns_CN≥300)"
    )
    result['specific_speed_optimal'] = ns_optimal_min <= ns <= ns_optimal_max
    result['specific_speed_unit_note'] = (
        f"SI: ns={ns:.4f} (rps, m³/s, m) | "
        f"中国: ns_CN={ns*NS_TO_NS_CN:.0f} (rpm, m³/h, m) | "
        f"US: Ns={ns*NS_TO_NS_US:.0f} (rpm, gpm, ft) | "
        f"换算: ns_CN = 3600×ns, Ns = 2730×ns"
    )
    if adjust_notes:
        result['adjust_notes'] = adjust_notes

    # 管径估算（若未提供，基于经济流速 1.5 m/s）
    if d_pipe is None:
        u_econ = 1.5  # m/s 经济流速
        d_pipe = math.sqrt(4 * Q_m3s / (math.pi * u_econ))
        result['d_pipe_estimated'] = True
    else:
        result['d_pipe_estimated'] = False
    result['d_pipe_m'] = round(d_pipe, 4)
    A_pipe = math.pi * d_pipe ** 2 / 4
    u_pipe = Q_m3s / A_pipe if A_pipe > 0 else 0
    result['u_pipe_m_s'] = round(u_pipe, 3)

    # 流体物性补全（用户给定值优先）
    fp = result.get('fluid_properties', {})
    if rho is not None:
        fp['rho'] = float(rho)
        conversion_notes.append(f"[用户给定] 流体密度 rho = {rho} kg/m³")
    else:
        fp.setdefault('rho', 1000.0)
    if mu is not None:
        fp['mu'] = float(mu)
        conversion_notes.append(f"[用户给定] 流体粘度 mu = {mu} Pa·s")
    else:
        fp.setdefault('mu', 0.001)
    if Pv is not None:
        fp['Pv'] = float(Pv)
        conversion_notes.append(f"[用户给定] 饱和蒸气压 Pv = {Pv} Pa")
    rho = fp.get('rho', 1000.0)
    mu = fp.get('mu', 0.001)
    result['fluid_properties']['kinematic_viscosity_m2_s'] = round(mu / rho if rho > 0 else 0, 8)

    # 功率详情补全
    power = result.get('power', {})
    if power:
        P_hyd = power.get('P_hydraulic', 0)
        result['power_details'] = {
            'P_hydraulic_kW': round(P_hyd / 1000, 3) if P_hyd else None,
            'P_shaft_kW': round(power.get('P_shaft', 0) / 1000, 3) if power.get('P_shaft') else None,
            'P_motor_kW': round(power.get('P_motor', 0) / 1000, 3) if power.get('P_motor') else None,
            'motor_efficiency': power.get('motor_efficiency'),
            'service_factor': power.get('service_factor'),
        }

    # NPSH 安全余量
    npsh = result.get('npsh', {})
    npshr = npsh.get('npshr', 0)
    result['npsh']['NPSHa_user'] = NPSHa
    result['npsh']['margin_user'] = round(NPSHa - npshr, 3) if npshr else None
    result['npsh']['safe_user'] = (NPSHa - npshr) >= 0.5 if npshr else None

    # 效率（用户给定值优先）
    if efficiency is not None:
        eff = float(efficiency)
        result['efficiency'] = eff
        conversion_notes.append(f"[用户给定] 泵效率 = {eff*100:.1f}%")
    else:
        eff = result.get('efficiency', 0)
    result['efficiency_percent'] = round(eff * 100, 1) if eff else None

    # 约束校验
    validation = validate_pump_selection(result)
    # 追加 NPSHa 用户值校验
    if npshr and NPSHa - npshr < 0.5:
        validation['warnings'].append(
            f"用户提供的 NPSHa={NPSHa:.1f}m 与 NPSHr={npshr:.1f}m 余量不足（<0.5m），"
            f"建议降低吸入高度或增大入口压力。"
        )
    result['validation'] = validation
    result['converged'] = validation.get('valid', True)
    result['warnings'] = validation.get('warnings', []) + conversion_notes + adjust_notes

    # 用户默认参数报告
    if _user_defaulted_params:
        _desc_map = {"n": "电机转速"}
        lines = ["\n---\n⚙️ **泵默认参数提醒**（以下参数用户未提供，已使用工程默认值）：\n"]
        for _k, _v in _user_defaulted_params.items():
            _label = _desc_map.get(_k, _k)
            lines.append(f"- **{_label}** = {_v}")
        result["user_defaulted_params_report"] = "\n".join(lines)

    return result

def register_pump_device(registry: ToolRegistry) -> ToolRegistry:
    registry.register_function(name="pump_selection_design",
        description="泵选型完整设计。封装完整设计流水线: 物性查询→计算→校验→报告。",
        func=_full_design,
        param_descriptions={
            "Q": "(必填) 流量 (m³/h)，无默认值", 
            "H": "(必填) 扬程 (m)，无默认值", 
            "fluid": "(必填) 流体名称，如 'water'，无默认值", 
            "T": "(必填) 温度 (K)，无默认值", 
            "NPSHa": "(必填) 有效汽蚀余量 (m)，无默认值",
            "pump_type": "(选填) 泵类型，默认自动推荐",
            "n": "(选填) 电机转速 rpm，默认 2900(4极)，可选: 2900/1450/960/730"
        },
        category="device", tags=["pump", "centrifugal", "泵选型"])
    return registry

def register_to_universal(registry: UniversalToolRegistry) -> None:
    capability = ToolCapability(
        name="pump_selection_design",
        description="泵选型完整设计。根据流量、扬程等参数选择合适的泵类型和型号，计算功率、效率、NPSH等参数。",
        category="device_design",
        parameters={
            "Q": {"type": "number", "required": True, "description": "(必填) 流量 m³/h，无默认值"},
            "H": {"type": "number", "required": True, "description": "(必填) 扬程 m，无默认值"},
            "fluid": {"type": "string", "required": True, "description": "(必填) 流体名称，如 'water'，无默认值"},
            "T": {"type": "number", "required": True, "description": "(必填) 温度，单位 K，无默认值"},
            "NPSHa": {"type": "number", "required": True, "description": "(必填) 有效汽蚀余量，单位 m，无默认值"},
            "pump_type": {"type": "string", "required": False, "description": "(选填) 泵类型，不填则系统自动匹配"},
            "n": {"type": "number", "required": False, "description": "(选填) 电机转速 rpm，默认 2900(4极电机)，可选: 2900/1450/960/730"},
            # ── 可选物性参数 ──
            "rho": {"type": "number", "required": False, "description": "(选填) 流体密度 kg/m³，给定后跳过物性查询"},
            "mu": {"type": "number", "required": False, "description": "(选填) 流体动力粘度 Pa·s，给定后跳过物性查询"},
            "Pv": {"type": "number", "required": False, "description": "(选填) 饱和蒸气压 Pa，给定后跳过物性查询"},
            "efficiency": {"type": "number", "required": False, "description": "(选填) 泵效率 (0-1)，给定后跳过效率估算"},
        },
        examples=[{"description": "水泵选型", "input": {"Q": 36, "H": 30, "fluid": "water", "T": 298.15}}],
        tags=["pump", "泵", "选型", "离心泵"],
    )
    registry.register(capability, _full_design)