"""process_engine 工具注册模块 — 将 process_engine 中的设备级 run_* 函数注册为独立工具。

每个设备作为一个可串联的独立工具，上游设备的出口参数可直接作为下游设备的入口参数，
实现多设备串联流程编排。

支持的设备（10 种）：
  1. run_pump                    — 泵（升压）
  2. run_pressure_reducing_valve — 减压阀（降压）
  3. run_heat_exchanger          — 换热器（加热/冷却）
  4. run_stoichiometric_reactor  — 化学计量反应器
  5. run_cstr_reactor            — CSTR 全混流反应器
  6. run_pfr_reactor             — PFR 平推流反应器
  7. run_distillation_column     — 板式精馏塔
  8. run_flash_drum              — 闪蒸罐（气液分离）
  9. run_mixer                   — 混合器（多股进料合并）
 10. run_flowsheet               — 流程求解器（带循环收敛的完整流程模拟）
"""

import os
import sys

# ─────────────────────────────────────────────────────────────────
# 路径设置：确保 process_engine 下的子模块可被导入
# ─────────────────────────────────────────────────────────────────
_PROCESS_ENGINE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "process_engine",
)
if _PROCESS_ENGINE_DIR not in sys.path:
    sys.path.insert(0, _PROCESS_ENGINE_DIR)

# 导入各设备级封装函数
from multi_equipment_series.funcs import (
    run_heat_exchanger as _run_heat_exchanger,
    run_stoichiometric_reactor as _run_stoichiometric_reactor,
    run_cstr_reactor as _run_cstr_reactor,
    run_pfr_reactor as _run_pfr_reactor,
    run_distillation_column as _run_distillation_column,
    run_pump as _run_pump,
    run_pressure_reducing_valve as _run_pressure_reducing_valve,
    run_flash_drum as _run_flash_drum,
    run_mixer as _run_mixer,
)

# 导入回流迭代引擎
from flowsheet.engine import FlowsheetEngine as _FlowsheetEngine


# ─────────────────────────────────────────────────────────────────
# 包装函数：human_ask — 人工参数补充请求
# ─────────────────────────────────────────────────────────────────
def human_ask(**kwargs):
    """人工参数请求工具：当核心参数缺失时，暂停流程并请求用户补充。

    参数:
        question (str): 向用户提出的问题
        missing_params (list): 缺失的参数名列表
        context (str): 补充说明上下文
    """
    question = kwargs.get("question", "请补充缺失的参数")
    missing = kwargs.get("missing_params", [])
    context = kwargs.get("context", "")
    prompt_parts = [question]
    if missing:
        prompt_parts.append(f"缺失参数: {', '.join(str(p) for p in missing)}")
    if context:
        prompt_parts.append(f"背景: {context}")
    prompt = "\n".join(prompt_parts)
    return {"needs_human": True, "human_prompt": prompt}


# ─────────────────────────────────────────────────────────────────
# 包装函数：流程求解器 — 带循环收敛的完整流程模拟
# ─────────────────────────────────────────────────────────────────
def run_flowsheet(**kwargs):
    """流程求解器包装函数：求解带循环的化工流程，返回收敛结果 + 后处理。

    参数:
        flowsheet_config (dict): 流程配置字典，包含 feeds/equipment/recycles
        max_iter (int): 最大迭代次数，默认 500
        tol (float): 收敛容差，默认 1e-4
        acceleration (str): 加速方法 "wegstein"/"successive"/"none"，默认 "wegstein"
        omega (float): 阻尼因子，默认 0.5
        do_post_process (bool): 是否执行后处理，默认 True
    """
    config = kwargs.get("flowsheet_config")
    if config is None:
        return {"error": "缺少必填参数 flowsheet_config（流程配置字典）"}

    # 容错：LLM 有时将 config 以 JSON 字符串形式传入
    if isinstance(config, str):
        import json
        try:
            config = json.loads(config)
            kwargs["flowsheet_config"] = config
        except (json.JSONDecodeError, TypeError):
            return {"error": f"flowsheet_config 是字符串但无法解析为 JSON: {config[:200]}"}

    max_iter = int(kwargs.get("max_iter", 500))
    tol = float(kwargs.get("tol", 1e-4))
    acceleration = kwargs.get("acceleration", "wegstein")
    omega = float(kwargs.get("omega", 0.5))
    do_post_process = kwargs.get("do_post_process", True)

    eng = _FlowsheetEngine(config)
    result = eng.solve_full(
        max_iter=max_iter,
        tol=tol,
        acceleration=acceleration,
        omega=omega,
        do_post_process=do_post_process,
    )
    return result


# ─────────────────────────────────────────────────────────────────
# 包装函数：为 CSTR/PFR 的催化剂和动力学参数提供工程默认值
# ─────────────────────────────────────────────────────────────────


def run_cstr_reactor(**kwargs):
    """CSTR 包装函数：自动填充缺失的催化剂/动力学参数为工程默认值。"""
    # CSTR 使用的默认参数（不含 catalyst_bed_void_fraction，那是 PFR 的）
    _cstr_defaults = {
        "catalyst_particle_diameter_m": 0.003,
        "catalyst_particle_density_kg_per_m3": 2500.0,
        "arrhenius_pre_exponential_factor": 1.0,
        "activation_energy_J_per_mol": 50000.0,
        "reaction_order": 1,
    }
    for key, default in _cstr_defaults.items():
        if key not in kwargs or kwargs[key] is None:
            kwargs[key] = default
    result = _run_cstr_reactor(**kwargs)
    defaulted = [k for k in _cstr_defaults
                 if k in kwargs and kwargs[k] == _cstr_defaults[k]
                 and k not in ("reaction_order",)]
    if defaulted:
        result["_defaulted_params"] = defaulted
        result["_note"] = "以下参数使用了工程默认值，建议根据实际工况调整: " + ", ".join(defaulted)
    return result


def run_pfr_reactor(**kwargs):
    """PFR 包装函数：自动填充缺失的催化剂/动力学参数为工程默认值。"""
    # PFR 使用的默认参数（不含 catalyst_particle_diameter_m，那是 CSTR 的）
    _pfr_defaults = {
        "catalyst_particle_density_kg_per_m3": 2500.0,
        "catalyst_bed_void_fraction": 0.4,
        "arrhenius_pre_exponential_factor": 1.0,
        "activation_energy_J_per_mol": 50000.0,
        "reaction_order": 1,
        "tube_inner_diameter_m": 0.05,
        "superficial_velocity_m_per_s": 1.0,
    }
    for key, default in _pfr_defaults.items():
        if key not in kwargs or kwargs[key] is None:
            kwargs[key] = default
    result = _run_pfr_reactor(**kwargs)
    defaulted = [k for k in _pfr_defaults
                 if k in kwargs and kwargs[k] == _pfr_defaults[k]
                 and k not in ("reaction_order",)]
    if defaulted:
        result["_defaulted_params"] = defaulted
        result["_note"] = "以下参数使用了工程默认值，建议根据实际工况调整: " + ", ".join(defaulted)
    return result


# 直接使用原始函数（无需包装的设备）
run_pump = _run_pump


# ─────────────────────────────────────────────────────────────────
# 包装函数：反应器 — 确保流量字典和系数为数值类型
# ─────────────────────────────────────────────────────────────────
def run_stoichiometric_reactor(**kwargs):
    """反应器包装函数：确保所有数值参数类型正确。"""
    for dict_key in ("inlet_molar_flows_mol_per_s", "stoichiometric_coefficients"):
        if dict_key in kwargs:
            kwargs[dict_key] = _ensure_numeric_dict_values(kwargs[dict_key])
    for key in ("key_component_conversion", "reactor_temp_K", "reactor_pressure_Pa"):
        if key in kwargs:
            kwargs[key] = _ensure_numeric(kwargs[key])
    return _run_stoichiometric_reactor(**kwargs)


# ─────────────────────────────────────────────────────────────────
# 包装函数：闪蒸罐 — 确保流量字典和回收率为数值类型
# ─────────────────────────────────────────────────────────────────
def run_flash_drum(**kwargs):
    """闪蒸罐包装函数：确保所有数值参数类型正确。"""
    if "inlet_molar_flows_mol_per_s" in kwargs:
        kwargs["inlet_molar_flows_mol_per_s"] = _ensure_numeric_dict_values(
            kwargs["inlet_molar_flows_mol_per_s"]
        )
    for key in ("flash_pressure_Pa", "key_component_recovery"):
        if key in kwargs:
            kwargs[key] = _ensure_numeric(kwargs[key])
    return _run_flash_drum(**kwargs)


# ─────────────────────────────────────────────────────────────────
# 类型安全辅助：递归将 dict/list 中的字符串数值转为 float
# ─────────────────────────────────────────────────────────────────
def _ensure_numeric_dict_values(d):
    """确保字典中所有值为数值类型，将字符串数值转为 float。"""
    if not isinstance(d, dict):
        return d
    out = {}
    for k, v in d.items():
        if isinstance(v, str):
            try:
                v = float(v)
            except (ValueError, TypeError):
                pass
        out[k] = v
    return out


def _ensure_numeric(value):
    """确保单个值为数值类型。"""
    if isinstance(value, str):
        try:
            return float(value)
        except (ValueError, TypeError):
            return value
    return value


# ─────────────────────────────────────────────────────────────────
# 包装函数：换热器 — 确保流量字典和温度/压力参数为数值类型
# 换热器不改变组分，透传入口流量为出口流量 outlet_molar_flows_mol_per_s
# ─────────────────────────────────────────────────────────────────
def run_heat_exchanger(**kwargs):
    """换热器包装函数：确保所有数值参数类型正确，并透传入口流量。"""
    # 确保流量字典中的值为数值
    inlet_flows = kwargs.get("process_fluid_molar_flows_mol_per_s")
    if inlet_flows is not None:
        kwargs["process_fluid_molar_flows_mol_per_s"] = _ensure_numeric_dict_values(inlet_flows)
    # 确保标量参数为数值
    for key in ("process_fluid_temp_in_K", "process_fluid_temp_out_K",
                "process_fluid_pressure_Pa", "utility_fluid_temp_in_K",
                "utility_fluid_temp_out_K"):
        if key in kwargs:
            kwargs[key] = _ensure_numeric(kwargs[key])
    result = _run_heat_exchanger(**kwargs)
    # 透传：换热器不改变组分，入口流量 = 出口流量
    if inlet_flows is not None:
        result["outlet_molar_flows_mol_per_s"] = _ensure_numeric_dict_values(inlet_flows)
        result["outlet_temp_K"] = _ensure_numeric(kwargs.get("process_fluid_temp_out_K"))
        result["outlet_pressure_Pa"] = _ensure_numeric(kwargs.get("process_fluid_pressure_Pa"))
    return result


# ─────────────────────────────────────────────────────────────────
# 包装函数：减压阀 — 接收并透传入口流量
# 减压阀不改变组分，透传入口流量为出口流量 outlet_molar_flows_mol_per_s
# ─────────────────────────────────────────────────────────────────
def run_pressure_reducing_valve(**kwargs):
    """减压阀包装函数：接收可选的入口流量参数并透传。"""
    # 提取可选的入口流量（原始函数不需要此参数）
    inlet_flows = kwargs.pop("inlet_molar_flows_mol_per_s", None)
    # 确保压力参数为数值
    for key in ("inlet_pressure_Pa", "target_pressure_Pa"):
        if key in kwargs:
            kwargs[key] = _ensure_numeric(kwargs[key])
    result = _run_pressure_reducing_valve(**kwargs)
    # 透传：减压阀不改变组分
    if inlet_flows is not None:
        result["outlet_molar_flows_mol_per_s"] = _ensure_numeric_dict_values(inlet_flows)
    return result


# ─────────────────────────────────────────────────────────────────
# 包装函数：精馏塔 — 确保流量字典、纯度、温度/压力参数为数值类型
# ─────────────────────────────────────────────────────────────────
def run_distillation_column(**kwargs):
    """精馏塔包装函数：确保所有数值参数类型正确。"""
    # 确保流量字典中的值为数值
    if "feed_molar_flows_mol_per_s" in kwargs:
        kwargs["feed_molar_flows_mol_per_s"] = _ensure_numeric_dict_values(
            kwargs["feed_molar_flows_mol_per_s"]
        )
    # 确保标量参数为数值
    for key in ("distillate_purity", "bottoms_purity",
                "feed_temp_K", "feed_pressure_Pa", "reflux_factor"):
        if key in kwargs:
            kwargs[key] = _ensure_numeric(kwargs[key])
    return _run_distillation_column(**kwargs)

# ─────────────────────────────────────────────────────────────────
# 包装函数：混合器 — 确保流量列表和温度/压力参数类型安全
# ─────────────────────────────────────────────────────────────────
def run_mixer(**kwargs):
    """混合器包装函数：确保多股进料参数类型正确。"""
    # 确保 inlet_molar_flows_list 中每个字典的值为数值
    if "inlet_molar_flows_list" in kwargs:
        flows_list = kwargs["inlet_molar_flows_list"]
        if isinstance(flows_list, list):
            kwargs["inlet_molar_flows_list"] = [
                _ensure_numeric_dict_values(f) if isinstance(f, dict) else f
                for f in flows_list
            ]
    # 确保温度/压力列表为数值
    for key in ("inlet_temperatures_K", "inlet_pressures_Pa"):
        if key in kwargs and isinstance(kwargs[key], list):
            kwargs[key] = [_ensure_numeric(v) for v in kwargs[key]]
    return _run_mixer(**kwargs)


# ─────────────────────────────────────────────────────────────────
# 延迟导入 UniversalToolRegistry / ToolCapability，避免循环依赖
# ─────────────────────────────────────────────────────────────────

def register_process_tools(registry) -> None:
    """将 process_engine 的 9 种设备级工具注册到 UniversalToolRegistry。

    Args:
        registry: UniversalToolRegistry 实例
    """
    from tools.registry import ToolCapability

    # ═════════════════════════════════════════════════════════════
    # 设备 1：泵
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="run_pump",
            description=(
                "泵升压设备。将流体从低压升至目标高压，返回出口压力和压升。"
                "只能升压（target > inlet）。"
                "出口参数 outlet_pressure_Pa 可直接作为下游设备的入口压力。"
            ),
            parameters={
                "inlet_pressure_Pa": {
                    "type": "number", "required": True,
                    "description": "泵入口压力 (Pa)，常压约 101325 Pa",
                },
                "target_pressure_Pa": {
                    "type": "number", "required": True,
                    "description": "目标出口压力 (Pa)，必须大于入口压力",
                },
            },
            returns={
                "outlet_pressure_Pa": "出口压力 (Pa)",
                "pressure_rise_Pa": "压升 ΔP (Pa)",
            },
            category="process_equipment",
            tags=["pump", "pressure", "升压", "泵"],
            examples=[
                {
                    "description": "常压进料加压至 2.8 MPa",
                    "parameters": {
                        "inlet_pressure_Pa": 101325,
                        "target_pressure_Pa": 2800000,
                    },
                }
            ],
        ),
        run_pump,
    )

    # ═════════════════════════════════════════════════════════════
    # 设备 2：减压阀
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="run_pressure_reducing_valve",
            description=(
                "减压阀降压设备。通过节流作用将流体从高压降至目标低压。"
                "只能降压（target < inlet）。"
                "出口参数 outlet_pressure_Pa 可直接作为下游设备的入口压力。"
                "可选传入 inlet_molar_flows_mol_per_s 透传流量，用于下游设备引用。"
            ),
            parameters={
                "inlet_pressure_Pa": {
                    "type": "number", "required": True,
                    "description": "入口压力 (Pa)",
                },
                "target_pressure_Pa": {
                    "type": "number", "required": True,
                    "description": "目标出口压力 (Pa)，必须小于入口压力",
                },
                "inlet_molar_flows_mol_per_s": {
                    "type": "object", "required": False,
                    "description": "入口各组分摩尔流量字典 (mol/s)，可选，透传至出口",
                },
            },
            returns={
                "outlet_pressure_Pa": "出口压力 (Pa)",
                "pressure_drop_Pa": "压降 ΔP (Pa)",
                "outlet_molar_flows_mol_per_s": "出口各组分流量字典 (mol/s)（透传，需传入 inlet）",
            },
            category="process_equipment",
            tags=["valve", "pressure", "降压", "减压阀"],
        ),
        run_pressure_reducing_valve,
    )

    # ═════════════════════════════════════════════════════════════
    # 设备 3：换热器
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="run_heat_exchanger",
            description=(
                "换热器设备，加热或冷却工艺流体。"
                "返回热负荷 Q、对数平均温差 LMTD、总传热系数 K、换热面积 A。"
                "加热时公用工程温度必须高于工艺流体温度；冷却时反之。"
                "utility_type 可选 'cooling_water'（循环水）或 'steam'（蒸汽）。"
                "换热器不改变组分，入口流量自动透传为出口流量 outlet_molar_flows_mol_per_s。"
            ),
            parameters={
                "process_fluid_molar_flows_mol_per_s": {
                    "type": "object", "required": True,
                    "description": "工艺流体各组分摩尔流量字典 (mol/s)，如 {'Benzene': 3000, 'Propylene': 1000}",
                },
                "process_fluid_temp_in_K": {
                    "type": "number", "required": True,
                    "description": "工艺流体入口温度 (K)",
                },
                "process_fluid_temp_out_K": {
                    "type": "number", "required": True,
                    "description": "工艺流体出口温度 (K)",
                },
                "process_fluid_pressure_Pa": {
                    "type": "number", "required": True,
                    "description": "工艺流体压力 (Pa)",
                },
                "utility_fluid_temp_in_K": {
                    "type": "number", "required": True,
                    "description": "公用工程入口温度 (K)",
                },
                "utility_fluid_temp_out_K": {
                    "type": "number", "required": True,
                    "description": "公用工程出口温度 (K)",
                },
                "utility_type": {
                    "type": "string", "required": False,
                    "description": "公用工程类型: 'cooling_water'（默认）或 'steam'",
                },
            },
            returns={
                "heat_duty_W": "热负荷 (W)",
                "log_mean_temp_difference_K": "对数平均温差 (K)",
                "overall_k_selected_W_per_m2K": "总传热系数 (W/(m²·K))",
                "heat_transfer_area_m2": "换热面积 (m²)",
                "outlet_molar_flows_mol_per_s": "出口各组分流量字典 (mol/s)（自动透传）",
                "outlet_temp_K": "出口温度 (K)",
                "outlet_pressure_Pa": "出口压力 (Pa)",
            },
            category="process_equipment",
            tags=["heat_exchanger", "换热器", "加热", "冷却"],
            examples=[
                {
                    "description": "苯+丙烯混合物从 25℃ 加热至 160℃，蒸汽加热",
                    "parameters": {
                        "process_fluid_molar_flows_mol_per_s": {"Benzene": 3000, "Propylene": 1000},
                        "process_fluid_temp_in_K": 298.15,
                        "process_fluid_temp_out_K": 433.15,
                        "process_fluid_pressure_Pa": 2800000,
                        "utility_fluid_temp_in_K": 523.15,
                        "utility_fluid_temp_out_K": 503.15,
                        "utility_type": "steam",
                    },
                }
            ],
        ),
        run_heat_exchanger,
    )

    # ═════════════════════════════════════════════════════════════
    # 设备 4：化学计量反应器
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="run_stoichiometric_reactor",
            description=(
                "化学计量反应器。按化学计量比进行反应，给定关键组分转化率。"
                "返回反应进度、出口各组分流量、出口总流量、出口体积流量、出口浓度。"
                "出口参数 outlet_molar_flows_mol_per_s 可直接作为下游设备的入口流量。"
                "入口流量必须包含所有组分（包括产物，初始为 0）。"
            ),
            parameters={
                "inlet_molar_flows_mol_per_s": {
                    "type": "object", "required": True,
                    "description": "入口各组分摩尔流量字典 (mol/s)，必须包含所有组分（含产物，流量为 0）",
                },
                "stoichiometric_coefficients": {
                    "type": "object", "required": True,
                    "description": "化学计量系数字典，反应物为负、产物为正",
                },
                "key_component": {
                    "type": "string", "required": True,
                    "description": "关键组分名称（必须在 inlet 和 stoichiometry 中）",
                },
                "key_component_conversion": {
                    "type": "number", "required": True,
                    "description": "关键组分转化率 (0~1)",
                },
                "reactor_temp_K": {
                    "type": "number", "required": True,
                    "description": "反应温度 (K)",
                },
                "reactor_pressure_Pa": {
                    "type": "number", "required": True,
                    "description": "反应压力 (Pa)",
                },
            },
            returns={
                "reaction_extent_mol_per_s": "反应进度 (mol/s)",
                "outlet_molar_flows_mol_per_s": "出口各组分流量字典 (mol/s)",
                "outlet_total_molar_flow_mol_per_s": "出口总摩尔流量 (mol/s)",
                "outlet_volumetric_flow_m3_per_s": "出口体积流量 (m³/s)",
                "outlet_concentrations_mol_per_m3": "出口各组分浓度字典 (mol/m³)",
            },
            category="process_equipment",
            tags=["reactor", "stoichiometric", "化学计量反应器"],
            aliases=["reactor_design", "stoichiometric_reactor", "run_reactor"],
        ),
        run_stoichiometric_reactor,
    )

    # ═════════════════════════════════════════════════════════════
    # 设备 5：CSTR 反应器
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="run_cstr_reactor",
            description=(
                "CSTR 全混流反应器。考虑催化剂和动力学参数，返回 12 个底层工具的全部结果。"
                "包括反应进度、出口流量、流体物性（密度/粘度）、流化参数、反应速率、"
                "反应器体积、停留时间。"
                "出口参数 outlet_molar_flows_mol_per_s 可直接作为下游设备的入口流量。"
            ),
            parameters={
                "inlet_molar_flows_mol_per_s": {
                    "type": "object", "required": True,
                    "description": "入口各组分摩尔流量字典 (mol/s)，含产物（流量为 0）",
                },
                "stoichiometric_coefficients": {
                    "type": "object", "required": True,
                    "description": "化学计量系数，反应物负、产物正",
                },
                "key_component": {
                    "type": "string", "required": True,
                    "description": "关键组分名称",
                },
                "key_component_conversion": {
                    "type": "number", "required": True,
                    "description": "关键组分转化率 (0~1)",
                },
                "reactor_temp_K": {
                    "type": "number", "required": True,
                    "description": "反应温度 (K)",
                },
                "reactor_pressure_Pa": {
                    "type": "number", "required": True,
                    "description": "反应压力 (Pa)",
                },
                "catalyst_particle_diameter_m": {
                    "type": "number", "required": False,
                    "description": "催化剂颗粒直径 (m)，默认 0.003 (3mm)",
                },
                "catalyst_particle_density_kg_per_m3": {
                    "type": "number", "required": False,
                    "description": "催化剂颗粒密度 (kg/m³)，默认 2500",
                },
                "arrhenius_pre_exponential_factor": {
                    "type": "number", "required": False,
                    "description": "Arrhenius 指前因子 A，默认 1.0（需根据实际反应调整）",
                },
                "activation_energy_J_per_mol": {
                    "type": "number", "required": False,
                    "description": "活化能 Ea (J/mol)，默认 50000 (50kJ/mol)",
                },
                "reaction_order": {
                    "type": "number", "required": False,
                    "description": "反应级数（默认 1）",
                },
            },
            returns={
                "reaction_extent_mol_per_s": "反应进度 (mol/s)",
                "outlet_molar_flows_mol_per_s": "出口各组分流量字典 (mol/s)",
                "outlet_total_molar_flow_mol_per_s": "出口总摩尔流量 (mol/s)",
                "outlet_volumetric_flow_m3_per_s": "出口体积流量 (m³/s)",
                "outlet_concentrations_mol_per_m3": "出口各组分浓度字典 (mol/m³)",
                "fluid_density_kg_per_m3": "流体密度 (kg/m³)",
                "fluid_viscosity_Pa_s": "流体粘度 (Pa·s)",
                "archimedes_number": "阿基米德数",
                "minimum_fluidization_velocity_m_per_s": "最小流化速度 (m/s)",
                "reaction_rate_mol_per_m3_s": "反应速率 (mol/(m³·s))",
                "cstr_volume_m3": "CSTR 体积 (m³)",
                "residence_time_s": "停留时间 (s)",
            },
            category="process_equipment",
            tags=["reactor", "cstr", "全混流反应器"],
            aliases=["cstr_reactor", "cstr"],
        ),
        run_cstr_reactor,
    )

    # ═════════════════════════════════════════════════════════════
    # 设备 6：PFR 反应器
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="run_pfr_reactor",
            description=(
                "PFR 平推流反应器（管式/固定床）。考虑催化剂床层和动力学参数，"
                "返回 14 个底层工具的全部结果。"
                "包括反应进度、出口流量、入口参数、床层体积、催化剂质量、"
                "管数、床层长度等。"
                "出口参数 outlet_molar_flows_mol_per_s 可直接作为下游设备的入口流量。"
            ),
            parameters={
                "inlet_molar_flows_mol_per_s": {
                    "type": "object", "required": True,
                    "description": "入口各组分摩尔流量字典 (mol/s)，含产物（流量为 0）",
                },
                "stoichiometric_coefficients": {
                    "type": "object", "required": True,
                    "description": "化学计量系数，反应物负、产物正",
                },
                "key_component": {
                    "type": "string", "required": True,
                    "description": "关键组分名称",
                },
                "key_component_conversion": {
                    "type": "number", "required": True,
                    "description": "关键组分转化率 (0~1)",
                },
                "reactor_temp_K": {
                    "type": "number", "required": True,
                    "description": "反应温度 (K)",
                },
                "reactor_pressure_Pa": {
                    "type": "number", "required": True,
                    "description": "反应压力 (Pa)",
                },
                "catalyst_particle_density_kg_per_m3": {
                    "type": "number", "required": False,
                    "description": "催化剂颗粒密度 (kg/m³)，默认 2500",
                },
                "catalyst_bed_void_fraction": {
                    "type": "number", "required": False,
                    "description": "催化剂床层空隙率 (0~1)，默认 0.4",
                },
                "arrhenius_pre_exponential_factor": {
                    "type": "number", "required": False,
                    "description": "Arrhenius 指前因子 A，默认 1.0（需根据实际反应调整）",
                },
                "activation_energy_J_per_mol": {
                    "type": "number", "required": False,
                    "description": "活化能 Ea (J/mol)，默认 50000 (50kJ/mol)",
                },
                "reaction_order": {
                    "type": "number", "required": False,
                    "description": "反应级数（默认 1）",
                },
                "tube_inner_diameter_m": {
                    "type": "number", "required": False,
                    "description": "反应管内径 (m)，默认 0.05",
                },
                "superficial_velocity_m_per_s": {
                    "type": "number", "required": False,
                    "description": "表观流速 (m/s)，默认 1.0",
                },
            },
            returns={
                "reaction_extent_mol_per_s": "反应进度 (mol/s)",
                "outlet_molar_flows_mol_per_s": "出口各组分流量字典 (mol/s)",
                "outlet_total_molar_flow_mol_per_s": "出口总摩尔流量 (mol/s)",
                "outlet_volumetric_flow_m3_per_s": "出口体积流量 (m³/s)",
                "outlet_concentrations_mol_per_m3": "出口各组分浓度字典 (mol/m³)",
                "bed_volume_m3": "床层体积 (m³)",
                "catalyst_mass_kg": "催化剂质量 (kg)",
                "tube_count": "管数",
                "bed_length_m": "床层长度 (m)",
            },
            category="process_equipment",
            tags=["reactor", "pfr", "平推流反应器", "固定床"],
            aliases=["pfr_reactor", "pfr", "fixed_bed_reactor"],
        ),
        run_pfr_reactor,
    )

    # ═════════════════════════════════════════════════════════════
    # 设备 7：精馏塔
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="run_distillation_column",
            description=(
                "板式精馏塔。完成物料衡算、操作条件、进料热状态、最小回流比（Underwood）、"
                "最小理论板（Fenske）、理论板数（Gilliland）、实际塔板数、进料板位置（Kirkbride）"
                "8 步完整计算。"
                "出口参数 distillate_flows_mol_per_s / bottoms_flows_mol_per_s 可作为下游设备入口。"
            ),
            parameters={
                "feed_molar_flows_mol_per_s": {
                    "type": "object", "required": True,
                    "description": "进料各组分摩尔流量字典 (mol/s)",
                },
                "distillate_purity": {
                    "type": "number", "required": True,
                    "description": "塔顶轻关键组分质量纯度 (0~1)",
                },
                "bottoms_purity": {
                    "type": "number", "required": True,
                    "description": "塔底重关键组分质量纯度 (0~1)",
                },
                "light_key_component": {
                    "type": "string", "required": True,
                    "description": "轻关键组分名称",
                },
                "heavy_key_component": {
                    "type": "string", "required": True,
                    "description": "重关键组分名称",
                },
                "light_components": {
                    "type": "array", "required": True,
                    "description": "轻组分列表（含轻关键组分）",
                },
                "heavy_components": {
                    "type": "array", "required": True,
                    "description": "重组分列表（含重关键组分）",
                },
                "feed_temp_K": {
                    "type": "number", "required": True,
                    "description": "进料温度 (K)",
                },
                "feed_pressure_Pa": {
                    "type": "number", "required": True,
                    "description": "进料压力 (Pa)",
                },
                "reflux_factor": {
                    "type": "number", "required": False,
                    "description": "回流比倍数 R/R_min（默认 1.2）",
                },
            },
            returns={
                "distillate_flows_mol_per_s": "塔顶各组分流量字典 (mol/s)",
                "bottoms_flows_mol_per_s": "塔底各组分流量字典 (mol/s)",
                "R_operating": "操作回流比",
                "N_actual": "实际塔板数",
                "feed_tray_position_from_top": "进料板位置（从塔顶计）",
                "column_top_temperature_K": "塔顶温度 (K)",
                "column_bottom_temperature_K": "塔底温度 (K)",
            },
            category="process_equipment",
            tags=["distillation", "精馏塔", "分离"],
        ),
        run_distillation_column,
    )

    # ═════════════════════════════════════════════════════════════
    # 设备 8：闪蒸罐
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="run_flash_drum",
            description=(
                "闪蒸罐气液分离设备。给定闪蒸压力和关键组分回收率，二分法反算闪蒸温度，"
                "返回气液两相各组分流量。"
                "单自由度系统：给定 P 和进料后，一个回收率约束唯一确定 T。"
                "出口参数 vapor_molar_flows_mol_per_s / liquid_molar_flows_mol_per_s "
                "可作为下游设备入口。"
            ),
            parameters={
                "flash_pressure_Pa": {
                    "type": "number", "required": True,
                    "description": "闪蒸压力 (Pa)",
                },
                "inlet_molar_flows_mol_per_s": {
                    "type": "object", "required": True,
                    "description": "进料各组分摩尔流量字典 (mol/s)，至少 2 组分",
                },
                "key_component": {
                    "type": "string", "required": True,
                    "description": "关键组分名称（必须在进料中）",
                },
                "key_component_recovery": {
                    "type": "number", "required": True,
                    "description": "关键组分目标回收率 (0~1，开区间)",
                },
                "recovery_phase": {
                    "type": "string", "required": False,
                    "description": "回收相: 'vapor'（默认）或 'liquid'",
                },
            },
            returns={
                "flash_temperature_K": "闪蒸温度 (K)",
                "vapor_molar_flows_mol_per_s": "气相各组分流量字典 (mol/s)",
                "liquid_molar_flows_mol_per_s": "液相各组分流量字典 (mol/s)",
            },
            category="process_equipment",
            tags=["flash", "闪蒸罐", "气液分离"],
        ),
        run_flash_drum,
    )

    # ═════════════════════════════════════════════════════════════
    # 设备 9：混合器
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="run_mixer",
            description=(
                "混合器设备。将 N 股进料合并为 1 股出料，计算绝热混合温度、出口压力、"
                "出口流量与组成。"
                "出口参数 outlet_molar_flows_mol_per_s 可直接作为下游设备的入口流量。"
            ),
            parameters={
                "inlet_molar_flows_list": {
                    "type": "array", "required": True,
                    "description": (
                        "多股进料的摩尔流量列表，每个元素为组分流量字典 (mol/s)。"
                        "例如 [{'H2O': 100, 'MeOH': 50}, {'H2O': 200, 'MeOH': 30}]"
                    ),
                },
                "inlet_temperatures_K": {
                    "type": "array", "required": True,
                    "description": "各股进料温度列表 (K)，长度须与 inlet_molar_flows_list 一致",
                },
                "inlet_pressures_Pa": {
                    "type": "array", "required": True,
                    "description": "各股进料压力列表 (Pa)，长度须与 inlet_molar_flows_list 一致",
                },
            },
            returns={
                "outlet_molar_flows_mol_per_s": "出口各组分流量字典 (mol/s)",
                "outlet_total_molar_flow_mol_per_s": "出口总摩尔流量 (mol/s)",
                "outlet_composition": "出口组成字典（摩尔分率）",
                "outlet_temperature_K": "出口绝热混合温度 (K)",
                "outlet_pressure_Pa": "出口压力 (Pa)（取各股最低值）",
                "outlet_volumetric_flow_m3_per_s": "出口体积流量 (m³/s)",
            },
            category="process_equipment",
            tags=["mixer", "混合器", "合并", "多股进料"],
            examples=[
                {
                    "description": "两股进料（水+甲醇）混合",
                    "parameters": {
                        "inlet_molar_flows_list": [
                            {"water": 100, "methanol": 50},
                            {"water": 200, "methanol": 30},
                        ],
                        "inlet_temperatures_K": [350.0, 320.0],
                        "inlet_pressures_Pa": [200000.0, 250000.0],
                    },
                }
            ],
        ),
        run_mixer,
    )

    # ═════════════════════════════════════════════════════════════
    # 设备 10：流程求解器（带循环收敛）
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="run_flowsheet",
            description=(
                "流程求解器 — 带循环收敛的完整化工流程模拟。"
                "基于配置字典描述流程拓扑（设备连接 + 循环流），自动识别 Tear Stream，"
                "使用 Wegstein 割线法加速迭代至收敛。"
                "支持设备类型: mixer/pump/valve/heater/cooler/reactor/distillation/flash。"
                "返回收敛状态、所有流股、设备摘要及完整后处理结果（换热面积、塔板数等）。"
            ),
            parameters={
                "flowsheet_config": {
                    "type": "object", "required": True,
                    "description": (
                        "流程配置字典，包含:\n"
                        "  - feeds: dict — 新鲜进料 {stream_name: {component: flow_mol_per_s, _T: K, _P: Pa}}\n"
                        "  - equipment: list[dict] — 设备列表（按执行顺序），每个设备包含 id/type/inlet/outlet/params\n"
                        "  - recycles: list[dict] — 循环流定义 [{from: 源流股, to: 目标流股名}]\n"
                        "  - defaults: dict(可选) — 默认温度/压力"
                    ),
                },
                "max_iter": {
                    "type": "number", "required": False,
                    "description": "最大迭代次数（默认 500）",
                },
                "tol": {
                    "type": "number", "required": False,
                    "description": "收敛容差 — 循环流总残差（默认 1e-4）",
                },
                "acceleration": {
                    "type": "string", "required": False,
                    "description": "加速方法: 'wegstein'（默认）/ 'successive' / 'none'",
                },
                "omega": {
                    "type": "number", "required": False,
                    "description": "阻尼因子（默认 0.5）",
                },
                "do_post_process": {
                    "type": "boolean", "required": False,
                    "description": "是否执行后处理获取完整设备设计参数（默认 true）",
                },
            },
            returns={
                "converged": "是否收敛 (bool)",
                "iterations": "迭代次数 (int)",
                "residual": "最终残差 (float)",
                "streams": "所有流股字典 {name: {component: flow, _T: K, _P: Pa}}",
                "feeds": "新鲜进料字典",
                "equipment_summary": "设备执行摘要列表",
                "post_process": "完整设备设计参数 {eq_id: {heat_duty, area, N_actual, ...}}",
            },
            category="process_equipment",
            tags=["flowsheet", "流程模拟", "循环收敛", "Wegstein", "带回流"],
        ),
        run_flowsheet,
    )

    # ═════════════════════════════════════════════════════════════
    # 特殊工具：human_ask — 请求人工补充缺失参数
    # ═════════════════════════════════════════════════════════════
    registry.register(
        ToolCapability(
            name="human_ask",
            description=(
                "人工参数请求工具。当核心工艺参数缺失且无法使用工程默认值时，"
                "调用此工具暂停流程并请求用户补充参数。"
                "注意：只有核心参数（如 reactor_type、X_target、k、T_in、components 等）缺失时才调用，"
                "辅助参数（如 rho_cat、particle_diameter、U 等）应使用默认值，不要调用此工具。"
            ),
            parameters={
                "question": {
                    "type": "string", "required": True,
                    "description": "向用户提出的问题，说明需要补充什么参数",
                },
                "missing_params": {
                    "type": "array", "required": False,
                    "description": "缺失的参数名列表",
                },
                "context": {
                    "type": "string", "required": False,
                    "description": "补充说明上下文，帮助用户理解需要提供什么信息",
                },
            },
            returns={
                "needs_human": "始终为 True，表示需要人工输入",
                "human_prompt": "展示给用户的问题文本",
            },
            category="system",
            tags=["human_ask", "人工输入", "参数补充", "human-in-the-loop"],
        ),
        human_ask,
    )


__all__ = ["register_process_tools"]
