from __future__ import annotations
import re
from typing import Dict, Any, List, TYPE_CHECKING

if TYPE_CHECKING:
    from graph.state import AgentState


class ExecutorAgent:
    """执行智能体 — 通用工具执行器"""

    def __init__(self, registry=None):
        self.registry = registry
        self.name = "Executor"

    def __call__(self, state): return self.execute(state)

    def _resolve_param_references(self, value: Any, execution_results: dict) -> Any:
        """解析参数中的引用表达式，替换为实际值
        
        支持的引用格式:
        - "$tool_name" - 引用工具返回的完整结果
        - "$tool_name.field" - 引用工具返回结果中的特定字段
        - "$tool_name.field.subfield" - 引用嵌套字典字段（如 $run_reactor.outlet_molar_flows_mol_per_s.Benzene）
        - "$tool_name.result" - 引用工具返回的 result 字段值
        
        当值为字典或列表时，递归解析其中包含的引用表达式。
        支持回退查找：当 tool_name 不存在时，自动搜索 tool_name#N 后缀的 key。
        支持字段回溯：当字段在当前工具结果中不存在时，反向搜索上游工具的结果。
        """
        if isinstance(value, dict):
            return {k: self._resolve_param_references(v, execution_results) for k, v in value.items()}
        if isinstance(value, list):
            return [self._resolve_param_references(item, execution_results) for item in value]
        if not isinstance(value, str):
            return value
        
        # 检查是否是引用表达式（以 $ 开头）
        if not value.startswith("$"):
            return value
        
        # 解析引用: $tool_name 或 $tool_name#N 或 $tool_name.field1.field2...
        # 支持 #N 后缀（同名工具多次调用时的唯一标识）
        ref_pattern = r'^\$([\w]+(?:#\d+)?)(?:\.(.+))?$'
        match = re.match(ref_pattern, value)
        if not match:
            return value
        
        tool_name = match.group(1)
        field_path = match.group(2)  # 可能是 "field" 或 "field.subfield" 或 None
        
        # 从执行结果中查找（支持回退查找 #N 后缀和 _N 后缀）
        tool_result = None
        if tool_name in execution_results:
            tool_result = execution_results[tool_name]
        else:
            # 回退查找 1：搜索 tool_name#N 后缀的 key，取最后一个（最新的）
            for key in sorted(execution_results.keys(), reverse=True):
                if key.startswith(f"{tool_name}#"):
                    tool_result = execution_results[key]
                    break
            # 回退查找 2：Planner 可能错误地使用 _N 后缀（如 run_heat_exchanger_1）
            # 尝试剥离尾部 _数字 后缀，查找基础工具名
            if tool_result is None:
                base_match = re.match(r'^(.+)_(\d+)$', tool_name)
                if base_match:
                    base_name = base_match.group(1)
                    # 优先查找 base_name#N（最新执行的同名工具）
                    for key in sorted(execution_results.keys(), reverse=True):
                        if key.startswith(f"{base_name}#"):
                            tool_result = execution_results[key]
                            break
                    # 如果没有 #N 变体，再查找基础名
                    if tool_result is None and base_name in execution_results:
                        tool_result = execution_results[base_name]
        
        if tool_result is None:
            return value  # 找不到引用，保持原值
        
        if field_path is None:
            # 没有指定字段，返回完整结果
            return tool_result
        
        # 按点号分割字段路径，支持嵌套访问
        fields = field_path.split(".")
        current = tool_result
        resolved = False
        
        for i, field in enumerate(fields):
            if isinstance(current, dict) and field in current:
                current = current[field]
                resolved = (i == len(fields) - 1)
            elif i == 0 and isinstance(current, dict) and "result" in current and isinstance(current["result"], dict):
                # 自动进入 result 子字典（处理 {"success": True, "result": {...}} 结构）
                inner = current["result"]
                if field in inner:
                    current = inner[field]
                    resolved = (i == len(fields) - 1)
                elif field == "result":
                    current = inner
                    resolved = False
                else:
                    # 字段不存在，触发回溯查找
                    current = self._traceback_field(field, fields[i+1:], tool_name, execution_results)
                    if current is not None:
                        resolved = True
                    else:
                        return value
            else:
                # 字段不存在，触发回溯查找
                if i == 0:
                    current = self._traceback_field(field, fields[i+1:], tool_name, execution_results)
                else:
                    current = None
                if current is not None:
                    resolved = True
                else:
                    return value
        
        return current

    def _traceback_field(self, field: str, remaining_fields: list, exclude_tool: str, execution_results: dict):
        """回溯查找：在所有已执行工具中搜索包含指定字段的最近工具。
        
        当 $tool_A.field 中 tool_A 不包含 field 时，自动搜索上游工具。
        """
        # 按执行顺序反向遍历（最新的优先）
        for key in reversed(list(execution_results.keys())):
            if key == exclude_tool:
                continue  # 跳过原始工具（已经找不到了）
            result = execution_results[key]
            # 检查 result 本身
            if isinstance(result, dict):
                if field in result:
                    current = result[field]
                elif "result" in result and isinstance(result["result"], dict) and field in result["result"]:
                    current = result["result"][field]
                else:
                    continue
                # 继续解析剩余字段
                for rf in remaining_fields:
                    if isinstance(current, dict) and rf in current:
                        current = current[rf]
                    else:
                        current = None
                        break
                if current is not None:
                    return current
        return None

    def _resolve_parameters(self, parameters: dict, execution_results: dict) -> dict:
        """解析参数字典中所有的引用表达式"""
        resolved = {}
        for key, value in parameters.items():
            resolved[key] = self._resolve_param_references(value, execution_results)
        return resolved

    def _convert_param_types(self, parameters: dict) -> dict:
        """将字符串类型的数值转换为数字，防止 LLM 输出格式错误。
        
        递归处理嵌套字典和列表中的字符串数值。
        """
        converted = {}
        for key, value in parameters.items():
            # 清理键名：去除 LLM 可能输出的多余引号
            if isinstance(key, str):
                key = key.strip().strip("'\"").strip()
            converted[key] = self._convert_value_recursive(value)
        return converted

    @staticmethod
    def _convert_value_recursive(value: Any) -> Any:
        """递归将字符串数值转换为数字，处理嵌套 dict/list。"""
        if isinstance(value, dict):
            return {k: ExecutorAgent._convert_value_recursive(v) for k, v in value.items()}
        if isinstance(value, list):
            return [ExecutorAgent._convert_value_recursive(item) for item in value]
        if not isinstance(value, str):
            return value
        # 尝试转换为整数
        try:
            return int(value)
        except ValueError:
            pass
        # 尝试转换为浮点数
        try:
            return float(value)
        except ValueError:
            pass
        return value

    def execute(self, state: AgentState) -> AgentState:
        """执行计划中的工具任务"""
        plan = state.get("plan", {})
        tasks = plan.get("tasks", [])
        current_index = plan.get("current_task_index", 0)

        if not tasks:
            state["execution_results"] = {"error": "没有可执行的任务"}
            state["current_phase"] = "checking"
            state["messages"].append({
                "role": "executor",
                "content": "没有可执行的任务",
            })
            return state

        # 获取当前要执行的任务
        if current_index >= len(tasks):
            # 所有任务已完成
            state["current_phase"] = "checking"
            state["messages"].append({
                "role": "executor",
                "content": "所有任务已执行完成",
            })
            return state

        current_task = tasks[current_index]
        tool_name = current_task.get("tool", "")
        parameters = current_task.get("parameters", {})
        purpose = current_task.get("purpose", "")

        # 规范化工具名：剥离 Planner 可能错误添加的 _N 后缀
        # 例如 run_heat_exchanger_1 → run_heat_exchanger
        normalized_tool_name = re.sub(r'_\d+$', '', tool_name)

        # 解析参数中的引用表达式（如 $get_fluid_density.result）
        execution_results = state.get("execution_results") or {}
        resolved_params = self._resolve_parameters(parameters, execution_results)

        # 类型转换保护：将字符串类型的数值转换为数字
        resolved_params = self._convert_param_types(resolved_params)

        # 执行工具（使用解析后的参数，传入原始 tool_name 以便 _execute_tool 也能处理）
        result = self._execute_tool(tool_name, resolved_params)

        # 收集结果（确保 execution_results 是字典，不是 None）
        if not state.get("execution_results"):
            state["execution_results"] = {}

        # 使用唯一 key 防止同名工具多次调用时结果被覆盖
        # 使用规范化后的工具名（已剥离 _N 后缀）
        # 格式: "tool_name" (首次) 或 "tool_name#N" (后续同名)
        if normalized_tool_name in state["execution_results"]:
            # 已存在同名工具结果，用带序号的 key
            unique_key = f"{normalized_tool_name}#{current_index}"
        else:
            unique_key = normalized_tool_name
        state["execution_results"][unique_key] = result
        state["validation_results"] = result.get("validation") or {}

        # 更新任务索引
        plan["current_task_index"] = current_index + 1
        state["plan"] = plan

        # 检查是否还有更多任务
        if current_index + 1 < len(tasks):
            state["current_phase"] = "executing"  # 继续执行下一个任务
        else:
            state["current_phase"] = "checking"  # 所有任务完成，交给 Planner 检查

        # 生成结果摘要消息
        summary = self._build_result_summary(tool_name, purpose, result)
        state["messages"].append({
            "role": "executor",
            "content": summary,
        })
        return state

    def _build_result_summary(self, tool_name: str, purpose: str, result: dict) -> str:
        """构建工具执行结果摘要，包含关键输出参数。"""
        if not result.get("success"):
            return f"工具 {tool_name} 执行失败: {result.get('error', 'unknown')}"

        r = result.get("result", {})
        if not isinstance(r, dict):
            return f"工具 {tool_name} 执行完成: {purpose}"

        parts = [f"工具 {tool_name} 执行完成: {purpose}"]

        # 根据不同工具类型提取关键结果
        # 流量类
        for flow_key in ("outlet_molar_flows_mol_per_s", "distillate_flows_mol_per_s",
                         "bottoms_flows_mol_per_s", "vapor_molar_flows_mol_per_s",
                         "liquid_molar_flows_mol_per_s"):
            if flow_key in r and isinstance(r[flow_key], dict):
                flows = r[flow_key]
                flow_str = ", ".join(f"{k}={v:.1f}" for k, v in flows.items() if isinstance(v, (int, float)))
                label = {
                    "outlet_molar_flows_mol_per_s": "出口流量",
                    "distillate_flows_mol_per_s": "塔顶流量",
                    "bottoms_flows_mol_per_s": "塔底流量",
                    "vapor_molar_flows_mol_per_s": "气相流量",
                    "liquid_molar_flows_mol_per_s": "液相流量",
                }.get(flow_key, flow_key)
                parts.append(f"  {label}(mol/s): {flow_str}")

        # 压力类
        for p_key, p_label in [("outlet_pressure_Pa", "出口压力"),
                                ("outlet_temp_K", "出口温度")]:
            if p_key in r and isinstance(r[p_key], (int, float)):
                unit = "Pa" if "pressure" in p_key else "K"
                parts.append(f"  {p_label}: {r[p_key]:.2f} {unit}")

        # 换热器
        if "heat_duty_W" in r:
            parts.append(f"  热负荷: {r['heat_duty_W']:.0f} W")
        if "heat_transfer_area_m2" in r:
            parts.append(f"  换热面积: {r['heat_transfer_area_m2']:.2f} m²")

        # 反应器
        if "reaction_extent_mol_per_s" in r:
            parts.append(f"  反应进度: {r['reaction_extent_mol_per_s']:.1f} mol/s")

        # 精馏塔
        if "N_actual" in r:
            parts.append(f"  实际塔板数: {r['N_actual']:.0f}")
        if "R_operating" in r:
            parts.append(f"  操作回流比: {r['R_operating']:.2f}")

        # 泵
        if "hydraulic_power_W" in r:
            parts.append(f"  水力功率: {r['hydraulic_power_W']:.0f} W")

        # 减压阀
        if "pressure_drop_Pa" in r:
            parts.append(f"  压降: {r['pressure_drop_Pa']:.0f} Pa")

        return "\n".join(parts)

    def _execute_tool(self, tool_name: str, parameters: dict) -> dict:
        """执行单个工具"""
        if not self.registry:
            return {"error": "工具注册中心未配置"}

        # 尝试从通用注册中心获取
        if hasattr(self.registry, "execute"):
            result = self.registry.execute(tool_name, parameters)
            # 如果工具未找到，尝试剥离 Planner 错误添加的 _N 后缀
            if result.get("error") and "未找到" in str(result.get("error", "")):
                base_match = re.match(r'^(.+)_(\d+)$', tool_name)
                if base_match:
                    base_name = base_match.group(1)
                    result = self.registry.execute(base_name, parameters)
            return result

        # 尝试从传统注册中心获取
        if hasattr(self.registry, "get"):
            handler = self.registry.get(tool_name)
            if handler:
                try:
                    # 清理参数，移除内部标记
                    clean_params = {k: v for k, v in parameters.items() if not k.startswith("_")}
                    result = handler.execute(clean_params)
                    return {"success": True, "result": result}
                except Exception as e:
                    return {"success": False, "error": str(e)}

        return {"error": f"工具 '{tool_name}' 未找到"}


class MultiToolExecutor:
    """多工具执行器 — 支持复杂的工具链"""

    def __init__(self, registry=None):
        self.registry = registry
        self.name = "MultiToolExecutor"

    def execute_chain(self, tasks: List[Dict[str, Any]]) -> Dict[str, Any]:
        """执行工具链"""
        results = {}
        context = {}  # 共享上下文，用于工具间传递数据

        for i, task in enumerate(tasks):
            tool_name = task.get("tool", "")
            parameters = task.get("parameters", {})

            # 将上下文中的数据合并到参数
            merged_params = {**context, **parameters}

            # 执行工具
            result = self._execute_single(tool_name, merged_params)
            results[tool_name] = result

            # 如果执行成功，将结果添加到上下文
            if result.get("success") and result.get("result"):
                # 将结果展平添加到上下文
                tool_result = result["result"]
                if isinstance(tool_result, dict):
                    for k, v in tool_result.items():
                        if isinstance(v, (int, float, str, bool)):
                            context[f"{tool_name}_{k}"] = v

        return results

    def _execute_single(self, tool_name: str, parameters: dict) -> dict:
        """执行单个工具"""
        if not self.registry:
            return {"error": "工具注册中心未配置"}

        if hasattr(self.registry, "execute"):
            return self.registry.execute(tool_name, parameters)

        return {"error": f"工具 '{tool_name}' 未找到"}
