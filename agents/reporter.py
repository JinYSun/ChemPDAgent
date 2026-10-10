"""输出智能体 — 生成通用自然语言报告"""
from __future__ import annotations
import json
import re
from datetime import datetime
from typing import Optional, TYPE_CHECKING

from langchain_core.messages import HumanMessage, SystemMessage

if TYPE_CHECKING:
    from graph.state import AgentState

# ============================================================
# 提示模板
# ============================================================

REPORTER_SYSTEM_PROMPT = """你是一个通用任务结果报告专家。你的职责是根据用户的原始请求和工具的执行结果，生成一份结构清晰、内容完整、语言自然的报告。

=== 报告原则 ===
1. **直接回应用户需求** — 开头必须用一两句话直接回答用户的问题或目标，不要用空洞的开场白
2. **忠实于数据** — 只陈述工具实际返回的结果，不要推断或编造未出现的数据，即使计算错误也要给出最后一次迭代后的结果，列出计算过程的所有公式
3. **解释关键结果** — 对所有的中间过程和结果输出数值或结论，简要说明其含义和意义
4. **指出问题与建议** — 如果执行过程中出现警告、错误或不理想的结果，明确指出并给出改进建议
5. **结构分明** — 使用 Markdown 标题和列表组织内容，层次清晰
6. **多设备流程必须逐设备输出** — 当涉及多个设备的工艺流程时，必须按流程顺序逐个设备列出计算结果，包括每个设备的：入口条件、出口条件、关键设计参数（如换热面积、塔板数、压降等）、物料平衡。不得省略任何设备的结果。

=== 报告结构参考（根据实际内容灵活调整）===
1. **结论摘要** — 一段话直接回答用户的问题
2. **工程分析** — 如果上下文中有"工程分析推理"字段，将其整理为简洁的工程分析段落，说明方案选择的物性依据和技术判断。如果无此字段，则跳过此节。
3. **执行概况** — 调用了哪些工具，整体是否成功
4. **结果详情** — 按工具或逻辑模块列出关键输出，要列出所有的中间变量计算过程和结果，每个数值都附上详细解释
5. **问题与警告** — 如有错误或警告，列出原因及影响
6. **建议与后续步骤** — 基于结果给出可操作的后续建议

=== 缺失参数与默认值提示（重要）===
工具结果中可能包含以下字段，你必须在报告中明确告知用户：

- **`user_defaulted_params_report`**：这是最重要的提醒字段，记录了用户未提供而工具自动使用工程默认值的参数。你必须在报告中用专门的章节展示这些内容，告知用户：
  - 哪些参数使用了默认值及具体数值
  - 每个参数的工程参考范围
  - 建议用户根据实际工况提供参数以提高精度

- **`defaulted_params`**：记录了所有使用默认值的参数。与 user_defaulted_params_report 配合使用。

- **`suggestions`**：工具生成的工程建议，包括参数分类提示。你必须在报告的"建议"部分完整呈现。

- **`notes_missing_params`**：记录了哪些参数是通过自动推算得到的（如反应焓由赫斯定律推算、C_A0 由混合物密度推算）。你应在报告中说明这些参数的推算来源和可信度。

- **`conversion_log`**：记录了参数推算的详细过程（如密度→C_A0、生成焓→delta_H_rxn）。如果推算过程有参考价值，可简要提及。

请用与用户提问相同的语言输出报告。"""


class ReporterAgent:
    """输出智能体 — 生成通用自然语言报告"""

    def __init__(self, llm=None):
        self.name = "Reporter"
        self._llm = llm

    def __call__(self, state): return self.report(state)

    def report(self, state: AgentState) -> AgentState:
        """生成最终报告；LLM 可用时生成自然语言，否则回退到结构化文本"""
        llm = self._llm or self._get_default_llm()
        if llm is not None:
            result = self._report_with_llm(llm, state)
            if result:
                state["final_report"] = result
                state["current_phase"] = "done"
                state["messages"].append({"role": "reporter", "content": "LLM 报告已生成"})
                return state

        state["final_report"] = self._report_structured(state)
        state["current_phase"] = "done"
        state["messages"].append({"role": "reporter", "content": "结构化报告已生成（LLM 不可用）"})
        return state

    def _get_default_llm(self):
        from agents.planner import get_llm
        return get_llm()

    # ------------------------------------------------------------------ #
    # LLM 驱动的报告生成
    # ------------------------------------------------------------------ #

    def _report_with_llm(self, llm, state: dict) -> Optional[str]:
        plan            = state.get("plan", {})
        raw_results     = state.get("execution_results") or {}
        validation      = state.get("validation_results") or {}
        user_request    = state.get("user_request", "")
        errors          = state.get("errors") or []
        iteration_count = state.get("iteration", 0)

        simplified_results = self._normalize_results(raw_results)

        context = {
            "用户请求":     user_request,
            "任务理解":     plan.get("understanding", ""),
            "工程分析推理": plan.get("reasoning", "") or None,
            "备注与假设":   plan.get("notes", []),
            "缺失参数":     plan.get("missing_params", []),
            "工具执行结果": simplified_results,
            "校验信息":     validation if validation else None,
            "系统错误":     errors if errors else None,
            "迭代修正次数": iteration_count if iteration_count > 0 else None,
            "执行日志": [
                m.get("content", "")
                for m in state.get("messages", [])
                if m.get("role") in ("planner", "executor", "system")
            ],
        }
        context = {k: v for k, v in context.items() if v is not None}
        context_str = json.dumps(context, ensure_ascii=False, indent=2, default=str)

        # 约束违反时，在 prompt 中强调诊断输出
        diagnostic_hint = ""
        if iteration_count > 0:
            diagnostic_hint = (
                "\n\n【重要】本次设计存在约束违反问题，经过多次迭代仍未解决。"
                "你必须在报告中用专门的章节详细说明："
                "(1) 哪些物理约束被违反 (2) 尝试了哪些修正 (3) 为什么无法解决 (4) 建议用户如何调整"
            )

        user_message = (
            f"用户请求：{user_request}\n\n"
            f"执行上下文（含工具结果）：\n{context_str}\n\n"
            f"请根据以上信息生成完整的报告。{diagnostic_hint}"
        )

        try:
            response = llm.invoke(
                [SystemMessage(content=REPORTER_SYSTEM_PROMPT), HumanMessage(content=user_message)],
                config={"timeout": 60},
            )
            text = re.sub(r'<think>.*?</think>', '', response.content, flags=re.DOTALL).strip()
            # 空响应检测：如果 LLM 返回空内容，回退到结构化报告
            if len(text) < 50:
                return None
            return text + f"\n\n---\n*报告生成时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*"
        except Exception:
            return None

    def _normalize_results(self, results: dict) -> dict:
        if not results or not isinstance(results, dict):
            return {}
        normalized = {}
        for tool_name, tool_result in results.items():
            if not isinstance(tool_result, dict):
                normalized[tool_name] = tool_result
            elif tool_result.get("success") and tool_result.get("result") is not None:
                normalized[tool_name] = tool_result["result"]
            elif tool_result.get("error"):
                normalized[tool_name] = {"error": tool_result["error"]}
            else:
                normalized[tool_name] = tool_result
        return normalized

    # ------------------------------------------------------------------ #
    # 结构化回退报告
    # ------------------------------------------------------------------ #

    def _report_structured(self, state: AgentState) -> str:
        plan          = state.get("plan", {})
        raw_results   = state.get("execution_results") or {}
        validation    = state.get("validation_results") or {}
        errors        = state.get("errors") or []
        understanding = plan.get("understanding", "")
        iteration     = state.get("iteration", 0)

        lines = ["=" * 70, "  任务执行报告", "=" * 70, ""]

        # ── 约束违反诊断摘要（置顶，最优先展示）──
        if iteration > 0:
            lines.append(f"【⚠️ 约束违反诊断（迭代 {iteration} 次）】")
            lines.append(f"  系统经过 {iteration} 次迭代修正，但物理约束仍被反复违反。")
            # 提取所有审核决策消息作为诊断依据
            for m in state.get("messages", []):
                content = m.get("content", "")
                role = m.get("role", "")
                if role == "planner" and "审核决策" in content:
                    lines.append(f"  ▸ {content}")
                elif role == "planner" and "迭代上限" in content:
                    lines.append(f"  ▸ {content}")
                elif role == "system" and "人工介入" in content:
                    lines.append(f"  ▸ {content}")
            # 从执行结果中提取关键违反信息
            for tool_name, tool_result in raw_results.items():
                if not isinstance(tool_result, dict):
                    continue
                actual = tool_result.get("result", tool_result)
                if not isinstance(actual, dict):
                    continue
                # 提取 suggestions 中的约束违反信息
                for s in actual.get("suggestions", []):
                    if isinstance(s, str) and any(kw in s for kw in ["违反", "超标", "不可行", "失败", "过大", "过高"]):
                        lines.append(f"  ▸ [{tool_name}] {s}")
            lines.append("")
            lines.append("【建议】")
            lines.append("  1. 检查输入参数是否在工程合理范围内")
            lines.append("  2. 考虑调整设计规格（如纯度要求、流量等）")
            lines.append("  3. 考虑更换设备类型或修改工艺路线")
            lines.append("")

        if state.get("user_request"):
            lines += ["【用户请求】", f"  {state['user_request']}", ""]
        if understanding:
            lines += ["【任务理解】", f"  {understanding}", ""]
        if plan.get("reasoning"):
            lines += ["【工程分析推理】", f"  {plan['reasoning']}", ""]
        if plan.get("notes"):
            lines.append("【备注与假设】")
            for note in plan["notes"]:
                lines.append(f"  • {note}")
            lines.append("")
        if plan.get("missing_params"):
            lines.append("【未提供的参数】")
            for p in plan["missing_params"]:
                lines.append(f"  ? {p}")
            lines.append("")

        normalized = self._normalize_results(raw_results)
        if normalized:
            lines.append("【工具执行结果】")
            # 按任务顺序输出每个设备的结果
            tasks = plan.get("tasks", [])
            reported_keys = set()
            for task in tasks:
                tn = task.get("tool", "")
                purpose = task.get("purpose", "")
                # 找到对应的结果 key
                result_key = None
                for k in normalized:
                    if k == tn or k.startswith(f"{tn}#"):
                        if k not in reported_keys:
                            result_key = k
                            break
                if result_key:
                    reported_keys.add(result_key)
                    lines.append(f"  ▸ {purpose or tn} ({result_key})")
                    self._format_value(normalized[result_key], lines, indent=4)
            # 输出未被任务覆盖的结果
            for tool_name, result in normalized.items():
                if tool_name not in reported_keys:
                    lines.append(f"  ▸ {tool_name}")
                    self._format_value(result, lines, indent=4)
            lines.append("")

        # ── 缺失参数与默认值提示 ──
        defaulted_all = {}
        suggestions_all = []
        notes_missing_all = []
        user_defaulted_reports = []
        for tool_name, tool_result in raw_results.items():
            if not isinstance(tool_result, dict):
                continue
            # 从 result 子对象或顶层提取
            actual = tool_result.get("result", tool_result)
            if not isinstance(actual, dict):
                continue
            if actual.get("user_defaulted_params_report"):
                user_defaulted_reports.append(actual["user_defaulted_params_report"])
            if actual.get("defaulted_params"):
                defaulted_all.update(actual["defaulted_params"])
            if actual.get("suggestions"):
                suggestions_all.extend(actual["suggestions"])
            if actual.get("notes_missing_params"):
                notes_missing_all.extend(actual["notes_missing_params"])

        if notes_missing_all:
            lines.append("【自动推算参数】")
            for note in notes_missing_all:
                lines.append(f"  ▸ {note}")
            lines.append("")

        if user_defaulted_reports:
            lines.append("【ℹ️ 使用工程默认值的参数（提醒）】")
            for report in user_defaulted_reports:
                for line in report.split("\n"):
                    lines.append(f"  {line}")
            lines.append("")
        elif defaulted_all:
            lines.append("【使用默认值的参数】")
            for pname, pval in defaulted_all.items():
                lines.append(f"  ⚠ {pname}: 当前默认值 {pval}（提供实际值后可完成更多计算）")
            lines.append("")

        if suggestions_all and iteration == 0:
            # 仅在无约束违反时输出常规建议；有违反时已在诊断摘要中展示
            lines.append("【工程建议】")
            for s in suggestions_all:
                if isinstance(s, str) and "\n" in s:
                    lines.append(f"  {s}")
                else:
                    lines.append(f"  • {s}")
            lines.append("")

        if validation:
            lines.append("【校验信息】")
            self._format_value(validation, lines, indent=2)
            lines.append("")

        if errors:
            lines.append("【错误信息】")
            for e in errors:
                lines.append(f"  ✗ {e}")
            lines.append("")

        # ── 完整执行日志（包含所有智能体消息）──
        agent_messages = [
            m for m in state.get("messages", [])
            if m.get("role") in ("planner", "executor", "system")
        ]
        if agent_messages:
            lines.append("【执行日志】")
            for m in agent_messages:
                lines.append(f"  [{m.get('role', '?')}] {m.get('content', '')}")
            lines.append("")

        lines += [
            "=" * 70,
            f"  报告时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "=" * 70,
        ]
        return "\n".join(lines)

    def _format_value(self, value, lines: list, indent: int = 0):
        prefix = " " * indent
        if isinstance(value, dict):
            for k, v in value.items():
                if isinstance(v, dict):
                    lines.append(f"{prefix}{k}:")
                    self._format_value(v, lines, indent + 2)
                elif isinstance(v, (list, tuple)):
                    lines.append(f"{prefix}{k}:")
                    for item in v:
                        if isinstance(item, (dict, list)):
                            self._format_value(item, lines, indent + 2)
                        else:
                            lines.append(f"{prefix}  - {item}")
                else:
                    lines.append(f"{prefix}{k}: {v}")
        elif isinstance(value, (list, tuple)):
            for item in value:
                self._format_value(item, lines, indent)
        else:
            lines.append(f"{prefix}{value}")