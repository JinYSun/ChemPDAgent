"""LangGraph 多智能体工作流编排 — 通用版本"""
from __future__ import annotations
import traceback
from typing import Dict, Any, Optional

from langgraph.graph import StateGraph, END

from .state import AgentState, create_initial_state
from agents.planner import PlannerAgent, get_llm, set_llm
from agents.executor import ExecutorAgent
from agents.reporter import ReporterAgent


def _ensure_reporter(state: AgentState) -> AgentState:
    """确保所有结果（包括错误）都经过 Reporter 生成最终报告。
    如果 final_report 为空，手动调用 ReporterAgent 生成报告。"""
    if state.get("final_report"):
        return state
    # 将错误信息记录到 state，Reporter 会读取
    if not state.get("errors"):
        state["errors"] = []
    try:
        llm = state.get("_llm") or get_llm()
        state = ReporterAgent(llm=llm).report(state)
    except Exception as e:
        # Reporter 本身也失败时，生成最简报告
        errors = state.get("errors") or []
        state["final_report"] = (
            f"报告生成失败: {e}\n\n"
            f"错误信息: {', '.join(str(err) for err in errors) if errors else '无'}"
        )
        state["current_phase"] = "done"
    return state

# ============================================================
# LangGraph 节点函数
# ============================================================

def _planning_node(state: AgentState) -> AgentState:
    phase = state.get("current_phase", "analyzing")
    phase_labels = {
        "analyzing": "正在分析任务，规划执行步骤...",
        "human_feedback": "正在处理用户补充信息，重新规划...",
        "checking": "正在审核执行结果...",
        "plan_reviewing": "正在审核执行计划...",
    }
    print(f"  [规划] {phase_labels.get(phase, f'规划中({phase})...')}")
    llm = state.get("_llm") or get_llm()
    try:
        return PlannerAgent(llm=llm).plan(state)
    except Exception as e:
        tb = traceback.format_exc()
        print(f"  [规划错误] {e}\n{tb}")
        state["errors"] = state.get("errors", []) + [f"规划器错误: {e}"]
        state["current_phase"] = "reporting"
        return state


def _execution_node(state: AgentState) -> AgentState:
    registry = state.get("_registry")
    if registry is None:
        state["errors"] = state.get("errors", []) + ["工具注册中心未注入到 state"]
        state["current_phase"] = "checking"
        return state

    # 显示当前执行进度
    plan = state.get("plan", {})
    tasks = plan.get("tasks", [])
    idx = plan.get("current_task_index", 0)
    if idx < len(tasks):
        tool_name = tasks[idx].get("tool", "")
        purpose = tasks[idx].get("purpose", "")
        print(f"  [执行] 正在调用工具: {tool_name} ({idx+1}/{len(tasks)}) — {purpose}")

    try:
        result_state = ExecutorAgent(registry).execute(state)
    except Exception as e:
        tb = traceback.format_exc()
        print(f"  [执行错误] {e}\n{tb}")
        result_state = state
        result_state["errors"] = result_state.get("errors", []) + [f"执行器错误: {e}"]
        result_state["current_phase"] = "checking"
        return result_state
    # 检测执行结果中是否有 needs_human 请求（如 human_ask 工具）
    exec_results = result_state.get("execution_results") or {}
    for tool_name, tool_result in exec_results.items():
        if isinstance(tool_result, dict):
            inner = tool_result.get("result", tool_result)
            if isinstance(inner, dict) and inner.get("needs_human"):
                result_state["human_prompt"] = inner.get("human_prompt", "请补充缺失的参数")
                result_state["current_phase"] = "human_feedback"
                result_state["messages"].append({
                    "role": "executor",
                    "content": f"工具 {tool_name} 请求人工输入：{inner.get('human_prompt', '')}",
                })
                return result_state
    return result_state


def _reporting_node(state: AgentState) -> AgentState:
    print("  [报告] 正在生成设计报告...")
    llm = state.get("_llm") or get_llm()
    return ReporterAgent(llm=llm).report(state)


def _human_feedback_node(state: AgentState) -> AgentState:
    """人机交互节点 — 暂停执行，等待用户补充必填参数。"""
    prompt = state.get('human_prompt', '请补充缺失参数')
    state["messages"].append({
        "role": "system",
        "content": f"[等待人工输入] {prompt}",
    })
    # human_feedback 阶段由外部交互循环处理，此处仅占位
    return state


# ============================================================
# 条件路由
# ============================================================

def _route_after_planning(state: AgentState) -> str:
    phase = state.get("current_phase", "analyzing")
    if phase == "executing":         return "execution"
    if phase == "human_feedback":    return "human_feedback"
    if phase == "plan_reviewing":    return "planning"  # 回到 planning 节点处理 plan_reviewing 阶段
    if phase == "checking":          return "planning"
    if phase == "reporting":         return "reporting"
    return "execution"


def _route_after_execution(state: AgentState) -> str:
    phase = state.get("current_phase", "checking")
    if phase == "human_feedback": return "human_feedback"
    if phase == "checking":       return "planning"
    if phase == "executing":      return "execution"
    return "reporting"


# ============================================================
# 工作流构建
# ============================================================

def create_workflow(registry, llm=None, interactive: bool = True):
    """构建并返回可调用的多智能体工作流

    Args:
        registry: 工具注册中心（需实现 execute / get_tools_prompt 等接口）
        llm:      可选的 LangChain LLM 实例（不传则使用环境变量配置）
        interactive: 是否启用人工交互循环（默认 True，测试时设 False）

    Returns:
        run(request, params=None) -> AgentState
    """
    workflow = StateGraph(AgentState)

    workflow.add_node("planning",  _planning_node)
    workflow.add_node("execution", _execution_node)
    workflow.add_node("reporting", _reporting_node)
    workflow.add_node("human_feedback", _human_feedback_node)

    workflow.set_entry_point("planning")

    workflow.add_conditional_edges(
        "planning",
        _route_after_planning,
        {"execution": "execution", "planning": "planning", "reporting": "reporting", "human_feedback": "human_feedback"},
    )
    workflow.add_conditional_edges(
        "execution",
        _route_after_execution,
        {"planning": "planning", "execution": "execution", "reporting": "reporting", "human_feedback": "human_feedback"},
    )
    # human_feedback 节点结束后终止本次执行。
    # 交互循环（run 函数）会检测 human_feedback 状态并提示用户，
    # 用户输入后重新调用 app.invoke，通过 Planner._handle_human_feedback 合并输入。
    workflow.add_edge("human_feedback", END)
    workflow.add_edge("reporting", END)

    app = workflow.compile()

    def run(request: str, params: Dict[str, Any] = None) -> AgentState:
        initial = create_initial_state(request)
        initial["_registry"] = registry
        initial["_llm"] = llm or get_llm()
        if params:
            initial["plan"]["parameters"] = params

        print("\n正在执行多智能体工作流...\n")
        try:
            result = app.invoke(initial)
        except Exception as e:
            tb = traceback.format_exc()
            print(f"\n[错误] 工作流执行失败: {e}")
            print(f"[调试] 完整堆栈:\n{tb}")
            initial["errors"] = initial.get("errors", []) + [f"工作流执行失败: {e}\n堆栈: {tb}"]
            initial["current_phase"] = "reporting"
            return _ensure_reporter(initial)

        # 检测是否需要人工介入（仅在交互模式下暂停）
        MAX_HUMAN_FEEDBACK_ROUNDS = 3
        if interactive:
            while result.get("current_phase") == "human_feedback" and result.get("human_prompt"):
                feedback_count = result.get("human_feedback_count", 0)
                if feedback_count >= MAX_HUMAN_FEEDBACK_ROUNDS:
                    print(f"\n[提示] 已达到最大补充轮次 ({MAX_HUMAN_FEEDBACK_ROUNDS})，将使用已有参数继续执行。")
                    # 强制跳过 human_feedback，让 planner 直接执行
                    result["current_phase"] = "analyzing"  # 触发重新规划
                    result["human_prompt"] = ""  # 清除旧 prompt
                    try:
                        result = app.invoke(result)
                    except Exception as e:
                        result["errors"] = result.get("errors", []) + [f"强制继续执行失败: {e}"]
                        result["current_phase"] = "reporting"
                        print(f"\n[错误] 强制继续执行失败: {e}")
                        return _ensure_reporter(result)
                    break

                print("\n" + "=" * 60)
                print(result["human_prompt"])
                print("=" * 60)
                user_reply = input("补充信息 >>> ").strip()
                if user_reply:
                    result["human_input"] = user_reply
                    result["current_phase"] = "human_feedback"  # 触发 Planner 重新处理
                    result["human_prompt"] = ""  # 清除旧 prompt，避免重复显示
                    print("\n正在重新规划并执行，请稍候...\n")
                    try:
                        result = app.invoke(result)
                    except Exception as e:
                        result["errors"] = result.get("errors", []) + [f"重新规划执行失败: {e}"]
                        result["current_phase"] = "reporting"
                        print(f"\n[错误] 重新规划后执行失败: {e}")
                        return _ensure_reporter(result)
                else:
                    # 用户未输入，强制继续执行
                    print("\n[提示] 未提供补充信息，将使用已有参数继续执行。")
                    result["current_phase"] = "analyzing"
                    result["human_prompt"] = ""
                    try:
                        result = app.invoke(result)
                    except Exception as e:
                        result["errors"] = result.get("errors", []) + [f"强制继续执行失败: {e}"]
                        result["current_phase"] = "reporting"
                        print(f"\n[错误] 强制继续执行失败: {e}")
                        return _ensure_reporter(result)
                    break

        # 确保所有结果都经过 Reporter
        return _ensure_reporter(result)

    return run


# ============================================================
# 便捷入口
# ============================================================

def run_task(request: str, registry, params: Dict[str, Any] = None, llm=None, interactive: bool = True) -> str:
    """运行完整多智能体流程，返回最终报告文本"""
    workflow = create_workflow(registry, llm=llm, interactive=interactive)
    result = workflow(request, params)
    # final_report 已由 _ensure_reporter 保证生成
    report = result.get("final_report", "")
    if not report:
        # 极端情况：Reporter 本身也失败
        errors = result.get("errors", [])
        report = f"任务执行未能生成最终报告。错误信息: {errors}" if errors else "任务执行未能生成最终报告。"
    return report


def create_universal_workflow(registry=None, tools_module=None, llm=None):
    """创建通用工作流，自动注册工具模块中的所有工具

    Args:
        registry:     已配置好的工具注册中心（优先使用）
        tools_module: 工具模块（需提供 register_all_tools 或 register_all_devices 方法）
        llm:          可选的 LLM 实例
    """
    if registry is None:
        from tools.registry import UniversalToolRegistry
        registry = UniversalToolRegistry()

    if tools_module is not None:
        if hasattr(tools_module, "register_all_tools"):
            tools_module.register_all_tools(registry)
        elif hasattr(tools_module, "register_all_devices"):
            tools_module.register_all_devices(registry)

    return create_workflow(registry, llm=llm)