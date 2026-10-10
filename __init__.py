"""MCP Agents — 基于 LangGraph 的化工设备设计多智能体系统

完全独立的库，所有物理计算代码内嵌在 physics_engine 中，无外部路径依赖。

通用架构:
  用户输入 → Planner(LLM分析) → Executor(工具执行) → Planner(LLM检查) → Reporter(报告生成)

特点:
  1. LLM 驱动 — 所有决策由 LLM 做出，不硬编码规则
  2. 动态工具 — 支持任意工具注册
  3. 零配置扩展 — 添加工具无需修改智能体代码

使用:
  from mcp_agents import create_universal_workflow
  workflow = create_universal_workflow()
  result = workflow("设计一台水泵，流量36m³/h，扬程30m")
"""
__all__ = ["registry", "run_design", "build_registry", "create_universal_workflow"]

from server import registry, build_registry, run_design


def create_universal_workflow(llm=None):
    """创建通用工作流 — 自动注册所有工具

    Args:
        llm: 可选的 LLM 实例

    Returns:
        可调用的工作流函数
    """
    from tools.registry import UniversalToolRegistry
    import device_tools
    from graph.workflow import create_workflow

    registry = UniversalToolRegistry()
    device_tools.register_all_devices(registry)

    return create_workflow(registry, llm=llm)
