"""LangGraph 多智能体状态定义 — 通用版本"""
from __future__ import annotations
from typing import TypedDict, Annotated, List, Dict, Any
import operator

from langgraph.graph import StateGraph  # noqa: F401  — 确认 langgraph 已安装


class AgentState(TypedDict):
    messages:          Annotated[list, operator.add]
    user_request:      str
    plan:              Dict[str, Any]
    execution_results: Dict[str, Any]
    validation_results: Dict[str, Any]
    final_report:      str
    current_phase:     str
    errors:            List[str]
    iteration:         int
    # 人机交互字段
    human_prompt:      str          # 向用户提问的内容，非空表示需要人工介入
    human_input:       str          # 用户补充的输入，注入后重新规划
    human_feedback_count: int       # 用户补充信息的轮次计数（防无限循环）
    plan_review_count:    int       # 计划审核轮次计数（防无限循环）
    # 内部注入字段（不参与 LLM 交互，仅用于节点间传递）
    _registry:           Any               # ToolRegistry 实例
    _llm:                Any               # LangChain LLM 实例
    _human_supplements:  List[str]         # 累积的用户补充信息列表


def create_initial_state(user_request: str) -> AgentState:
    return {
        "messages":          [],
        "user_request":      user_request,
        "plan":              {},
        "execution_results": {},
        "validation_results": {},
        "final_report":      "",
        "current_phase":     "analyzing",
        "errors":            [],
        "iteration":         0,
        "human_prompt":      "",
        "human_input":       "",
        "human_feedback_count": 0,
        "plan_review_count":  0,
        "_registry":         None,
        "_llm":              None,
        "_human_supplements": [],
    }