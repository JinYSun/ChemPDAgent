"""工具基础类"""
from typing import Any, Callable, Optional

class ToolRegistry:
    """工具注册中心基类"""
    def __init__(self):
        self._tools = {}
    
    def register_function(self, name: str, func: Callable, description: str = "", **kwargs):
        self._tools[name] = {"func": func, "description": description, **kwargs}
        return self
    
    def get_tool(self, name: str) -> Optional[Callable]:
        tool = self._tools.get(name)
        return tool["func"] if tool else None
