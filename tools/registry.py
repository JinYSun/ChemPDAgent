"""通用工具注册中心"""
from typing import Any, Callable, Optional
from dataclasses import dataclass, field

@dataclass
class ToolCapability:
    """工具能力描述"""
    name: str
    description: str
    parameters: dict = field(default_factory=dict)
    returns: dict = field(default_factory=dict)
    category: str = ""
    tags: list = field(default_factory=list)
    examples: list = field(default_factory=list)
    aliases: list = field(default_factory=list)

class UniversalToolRegistry:
    """通用工具注册中心"""
    def __init__(self):
        self._tools = {}
        self._aliases = {}  # alias -> canonical name
    
    def register(self, capability: ToolCapability, func: Callable):
        self._tools[capability.name] = {"capability": capability, "func": func}
        for alias in capability.aliases:
            self._aliases[alias] = capability.name
        return self
    
    def get_tool(self, name: str) -> Optional[Callable]:
        # 1. 精确匹配（工具名）
        tool = self._tools.get(name)
        if tool:
            return tool["func"]
        # 2. 精确匹配（别名）
        canonical = self._aliases.get(name)
        if canonical:
            tool = self._tools.get(canonical)
            if tool:
                return tool["func"]
        # 3. 大小写不敏感匹配（工具名 + 别名）
        name_lower = name.lower()
        for key, entry in self._tools.items():
            if key.lower() == name_lower:
                return entry["func"]
        for alias, canonical in self._aliases.items():
            if alias.lower() == name_lower:
                tool = self._tools.get(canonical)
                if tool:
                    return tool["func"]
        return None

    @property
    def tool_count(self) -> int:
        return len(self._tools)

    @property
    def categories(self) -> list:
        cats = set()
        for entry in self._tools.values():
            cat = entry["capability"].category
            if cat:
                cats.add(cat)
        return sorted(cats)

    def list_tools(self) -> list:
        """返回所有已注册工具的 ToolCapability 列表"""
        return [entry["capability"] for entry in self._tools.values()]

    def get_tools_prompt(self) -> str:
        """生成供 LLM Planner 阅读的工具列表文本"""
        if not self._tools:
            return "（无可用工具）"
        
        lines = []
        # 按类别分组
        by_category = {}
        for entry in self._tools.values():
            cap = entry["capability"]
            cat = cap.category or "other"
            by_category.setdefault(cat, []).append(cap)
        
        for cat in sorted(by_category.keys()):
            lines.append(f"\n--- 类别: {cat} ---")
            for cap in by_category[cat]:
                lines.append(f"工具名: {cap.name}")
                lines.append(f"  描述: {cap.description}")
                if cap.aliases:
                    lines.append(f"  别名: {', '.join(cap.aliases[:5])}")
                if cap.parameters:
                    lines.append("  参数:")
                    for pname, pinfo in cap.parameters.items():
                        req = "必填" if pinfo.get("required") else "可选"
                        desc = pinfo.get("description", "")
                        lines.append(f"    - {pname} ({req}): {desc}")
                lines.append("")
        
        return "\n".join(lines)

    def execute(self, tool_name: str, parameters: dict) -> dict:
        """执行指定工具，返回结果字典
        
        Args:
            tool_name: 工具名称（支持别名）
            parameters: 工具参数字典
            
        Returns:
            {"success": True, "result": ...} 或 {"error": "..."}
        """
        func = self.get_tool(tool_name)
        if func is None:
            return {"error": f"工具 '{tool_name}' 未找到"}
        try:
            # 清理参数键名（去除多余引号）和值（移除内部标记 + 类型转换）
            clean_params = {k: v for k, v in parameters.items() if not k.startswith("_")}
            clean_params = self._sanitize_param_keys(clean_params)
            clean_params = self._convert_string_numbers(clean_params)
            result = func(**clean_params)
            return {"success": True, "result": result}
        except TypeError as e:
            # 参数不匹配时尝试位置参数或简化调用
            try:
                result = func(**clean_params)
                return {"success": True, "result": result}
            except Exception as e2:
                return {"success": False, "error": f"参数错误: {e2}"}
        except Exception as e:
            return {"success": False, "error": str(e)}

    @staticmethod
    def _convert_string_numbers(params: dict) -> dict:
        """递归将字符串类型的数值转换为数字，防止 LLM 输出格式错误。
        
        处理嵌套 dict/list 中的字符串数值（如流量字典 {"Benzene": "3000"}）。
        """
        converted = {}
        for k, v in params.items():
            converted[k] = UniversalToolRegistry._convert_value_recursive(v)
        return converted

    @staticmethod
    def _convert_value_recursive(value) -> any:
        """递归将字符串数值转换为 int/float，处理嵌套 dict/list。"""
        if isinstance(value, dict):
            return {kk: UniversalToolRegistry._convert_value_recursive(vv) for kk, vv in value.items()}
        if isinstance(value, (list, tuple)):
            return [UniversalToolRegistry._convert_value_recursive(item) for item in value]
        if isinstance(value, str):
            try:
                return int(value)
            except ValueError:
                pass
            try:
                return float(value)
            except ValueError:
                pass
        return value

    @staticmethod
    def _sanitize_param_keys(params: dict) -> dict:
        """清理参数键名，去除 LLM 可能输出的多余引号或空白
        
        处理场景：LLM 输出 '"T_hot_in": 400' 导致键名为 '"T_hot_in"'
        """
        sanitized = {}
        for k, v in params.items():
            if isinstance(k, str):
                # 去除首尾的引号（单引号、双引号）和空白
                k = k.strip().strip("'\"").strip()
            sanitized[k] = v
        return sanitized
