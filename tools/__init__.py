"""工具模块"""
from .base import ToolRegistry
from .registry import UniversalToolRegistry, ToolCapability

__all__ = ["ToolRegistry", "UniversalToolRegistry", "ToolCapability", "register_all_tools_to_universal"]


def register_all_tools_to_universal(registry: UniversalToolRegistry) -> UniversalToolRegistry:
    """将所有工具注册到通用注册中心（设备级 + 物理计算级 + 扩展物性）"""
    # 1. 设备级工具 (device_design)
    from device_tools import register_all_devices
    register_all_devices(registry)
    # 2. 物理计算工具 (thermo/pipe/pump/valve/...)
    from .physics_tools import register_physics_tools
    register_physics_tools(registry)
    # 3. 扩展物性计算 + 设备参数计算工具
    from .extended_physics_tools import register_extended_physics_tools
    register_extended_physics_tools(registry)
    return registry
