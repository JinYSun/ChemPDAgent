"""设备级 MCP 工具集"""
from device_tools.distillation_device import register_distillation_device
from device_tools.flash_drum_device import register_flash_drum_device
from device_tools.heatexchanger_device import register_heatexchanger_device
from device_tools.reactor_device import register_reactor_device
from device_tools.pump_device import register_pump_device
from device_tools.storage_tank_device import register_storage_tank_device
from tools.base import ToolRegistry
from tools.registry import UniversalToolRegistry

def register_all_device_tools(registry: ToolRegistry) -> ToolRegistry:
    """注册所有设备工具到传统注册中心"""
    register_distillation_device(registry)
    register_flash_drum_device(registry)
    register_heatexchanger_device(registry)
    register_reactor_device(registry)
    register_pump_device(registry)
    register_storage_tank_device(registry)
    return registry


def register_all_devices(registry: UniversalToolRegistry) -> UniversalToolRegistry:
    """注册所有设备到通用注册中心 — 自动提取能力描述"""
    from . import distillation_device, flash_drum_device, heatexchanger_device
    from . import reactor_device, pump_device, storage_tank_device

    # 注册每个设备模块
    for module in [distillation_device, flash_drum_device, heatexchanger_device,
                   reactor_device, pump_device, storage_tank_device]:
        if hasattr(module, "register_to_universal"):
            module.register_to_universal(registry)

    return registry
