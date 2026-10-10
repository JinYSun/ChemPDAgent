"""减压阀工具函数。

设计原则：
- 每个 public calc_* 只接收原始题目参数
- public calc_* 之间禁止互相调用
- 每个函数只返回自己的最终结果
"""

import os

_VERBOSE = os.getenv("FUNCS_VERBOSE", "0") == "1"


def _vprint(*args, **kwargs):
    if _VERBOSE:
        print(*args, **kwargs)


def calc_outlet_pressure(
    inlet_pressure_Pa: float,
    target_pressure_Pa: float,
) -> dict:
    """计算减压阀出口压力（单位 Pa）。

    减压阀通过节流作用，将流体从高压降低到目标压力。本函数记录流程中的压力变化：
    验证目标压力必须低于入口压力（减压阀只能降压），返回出口压力和压降。

    参数:
        inlet_pressure_Pa  — 减压阀入口压力 (Pa)
        target_pressure_Pa — 目标出口压力 (Pa)

    返回:
        {"outlet_pressure_Pa": 出口压力 (Pa),
         "pressure_drop_Pa": 压降 ΔP (Pa)}
    """
    if target_pressure_Pa >= inlet_pressure_Pa:
        raise ValueError(
            f"目标压力 {target_pressure_Pa:.0f} Pa 必须小于入口压力 {inlet_pressure_Pa:.0f} Pa，"
            f"减压阀只能降压。"
        )
    pressure_drop_Pa = inlet_pressure_Pa - target_pressure_Pa
    _vprint(f"\n【calc_outlet_pressure (减压阀)】")
    _vprint(f"  入口压力: {inlet_pressure_Pa:.0f} Pa ({inlet_pressure_Pa/1e6:.4f} MPa)")
    _vprint(f"  出口压力: {target_pressure_Pa:.0f} Pa ({target_pressure_Pa/1e6:.4f} MPa)")
    _vprint(f"  压降 ΔP: {pressure_drop_Pa:.0f} Pa ({pressure_drop_Pa/1e6:.4f} MPa)")
    return {
        "outlet_pressure_Pa": target_pressure_Pa,
        "pressure_drop_Pa": pressure_drop_Pa,
    }
