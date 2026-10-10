"""
通用守恒方程与热力学基础工具
"""
from __future__ import annotations
import numpy as np
from typing import Optional, Dict, Any, Union

# ============================================================
# 类型安全转换 — 防止 LLM 输出字符串导致数学运算报错
# ============================================================

def safe_float(value: Any, default: float = None, name: str = None) -> Optional[float]:
    """安全地将值转换为 float，防止字符串类型导致 TypeError。
    
    Args:
        value: 待转换的值（可以是 str/int/float/None）
        default: 转换失败或值为 None 时的默认值
        name: 参数名（用于错误提示）
    
    Returns:
        float 值，或 default
    """
    if value is None:
        return default
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            if default is not None:
                return default
            raise ValueError(f"参数 '{name}' 无法转换为 float: {value!r}")
    # numpy 类型等
    try:
        return float(value)
    except (TypeError, ValueError):
        if default is not None:
            return default
        raise ValueError(f"参数 '{name}' 类型不支持: {type(value)}")


def safe_int(value: Any, default: int = None, name: str = None) -> Optional[int]:
    """安全地将值转换为 int，防止字符串类型导致 TypeError。"""
    if value is None:
        return default
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError:
            # 尝试先转 float 再转 int（如 "3.0" → 3）
            try:
                return int(float(value))
            except ValueError:
                if default is not None:
                    return default
                raise ValueError(f"参数 '{name}' 无法转换为 int: {value!r}")
    try:
        return int(value)
    except (TypeError, ValueError):
        if default is not None:
            return default
        raise ValueError(f"参数 '{name}' 类型不支持: {type(value)}")


# ============================================================
# 必填参数哨兵 — 取消硬编码默认值，缺失时触发 human-in-the-loop
# ============================================================

class _RequiredParam:
    """必填参数哨兵：当参数为此值时，表示用户必须提供，不可使用默认值"""
    def __repr__(self): return "<REQUIRED>"
    def __bool__(self): return False
    def __eq__(self, other): return isinstance(other, _RequiredParam)

_REQUIRED = _RequiredParam()


def check_required_params(locals_dict: Dict[str, Any], param_descriptions: Dict[str, str]) -> Dict[str, Any]:
    """检查必填参数是否缺失。
    
    用法：在 _full_design 等函数开头调用：
        missing_result = check_required_params(locals(), {
            'F': '进料流量 (kmol/h)',
            'P': '操作压力 (Pa)',
        })
        if missing_result:
            return missing_result
    
    Returns:
        若所有必填参数均已提供 → 返回 None（继续执行）
        若有必填参数缺失 → 返回 {"needs_human": True, "human_prompt": ..., "missing_params": [...]}
    """
    missing = []
    for param_name, desc in param_descriptions.items():
        value = locals_dict.get(param_name, _REQUIRED)
        if isinstance(value, _RequiredParam):
            missing.append(f"  • {param_name}：{desc}")
    
    if not missing:
        return None
    
    return {
        "needs_human": True,
        "human_prompt": (
            "缺少以下必填参数，请补充：\n\n" +
            "\n".join(missing) +
            "\n\n请直接输入补充信息（如：进料流量 100 kmol/h，操作压力 101325 Pa）"
        ),
        "error": "missing_required_params",
        "missing_params": [p for p in param_descriptions if isinstance(locals_dict.get(p, _REQUIRED), _RequiredParam)],
    }


def material_balance(
    in_flows: list[float],
    out_flows: list[float],
    generation: float = 0.0,
    consumption: float = 0.0,
    tolerance: float = 0.01
) -> dict:
    """
    物料平衡校验
    Σ(F_in) + generation = Σ(F_out) + consumption

    Returns:
        {"balanced": bool, "in_total": float, "out_total": float, "deviation": float}
    """
    in_total = sum(in_flows)
    out_total = sum(out_flows)
    net = (in_total + generation) - (out_total + consumption)
    deviation = abs(net) / max(abs(in_total + generation), 1e-10)

    return {
        "balanced": deviation <= tolerance,
        "in_total": in_total,
        "out_total": out_total,
        "net_accumulation": net,
        "deviation": deviation
    }


def energy_balance(
    in_enthalpies: list[float],
    out_enthalpies: list[float],
    Q_in: float = 0.0,
    Q_out: float = 0.0,
    W_in: float = 0.0,
    W_out: float = 0.0,
    tolerance: float = 0.01
) -> dict:
    """
    能量平衡校验
    Σ(H_in) + Q_in + W_in = Σ(H_out) + Q_out + W_out

    Returns:
        {"balanced": bool, "deviation": float, ...}
    """
    lhs = sum(in_enthalpies) + Q_in + W_in
    rhs = sum(out_enthalpies) + Q_out + W_out
    net = lhs - rhs
    deviation = abs(net) / max(abs(lhs), 1e-10)

    return {
        "balanced": deviation <= tolerance,
        "lhs": lhs,
        "rhs": rhs,
        "net_energy": net,
        "deviation": deviation
    }


def clausius_clapeyron(
    T1: float, P1: float, T2: float, P2: float, delta_Hvap: float
) -> dict:
    """
    Clausius-Clapeyron方程
    ln(P2/P1) = (ΔHvap/R) * (1/T1 - 1/T2)

    可用于已知两点蒸气压反推蒸发焓
    """
    R = 8.314  # J/(mol·K)
    if T1 <= 0 or T2 <= 0 or P1 <= 0 or P2 <= 0:
        return {"valid": False, "error": "温度和压力必须为正数"}

    lhs = np.log(P2 / P1)
    rhs = (delta_Hvap / R) * (1.0 / T1 - 1.0 / T2)

    return {
        "valid": True,
        "ln_P_ratio": lhs,
        "calculated_rhs": rhs,
        "consistent": abs(lhs - rhs) / max(abs(lhs), 1e-10) < 0.05
    }


def antoine_equation(
    A: float, B: float, C: float, T: float,
    base: float = 10.0, unit_P: str = "Pa"
) -> float:
    """
    Antoine方程计算蒸气压
    log10(P) = A - B/(T + C)   (P in mmHg)
    或       ln(P) = A - B/(T + C) (自然对数形式)

    Args:
        A, B, C: Antoine常数
        T: 温度 K (或 °C，取决于常数C的量纲)
        base: 对数底数 (10 或 e)
        unit_P: 输出压力单位

    Returns:
        蒸气压 (对应 unit_P 指定的单位，默认为 Pa)
    """
    # 计算出基准的 Pa 压力
    if base == 10:
        logP = A - B / (T + C)
        P_mmHg = 10 ** logP
        P_Pa = P_mmHg * 133.3224  # mmHg → Pa
    else:
        lnP = A - B / (T + C)
        P_Pa = np.exp(lnP)        # 假设 base 为 e 时，常数通常默认基于 Pa 或自然单位

    # 根据指定的输出单位进行转换
    unit_lower = unit_P.lower()
    if unit_lower == "pa":
        return P_Pa
    elif unit_lower == "kpa":
        return P_Pa / 1000.0
    elif unit_lower == "mpa":
        return P_Pa / 1e6
    elif unit_lower == "bar":
        return P_Pa / 100000.0
    elif unit_lower == "atm":
        return P_Pa / 101325.0
    elif unit_lower == "mmhg":
        return P_Pa / 133.3224
    else:
        return P_Pa  # 兜底返回 Pa


def raoults_law(
    x: list[float],
    P_sat: list[float],
    gamma: Optional[list[float]] = None
) -> float:
    """
    Raoult定律 (修正)
    P_total = Σ(xi * γi * P_sat_i)

    Args:
        x: 液相摩尔分数列表
        P_sat: 各组分蒸气压列表 (Pa)
        gamma: 活度系数列表 (None = 理想溶液，全为1)

    Returns:
        总压 Pa
    """
    if gamma is None:
        gamma = [1.0] * len(x)

    P_total = sum(xi * gi * ps for xi, gi, ps in zip(x, gamma, P_sat))
    return P_total