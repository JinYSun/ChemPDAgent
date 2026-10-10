"""
调节阀流量系数计算与口径选型示例
==================================

利用 physics_engine.valve 物理引擎，完成：
  1. 水的物性计算 (密度、蒸气压)
  2. 流量系数 Kv/Cv 计算 (IEC 60534-2-1)
  3. 阻塞流判别与空化校核
  4. 阀门口径推荐

已知条件:
  最大流量 Q = 15 m3/h
  阀前后压差 ΔP = 150 kPa
  介质: 水
  操作压力: 1 bar (表压), 即阀前约 2 bar abs
  温度: 20°C (常温)
"""

import sys
import os
import io
import math

if sys.stdout.encoding and sys.stdout.encoding.lower() in ('gbk', 'cp936'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

_project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from physics_engine.valve import size_control_valve_liquid
from physics_engine.thermo_helper import (
    get_fluid_density,
    get_fluid_vapor_pressure,
)


def sep(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def main():
    P_atm = 101325.0  # 标准大气压 Pa

    # ================================================================
    # 1. 已知条件
    # ================================================================
    sep("1. 已知条件")

    Q = 15.0                     # 最大流量 (m3/h)
    delta_P = 150e3              # 阀前后压差 (Pa) = 150 kPa
    T_C = 20.0                   # 操作温度 (°C)
    T_K = T_C + 273.15           # 操作温度 (K)
    P_gauge = 1e5                # 操作表压 1 bar = 100 kPa

    # 阀前绝对压力 = 大气压 + 表压
    P1 = P_atm + P_gauge         # = 201325 Pa ≈ 2 bar abs
    # 阀后绝对压力 = 阀前 - 压差
    P2 = P1 - delta_P            # = 51325 Pa ≈ 0.5 bar abs

    print(f"  工艺参数:")
    print(f"    最大流量 Q        = {Q} m3/h")
    print(f"    阀前后压差 ΔP     = {delta_P/1000:.0f} kPa = {delta_P} Pa")
    print(f"    操作温度          = {T_C} °C = {T_K:.2f} K")
    print(f"    操作压力 (表压)   = {P_gauge/1e5:.0f} bar = {P_gauge:.0f} Pa")
    print(f"")
    print(f"  压力分析:")
    print(f"    阀前绝对压力 P1   = P_atm + P_gauge = {P_atm} + {P_gauge} = {P1:.0f} Pa")
    print(f"    阀后绝对压力 P2   = P1 - ΔP = {P1:.0f} - {delta_P:.0f} = {P2:.0f} Pa")
    print(f"    P1 = {P1/1e5:.3f} bar abs, P2 = {P2/1e5:.3f} bar abs")

    # ================================================================
    # 2. 水的物性 (引擎计算)
    # ================================================================
    sep("2. 水的物性计算 (thermo_helper)")

    rho = get_fluid_density("water", T_K, P1)
    Pv = get_fluid_vapor_pressure("water", T_K)
    Pc = 22.064e6  # 水的临界压力 (Pa)

    print(f"\n  调用 get_fluid_density('water', {T_K:.2f}K, {P1:.0f}Pa):")
    print(f"    ρ = {rho:.2f} kg/m3")
    print(f"\n  调用 get_fluid_vapor_pressure('water', {T_K:.2f}K):")
    print(f"    Pv = {Pv:.0f} Pa = {Pv/1000:.3f} kPa")
    print(f"\n  水的临界压力 Pc = {Pc/1e6:.3f} MPa")

    SG = rho / 1000.0
    print(f"\n  比重 SG = ρ/1000 = {SG:.4f}")

    # ================================================================
    # 3. 流量系数计算 (IEC 60534-2-1)
    # ================================================================
    sep("3. 流量系数 Kv/Cv 计算 (IEC 60534-2-1)")

    # 3.1 手动计算展示过程
    print(f"\n  3.1 基本公式:")
    print(f"    Kv = Q × √(SG / ΔP_bar)")
    print(f"    Cv = 1.156 × Kv")
    print(f"")

    delta_P_bar = delta_P / 1e5
    Kv_basic = Q * math.sqrt(SG / delta_P_bar)
    Cv_basic = 1.156 * Kv_basic

    print(f"  3.2 代入数值 (非阻塞流):")
    print(f"    ΔP = {delta_P_bar:.2f} bar")
    print(f"    Kv = {Q} × √({SG:.4f} / {delta_P_bar:.2f})")
    print(f"       = {Q} × {math.sqrt(SG / delta_P_bar):.4f}")
    print(f"       = {Kv_basic:.2f}")
    print(f"    Cv = 1.156 × {Kv_basic:.2f} = {Cv_basic:.2f}")

    # 3.3 引擎计算 (多种阀型对比)
    print(f"\n  3.3 引擎计算 (size_control_valve_liquid) — 多阀型对比:")
    valve_types = ["globe", "ball", "butterfly"]
    valve_names = {"globe": "截止阀", "ball": "球阀", "butterfly": "蝶阀"}

    print(f"\n  {'阀型':>10s} | {'Kv':>6s} | {'Cv':>6s} | {'Kv_req(×1.3)':>12s} | {'Cv_req(×1.3)':>12s} | {'流态':>10s} | {'空化':>6s}")
    print(f"  {'-'*10}-+-{'-'*6}-+-{'-'*6}-+-{'-'*12}-+-{'-'*12}-+-{'-'*10}-+-{'-'*6}")

    results = {}
    for vt in valve_types:
        res = size_control_valve_liquid(
            Q=Q, rho=rho, P1=P1, P2=P2, Pv=Pv, Pc=Pc,
            valve_type=vt
        )
        results[vt] = res
        print(f"  {valve_names[vt]:>10s} | {res['Kv']:6.2f} | {res['Cv']:6.2f} | {res['Kv_required']:12.2f} | {res['Cv_required']:12.2f} | {res['flow_condition']:>10s} | {res['cavitation_risk']:>6s}")

    # ================================================================
    # 4. 阻塞流与空化分析 (以截止阀为例)
    # ================================================================
    sep("4. 阻塞流与空化分析 (截止阀)")

    res_globe = results["globe"]

    print(f"\n  压力恢复系数:")
    print(f"    FL (截止阀) = {res_globe['FL']:.2f}")
    print(f"    FF (饱和比修正) = {res_globe['FF']:.4f}")
    print(f"")
    print(f"  阻塞流判别:")
    print(f"    ΔP_choked = FL² × (P1 - FF × Pv)")
    delta_P_choked = res_globe['FL']**2 * (P1 - res_globe['FF'] * Pv)
    print(f"              = {res_globe['FL']:.2f}² × ({P1:.0f} - {res_globe['FF']:.4f} × {Pv:.0f})")
    print(f"              = {delta_P_choked:.0f} Pa = {delta_P_choked/1000:.1f} kPa")
    print(f"    实际 ΔP    = {delta_P:.0f} Pa = {delta_P/1000:.1f} kPa")
    if delta_P > delta_P_choked:
        print(f"    结论: ΔP > ΔP_choked → 阻塞流!")
    else:
        print(f"    结论: ΔP < ΔP_choked → 非阻塞流")

    print(f"\n  空化校核:")
    print(f"    σ_actual = (P2 - Pv) / ΔP = ({P2:.0f} - {Pv:.0f}) / {delta_P:.0f} = {res_globe['sigma_actual']:.4f}")
    print(f"    σ_allowable (截止阀) = {res_globe['sigma_allowable']:.2f}")
    if res_globe['sigma_actual'] < res_globe['sigma_allowable']:
        print(f"    结论: σ_actual < σ_allowable → 高空化风险!")
    elif res_globe['sigma_actual'] < res_globe['sigma_allowable'] * 1.5:
        print(f"    结论: σ_actual 接近 σ_allowable → 中等空化风险")
    else:
        print(f"    结论: σ_actual > 1.5×σ_allowable → 空化风险低")

    # 闪蒸判别
    print(f"\n  闪蒸判别:")
    print(f"    P2 = {P2:.0f} Pa, Pv = {Pv:.0f} Pa")
    if P2 < Pv:
        print(f"    P2 < Pv → 阀后发生闪蒸!")
    else:
        print(f"    P2 > Pv → 无闪蒸")

    # ================================================================
    # 5. 阀门口径推荐
    # ================================================================
    sep("5. 阀门口径推荐")

    # 标准调节阀口径系列 (DN / 英寸 / 典型 Kv 范围)
    # 基于等百分比阀芯, 全开时 Kv
    valve_sizes = [
        # DN, inch_str, Kv_range (全行程)
        (15,  '1/2',   1.0,  6.3),
        (20,  '3/4',   2.0,  12),
        (25,  '1',     3.5,  20),
        (32,  '1 1/4', 6.0,  32),
        (40,  '1 1/2', 9.0,  50),
        (50,  '2',     15,   80),
        (65,  '2 1/2', 25,   125),
        (80,  '3',     40,   200),
        (100, '4',     65,   320),
        (125, '5',     100,  500),
        (150, '6',     160,  800),
        (200, '8',     280,  1400),
    ]

    Kv_req = res_globe['Kv_required']  # 含30%裕量的Kv
    Cv_req = res_globe['Cv_required']

    print(f"\n  所需流量系数 (含30%裕量):")
    print(f"    Kv_required = {Kv_req:.2f}")
    print(f"    Cv_required = {Cv_req:.2f}")
    print(f"")
    print(f"  调节阀口径选型 (等百分比阀芯, 推荐开度 60~80%):")
    print(f"")
    print(f"  {'DN':>4s} | {'inch':>5s} | {'Kv_max':>7s} | {'Kv/Kv_req':>10s} | {'推荐开度':>8s} | {'评价':>8s}")
    print(f"  {'-'*4}-+-{'-'*5}-+-{'-'*7}-+-{'-'*10}-+-{'-'*8}-+-{'-'*8}")

    best_dn = None
    for dn, inch, kv_min, kv_max in valve_sizes:
        if kv_max >= Kv_req:
            # 推荐工作开度: 假设等百分比特性 R=30
            # Kv_actual = Kv_max * R^(x-1), where x = travel fraction
            # For x = 0.6~0.8: Kv_actual = Kv_max * 30^(x-1)
            # We want Kv_actual >= Kv_basic (不含裕量)
            # Solve: x = 1 + ln(Kv_basic/Kv_max) / ln(30)
            if kv_max > 0 and Kv_basic <= kv_max:
                R = 30  # 等百分比范围度
                x_ideal = 1 + math.log(Kv_basic / kv_max) / math.log(R)
                open_pct = f"{x_ideal*100:.0f}%"
            else:
                open_pct = "N/A"

            ratio = kv_max / Kv_req
            if 1.0 <= ratio <= 2.0:
                eval_str = "合适"
                if best_dn is None:
                    best_dn = dn
                    best_inch = inch
                    best_kv_max = kv_max
            elif ratio > 2.0:
                eval_str = "偏大"
            else:
                eval_str = "偏小"

            print(f"  DN{dn:<3d} | {inch:>5s} | {kv_max:>7.1f} | {ratio:>10.2f} | {open_pct:>8s} | {eval_str:>8s}")

    # ================================================================
    # 6. 设计结论
    # ================================================================
    sep("6. 设计结论")

    print(f"""
  工艺条件:
    最大流量 Q     = {Q} m3/h
    阀前后压差 ΔP  = {delta_P/1000:.0f} kPa
    介质           = 水 ({T_C}°C)
    阀前压力 P1    = {P1/1e5:.3f} bar abs
    阀后压力 P2    = {P2/1e5:.3f} bar abs

  物性参数:
    密度 ρ         = {rho:.2f} kg/m3
    蒸气压 Pv      = {Pv:.0f} Pa
    比重 SG        = {SG:.4f}

  流量系数:
    Kv (计算值)    = {Kv_basic:.2f}
    Cv (计算值)    = {Cv_basic:.2f}
    Kv (含30%裕量) = {Kv_req:.2f}
    Cv (含30%裕量) = {Cv_req:.2f}

  流态分析 (截止阀):
    阻塞流判别     = {res_globe['flow_condition']}
    空化风险       = {res_globe['cavitation_risk']}
    σ_actual       = {res_globe['sigma_actual']:.4f}

  选型建议:""")

    if best_dn:
        print(f"    推荐口径: DN{best_dn} (NPS {best_inch})")
        print(f"    阀门类型: 截止阀 (globe valve)")
        print(f"    Kv_max   = {best_kv_max:.1f}")
        print(f"    阀门口径: DN{best_dn} 等百分比特性调节阀")
        print(f"    额定行程下 Kv = {best_kv_max:.1f} > Kv_required = {Kv_req:.2f} ✓")
    else:
        print(f"    超出标准系列范围，需定制")

    print(f"\n{'='*70}")
    print(f"  计算完成!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
