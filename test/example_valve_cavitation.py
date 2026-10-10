"""
调节阀空化与闪蒸判别示例
==========================

利用 physics_engine.valve + thermo_helper 物理引擎，完成：
  1. 65°C 水的物性计算 (密度、蒸气压)
  2. 阻塞流判别 (ΔP vs ΔP_choked)
  3. 空化系数校核 (σ_actual vs σ_allowable)
  4. 闪蒸判别 (P2 vs Pv)

已知条件:
  阀前绝压 P1 = 500 kPa
  阀后绝压 P2 = 120 kPa
  介质: 水, 65°C
  阀门压力恢复系数 FL = 0.9
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
    g = 9.81

    # ================================================================
    # 1. 已知条件
    # ================================================================
    sep("1. 已知条件")

    P1 = 500e3                   # 阀前绝对压力 (Pa) = 500 kPa
    P2 = 120e3                   # 阀后绝对压力 (Pa) = 120 kPa
    T_C = 65.0                   # 操作温度 (°C)
    T_K = T_C + 273.15           # 操作温度 (K)
    FL_given = 0.9               # 阀门压力恢复系数
    Pc = 22.064e6                # 水的临界压力 (Pa)

    delta_P = P1 - P2            # 阀前后压差 (Pa)

    print(f"  工艺参数:")
    print(f"    阀前绝压 P1       = {P1/1000:.0f} kPa = {P1} Pa")
    print(f"    阀后绝压 P2       = {P2/1000:.0f} kPa = {P2} Pa")
    print(f"    压差 ΔP           = {delta_P/1000:.0f} kPa = {delta_P} Pa")
    print(f"    操作温度          = {T_C} °C = {T_K:.2f} K")
    print(f"    压力恢复系数 FL   = {FL_given}")

    # ================================================================
    # 2. 65°C 水的物性 (引擎计算)
    # ================================================================
    sep("2. 65°C 水的物性计算 (thermo_helper)")

    rho = get_fluid_density("water", T_K, P1)
    Pv = get_fluid_vapor_pressure("water", T_K)

    print(f"\n  调用 get_fluid_density('water', {T_K:.2f}K, {P1:.0f}Pa):")
    print(f"    ρ = {rho:.2f} kg/m3")
    print(f"\n  调用 get_fluid_vapor_pressure('water', {T_K:.2f}K):")
    print(f"    Pv = {Pv:.0f} Pa = {Pv/1000:.3f} kPa")
    print(f"\n  水的临界压力 Pc = {Pc/1e6:.3f} MPa")

    SG = rho / 1000.0

    # ================================================================
    # 3. 阻塞流判别
    # ================================================================
    sep("3. 阻塞流判别 (IEC 60534-2-1)")

    # 饱和比修正系数 FF
    FF = 0.96 - 0.28 * math.sqrt(Pv / Pc)

    # 阻塞流临界压差
    delta_P_choked = FL_given**2 * (P1 - FF * Pv)

    print(f"\n  3.1 饱和比修正系数 FF:")
    print(f"    FF = 0.96 - 0.28 × √(Pv/Pc)")
    print(f"       = 0.96 - 0.28 × √({Pv:.0f}/{Pc:.0f})")
    print(f"       = 0.96 - 0.28 × {math.sqrt(Pv/Pc):.6f}")
    print(f"       = {FF:.4f}")

    print(f"\n  3.2 阻塞流临界压差:")
    print(f"    ΔP_choked = FL² × (P1 - FF × Pv)")
    print(f"              = {FL_given:.1f}² × ({P1:.0f} - {FF:.4f} × {Pv:.0f})")
    print(f"              = {FL_given**2:.4f} × ({P1:.0f} - {FF*Pv:.0f})")
    print(f"              = {FL_given**2:.4f} × {P1 - FF*Pv:.0f}")
    print(f"              = {delta_P_choked:.0f} Pa = {delta_P_choked/1000:.1f} kPa")

    print(f"\n  3.3 判别:")
    print(f"    实际 ΔP    = {delta_P/1000:.0f} kPa")
    print(f"    临界 ΔP    = {delta_P_choked/1000:.1f} kPa")
    if delta_P > delta_P_choked:
        print(f"    结论: ΔP ({delta_P/1000:.0f}) > ΔP_choked ({delta_P_choked/1000:.1f})")
        print(f"    → 阻塞流 (Choked Flow)! 阀门达到流通能力极限")
        is_choked = True
    else:
        print(f"    结论: ΔP ({delta_P/1000:.0f}) < ΔP_choked ({delta_P_choked/1000:.1f})")
        print(f"    → 非阻塞流 (Non-choked Flow)")
        is_choked = False

    # ================================================================
    # 4. 闪蒸判别
    # ================================================================
    sep("4. 闪蒸判别")

    print(f"\n  阀后压力 P2 = {P2/1000:.0f} kPa = {P2:.0f} Pa")
    print(f"  饱和蒸气压 Pv = {Pv/1000:.3f} kPa = {Pv:.0f} Pa")

    if P2 < Pv:
        print(f"\n  结论: P2 ({P2/1000:.0f} kPa) < Pv ({Pv/1000:.3f} kPa)")
        print(f"  → 阀后压力低于饱和蒸气压，发生闪蒸 (Flashing)!")
        print(f"    阀后将有气泡持续存在，液相部分汽化")
        is_flashing = True
    else:
        print(f"\n  结论: P2 ({P2/1000:.0f} kPa) > Pv ({Pv/1000:.3f} kPa)")
        print(f"  → 阀后压力高于饱和蒸气压，无闪蒸")
        is_flashing = False

    # ================================================================
    # 5. 空化判别 (FL-based 缩流断面压力法)
    # ================================================================
    sep("5. 空化判别 (FL-based 缩流断面压力法)")

    # 缩流断面压力: P_vc = P1 - ΔP / FL²
    P_vc = P1 - delta_P / (FL_given ** 2)
    sigma_FL = (P_vc - Pv) / delta_P  # FL-based 空化安全裕度

    print(f"\n  5.1 缩流断面压力 (Vena Contracta):")
    print(f"    P_vc = P1 - ΔP / FL²")
    print(f"         = {P1:.0f} - {delta_P:.0f} / {FL_given:.1f}²")
    print(f"         = {P1:.0f} - {delta_P:.0f} / {FL_given**2:.4f}")
    print(f"         = {P1:.0f} - {delta_P/(FL_given**2):.0f}")
    print(f"         = {P_vc:.0f} Pa = {P_vc/1000:.1f} kPa")

    print(f"\n  5.2 空化判据:")
    print(f"    缩流断面压力 P_vc = {P_vc/1000:.1f} kPa")
    print(f"    饱和蒸气压   Pv  = {Pv/1000:.3f} kPa")

    if P_vc < Pv:
        print(f"    结论: P_vc ({P_vc/1000:.1f} kPa) < Pv ({Pv/1000:.3f} kPa)")
        print(f"    → 缩流断面处压力低于蒸气压，发生空化!")
        is_cavitation = True
    else:
        print(f"    结论: P_vc ({P_vc/1000:.1f} kPa) > Pv ({Pv/1000:.3f} kPa)")
        print(f"    → 缩流断面处压力高于蒸气压，无空化")
        is_cavitation = False

    print(f"\n  5.3 空化安全裕度 σ_FL:")
    print(f"    σ_FL = (P_vc - Pv) / ΔP")
    print(f"         = ({P_vc:.0f} - {Pv:.0f}) / {delta_P:.0f}")
    print(f"         = {sigma_FL:.4f}")
    if sigma_FL > 0:
        print(f"    → σ_FL > 0，缩流断面压力高于蒸气压，安全")
    else:
        print(f"    → σ_FL < 0，缩流断面压力低于蒸气压，空化风险")

    # ================================================================
    # 6. 引擎综合计算 (以截止阀为例, 假设一个流量)
    # ================================================================
    sep("6. 引擎综合验证 (size_control_valve_liquid)")

    # 假设一个流量用于引擎调用 (此处取典型值, 不影响空化/闪蒸判别)
    Q_assume = 10.0  # m3/h (假设流量)
    res = size_control_valve_liquid(
        Q=Q_assume, rho=rho, P1=P1, P2=P2, Pv=Pv, Pc=Pc,
        valve_type="globe", FL=FL_given
    )

    print(f"\n  假设流量 Q = {Q_assume} m3/h")
    print(f"  引擎计算结果:")
    print(f"    Kv            = {res['Kv']:.2f}")
    print(f"    Cv            = {res['Cv']:.2f}")
    print(f"    流态          = {res['flow_condition']}")
    print(f"    空化风险      = {res['cavitation_risk']}")
    print(f"    P_vc          = {res['P_vc_Pa']:.0f} Pa = {res['P_vc_Pa']/1000:.1f} kPa")
    print(f"    σ_FL          = {res['sigma_FL']:.4f}")
    print(f"    ΔP_choked     = {res['delta_P_choked_Pa']:.0f} Pa = {res['delta_P_choked_Pa']/1000:.1f} kPa")
    print(f"    FL            = {res['FL']:.2f}")
    print(f"    FF            = {res['FF']:.4f}")

    # ================================================================
    # 7. 综合结论
    # ================================================================
    sep("7. 综合结论")

    print(f"""
  工况条件:
    阀前绝压 P1    = {P1/1000:.0f} kPa
    阀后绝压 P2    = {P2/1000:.0f} kPa
    压差 ΔP        = {delta_P/1000:.0f} kPa
    介质           = 水 ({T_C}°C)
    水蒸气压 Pv    = {Pv/1000:.3f} kPa
    压力恢复系数 FL = {FL_given}

  判别结果:""")

    # 闪蒸
    print(f"    闪蒸: ", end="")
    if is_flashing:
        print(f"发生! P2={P2/1000:.0f}kPa < Pv={Pv/1000:.3f}kPa")
        print(f"          → 阀后液体汽化，气泡持续存在")
    else:
        print(f"不发生。P2={P2/1000:.0f}kPa > Pv={Pv/1000:.3f}kPa")

    # 阻塞流
    print(f"    阻塞流: ", end="")
    if is_choked:
        print(f"发生! ΔP={delta_P/1000:.0f}kPa > ΔP_choked={delta_P_choked/1000:.1f}kPa")
        print(f"          → 阀门达到流通能力极限，流量不再随压差增大")
    else:
        print(f"不发生。ΔP={delta_P/1000:.0f}kPa < ΔP_choked={delta_P_choked/1000:.1f}kPa")

    # 空化 (FL-based 判据)
    print(f"    空化 (FL-based): ", end="")
    if is_cavitation:
        print(f"发生! P_vc={P_vc/1000:.1f}kPa < Pv={Pv/1000:.3f}kPa")
        print(f"          → 缩流断面处压力低于蒸气压，产生气泡并在下游溃灭")
        print(f"          → 造成噪声、振动、阀芯阀座冲蚀")
    else:
        print(f"不发生。P_vc={P_vc/1000:.1f}kPa > Pv={Pv/1000:.3f}kPa")

    # 综合建议
    print(f"\n  工程建议:")
    if is_flashing or is_choked or is_cavitation:
        print(f"    1. 选用抗空化阀内件 (多级降压笼、迷宫式阀芯)")
        print(f"    2. 采用压力恢复系数 FL 更低的阀门 (如球阀 FL=0.65)")
        print(f"    3. 提高阀后背压 (如增设节流孔板)")
        print(f"    4. 阀后管道采用抗冲蚀材料或加厚壁厚")
        print(f"    5. 考虑串联两级阀门分担压降")
    else:
        print(f"    该工况阀门可安全运行")

    print(f"\n{'='*70}")
    print(f"  计算完成!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
