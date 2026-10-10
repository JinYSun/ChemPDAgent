"""
管道校核计算示例 — DN100 钢管输送水
=====================================

利用 physics_engine.pipe + thermo_helper 物理引擎，完成：
  1. 水的物性计算（密度、粘度）
  2. 流速校核（是否在推荐范围内）
  3. 雷诺数及流态判断
  4. Churchill 摩擦因子计算
  5. 直管段摩擦压降 + 高程压降
  6. 系统总压降

已知条件:
  管道: DN100 钢管, 内径 D = 102.3 mm, 绝对粗糙度 ε = 0.045 mm
  流量: Q = 60 m3/h
  管长: L = 200 m
  标高差: Δz = 8 m (提升)
  流体: 水, 默认 20 degC, 1 atm
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

from physics_engine.thermo_helper import (
    get_fluid_density,
    get_fluid_viscosity,
)
from physics_engine.pipe import (
    calculate_pipe_pressure_drop,
    calculate_pipe_fitting_pressure_drop,
    calculate_total_pipe_pressure_drop,
    churchill_friction,
    FITTING_K_VALUES,
    _RECOMMENDED_VELOCITY,
)


def sep(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def main():
    # ================================================================
    # 1. 已知条件汇总
    # ================================================================
    sep("1. 已知条件")

    D_mm = 102.3              # 管道内径 (mm)
    D = D_mm / 1000.0         # 管道内径 (m)
    roughness_mm = 0.045      # 绝对粗糙度 (mm)
    roughness = roughness_mm / 1000.0  # 绝对粗糙度 (m)
    Q_m3h = 60.0              # 体积流量 (m3/h)
    Q = Q_m3h / 3600.0        # 体积流量 (m3/s)
    L = 200.0                 # 管长 (m)
    dz = 8.0                  # 标高差 (m), 提升为正
    T_C = 20.0                # 水温 (degC)
    T_K = T_C + 273.15        # 水温 (K)
    P = 101325.0              # 压力 (Pa), 1 atm

    print(f"  管道参数:")
    print(f"    公称直径       DN100")
    print(f"    内径 D         = {D_mm} mm = {D} m")
    print(f"    绝对粗糙度 ε   = {roughness_mm} mm = {roughness:.5e} m")
    print(f"    相对粗糙度 ε/D = {roughness/D:.6e}")
    print(f"    管长 L         = {L} m")
    print(f"    标高差 Δz      = {dz} m (提升)")
    print(f"  流体参数:")
    print(f"    流体           = 水")
    print(f"    温度 T         = {T_C} degC = {T_K:.2f} K")
    print(f"    压力 P         = {P:.0f} Pa = {P/101325:.2f} atm")
    print(f"  流量参数:")
    print(f"    体积流量 Q     = {Q_m3h} m3/h = {Q:.6f} m3/s")

    # ================================================================
    # 2. 水的物性计算 (thermo_helper)
    # ================================================================
    sep("2. 水的物性计算 (thermo_helper 引擎)")

    rho = get_fluid_density("water", T_K, P, phase="liquid")
    mu = get_fluid_viscosity("water", T_K, P, phase="liquid")

    print(f"  调用 get_fluid_density('water', T={T_K:.2f} K, P={P:.0f} Pa, phase='liquid')")
    print(f"    => 密度 rho = {rho:.4f} kg/m3")
    print(f"  调用 get_fluid_viscosity('water', T={T_K:.2f} K, P={P:.0f} Pa, phase='liquid')")
    print(f"    => 动力粘度 mu = {mu:.6e} Pa.s")

    nu = mu / rho  # 运动粘度
    print(f"  运动粘度 nu = mu/rho = {mu:.6e}/{rho:.4f} = {nu:.6e} m2/s")

    # ================================================================
    # 3. 流速计算与校核
    # ================================================================
    sep("3. 流速计算与校核")

    A_pipe = math.pi * (D / 2) ** 2   # 管道截面积 (m2)
    v = Q / A_pipe                      # 实际流速 (m/s)

    v_range = _RECOMMENDED_VELOCITY["liquid"]
    v_min, v_max = v_range

    print(f"  管道截面积 A = pi*(D/2)^2 = pi*({D}/2)^2 = {A_pipe:.6f} m2")
    print(f"  实际流速 v = Q/A = {Q:.6f}/{A_pipe:.6f} = {v:.4f} m/s")
    print(f"  液体推荐流速范围: {v_min}~{v_max} m/s")

    if v_min <= v <= v_max:
        verdict = "合理"
        detail = f"{v_min} <= {v:.4f} <= {v_max}, 流速在推荐范围内"
    elif v < v_min:
        verdict = "偏低"
        detail = f"{v:.4f} < {v_min}, 流速低于推荐下限, 可能引起沉积"
    else:
        verdict = "偏高"
        detail = f"{v:.4f} > {v_max}, 流速超出推荐上限, 冲刷/噪声风险增大"

    print(f"\n  >>> 流速校核结论: {verdict}")
    print(f"      {detail}")

    # ================================================================
    # 4. 雷诺数与流态判断
    # ================================================================
    sep("4. 雷诺数与流态判断")

    Re = rho * v * D / mu

    if Re < 2300:
        regime = "层流 (Laminar)"
    elif Re < 4000:
        regime = "过渡区 (Transitional)"
    else:
        regime = "湍流 (Turbulent)"

    print(f"  Re = rho*v*D/mu = {rho:.4f}*{v:.4f}*{D}/{mu:.6e}")
    print(f"     = {Re:.0f}")
    print(f"  流态判断: Re = {Re:.0f} => {regime}")

    # ================================================================
    # 5. Churchill 摩擦因子计算
    # ================================================================
    sep("5. Churchill 摩擦因子计算")

    f = churchill_friction(Re, roughness, D)

    rel_rough = roughness / D
    print(f"  相对粗糙度 ε/D = {roughness:.5e}/{D} = {rel_rough:.6e}")
    print(f"  调用 churchill_friction(Re={Re:.0f}, roughness={roughness:.5e}, D={D})")
    print(f"    => Darcy 摩擦因子 f = {f:.6e}")

    # ================================================================
    # 6. 直管段压降计算
    # ================================================================
    sep("6. 直管段压降计算 (Darcy-Weisbach)")

    # 使用引擎函数 (含高程)
    straight_result = calculate_pipe_pressure_drop(
        Q=Q, D=D, L=L, rho=rho, mu=mu,
        roughness=roughness, elevation_change=dz,
        method="churchill"
    )

    dP_friction = straight_result["friction_drop_Pa"]
    dP_elevation = straight_result["elevation_drop_Pa"]
    dP_straight_total = straight_result["pressure_drop_Pa"]

    # 手动验证
    dP_friction_manual = f * (L / D) * (rho * v ** 2 / 2)
    dP_elevation_manual = rho * 9.81 * dz

    print(f"  摩擦压降:")
    print(f"    dP_f = f*(L/D)*(rho*v^2/2)")
    print(f"         = {f:.6e}*({L}/{D})*({rho:.4f}*{v:.4f}^2/2)")
    print(f"         = {dP_friction_manual:.2f} Pa = {dP_friction_manual/1000:.4f} kPa")
    print(f"    引擎计算值: {dP_friction:.2f} Pa = {dP_friction/1000:.4f} kPa")
    print(f"  高程压降:")
    print(f"    dP_z = rho*g*dz = {rho:.4f}*9.81*{dz}")
    print(f"         = {dP_elevation_manual:.2f} Pa = {dP_elevation_manual/1000:.4f} kPa")
    print(f"    引擎计算值: {dP_elevation:.2f} Pa = {dP_elevation/1000:.4f} kPa")
    print(f"  直管段总压降:")
    print(f"    dP_straight = dP_f + dP_z = {dP_straight_total:.2f} Pa = {dP_straight_total/1000:.4f} kPa")
    print(f"  单位管长压降:")
    print(f"    dP/L = {straight_result['pressure_drop_per_m']:.4f} Pa/m")
    print(f"  流态: {straight_result['flow_regime']}")

    # ================================================================
    # 7. 系统总压降 (本题无管件, 总压降 = 直管段压降)
    # ================================================================
    sep("7. 系统总压降汇总")

    # 本题未指定管件, 但展示若加入典型管件后的计算方式
    print(f"  本题未指定管件, 系统总压降 = 直管摩擦压降 + 高程压降")
    print(f"")
    print(f"  {'项目':<20s} | {'压降 (Pa)':>12s} | {'压降 (kPa)':>12s}")
    print(f"  {'-'*20}-+-{'-'*12}-+-{'-'*12}")
    print(f"  {'直管摩擦压降':<18s} | {dP_friction:12.2f} | {dP_friction/1000:12.4f}")
    print(f"  {'高程压降':<18s} | {dP_elevation:12.2f} | {dP_elevation/1000:12.4f}")
    print(f"  {'-'*20}-+-{'-'*12}-+-{'-'*12}")
    print(f"  {'系统总压降':<18s} | {dP_straight_total:12.2f} | {dP_straight_total/1000:12.4f}")

    # 补充: 若加入常见管件的示例
    print(f"\n  [补充] 若管路包含典型管件的额外压降示例:")
    sample_fittings = [
        {"type": "elbow_90_long", "count": 5},
        {"type": "gate_valve", "count": 2},
    ]
    fitting_result = calculate_pipe_fitting_pressure_drop(Q, D, rho, sample_fittings)
    dP_fitting = fitting_result["pressure_drop_Pa"]
    v_head = fitting_result["velocity_head_Pa"]

    print(f"    速度头 rho*v^2/2 = {v_head:.4f} Pa")
    for item in fitting_result["fittings_breakdown"]:
        print(f"    {item['type']} x{item['count']}: K_total={item['K_total']:.2f}, "
              f"dP={item['pressure_drop_Pa']:.2f} Pa")
    print(f"    管件总压降 = {dP_fitting:.2f} Pa = {dP_fitting/1000:.4f} kPa")

    dP_grand_total = dP_straight_total + dP_fitting
    print(f"\n    含管件的系统总压降:")
    print(f"      = {dP_straight_total:.2f} + {dP_fitting:.2f} = {dP_grand_total:.2f} Pa = {dP_grand_total/1000:.4f} kPa")

    # ================================================================
    # 8. 使用 calculate_total_pipe_pressure_drop 一站式计算
    # ================================================================
    sep("8. 一站式总压降计算 (calculate_total_pipe_pressure_drop)")

    total_result = calculate_total_pipe_pressure_drop(
        Q=Q, D=D, L=L, rho=rho, mu=mu,
        roughness=roughness, elevation_change=dz,
        fittings=sample_fittings,
    )

    print(f"  调用 calculate_total_pipe_pressure_drop:")
    print(f"    Q={Q:.6f} m3/s, D={D} m, L={L} m")
    print(f"    rho={rho:.4f} kg/m3, mu={mu:.6e} Pa.s")
    print(f"    roughness={roughness:.5e} m, elevation_change={dz} m")
    print(f"    fittings=5x elbow_90_long + 2x gate_valve")
    print(f"")
    print(f"  结果:")
    print(f"    总压降         = {total_result['total_pressure_drop_Pa']:.2f} Pa "
          f"= {total_result['total_pressure_drop_Pa']/1000:.4f} kPa")
    print(f"    直管摩擦压降   = {total_result['straight_pipe_drop_Pa']:.2f} Pa "
          f"= {total_result['straight_pipe_drop_Pa']/1000:.4f} kPa")
    print(f"    管件压降       = {total_result['fitting_drop_Pa']:.2f} Pa "
          f"= {total_result['fitting_drop_Pa']/1000:.4f} kPa")
    print(f"    高程压降       = {total_result['elevation_drop_Pa']:.2f} Pa "
          f"= {total_result['elevation_drop_Pa']/1000:.4f} kPa")
    print(f"    流速           = {total_result['velocity_m_s']:.4f} m/s")
    print(f"    雷诺数         = {total_result['Re']:.0f}")
    print(f"    摩擦因子       = {total_result['friction_factor']:.6e}")
    print(f"    流态           = {total_result['flow_regime']}")

    # ================================================================
    # 9. 设计结论
    # ================================================================
    sep("9. 设计结论")

    print(f"""
  管道: DN100 钢管 (内径 {D_mm} mm, ε = {roughness_mm} mm)
  流体: 水 @ {T_C} degC, 1 atm
  流量: {Q_m3h} m3/h = {Q:.6f} m3/s

  物性参数 (thermo_helper 计算):
    密度 rho = {rho:.4f} kg/m3
    粘度 mu = {mu:.6e} Pa.s

  流速校核:
    实际流速 v = {v:.4f} m/s
    推荐范围   = {v_min}~{v_max} m/s (一般液体)
    结论: 流速{verdict} ({detail})

  水力计算:
    雷诺数 Re = {Re:.0f} ({regime})
    摩擦因子 f = {f:.6e}

  系统总压降 (仅直管 + 高程):
    直管摩擦压降 = {dP_friction:.2f} Pa = {dP_friction/1000:.4f} kPa
    高程压降     = {dP_elevation:.2f} Pa = {dP_elevation/1000:.4f} kPa
    总压降       = {dP_straight_total:.2f} Pa = {dP_straight_total/1000:.4f} kPa

  若含管件 (5x90deg弯头 + 2x闸阀):
    管件压降     = {dP_fitting:.2f} Pa = {dP_fitting/1000:.4f} kPa
    系统总压降   = {dP_grand_total:.2f} Pa = {dP_grand_total/1000:.4f} kPa
""")

    print(f"{'='*70}")
    print(f"  计算完成!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
