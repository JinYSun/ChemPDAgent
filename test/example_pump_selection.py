"""
泵选型计算示例 — 总扬程与轴功率计算
======================================

利用 physics_engine.pump 物理引擎，完成：
  1. 系统总扬程计算 (位差 + 压差 + 管路阻力)
  2. 水力功率与轴功率计算
  3. 电机功率选型
  4. 标准离心泵系列匹配
  5. 泵效率估算

已知条件:
  吸入液面与排出液面高差 Δz = 15 m
  管路阻力损失合计 hf = 6.5 m 液柱
  排出口表压 P_gauge = 0.3 MPa
  设计流量 Q = 50 m3/h
  介质密度 ρ = 1000 kg/m3
  泵效率 η = 70%
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

from physics_engine.pump import (
    pump_hydraulic_power,
    pump_shaft_power,
    motor_power,
    pump_power_full,
    size_centrifugal_pump_power,
    estimate_pump_efficiency,
    select_centrifugal_pump,
    _select_pump_from_series,
    _CENTRIFUGAL_PUMP_SERIES,
)


def sep(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def main():
    g = 9.81  # 重力加速度 m/s2

    # ================================================================
    # 1. 已知条件
    # ================================================================
    sep("1. 已知条件")

    dz = 15.0               # 吸入液面与排出液面高差 (m)
    hf = 6.5                # 管路阻力损失 (m 液柱)
    P_gauge = 0.3e6         # 排出口表压 (Pa) = 0.3 MPa
    Q_m3h = 50.0            # 设计流量 (m3/h)
    Q = Q_m3h / 3600.0      # 设计流量 (m3/s)
    rho = 1000.0            # 介质密度 (kg/m3)
    eta_pump = 0.70         # 泵效率
    eta_motor = 0.90        # 电机效率
    K_service = 1.15        # 安全系数

    print(f"  系统参数:")
    print(f"    高差 Δz          = {dz} m")
    print(f"    管路阻力损失 hf   = {hf} m 液柱")
    print(f"    排出口表压 P      = {P_gauge/1e6} MPa = {P_gauge} Pa")
    print(f"    吸入侧压力        = 大气压 (表压 = 0)")
    print(f"  流体参数:")
    print(f"    介质密度 ρ        = {rho} kg/m3")
    print(f"    设计流量 Q        = {Q_m3h} m3/h = {Q:.6f} m3/s")
    print(f"  效率参数:")
    print(f"    泵效率 η_pump     = {eta_pump*100:.0f}%")
    print(f"    电机效率 η_motor  = {eta_motor*100:.0f}%")
    print(f"    安全系数 K        = {K_service}")

    # ================================================================
    # 2. 总扬程计算
    # ================================================================
    sep("2. 总扬程计算")

    # 各分项扬程
    H_elevation = dz                                    # 位差扬程 (m)
    H_pressure = P_gauge / (rho * g)                    # 压差扬程 (m)
    H_friction = hf                                     # 阻力损失扬程 (m)

    # 总扬程
    H_total = H_elevation + H_pressure + H_friction

    print(f"\n  伯努利方程 (以吸入液面为基准):")
    print(f"    H = Δz + ΔP/(ρg) + hf")
    print(f"")
    print(f"  各位分项:")
    print(f"    位差扬程   H_z = Δz = {H_elevation} m")
    print(f"    压差扬程   H_P = P/(ρg) = {P_gauge}/({rho}*{g}) = {H_pressure:.4f} m")
    print(f"    阻力扬程   H_f = hf = {H_friction} m")
    print(f"  ─────────────────────────────────")
    print(f"    总扬程     H   = {H_elevation} + {H_pressure:.4f} + {H_friction}")
    print(f"                  = {H_total:.4f} m")
    print(f"")
    print(f"  取设计裕量 10%: H_design = {H_total * 1.1:.2f} m")

    H_design = H_total * 1.1  # 10% 设计裕量

    # ================================================================
    # 3. 功率计算 (引擎函数)
    # ================================================================
    sep("3. 功率计算")

    # 3.1 水力功率
    P_hyd = pump_hydraulic_power(Q, H_total, rho)
    print(f"\n  3.1 水力功率 (pump_hydraulic_power)")
    print(f"    P_hyd = ρ*g*Q*H = {rho}*{g}*{Q:.6f}*{H_total:.4f}")
    print(f"          = {P_hyd:.2f} W = {P_hyd/1000:.4f} kW")

    # 3.2 轴功率
    P_shaft = pump_shaft_power(P_hyd, eta_pump)
    print(f"\n  3.2 轴功率 (pump_shaft_power)")
    print(f"    P_shaft = P_hyd / η_pump = {P_hyd:.2f} / {eta_pump}")
    print(f"            = {P_shaft:.2f} W = {P_shaft/1000:.4f} kW")

    # 3.3 电机功率
    P_motor = motor_power(P_shaft, eta_motor, K_service)
    print(f"\n  3.3 电机功率 (motor_power)")
    print(f"    P_motor = P_shaft / η_motor * K = {P_shaft:.2f} / {eta_motor} * {K_service}")
    print(f"            = {P_motor:.2f} W = {P_motor/1000:.4f} kW")

    # 3.4 一站式计算 (pump_power_full)
    print(f"\n  3.4 一站式计算 (pump_power_full)")
    power_result = pump_power_full(Q, H_total, rho, eta_pump, eta_motor, K_service)
    print(f"    水力功率  = {power_result['P_hydraulic']:.2f} W = {power_result['P_hydraulic']/1000:.4f} kW")
    print(f"    轴功率    = {power_result['P_shaft']:.2f} W = {power_result['P_shaft']/1000:.4f} kW")
    print(f"    电机功率  = {power_result['P_motor']:.2f} W = {power_result['P_motor']/1000:.4f} kW")

    # ================================================================
    # 4. 效率估算与选型
    # ================================================================
    sep("4. 泵效率估算与选型")

    # 引擎效率估算
    eff_est = estimate_pump_efficiency(Q, H_total)
    print(f"\n  引擎效率估算 (estimate_pump_efficiency):")
    print(f"    η_est = {eff_est*100:.1f}%")
    print(f"    (基于比转速经验关联式, 与用户给定 {eta_pump*100:.0f}% 对比)")

    # 离心泵选型
    print(f"\n  离心泵选型 (select_centrifugal_pump):")
    pump_match = select_centrifugal_pump(Q, H_design)
    if pump_match:
        print(f"    推荐型号: {pump_match.get('model', 'N/A')}")
        print(f"    额定流量: {pump_match.get('Q_rated', 0):.1f} m3/h")
        print(f"    额定扬程: {pump_match.get('H_rated', 0):.1f} m")
        print(f"    转速:     {pump_match.get('n', 0)} rpm")
        print(f"    效率:     {pump_match.get('efficiency', 0)*100:.0f}%")
        print(f"    匹配度:   {'合适' if pump_match.get('suitable') else '需复核'}")
        if pump_match.get('warning'):
            print(f"    警告:     {pump_match['warning']}")

    # 从系列中精确匹配
    print(f"\n  IS 系列离心泵精确匹配 (_select_pump_from_series):")
    series_match = _select_pump_from_series(
        _CENTRIFUGAL_PUMP_SERIES, Q_m3h, H_design, "IS清水离心泵"
    )
    print(f"    推荐型号: {series_match.get('model', 'N/A')}")
    print(f"    流量范围: {series_match.get('Q_min', 0)}~{series_match.get('Q_max', 0)} m3/h")
    print(f"    扬程范围: {series_match.get('H_min', 0)}~{series_match.get('H_max', 0)} m")
    print(f"    转速:     {series_match.get('n', 0)} rpm")
    print(f"    效率:     {series_match.get('efficiency', 0)*100:.0f}%")

    # ================================================================
    # 5. 不同效率下的功率对比
    # ================================================================
    sep("5. 不同泵效率下的功率对比")

    efficiencies = [0.50, 0.60, 0.70, 0.75, 0.80, 0.85]
    print(f"\n  {'η_pump':>7s} | {'P_hyd(kW)':>10s} | {'P_shaft(kW)':>12s} | {'P_motor(kW)':>12s}")
    print(f"  {'-'*7}-+-{'-'*10}-+-{'-'*12}-+-{'-'*12}")

    for eff in efficiencies:
        P_h = pump_hydraulic_power(Q, H_total, rho)
        P_s = pump_shaft_power(P_h, eff)
        P_m = motor_power(P_s, eta_motor, K_service)
        marker = " <-- 本题" if abs(eff - eta_pump) < 0.001 else ""
        print(f"  {eff*100:7.0f}% | {P_h/1000:10.4f} | {P_s/1000:12.4f} | {P_m/1000:12.4f}{marker}")

    # ================================================================
    # 6. 设计结论
    # ================================================================
    sep("6. 设计结论")

    # 标准电机功率档位
    motor_sizes = [0.75, 1.1, 1.5, 2.2, 3.0, 4.0, 5.5, 7.5, 11, 15, 18.5, 22, 30, 37, 45, 55, 75, 90, 110]
    P_motor_kw = P_motor / 1000
    selected_motor = None
    for size in motor_sizes:
        if size >= P_motor_kw:
            selected_motor = size
            break

    print(f"""
  系统参数:
    总流量 Q = {Q_m3h} m3/h
    总扬程 H = {H_total:.2f} m (含位差 {H_elevation}m + 压差 {H_pressure:.2f}m + 阻力 {H_friction}m)
    设计扬程 = {H_design:.2f} m (含10%裕量)

  功率计算 (η_pump = {eta_pump*100:.0f}%):
    水力功率 P_hyd   = {P_hyd:.2f} W = {P_hyd/1000:.4f} kW
    轴功率 P_shaft   = {P_shaft:.2f} W = {P_shaft/1000:.4f} kW
    电机功率 P_motor = {P_motor:.2f} W = {P_motor/1000:.4f} kW

  选型建议:
    推荐泵型: {series_match.get('model', 'N/A')} (IS系列清水离心泵)
    额定转速: {series_match.get('n', 0)} rpm
    泵效率:   {series_match.get('efficiency', 0)*100:.0f}%
    配套电机: {selected_motor} kW (标准电机功率档位)
""")

    print(f"{'='*70}")
    print(f"  计算完成!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
