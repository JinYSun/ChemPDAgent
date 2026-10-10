# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
"""
DN150 管道孔板流量计设计计算 (ISO 5167-2)
============================================
已知:
  管道 DN150, 设计流量 20~80 m3/h (水)
  差压变送器满量程 dP_max = 40 kPa
  水物性: rho = 998 kg/m3, mu = 0.001 Pa.s (20 degC)

求解:
  1. 孔板开孔直径 d 及直径比 beta
  2. 最大流量下雷诺数校核
  3. 流量范围下限雷诺数校核
"""
sys.path.insert(0, r"c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents")

import math
from physics_engine.instrument import orifice_sizing, calculate_orifice_flow

print("=" * 70)
print("       DN150 孔板流量计设计计算 (ISO 5167-2)")
print("=" * 70)

# ============================================================
# 1. 已知参数
# ============================================================
DN = 150                   # 公称直径
D_pipe = 0.150             # 管道内径 m (DN150 近似取内径=公称直径)
Q_min = 20.0 / 3600        # m3/s
Q_max = 80.0 / 3600        # m3/s
dP_max = 40000.0           # Pa (40 kPa)
rho = 998.0                # kg/m3 (20 degC 水)
mu = 0.001                 # Pa.s (20 degC 水)

A_pipe = math.pi / 4 * D_pipe ** 2

print(f"\n[已知条件]")
print(f"  管道: DN{DN}, D = {D_pipe*1000:.0f} mm")
print(f"  流量范围: {20:.0f} ~ {80:.0f} m3/h")
print(f"  Q_min = {Q_min:.6f} m3/s")
print(f"  Q_max = {Q_max:.6f} m3/s")
print(f"  dP_max = {dP_max:.0f} Pa (40 kPa)")
print(f"  水物性 (20 degC): rho = {rho} kg/m3, mu = {mu} Pa.s")
print(f"  管道截面积 A = {A_pipe:.6f} m2")

# 管道流速范围
v_min = Q_min / A_pipe
v_max = Q_max / A_pipe
print(f"  管道流速: {v_min:.3f} ~ {v_max:.3f} m/s")

# ============================================================
# 2. 孔板孔径设计 (以最大流量为基准)
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤1] 孔板孔径设计 (以 Q_max = {80:.0f} m3/h 为基准)")

sizing = orifice_sizing(
    Q_target=Q_max,
    D_pipe=D_pipe,
    rho=rho,
    mu=mu,
    dP_max=dP_max,
)

d_orifice = sizing["d_orifice_m"]
d_orifice_mm = sizing["d_orifice_mm"]
beta = sizing["beta"]
C_discharge = sizing["C_discharge"]
E_velocity = sizing["E_velocity"]
dP_design = sizing["dP_design_Pa"]
Re_D_max = sizing["Re_D"]

print(f"\n  二分法求解结果:")
print(f"    直径比 beta = d/D = {beta:.6f}")
print(f"    孔板开孔直径 d = {d_orifice_mm:.2f} mm")
print(f"    流出系数 C = {C_discharge:.6f}")
print(f"    渐进速度系数 E = {E_velocity:.6f}")
print(f"    设计差压 dP = {dP_design:.1f} Pa (目标 {dP_max:.0f} Pa)")
print(f"    管道雷诺数 Re_D = {Re_D_max:.0f}")

# ============================================================
# 3. 反算校核 (用 calculate_orifice_flow 验证)
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤2] 反算校核 (calculate_orifice_flow)")

verify = calculate_orifice_flow(
    dP=dP_max,
    D_pipe=D_pipe,
    d_orifice=d_orifice,
    rho=rho,
    mu=mu,
)

Q_verify = verify["volume_flow_m3_s"] * 3600  # m3/h
print(f"\n  在 dP = {dP_max:.0f} Pa 下:")
print(f"    计算流量 = {Q_verify:.2f} m3/h (目标 {80:.0f} m3/h)")
print(f"    beta = {verify['beta']:.6f}")
print(f"    C = {verify['C']:.6f}")
print(f"    Re_D = {verify['Re_D']:.0f}")
print(f"    永久压损 = {verify['permanent_pressure_loss_Pa']:.0f} Pa")

# ============================================================
# 4. 各流量点差压与雷诺数
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤3] 流量-差压-雷诺数特性")

# 差压与流量平方成正比: dP ~ Q^2
# 所以 dP(Q) = dP_max * (Q/Q_max)^2

print(f"\n  {'Q(m3/h)':<10} {'dP(Pa)':<12} {'Re_D':<14} {'v_pipe(m/s)':<12} {'v_orifice(m/s)':<14}")
print(f"  {'-'*64}")

A_orifice = math.pi / 4 * d_orifice ** 2

for Q_m3h in [20, 30, 40, 50, 60, 70, 80]:
    Q_s = Q_m3h / 3600.0
    # 差压: dP = dP_max * (Q/Q_max)^2
    dP = dP_max * (Q_m3h / 80.0) ** 2
    # 管道流速
    v_pipe = Q_s / A_pipe
    # 孔口流速 (近似)
    v_orifice = Q_s / A_orifice
    # 雷诺数
    Re = rho * v_pipe * D_pipe / mu
    print(f"  {Q_m3h:<10} {dP:<12.0f} {Re:<14.0f} {v_pipe:<12.3f} {v_orifice:<14.2f}")

# ============================================================
# 5. 雷诺数校核
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤4] 雷诺数校核 (ISO 5167 要求)")

Re_min = rho * v_min * D_pipe / mu
Re_max_calc = Re_D_max

print(f"\n  ISO 5167-2 要求: Re_D >= 5000 (湍流, 流出系数稳定)")
print(f"  推荐范围: Re_D >= 10000 (测量精度最佳)")
print(f"")
print(f"  Q = 20 m3/h (最小流量):")
print(f"    v_pipe = {v_min:.3f} m/s")
print(f"    Re_D = {Re_min:.0f}")
if Re_min >= 10000:
    print(f"    判定: OK (Re > 10000, 精度最佳范围)")
elif Re_min >= 5000:
    print(f"    判定: WARN (5000 < Re < 10000, 可用但精度略降)")
else:
    print(f"    判定: NG (Re < 5000, 流出系数不稳定, 测量误差大)")

print(f"\n  Q = 80 m3/h (最大流量):")
print(f"    v_pipe = {v_max:.3f} m/s")
print(f"    Re_D = {Re_max_calc:.0f}")
if Re_max_calc >= 10000:
    print(f"    判定: OK (Re > 10000, 精度最佳范围)")
elif Re_max_calc >= 5000:
    print(f"    判定: WARN (5000 < Re < 10000, 可用但精度略降)")
else:
    print(f"    判定: NG (Re < 5000, 流出系数不稳定)")

# 量程比
rangeability = Q_max / Q_min
print(f"\n  量程比 = Q_max / Q_min = {80:.0f}/{20:.0f} = {rangeability:.1f}")
# 孔板流量计典型量程比 3:1 ~ 4:1
print(f"  孔板流量计典型量程比: 3:1 ~ 4:1")
print(f"  当前量程比 {rangeability:.1f}:1 {'在推荐范围内' if rangeability <= 4 else '偏大, 建议考虑其他方案'}")

# ============================================================
# 6. 最终结论
# ============================================================
print(f"\n{'='*70}")
print(f"[设计结论]")
print(f"  管道: DN{DN}, D = {D_pipe*1000:.0f} mm")
print(f"  流量范围: 20 ~ 80 m3/h (水, 20 degC)")
print(f"  差压变送器满量程: {dP_max/1000:.0f} kPa")
print(f"")
print(f"  孔板设计:")
print(f"    直径比 beta = {beta:.4f}")
print(f"    开孔直径 d = {d_orifice_mm:.2f} mm")
print(f"    流出系数 C = {C_discharge:.4f}")
print(f"    渐进速度系数 E = {E_velocity:.4f}")
print(f"")
print(f"  雷诺数校核:")
print(f"    Q = 80 m3/h: Re_D = {Re_max_calc:.0f} {'OK' if Re_max_calc >= 5000 else 'NG'}")
print(f"    Q = 20 m3/h: Re_D = {Re_min:.0f} {'OK' if Re_min >= 5000 else 'WARN/NG'}")
print(f"")
print(f"  永久压损 (满量程): {verify['permanent_pressure_loss_Pa']:.0f} Pa")
print(f"{'='*70}")
