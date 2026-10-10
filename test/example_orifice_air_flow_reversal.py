# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
"""
已安装孔板流量计反算实际流量 (ISO 5167-2)
============================================
已知:
  孔径 d = 60 mm, 管道内径 D = 150 mm
  介质: 20 degC 空气 (rho=1.2 kg/m3, mu=1.8e-5 Pa.s)
  实测差压 dP = 15 kPa

求解:
  1. 反算实际质量流量与体积流量
  2. 判断读数是否处于有效测量范围
"""
sys.path.insert(0, r"c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents")

import math
from physics_engine.instrument import calculate_orifice_flow, _iso5167_discharge_coefficient

print("=" * 70)
print("       已安装孔板流量计 -- 反算实际流量")
print("=" * 70)

# ============================================================
# 1. 已知参数
# ============================================================
d_orifice = 0.060        # m
D_pipe = 0.150           # m
rho = 1.2                # kg/m3
mu = 1.8e-5              # Pa.s
dP = 15000.0             # Pa (15 kPa)
kappa = 1.4              # 空气绝热指数

beta = d_orifice / D_pipe
A_pipe = math.pi / 4 * D_pipe ** 2
A_orifice = math.pi / 4 * d_orifice ** 2

print(f"\n[已知条件]")
print(f"  孔径 d = {d_orifice*1000:.0f} mm")
print(f"  管道内径 D = {D_pipe*1000:.0f} mm")
print(f"  直径比 beta = d/D = {beta:.4f}")
print(f"  介质: 20 degC 空气")
print(f"  rho = {rho} kg/m3, mu = {mu} Pa.s")
print(f"  kappa = {kappa}")
print(f"  实测差压 dP = {dP/1000:.0f} kPa")
print(f"  管道截面积 A_pipe = {A_pipe:.6f} m2")
print(f"  孔口面积 A_orifice = {A_orifice:.6f} m2")

# ============================================================
# 2. 反算流量 (液体模式, epsilon=1.0)
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤1] 反算实际流量 (液体模式)")

result_liq = calculate_orifice_flow(
    dP=dP,
    D_pipe=D_pipe,
    d_orifice=d_orifice,
    rho=rho,
    mu=mu,
)

q_m_liq = result_liq["mass_flow_kg_s"]
q_v_liq = result_liq["volume_flow_m3_s"]
v_liq = result_liq["velocity_m_s"]
Re_liq = result_liq["Re_D"]
C_liq = result_liq["C"]

print(f"\n  计算结果 (epsilon = 1.0, 不可压缩近似):")
print(f"    质量流量 q_m = {q_m_liq:.4f} kg/s")
print(f"    体积流量 q_v = {q_v_liq:.6f} m3/s = {q_v_liq*3600:.2f} m3/h")
print(f"    管道流速 v = {v_liq:.4f} m/s")
print(f"    雷诺数 Re_D = {Re_liq:.0f}")
print(f"    流出系数 C = {C_liq:.6f}")
print(f"    永久压损 = {result_liq['permanent_pressure_loss_Pa']:.0f} Pa")

# ============================================================
# 3. 气体修正 (假设上游表压 100 kPaG, 绝对压力 201.325 kPa)
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤2] 气体可压缩性修正 (epsilon)")

# 典型工业空气管道, 假设上游表压 ~100 kPaG
P1_abs_options = [101325, 201325, 301325, 501325]
print(f"\n  空气为可压缩气体, 需上游绝对压力计算膨胀系数 epsilon")
print(f"  epsilon = 1 - (0.41 + 0.35*beta^4) * dP / (kappa * P1_abs)")
print(f"\n  以下列出不同上游压力下的修正结果:")
print(f"\n  {'P1_abs(kPa)':<14} {'epsilon':<10} {'q_m(kg/s)':<12} {'q_v(m3/h)':<12} {'Re_D':<12}")
print(f"  {'-'*62}")

for P1_abs in P1_abs_options:
    eps = 1.0 - (0.41 + 0.35 * beta**4) * dP / (kappa * P1_abs)
    eps = max(0.66, min(1.0, eps))
    # 质量流量正比于 epsilon
    q_m_gas = q_m_liq * eps  # 近似: q_m ~ epsilon (简化)
    # 更精确: 重新计算
    E = 1.0 / math.sqrt(1 - beta**4)
    C_gas = C_liq  # 近似不变
    q_m_exact = C_gas * E * eps * (math.pi/4) * d_orifice**2 * math.sqrt(2 * rho * dP)
    q_v_gas = q_m_exact / rho
    Re_gas = rho * (q_v_gas / A_pipe) * D_pipe / mu
    print(f"  {P1_abs/1000:<14.1f} {eps:<10.6f} {q_m_exact:<12.4f} {q_v_gas*3600:<12.2f} {Re_gas:<12.0f}")

# 采用大气压 (最保守) 作为基准
P1_abs_ref = 101325.0  # 1 atm
eps_ref = 1.0 - (0.41 + 0.35 * beta**4) * dP / (kappa * P1_abs_ref)
eps_ref = max(0.66, min(1.0, eps_ref))

result_gas = calculate_orifice_flow(
    dP=dP,
    D_pipe=D_pipe,
    d_orifice=d_orifice,
    rho=rho,
    mu=mu,
    P1_absolute=P1_abs_ref,
    kappa=kappa,
)

q_m_gas = result_gas["mass_flow_kg_s"]
q_v_gas = result_gas["volume_flow_m3_s"]
Re_gas = result_gas["Re_D"]

print(f"\n  基准工况 (P1_abs = {P1_abs_ref/1000:.1f} kPa, 即 1 atm):")
print(f"    epsilon = {result_gas['epsilon']:.6f}")
print(f"    q_m = {q_m_gas:.4f} kg/s")
print(f"    q_v = {q_v_gas*3600:.2f} m3/h")
print(f"    Re_D = {Re_gas:.0f}")

# ============================================================
# 4. 有效测量范围判断
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤3] 有效测量范围判断")

# ISO 5167 要求
print(f"\n  ISO 5167-2 有效测量条件:")

# (1) beta 范围
beta_ok = 0.1 <= beta <= 0.75
print(f"\n  (1) 直径比 beta = {beta:.4f}")
print(f"      要求: 0.1 <= beta <= 0.75")
print(f"      判定: {'OK' if beta_ok else 'NG'}")

# (2) 雷诺数范围
Re_min_ISO = 5000
# 更精确: ISO 5167 对角接取压有 Re 下限公式
Re_limit = 5000  # 简化
Re_ok = Re_gas >= Re_min_ISO
print(f"\n  (2) 管道雷诺数 Re_D = {Re_gas:.0f}")
print(f"      要求: Re_D >= {Re_min_ISO} (湍流, C 稳定)")
print(f"      判定: {'OK' if Re_ok else 'NG'}")

# (3) 差压与满量程比 (典型要求 差压在满量程 20%~80% 之间最佳)
# 假设差压变送器满量程 40 kPa (常见)
dP_fullscale = 40000.0
ratio = dP / dP_fullscale * 100
print(f"\n  (3) 差压占满量程比 (假设变送器满量程 40 kPa):")
print(f"      dP/dP_FS = {dP/1000:.0f}/{dP_fullscale/1000:.0f} = {ratio:.1f}%")
if 20 <= ratio <= 80:
    print(f"      判定: OK (20%~80%, 最佳测量区间)")
elif 10 <= ratio <= 90:
    print(f"      判定: WARN (10%~90%, 可用但精度略降)")
else:
    print(f"      判定: NG (偏离最佳区间)")

# (4) epsilon 有效性
eps_ok = result_gas["epsilon"] >= 0.66
print(f"\n  (4) 可膨胀性系数 epsilon = {result_gas['epsilon']:.6f}")
print(f"      要求: epsilon >= 0.66 (ISO 5167 下限)")
print(f"      判定: {'OK' if eps_ok else 'NG'}")

# (5) 流量测量不确定度分析
# 差压测量误差 -> 流量误差: dq/q = 0.5 * d(dP)/dP
# 典型差压变送器精度: +-0.1% FS
dP_error = 0.001 * dP_fullscale  # 0.1% of FS
q_error_pct = 0.5 * dP_error / dP * 100
print(f"\n  (5) 流量测量不确定度:")
print(f"      差压变送器精度: +-0.1% FS = +-{-dP_error:.0f} Pa")
print(f"      流量相对误差: dq/q = 0.5 * d(dP)/dP = +-{q_error_pct:.2f}%")

# ============================================================
# 5. 最终结论
# ============================================================
print(f"\n{'='*70}")
print(f"[设计结论]")
print(f"  已安装孔板: d = {d_orifice*1000:.0f} mm, D = {D_pipe*1000:.0f} mm, beta = {beta:.4f}")
print(f"  介质: 20 degC 空气 (rho={rho}, mu={mu})")
print(f"  实测差压: {dP/1000:.0f} kPa")
print(f"")
print(f"  反算流量:")
print(f"    不可压缩近似: q_v = {q_v_liq*3600:.2f} m3/h, q_m = {q_m_liq:.4f} kg/s")
print(f"    气体修正后:   q_v = {q_v_gas*3600:.2f} m3/h, q_m = {q_m_gas:.4f} kg/s")
print(f"    (气体修正基于 P1_abs = {P1_abs_ref/1000:.1f} kPa)")
print(f"")
print(f"  有效测量范围判定:")
print(f"    beta = {beta:.4f}: {'OK' if beta_ok else 'NG'}")
print(f"    Re_D = {Re_gas:.0f}: {'OK' if Re_ok else 'NG'}")
print(f"    epsilon = {result_gas['epsilon']:.4f}: {'OK' if eps_ok else 'NG'}")
print(f"    差压占比 {ratio:.1f}%: {'OK' if 20<=ratio<=80 else 'WARN'}")
print(f"    流量不确定度: +-{q_error_pct:.2f}%")
print(f"{'='*70}")
