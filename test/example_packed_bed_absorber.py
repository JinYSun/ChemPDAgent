# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
"""
吸收塔填料层压降计算
======================
已知:
  气相处理量: 2000 Nm3/h (标准状态 0 degC, 1 atm)
  气相密度: 1.3 kg/m3 (操作条件)
  液相流量: 8 m3/h
  液相密度: 1000 kg/m3
  填料: 50 mm 聚丙烯阶梯环 (CMR)

求解:
  1. 塔径选择与操作气速
  2. 填料层压降 (Ergun 干床 +  irrigated 修正)
"""
sys.path.insert(0, r"c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents")

import math
from physics_engine.packed_bed import (
    select_packed_tower_diameter,
    calculate_packed_bed_pressure_drop,
    calculate_packed_bed_pressure_drop_bs,
    _PACKING_DATA,
)

print("=" * 70)
print("       吸收塔填料层压降计算")
print("=" * 70)

# ============================================================
# 1. 已知参数
# ============================================================
Q_gas_Nm3h = 2000.0     # Nm3/h (标准状态: 0 degC, 101325 Pa)
rho_gas_op = 1.3         # 操作条件下气相密度 kg/m3
Q_liquid_m3h = 8.0       # 液相流量 m3/h
rho_liquid = 1000.0      # 液相密度 kg/m3
mu_gas = 1.8e-5          # 空气黏度 Pa.s (20 degC 近似)
mu_liquid = 0.001        # 水黏度 Pa.s (20 degC)
sigma = 0.0725           # 气液表面张力 N/m (水-空气, 20 degC)
packing_name = "step_cascade_ring_50"

# 标准状态换算
rho_gas_N = 1.293        # 空气标准状态密度 kg/m3 (0 degC, 1 atm)
# 实际操作体积流量: Q_actual = Q_N * rho_N / rho_op
Q_gas_actual_m3h = Q_gas_Nm3h * rho_gas_N / rho_gas_op
Q_gas_actual_m3s = Q_gas_actual_m3h / 3600.0
Q_liquid_m3s = Q_liquid_m3h / 3600.0

# 质量流量
m_gas = Q_gas_actual_m3s * rho_gas_op     # kg/s
m_liquid = Q_liquid_m3s * rho_liquid       # kg/s

print(f"\n[已知条件]")
print(f"  气相: {Q_gas_Nm3h} Nm3/h (标准状态)")
print(f"  操作条件密度: {rho_gas_op} kg/m3")
print(f"  实际体积流量: {Q_gas_Nm3h} * {rho_gas_N}/{rho_gas_op} = {Q_gas_actual_m3h:.1f} m3/h = {Q_gas_actual_m3s:.4f} m3/s")
print(f"  气相质量流量: {m_gas:.4f} kg/s")
print(f"  液相: {Q_liquid_m3h} m3/h, rho = {rho_liquid} kg/m3")
print(f"  液相质量流量: {m_liquid:.4f} kg/s")
print(f"  填料: 50mm 聚丙烯阶梯环 ({packing_name})")

# 填料参数
packing = _PACKING_DATA[packing_name]
dp = packing["dp"]
epsilon = packing["epsilon"]
Fp = packing["Fp"]
a = packing["a"]
print(f"\n  填料参数:")
print(f"    公称直径 dp = {dp*1000:.0f} mm")
print(f"    空隙率 eps = {epsilon}")
print(f"    填料因子 Fp = {Fp}")
print(f"    比表面积 a = {a} m2/m3")

# ============================================================
# 2. 塔径选择 (GPDC 泛点气速法)
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤1] 塔径选择 (Sherwood-Fair GPDC 泛点法)")

tower = select_packed_tower_diameter(
    Q_gas=Q_gas_actual_m3s,
    Q_liquid=Q_liquid_m3s,
    rho_gas=rho_gas_op,
    rho_liquid=rho_liquid,
    mu_liquid=mu_liquid,
    packing_name=packing_name,
    flooding_factor=0.7,    # 设计取 70% 液泛
)

D_tower = tower["D_selected_m"]
A_tower = tower["A_cross_section_m2"]
u_flood = tower["u_flood_m_s"]
u_actual = tower["u_actual_m_s"]
X_flow = tower["flow_parameter_X"]
Y_cap = tower["capacity_parameter_Y"]
flood_ratio = tower["flooding_ratio"]

print(f"\n  流动参数 X = (L/G)*sqrt(rho_g/rho_l) = {X_flow:.4f}")
print(f"  容量参数 Y = {Y_cap:.6f}")
print(f"  泛点气速 u_f = {u_flood:.4f} m/s")
print(f"  设计操作气速 (70%液泛) = {u_flood*0.7:.4f} m/s")
print(f"  计算塔径 = {tower['D_calculated_m']:.4f} m")
print(f"  选定标准塔径 D = {D_tower} m")
print(f"  塔截面积 A = {A_tower:.4f} m2")
print(f"  实际操作气速 u = {u_actual:.4f} m/s")
print(f"  实际液泛率 = {flood_ratio*100:.1f}%")

# ============================================================
# 3. 填料层压降计算 — 三种方法对比
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤2] 填料层压降计算 (三种方法对比)")

# 气相/液相表观流速
u_gas = Q_gas_actual_m3s / A_tower
u_L = Q_liquid_m3s / A_tower

# 当量直径
dp_equiv = 6 * (1 - epsilon) / a
Re_p = rho_gas_op * u_gas * dp_equiv / mu_gas

print(f"\n  气相表观流速 u_G = {u_gas:.4f} m/s")
print(f"  液相表观流速 u_L = {u_L:.6f} m/s")
print(f"  填料当量直径 d_p,equiv = {dp_equiv*1000:.1f} mm")
print(f"  颗粒雷诺数 Re_p = {Re_p:.0f}")

# --- 方法 A: 当量直径 Ergun + Leva 持液修正 ---
print(f"\n  --- 方法A: 当量直径 Ergun + Leva 持液修正 ---")
L_bed = 1.0
r_ergun = calculate_packed_bed_pressure_drop(
    Q=Q_gas_actual_m3s, D_bed=D_tower, L_bed=L_bed,
    rho=rho_gas_op, mu=mu_gas, d_particle=dp_equiv,
    epsilon=epsilon, method="ergun",
)
dP_ergun_dry = r_ergun["pressure_drop_per_m"]

# 简化持液量模型
h_static_simple = 0.03
h_dynamic_simple = 0.12 * u_L ** 0.4
h_total_simple = min(h_static_simple + h_dynamic_simple, epsilon * 0.5)
eps_eff_simple = epsilon - h_total_simple
corr_simple = (epsilon / eps_eff_simple) ** 3
dP_ergun_wet = dP_ergun_dry * corr_simple

print(f"    干压降 = {dP_ergun_dry:.1f} Pa/m")
print(f"    持液量 = {h_total_simple:.4f} m3/m3, eps_eff = {eps_eff_simple:.4f}")
print(f"    湿压降 = {dP_ergun_wet:.1f} Pa/m (修正因子 {corr_simple:.2f})")

# --- 方法 B: Billet-Schultes 严格模型 ---
print(f"\n  --- 方法B: Billet-Schultes 严格模型 ---")
bs = calculate_packed_bed_pressure_drop_bs(
    u_G=u_gas, u_L=u_L,
    rho_G=rho_gas_op, rho_L=rho_liquid,
    mu_G=mu_gas, mu_L=mu_liquid,
    sigma=sigma,
    epsilon=epsilon, a=a, Fp=Fp,
    dp_nominal=dp,
)

dP_bs_dry = bs["pressure_drop_dry_Pa_m"]
dP_bs_wet = bs["pressure_drop_wet_Pa_m"]
corr_bs = bs["correction_factor"]

print(f"    模型参数:")
print(f"      修正雷诺数 Re_G,mod = {bs['Re_G_modified']:.0f} ({bs['flow_regime']})")
print(f"      阻力系数 Psi0 = {bs['Psi0']:.4f} (B0={bs['B0']}, C0={bs['C0']})")
print(f"      液相修正雷诺数 Re_L,mod = {bs['Re_L_modified']:.2f}")
print(f"      Froude 数 Fr_L = {bs['Fr_L']:.6f}")
print(f"")
print(f"    持液量模型:")
print(f"      静态持液量 h_static = {bs['h_static']:.4f} m3/m3")
print(f"      动态持液量 h_dynamic = {bs['h_dynamic']:.4f} m3/m3")
print(f"      总持液量 h_total = {bs['h_total']:.4f} m3/m3")
print(f"      有效空隙率 eps_eff = {bs['epsilon_eff']:.4f}")
print(f"")
print(f"    干压降 = {dP_bs_dry:.1f} Pa/m")
print(f"    湿压降 = {dP_bs_wet:.1f} Pa/m (修正因子 {corr_bs:.2f})")
if bs["warning"]:
    print(f"    !! 警告: {bs['warning']}")

# --- 方法 C: GPDC 反推法 ---
g_const = 9.81
Y_actual = Y_cap * flood_ratio ** 2
dP_Fp = Y_actual * g_const * rho_gas_op * Fp
print(f"\n  --- 方法C: GPDC 反推法 (参考) ---")
print(f"    Y_actual = {Y_actual:.6f}")
print(f"    dP/dZ = {dP_Fp:.1f} Pa/m")

# --- 汇总对比 ---
print(f"\n  {'='*50}")
print(f"  三种方法压降对比:")
print(f"  {'方法':<30s} {'干压降(Pa/m)':<14s} {'湿压降(Pa/m)':<14s}")
print(f"  {'-'*50}")
print(f"  {'A. Ergun(当量直径)+Leva':<30s} {dP_ergun_dry:<14.1f} {dP_ergun_wet:<14.1f}")
print(f"  {'B. Billet-Schultes 严格模型':<30s} {dP_bs_dry:<14.1f} {dP_bs_wet:<14.1f}")
print(f"  {'C. GPDC 反推法(仅湿压降)':<30s} {'---':<14s} {dP_Fp:<14.1f}")
print(f"  {'='*50}")

# 以 Billet-Schultes 为设计基准
dP_design_dry = dP_bs_dry
dP_design_wet = dP_bs_wet

# ============================================================
# 4. 典型填料高度下的总压降 (Billet-Schultes)
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤3] 典型填料高度总压降 (Billet-Schultes 基准)")

heights = [3, 5, 6, 8, 10]
print(f"\n  {'填料高度(m)':<14} {'干压降(Pa)':<14} {'湿压降(Pa)':<14} {'湿压降(mmH2O)':<16}")
print(f"  {'-'*60}")
for H in heights:
    dP_dry_total = dP_design_dry * H
    dP_wet_total = dP_design_wet * H
    dP_mmH2O = dP_wet_total / 9.81
    print(f"  {H:<14} {dP_dry_total:<14.1f} {dP_wet_total:<14.1f} {dP_mmH2O:<16.1f}")

# ============================================================
# 5. 最终结论
# ============================================================
print(f"\n{'='*70}")
print(f"[设计结论]")
print(f"  吸收塔参数:")
print(f"    塔径 D = {D_tower} m")
print(f"    填料: 50mm PP 阶梯环 (Fp={Fp}, eps={epsilon}, a={a} m2/m3)")
print(f"    操作气速 u_G = {u_gas:.4f} m/s, 液泛率 = {flood_ratio*100:.1f}%")
print(f"")
print(f"  填料层压降 (每米, Billet-Schultes 严格模型):")
print(f"    干填料 = {dP_bs_dry:.1f} Pa/m")
print(f"    湿填料 = {dP_bs_wet:.1f} Pa/m")
print(f"")
print(f"  方法对比 (湿压降, Pa/m):")
print(f"    Ergun+Leva:           {dP_ergun_wet:.1f}")
print(f"    Billet-Schultes:      {dP_bs_wet:.1f}")
print(f"    GPDC 反推:            {dP_Fp:.1f}")
print(f"")
print(f"  持液量对比:")
print(f"    Ergun+Leva 简化:      h = {h_total_simple:.4f} m3/m3")
print(f"    Billet-Schultes:      h = {bs['h_total']:.4f} m3/m3 (静{bs['h_static']:.4f}+动{bs['h_dynamic']:.4f})")
print(f"")
print(f"  工程建议:")
if dP_design_wet < 200:
    print(f"    湿填料压降 < 200 Pa/m, 低负荷区, 操作稳定")
elif dP_design_wet < 600:
    print(f"    湿填料压降 200~600 Pa/m, 正常操作范围")
else:
    print(f"    湿填料压降 > 600 Pa/m, 偏高, 需关注液泛风险")
print(f"{'='*70}")
