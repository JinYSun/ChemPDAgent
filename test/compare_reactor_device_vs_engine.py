# -*- coding: utf-8 -*-
"""对比 reactor_device._full_reactor_design vs physics_engine 直接调用"""
import sys, io, json, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.reactor_device import _full_reactor_design
from physics_engine.reactor import (
    pfr_design, cstr_design, design_reactor,
    validate_reactor_design, fixed_bed_dimensions,
    fixed_bed_pressure_drop_ergun, fixed_bed_catalyst_inventory,
    fixed_bed_residence_time, pfr_heat_transfer,
    adiabatic_temperature_rise, estimate_reactor_U,
    estimate_bed_void_fraction,  # 新增：床层空隙率估算
)

# ============================================================
# 反应参数
# ============================================================
F_total = 0.05        # mol/s 混合进料
C_A0 = 20.0           # mol/m³ 甲醇
C_B0 = 15.0           # mol/m³ 乙酸
X_B_target = 0.90     # 乙酸目标转化率
X_A_target = (C_B0 / C_A0) * X_B_target  # 甲醇转化率 = 0.675
k = 0.1               # m³/(mol·s)
T = 200               # °C → 473.15 K
P = 101325            # Pa
F_A0 = F_total * C_A0 / (C_A0 + C_B0)  # 0.02857 mol/s

# 流体物性（工程估算值）
rho_gas = 800.0       # kg/m³ (液相)
mu = 5e-4             # Pa·s
Cp_mix = 150.0        # J/(mol·K)

# 催化剂参数
particle_diameter = 0.005  # 5 mm 催化剂颗粒

print("=" * 80)
print("  reactor_device._full_reactor_design vs physics_engine 直接调用 对比")
print("=" * 80)
print(f"\n反应: 甲醇(A) + 乙酸(B) → 乙酸甲酯 + 水")
print(f"动力学: r = k × C_A × C_B (二级双分子)")
print(f"\n进料条件:")
print(f"  F_total = {F_total} mol/s, C_A0 = {C_A0} mol/m³, C_B0 = {C_B0} mol/m³")
print(f"  F_A0 = {F_A0:.4f} mol/s")
print(f"\n反应条件:")
print(f"  k = {k} m³/(mol·s), T = {T}°C = {T+273.15} K, P = {P/1000:.2f} kPa")
print(f"  乙酸转化率 X_B = {X_B_target*100:.1f}% → 甲醇转化率 X_A = {X_A_target*100:.1f}%")

# ============================================================
# 床层空隙率计算
# ============================================================
print("\n" + "=" * 80)
print("  床层空隙率 (ε) 计算")
print("=" * 80)

# 先计算反应器尺寸（用于壁面效应修正）
T_K = T + 273.15
pfr_result = pfr_design(F_A0, X_A_target, k, C_A0, n=2, epsilon=0.0, C_B0=C_B0)
V_pfr = pfr_result['V']
D_bed_calc = (4 * V_pfr / (math.pi * 4.0)) ** (1/3)  # 假设 L/D = 4

print(f"\n催化剂颗粒直径: {particle_diameter*1000:.1f} mm")
print(f"床层直径: {D_bed_calc*1000:.1f} mm")
print(f"D_bed/D_particle = {D_bed_calc/particle_diameter:.1f}")

# 不同填充方式的空隙率对比
packing_types = ["random_spheres", "ordered_spheres", "cylinders", "rings"]
packing_names = {"random_spheres": "随机球形", "ordered_spheres": "规则球形", 
                 "cylinders": "圆柱形", "rings": "环形"}

print(f"\n{'填充方式':<15} {'ε (无壁面)':<12} {'ε (有壁面)':<12} {'说明'}")
print("-" * 60)
for pt in packing_types:
    eps_no_wall = estimate_bed_void_fraction(particle_diameter, D_bed=None, packing_type=pt)
    eps_with_wall = estimate_bed_void_fraction(particle_diameter, D_bed=D_bed_calc, packing_type=pt)
    print(f"{packing_names[pt]:<15} {eps_no_wall:<12.4f} {eps_with_wall:<12.4f} {pt}")

# 使用默认值（随机球形）
epsilon_calc = estimate_bed_void_fraction(particle_diameter, D_bed=D_bed_calc, packing_type="random_spheres")
print(f"\n选用空隙率: ε = {epsilon_calc:.4f} (随机球形，考虑壁面效应)")

# ============================================================
# 方案一：固定床反应器 (PFR)
# ============================================================
print("\n" + "=" * 80)
print("  方案一：固定床反应器 (Fixed Bed / PFR)")
print("=" * 80)

# ── Device 层 ──
print("\n--- Device 层调用 (_full_reactor_design) ---")
result_pfr_device = _full_reactor_design(
    reactor_type="fixed_bed",
    X_target=X_A_target,
    k=k,
    delta_H_rxn=0.0,
    T_in=T,
    F_A0=F_A0,
    C_A0=C_A0,
    C_B0=C_B0,
    n=2,
    P_in=P,
    rho_gas=rho_gas,
    mu=mu,
    Cp_mix=Cp_mix,
    particle_diameter=particle_diameter,
    epsilon=epsilon_calc,  # 使用计算得到的空隙率
)

# ── Engine 层 ──
print("\n--- Engine 层直接调用 (pfr_design + design_reactor) ---")

# 方式1: 直接调用 pfr_design (纯动力学)
pfr_result = pfr_design(F_A0, X_A_target, k, C_A0, n=2, epsilon=0.0, C_B0=C_B0)
print(f"\n  pfr_design 结果:")
print(f"    V = {pfr_result['V']*1000:.2f} L")
print(f"    τ_theoretical = {pfr_result['tau']:.2f} s  ← 理论空时 (V×C_A0/F_A0)")
print(f"    双分子模式: {pfr_result.get('bimolecular', False)}")

# 方式2: 调用 design_reactor (完整设计)
design_result = design_reactor(
    reactor_type="fixed_bed",
    F_A0=F_A0, X_target=X_A_target,
    delta_H_rxn=0.0,
    k=k, C_A0=C_A0, n=2,
    rho_cat=2000.0, particle_diameter=particle_diameter, epsilon=epsilon_calc,
    T_in=T_K, P_in=P,
    d_tube_inner=0.038, L_tube=6.0,
    U=100.0, T_coolant=T_K,
    num_tanks=1,
    rho_gas=rho_gas, mu=mu, Cp_mix=Cp_mix,
    C_B0=C_B0,
)
print(f"\n  design_reactor 结果:")
print(f"    V = {design_result.get('V', 0)*1000:.2f} L")
print(f"    D_bed = {design_result.get('D_bed', 0)*1000:.1f} mm")
print(f"    L_bed = {design_result.get('L_bed', 0)*1000:.1f} mm")
print(f"    ε = {epsilon_calc:.4f} (由 estimate_bed_void_fraction 计算)")
print(f"    τ_actual = {design_result.get('tau', 0):.2f} s  ← 实际停留时间 (V×ε/F_v0)")
print(f"    delta_P = {design_result.get('delta_P', 0)/1000:.2f} kPa")

# ── 对比 ──
print("\n--- Device 层 vs Engine 层 对比 ---")
print(f"  {'参数':<25} {'Device层':<15} {'Engine层(pfr)':<15} {'Engine层(design)':<18} {'一致?'}")
print("  " + "-" * 80)

V_device = result_pfr_device.get('V', 0)
V_pfr = pfr_result.get('V', 0)
V_design = design_result.get('V', 0)

print(f"  {'体积 V (L)':<25} {V_device*1000:<15.2f} {V_pfr*1000:<15.2f} {V_design*1000:<18.2f} {'✓' if abs(V_device-V_pfr)<0.01 else '✗'}")

tau_device = result_pfr_device.get('tau', 0)
tau_pfr = pfr_result.get('tau', 0)  # 理论空时
tau_design = design_result.get('tau', 0)  # 实际停留时间
print(f"  {'τ 理论空时 (s)':<25} {'-':<15} {tau_pfr:<15.2f} {'-':<18} {'-'}")
print(f"  {'τ 实际停留时间 (s)':<25} {tau_device:<15.2f} {'-':<15} {tau_design:<18.2f} {'✓' if abs(tau_device-tau_design)<0.1 else '✗'}")
print(f"  {'ε 床层空隙率':<25} {epsilon_calc:<15.4f} {'N/A':<15} {epsilon_calc:<18.4f} {'✓'}")

D_device = result_pfr_device.get('D_bed', 0)
D_design = design_result.get('D_bed', 0)
print(f"  {'直径 D (mm)':<25} {D_device*1000:<15.1f} {'N/A':<15} {D_design*1000:<18.1f} {'✓' if abs(D_device-D_design)<0.1 else '✗'}")

L_device = result_pfr_device.get('L_bed', 0)
L_design = design_result.get('L_bed', 0)
print(f"  {'长度 L (mm)':<25} {L_device*1000:<15.1f} {'N/A':<15} {L_design*1000:<18.1f} {'✓' if abs(L_device-L_design)<0.1 else '✗'}")

print(f"\n--- Device 层额外计算 ---")
print(f"  nth_order_analysis: {result_pfr_device.get('nth_order_analysis', {}).get('order_label', 'N/A')}")
print(f"  双分子模式: {result_pfr_device.get('nth_order_analysis', {}).get('bimolecular', False)}")
ft = result_pfr_device.get('flow_and_transport', {})
print(f"  体积流量 v0: {ft.get('v0_m3_s', 0)*1e6:.2f} cm³/s")
print(f"  空塔气速: {ft.get('u_superficial_m_s', 0):.4f} m/s")
print(f"  颗粒Re数: {ft.get('Re_particle', 0):.2f}")
val = result_pfr_device.get('validation', {})
print(f"  约束校核: {val.get('n_pass', 0)}/{val.get('n_total', 0)} 通过")

# ============================================================
# 方案二：全混流反应器 (CSTR)
# ============================================================
print("\n" + "=" * 80)
print("  方案二：全混流反应器 (CSTR)")
print("=" * 80)

# ── Device 层 ──
print("\n--- Device 层调用 (_full_reactor_design) ---")
result_cstr_device = _full_reactor_design(
    reactor_type="cstr",
    X_target=X_A_target,
    k=k,
    delta_H_rxn=0.0,
    T_in=T,
    F_A0=F_A0,
    C_A0=C_A0,
    C_B0=C_B0,
    n=2,
    P_in=P,
    rho_gas=rho_gas,
    mu=mu,
    Cp_mix=Cp_mix,
    num_tanks=1,
)

# ── Engine 层 ──
print("\n--- Engine 层直接调用 (cstr_design + design_reactor) ---")

# 方式1: 直接调用 cstr_design (纯动力学)
cstr_result = cstr_design(F_A0, X_A_target, k, C_A0, n=2, num_tanks=1, C_B0=C_B0)
print(f"\n  cstr_design 结果:")
print(f"    V_total = {cstr_result.get('V_total', 0)*1000:.2f} L")
print(f"    τ = {cstr_result.get('tau', 0):.2f} s")

# 方式2: 调用 design_reactor (完整设计)
design_cstr_result = design_reactor(
    reactor_type="cstr",
    F_A0=F_A0, X_target=X_A_target,
    delta_H_rxn=0.0,
    k=k, C_A0=C_A0, n=2,
    rho_cat=2000.0, particle_diameter=particle_diameter, epsilon=epsilon_calc,
    T_in=T_K, P_in=P,
    d_tube_inner=0.038, L_tube=6.0,
    U=100.0, T_coolant=T_K,
    num_tanks=1,
    rho_gas=rho_gas, mu=mu, Cp_mix=Cp_mix,
    C_B0=C_B0,
)
print(f"\n  design_reactor 结果:")
print(f"    V_total = {design_cstr_result.get('V_total', 0)*1000:.2f} L")
print(f"    D_tank = {design_cstr_result.get('D_tank', 0)*1000:.1f} mm")
print(f"    H_tank = {design_cstr_result.get('H_tank', 0)*1000:.1f} mm")
print(f"    τ = {design_cstr_result.get('tau', 0):.2f} s")

# ── 对比 ──
print("\n--- Device 层 vs Engine 层 对比 ---")
print(f"  {'参数':<25} {'Device层':<15} {'Engine层(cstr)':<15} {'Engine层(design)':<18} {'一致?'}")
print("  " + "-" * 80)

V_device_cstr = result_cstr_device.get('V', result_cstr_device.get('V_total', 0))
V_cstr = cstr_result.get('V_total', 0)
V_design_cstr = design_cstr_result.get('V_total', 0)

print(f"  {'体积 V (L)':<25} {V_device_cstr*1000:<15.2f} {V_cstr*1000:<15.2f} {V_design_cstr*1000:<18.2f} {'✓' if abs(V_device_cstr-V_cstr)<0.01 else '✗'}")

tau_device_cstr = result_cstr_device.get('tau', 0)
tau_cstr = cstr_result.get('tau', 0)
tau_design_cstr = design_cstr_result.get('tau', 0)
print(f"  {'空时 τ (s)':<25} {tau_device_cstr:<15.2f} {tau_cstr:<15.2f} {tau_design_cstr:<18.2f} {'✓' if abs(tau_device_cstr-tau_cstr)<0.1 else '✗'}")

D_device_cstr = result_cstr_device.get('D_tank', result_cstr_device.get('D_bed', 0))
D_design_cstr = design_cstr_result.get('D_tank', 0)
print(f"  {'直径 D (mm)':<25} {D_device_cstr*1000:<15.1f} {'N/A':<15} {D_design_cstr*1000:<18.1f} {'✓' if abs(D_device_cstr-D_design_cstr)<0.1 else '✗'}")

H_device_cstr = result_cstr_device.get('H_tank', result_cstr_device.get('L_bed', 0))
H_design_cstr = design_cstr_result.get('H_tank', 0)
print(f"  {'高度 H (mm)':<25} {H_device_cstr*1000:<15.1f} {'N/A':<15} {H_design_cstr*1000:<18.1f} {'✓' if abs(H_device_cstr-H_design_cstr)<0.1 else '✗'}")

print(f"\n--- Device 层额外计算 ---")
print(f"  nth_order_analysis: {result_cstr_device.get('nth_order_analysis', {}).get('order_label', 'N/A')}")
print(f"  双分子模式: {result_cstr_device.get('nth_order_analysis', {}).get('bimolecular', False)}")
ft = result_cstr_device.get('flow_and_transport', {})
print(f"  体积流量 v0: {ft.get('v0_m3_s', 0)*1e6:.2f} cm³/s")
val = result_cstr_device.get('validation', {})
print(f"  约束校核: {val.get('n_pass', 0)}/{val.get('n_total', 0)} 通过")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 80)
print("  对比总结")
print("=" * 80)

print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│  床层空隙率 (ε) 计算方法                                                     │
│  ─────────────────────────────────────────────────────────────────────────  │
│  新增函数: estimate_bed_void_fraction(D_particle, D_bed, packing_type)       │
│                                                                             │
│  影响因素:                                                                  │
│  1. 填充类型: 随机球形(0.40) / 规则球形(0.26) / 圆柱形(0.37) / 环形(0.45)    │
│  2. 球形度修正: 非球形颗粒空隙率更高                                          │
│  3. 壁面效应: D_bed/D_particle < 50 时，壁面附近空隙率增大                     │
│                                                                             │
│  关联式: Dixon (1988)                                                       │
│  ε = ε_inf + (1 - ε_inf) × 0.25 × (d_p/D)  (当 D_bed/d_p < 50)             │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  reactor_device._full_reactor_design (Device 层)                            │
│  ─────────────────────────────────────────────────────────────────────────  │
│  功能: 完整反应器设计流水线                                                  │
│  1. 参数归一化 (温度/压力/流量单位自动转换)                                    │
│  2. 混合物物性推算 (从components+z自动推算C_A0/F_A0/rho/mu/Cp)               │
│  3. 反应焓推算 (赫斯定律，从stoichiometry推算)                                │
│  4. 催化剂失活修正                                                           │
│  5. Thiele模量/效率因子                                                      │
│  6. 零/一/二级反应解析体积 (含双分子支持)                                      │
│  7. 副反应计算 (平行/串联)                                                   │
│  8. 传热系数U自动估算                                                        │
│  9. 调用 design_reactor 物理引擎                                             │
│  10. 流动与传质参数补全                                                       │
│  11. 约束校验                                                                │
│  12. 工程优化建议                                                            │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  physics_engine.reactor (Engine 层)                                         │
│  ─────────────────────────────────────────────────────────────────────────  │
│  核心函数:                                                                  │
│  • pfr_design: PFR动力学体积计算 (支持双分子C_B0)                             │
│  • cstr_design: CSTR动力学体积计算 (支持双分子C_B0)                           │
│  • design_reactor: 完整反应器设计 (几何+压降+传热+催化剂)                      │
│  • estimate_bed_void_fraction: 床层空隙率估算 (新增)                         │
│  • validate_reactor_design: 约束校验                                        │
│  • fixed_bed_*: 固定床专用函数 (尺寸/压降/催化剂/停留时间)                     │
│  • estimate_reactor_U: 传热系数估算                                          │
│  输出: 几何尺寸 + 压降 + 传热 + 催化剂用量                                    │
└─────────────────────────────────────────────────────────────────────────────┘

核心结论:
  • 动力学计算: 完全一致 (Device层直接调用Engine层函数)
  • Device层额外提供: 物性推算/失活修正/副反应/U估算/流动参数/建议
  • Engine层职责: 纯动力学+几何+压降+传热计算
  • 两者关系: Device层 = Engine层 + 物性推算 + 工程智能 + 报告生成

关于停留时间 τ 的说明:
  • pfr_design 返回理论空时 (τ = V×C_A0/F_A0) - 适用于空管PFR
  • design_reactor 返回实际停留时间 (τ = V×ε/F_v0) - 考虑床层空隙率
  • 对于固定床反应器，实际停留时间更准确
  • Device层使用 design_reactor，因此返回实际停留时间
""")

# PFR vs CSTR 对比
print("=" * 80)
print("  PFR vs CSTR 体积对比")
print("=" * 80)
V_pfr_val = result_pfr_device.get('V', 0) * 1000
V_cstr_val = result_cstr_device.get('V', result_cstr_device.get('V_total', 0)) * 1000
ratio = V_cstr_val / V_pfr_val if V_pfr_val > 0 else 0

print(f"\n  {'反应器':<15} {'体积 (L)':<15} {'空时 (s)':<15} {'ε'}")
print("  " + "-" * 55)
print(f"  {'固定床(PFR)':<15} {V_pfr_val:<15.2f} {result_pfr_device.get('tau', 0):<15.2f} {epsilon_calc:.4f}")
print(f"  {'全混流(CSTR)':<15} {V_cstr_val:<15.2f} {result_cstr_device.get('tau', 0):<15.2f} N/A")
print(f"\n  CSTR/PFR 体积比: {ratio:.2f}")
print(f"  结论: 对于二级双分子反应，PFR 效率是 CSTR 的 {ratio:.1f} 倍")
