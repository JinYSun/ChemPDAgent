# -*- coding: utf-8 -*-
"""管壳式换热器设计 - 液-液换热 (无相变)"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.heatexchanger_device import _full_he_design

print("=" * 80)
print("  管壳式换热器设计 (_full_he_design)")
print("=" * 80)

print("\n设计条件:")
print("  热侧: 热水 350 K (77°C) → 330 K (57°C), m = 1.0 kg/s")
print("  冷侧: 冷水 300 K (27°C) → 320 K (47°C)")
print("  压力: 常压 (1 atm)")

# ============================================================
# 调用换热器设计函数
# ============================================================
result = _full_he_design(
    hot_fluid="water",
    cold_fluid="water",
    T_hot_in=350.0,
    T_hot_out=330.0,
    m_hot=1.0,
    P_hot=101325.0,
    T_cold_in=300.0,
    T_cold_out=320.0,
    P_cold=101325.0,
)

# 从嵌套字典中提取数据
ds = result.get('design_summary', {})
hb = result.get('heat_balance', {})
entu = result.get('epsilon_ntu', {})
udec = result.get('U_decomposition', {})
pd = result.get('pressure_drop', {})
ht = result.get('heat_transfer', {})
prop = result.get('properties', {})
val = result.get('validation', {})

print("\n" + "=" * 80)
print("  设计结果")
print("=" * 80)

print(f"\n--- 热负荷 ---")
print(f"  热负荷 Q = {ds.get('Q_kW', 0):.2f} kW")
print(f"  热侧放热 = {hb.get('Q_hot_kW', 0):.2f} kW")
print(f"  冷侧吸热 = {hb.get('Q_cold_kW', 0):.2f} kW")
print(f"  热平衡偏差 = {hb.get('imbalance_pct', 0):.2f}%")

print(f"\n--- 传热参数 ---")
print(f"  LMTD (raw) = {ds.get('LMTD_raw_K', 0):.2f} K")
print(f"  LMTD (修正) = {ds.get('LMTD_K', 0):.2f} K")
print(f"  F 修正系数 = {ds.get('F_factor', 0):.4f} ({ds.get('F_quality', 'N/A')})")
print(f"  总传热系数 U = {ds.get('U_W_m2K', 0):.1f} W/(m²·K)")

print(f"\n--- 换热面积 ---")
print(f"  所需面积 A_req = {ds.get('A_required_m2', 0):.2f} m²")
print(f"  设计面积 A_design = {ds.get('A_design_m2', 0):.2f} m²")
print(f"  裕量 = {ds.get('overspec_percent', 0):.1f}%")

print(f"\n--- 结构参数 ---")
print(f"  管数 N_tubes = {ds.get('N_tubes', 0)}")
print(f"  管程数 N_pass = {ds.get('N_tube_passes', 0)}")
print(f"  壳径 D_shell = {ds.get('D_shell_m', 0):.3f} m = {ds.get('D_shell_m', 0)*1000:.0f} mm")
print(f"  管长 L_tube = {ds.get('L_tube_m', 0):.2f} m")
print(f"  管外径 d_o = {ds.get('d_o_m', 0)*1000:.0f} mm")
print(f"  管内径 d_i = {ds.get('d_i_m', 0)*1000:.3f} mm")
print(f"  管间距 pitch = {ds.get('pitch_m', 0)*1000:.1f} mm")
print(f"  长径比 L/D = {ds.get('L_D_ratio', 0):.1f}")
print(f"  管排列 = {ds.get('tube_layout', 'N/A')}")
print(f"  折流板数 = {ds.get('N_baffles', 0)}")
print(f"  折流板间距 = {ds.get('baffle_spacing_m', 0):.3f} m")
print(f"  折流板切口 = {ds.get('baffle_cut', 0)*100:.0f}%")

print(f"\n--- 压降 ---")
print(f"  管程压降 ΔP_tube = {pd.get('tube_Pa', 0)/1000:.2f} kPa")
print(f"  壳程压降 ΔP_shell = {pd.get('shell_Pa', 0)/1000:.2f} kPa")

print(f"\n--- 流速 ---")
print(f"  管程流速 = {ht.get('u_tube_m_s', 0):.3f} m/s")
print(f"  壳程流速 = {ht.get('u_shell_m_s', 0):.3f} m/s")
print(f"  管程 Re = {ht.get('Re_tube', 0):.0f} ({ht.get('tube_flow_regime', 'N/A')})")
print(f"  壳程 Re = {ht.get('Re_shell', 0):.0f} ({ht.get('shell_flow_regime', 'N/A')})")

print(f"\n--- ε-NTU 校核 ---")
print(f"  ε = {entu.get('epsilon', 0):.4f}")
print(f"  NTU = {entu.get('NTU', 0):.3f}")
print(f"  C_r = {entu.get('C_r', 0):.3f}")
print(f"  Q_max = {entu.get('Q_max_kW', 0):.1f} kW")
print(f"  Q_actual = {entu.get('Q_actual_kW', 0):.1f} kW")

print(f"\n--- U 分解 ---")
print(f"  U_overall = {udec.get('U_overall', 0):.1f} W/(m²·K)")
print(f"  R_total = {udec.get('R_total', 0):.6f} m²·K/W")
print(f"  R_conv_tube = {udec.get('R_conv_tube', 0):.6f} m²·K/W")
print(f"  R_conv_shell = {udec.get('R_conv_shell', 0):.6f} m²·K/W")
print(f"  R_wall = {udec.get('R_wall', 0):.6f} m²·K/W")
print(f"  R_fouling_tube = {udec.get('R_fouling_tube', 0):.6f} m²·K/W")
print(f"  R_fouling_shell = {udec.get('R_fouling_shell', 0):.6f} m²·K/W")
print(f"  污垢占比 = {udec.get('fouling_fraction_pct', 0):.1f}%")

print(f"\n--- 流体物性 ---")
print(f"  热侧: ρ={prop.get('rho_hot_kg_m3',0):.1f} kg/m³, μ={prop.get('mu_hot_Pas',0):.6f} Pa·s, Cp={prop.get('Cp_hot_J_kgK',0):.1f} J/(kg·K), k={prop.get('k_hot_W_mK',0):.4f} W/(m·K)")
print(f"  冷侧: ρ={prop.get('rho_cold_kg_m3',0):.1f} kg/m³, μ={prop.get('mu_cold_Pas',0):.6f} Pa·s, Cp={prop.get('Cp_cold_J_kgK',0):.1f} J/(kg·K), k={prop.get('k_cold_W_mK',0):.4f} W/(m·K)")

print(f"\n--- 校核结果 ---")
checks = val.get('checks', {})
for name, check in checks.items():
    status = "✓" if check.get('pass', False) else "✗"
    print(f"  {status} {name}: {check.get('value', 0):.2f} ({check.get('criterion', 'N/A')})")
print(f"  通过率: {val.get('n_pass', 0)}/{val.get('n_total', 0)} = {val.get('pass_rate', 0)*100:.0f}%")

print(f"\n--- 警告 ---")
for w in result.get('warnings', []):
    print(f"  ⚠ {w}")

# 行业标准校核
print(f"\n--- 行业标准校核 (GB/T 151, SH/T 3121, TEMA) ---")
isc = result.get('industry_standard_checks', {})

# 管程流速
tv = isc.get('tube_velocity', {})
status = "✓" if tv.get('pass') else "✗"
print(f"  {status} 管程流速: {tv.get('value_m_s', 0):.3f} m/s (范围: {tv.get('min_m_s', 0)}-{tv.get('max_m_s', 0)} m/s, 最优: {tv.get('optimal_m_s', 0)} m/s)")
if tv.get('warning'):
    print(f"      ⚠ {tv['warning']}")

# 壳程流速
sv = isc.get('shell_velocity', {})
status = "✓" if sv.get('pass') else "✗"
print(f"  {status} 壳程流速: {sv.get('value_m_s', 0):.3f} m/s (范围: {sv.get('min_m_s', 0)}-{sv.get('max_m_s', 0)} m/s, 最优: {sv.get('optimal_m_s', 0)} m/s)")
if sv.get('warning'):
    print(f"      ⚠ {sv['warning']}")

# 折流板间距比
br = isc.get('baffle_spacing_ratio', {})
if br:
    status = "✓" if br.get('pass') else "✗"
    print(f"  {status} 折流板间距比 B/D: {br.get('B_D_ratio', 0):.3f} (范围: {br.get('min', 0)}-{br.get('max', 0)}, 最优: {br.get('optimal_range', 'N/A')})")
    if br.get('warning'):
        print(f"      ⚠ {br['warning']}")

# 管束振动
vb = isc.get('tube_bundle_vibration', {})
if vb:
    status = "✓" if vb.get('pass') else "✗"
    print(f"  {status} 管束振动风险: {vb.get('risk_level', 'N/A')} (壳程流速: {vb.get('u_shell_m_s', 0):.3f} m/s)")
    if vb.get('warning'):
        print(f"      ⚠ {vb['warning']}")

# 热膨胀
te = isc.get('thermal_expansion', {})
if te:
    status = "✓" if te.get('pass') else "⚠"
    print(f"  {status} 热膨胀差: ΔL = {te.get('delta_L_mm', 0):.2f} mm (允许: {te.get('expansion_joint_max_mm', 0)} mm)")
    print(f"      管壳温差: {te.get('delta_T_K', 0):.1f} K, 管长: {te.get('L_tube_m', 0):.2f} m, α = {te.get('alpha_1_K', 0):.2e} 1/K")
    if te.get('warning'):
        print(f"      ⚠ {te['warning']}")

print("\n" + "=" * 80)
print("  设计总结")
print("=" * 80)
print(f"""
┌─────────────────────────────────────────────────────────────────────────────┐
│  管壳式换热器设计结果                                                       │
│  ─────────────────────────────────────────────────────────────────────────  │
│  操作条件:                                                                  │
│    - 热侧: 热水 77°C → 57°C (常压), m = 1.0 kg/s                          │
│    - 冷侧: 冷水 27°C → 47°C (常压), m = {hb.get('Q_cold_kW',0)/4180/20:.2f} kg/s                          │
│                                                                             │
│  热负荷: {ds.get('Q_kW', 0):.1f} kW                                                         │
│  换热面积: {ds.get('A_design_m2', 0):.2f} m² (裕量 {ds.get('overspec_percent', 0):.0f}%)                           │
│  结构: {ds.get('N_tubes', 0)}管 × {ds.get('N_tube_passes', 0)}管程, 壳径 {ds.get('D_shell_m', 0)*1000:.0f}mm, 管长 {ds.get('L_tube_m', 0):.0f}m              │
│  校核: {val.get('n_pass', 0)}/{val.get('n_total', 0)} 通过 ({val.get('pass_rate', 0)*100:.0f}%)                                                    │
│                                                                             │
│  设计特点:                                                                  │
│    - 液-液换热，无相变                                                      │
│    - 温度均低于 100°C，常压操作安全                                         │
│    - 三角排列，换热效率高                                                   │
└─────────────────────────────────────────────────────────────────────────────┘
""")
