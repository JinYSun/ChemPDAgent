# -*- coding: utf-8 -*-
"""
换热器 Device vs Engine 对比
热侧: 热水 360K -> 350K, m=1 kg/s
冷侧: 冷水 300K -> 320K
"""
import sys, json
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, r"c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents")

from device_tools.heatexchanger_device import _full_he_design
from physics_engine.heatexchanger import design_heat_exchanger
from physics_engine.thermo_helper import get_fluid_Cp, get_fluid_density, get_fluid_viscosity, get_fluid_thermal_conductivity, get_molecular_weight

print("=" * 80)
print("  管壳式换热器设计对比: Device层 vs Engine层")
print("=" * 80)
print("\n设计条件:")
print("  热侧: 热水 360 K (87 C) -> 350 K (77 C), m = 1.0 kg/s")
print("  冷侧: 冷水 300 K (27 C) -> 320 K (47 C)")
print("  压力: 常压 (1 atm)")

# ============================================================
# Device层设计
# ============================================================
print("\n" + "=" * 80)
print("  [1] Device层 (_full_he_design)")
print("=" * 80)

device_result = _full_he_design(
    hot_fluid="water", cold_fluid="water",
    T_hot_in=360.0, T_hot_out=350.0, m_hot=1.0, P_hot=101325.0,
    T_cold_in=300.0, T_cold_out=320.0, P_cold=101325.0,
)

ds = device_result.get('design_summary', {})
hb = device_result.get('heat_balance', {})
entu = device_result.get('epsilon_ntu', {})
udec = device_result.get('U_decomposition', {})
pd = device_result.get('pressure_drop', {})
ht = device_result.get('heat_transfer', {})
prop = device_result.get('properties', {})
val = device_result.get('validation', {})
isc = device_result.get('industry_standard_checks', {})

print("\n--- 热负荷 ---")
print(f"  热负荷 Q = {ds.get('Q_kW', 0):.2f} kW")
print(f"  热侧放热 = {hb.get('Q_hot_kW', 0):.2f} kW")
print(f"  冷侧吸热 = {hb.get('Q_cold_kW', 0):.2f} kW")
print(f"  热平衡偏差 = {hb.get('imbalance_percent', 0):.2f}%")

print("\n--- 传热参数 ---")
print(f"  LMTD (raw) = {ds.get('LMTD_raw_K', 0):.2f} K")
print(f"  LMTD (修正) = {ds.get('LMTD_K', 0):.2f} K")
print(f"  F 修正系数 = {ds.get('F_factor', 0):.4f} ({ds.get('F_quality', '')})")
print(f"  总传热系数 U = {ds.get('U_W_m2K', 0):.1f} W/(m2.K)")

print("\n--- 换热面积 ---")
print(f"  所需面积 A_req = {ds.get('A_required_m2', 0):.2f} m2")
print(f"  设计面积 A_design = {ds.get('A_design_m2', 0):.2f} m2")
print(f"  裕量 = {ds.get('overspec_percent', 0):.1f}%")

print("\n--- 结构参数 ---")
print(f"  管数 N_tubes = {ds.get('N_tubes', 0)}")
print(f"  管程数 N_pass = {ds.get('N_tube_passes', 0)}")
print(f"  壳径 D_shell = {ds.get('D_shell_m', 0):.3f} m = {ds.get('D_shell_m', 0)*1000:.0f} mm")
print(f"  管长 L_tube = {ds.get('L_tube_m', 0):.2f} m")
print(f"  管外径 d_o = {ds.get('d_o_m', 0)*1000:.0f} mm")
print(f"  管内径 d_i = {ds.get('d_i_m', 0)*1000:.3f} mm")
print(f"  管间距 pitch = {ds.get('pitch_m', 0)*1000:.1f} mm")
print(f"  长径比 L/D = {ds.get('L_D_ratio', 0):.1f}")
print(f"  管排列 = {ds.get('tube_layout', '')}")
print(f"  折流板数 = {ds.get('N_baffles', 0)}")
print(f"  折流板间距 = {ds.get('baffle_spacing_m', 0):.3f} m")
print(f"  折流板切口 = {ds.get('baffle_cut', 0)*100:.0f}%")

print("\n--- 压降 ---")
print(f"  管程压降 DP_tube = {pd.get('tube_Pa', 0)/1000:.2f} kPa")
print(f"  壳程压降 DP_shell = {pd.get('shell_Pa', 0)/1000:.2f} kPa")

print("\n--- 流速 ---")
print(f"  管程流速 = {ht.get('u_tube_m_s', 0):.3f} m/s")
print(f"  壳程流速 = {ht.get('u_shell_m_s', 0):.3f} m/s")
print(f"  管程 Re = {ht.get('Re_tube', 0):.0f} ({ht.get('regime_tube', '')})")
print(f"  壳程 Re = {ht.get('Re_shell', 0):.0f} ({ht.get('regime_shell', '')})")

print("\n--- e-NTU 校核 ---")
print(f"  epsilon = {entu.get('epsilon', 0):.4f}")
print(f"  NTU = {entu.get('NTU', 0):.3f}")
print(f"  C_r = {entu.get('C_r', 0):.3f}")
print(f"  Q_max = {entu.get('Q_max_kW', 0):.1f} kW")
print(f"  Q_actual = {entu.get('Q_actual_kW', 0):.1f} kW")

print("\n--- U 分解 ---")
print(f"  U_overall = {udec.get('U_overall', 0):.1f} W/(m2.K)")
print(f"  R_total = {udec.get('R_total', 0):.6f} m2.K/W")
print(f"  R_conv_tube = {udec.get('R_conv_tube', 0):.6f} m2.K/W")
print(f"  R_conv_shell = {udec.get('R_conv_shell', 0):.6f} m2.K/W")
print(f"  R_wall = {udec.get('R_wall', 0):.6f} m2.K/W")
print(f"  R_fouling_tube = {udec.get('R_fouling_tube', 0):.6f} m2.K/W")
print(f"  R_fouling_shell = {udec.get('R_fouling_shell', 0):.6f} m2.K/W")
print(f"  污垢占比 = {udec.get('fouling_fraction_pct', ds.get('fouling_fraction', 0)):.1f}%")

print("\n--- 流体物性 ---")
print(f"  热侧: rho={prop.get('rho_hot_kg_m3', 0):.1f} kg/m3, mu={prop.get('mu_hot_Pa_s', 0):.6f} Pa.s, Cp={prop.get('Cp_hot_J_kgK', 0):.1f} J/(kg.K), k={prop.get('k_hot_W_mK', 0):.4f} W/(m.K)")
print(f"  冷侧: rho={prop.get('rho_cold_kg_m3', 0):.1f} kg/m3, mu={prop.get('mu_cold_Pa_s', 0):.6f} Pa.s, Cp={prop.get('Cp_cold_J_kgK', 0):.1f} J/(kg.K), k={prop.get('k_cold_W_mK', 0):.4f} W/(m.K)")

print("\n--- 基础校核结果 ---")
checks = val.get('checks', {})
pass_count = val.get('n_pass', 0)
total_count = val.get('n_total', 0)
pass_rate = val.get('pass_rate', 0)
for name, c in checks.items():
    status = "PASS" if c.get('pass') else "FAIL"
    print(f"  [{status}] {name}: {c.get('value', 0):.2f} ({c.get('criterion', '')})")
print(f"  通过率: {pass_count}/{total_count} = {pass_rate*100:.0f}%")

print("\n--- 行业标准校核 (GB/T 151, SH/T 3121, TEMA) ---")
if isc:
    tv = isc.get('tube_velocity', {})
    s = "PASS" if tv.get('pass') else "FAIL"
    print(f"  [{s}] 管程流速: {tv.get('value_m_s', 0):.3f} m/s (范围: {tv.get('min_m_s', 0)}-{tv.get('max_m_s', 0)} m/s)")
    if tv.get('warning'):
        print(f"      ! {tv['warning']}")
    
    sv = isc.get('shell_velocity', {})
    s = "PASS" if sv.get('pass') else "FAIL"
    print(f"  [{s}] 壳程流速: {sv.get('value_m_s', 0):.3f} m/s (范围: {sv.get('min_m_s', 0)}-{sv.get('max_m_s', 0)} m/s)")
    if sv.get('warning'):
        print(f"      ! {sv['warning']}")
    
    vib = isc.get('vibration_risk', {})
    s = "PASS" if vib.get('pass') else "FAIL"
    print(f"  [{s}] 管束振动风险: {vib.get('risk_level', '')} (壳程流速: {vib.get('shell_velocity_m_s', 0):.3f} m/s)")
    
    te = isc.get('thermal_expansion', {})
    s = "PASS" if te.get('pass') else "FAIL"
    print(f"  [{s}] 热膨胀差: DL = {te.get('delta_L_mm', 0):.2f} mm (允许: {te.get('allowable_mm', 0)} mm)")

warnings = device_result.get('warnings', [])
if warnings:
    print("\n--- 迭代日志/警告 ---")
    for w in warnings:
        print(f"  * {w}")

# ============================================================
# Engine层设计
# ============================================================
print("\n" + "=" * 80)
print("  [2] Engine层 (design_heat_exchanger)")
print("=" * 80)

T_hot_avg = (360 + 350) / 2
T_cold_avg = (300 + 320) / 2

Cp_h_molar = get_fluid_Cp("water", T=T_hot_avg)
Cp_c_molar = get_fluid_Cp("water", T=T_cold_avg)
MW_h = get_molecular_weight("water")
MW_c = get_molecular_weight("water")
Cp_h = Cp_h_molar / (MW_h / 1000.0)
Cp_c = Cp_c_molar / (MW_c / 1000.0)

rho_h = get_fluid_density("water", T=T_hot_avg)
rho_c = get_fluid_density("water", T=T_cold_avg)
mu_h = get_fluid_viscosity("water", T=T_hot_avg)
mu_c = get_fluid_viscosity("water", T=T_cold_avg)
k_h = get_fluid_thermal_conductivity("water", T=T_hot_avg)
k_c = get_fluid_thermal_conductivity("water", T=T_cold_avg)

Q_hot = abs(Cp_h * 1.0 * (360 - 350))
delta_T_cold = abs(320 - 300)
m_cold = Q_hot / (Cp_c * delta_T_cold)

print(f"\n物性参数 (Engine层手动获取):")
print(f"  热侧 (T_avg={T_hot_avg}K): Cp={Cp_h:.1f} J/(kg.K), rho={rho_h:.1f} kg/m3, mu={mu_h:.6f} Pa.s, k={k_h:.4f} W/(m.K)")
print(f"  冷侧 (T_avg={T_cold_avg}K): Cp={Cp_c:.1f} J/(kg.K), rho={rho_c:.1f} kg/m3, mu={mu_c:.6f} Pa.s, k={k_c:.4f} W/(m.K)")
print(f"  冷侧估算流量: m_cold = {m_cold:.4f} kg/s")

engine_result = design_heat_exchanger(
    T_hot_in=360.0, T_hot_out=350.0,
    T_cold_in=300.0, T_cold_out=320.0,
    m_hot=1.0, m_cold=m_cold,
    Cp_hot=Cp_h, Cp_cold=Cp_c,
    rho_hot=rho_h, rho_cold=rho_c,
    mu_hot=mu_h, mu_cold=mu_c,
    k_hot=k_h, k_cold=k_c,
)

print("\n--- Engine层设计结果 ---")
print(f"  Q = {engine_result.get('Q', 0)/1000:.2f} kW")
print(f"  LMTD = {engine_result.get('LMTD', 0):.2f} K")
print(f"  F = {engine_result.get('F', 0):.4f}")
print(f"  U = {engine_result.get('U', 0):.1f} W/(m2.K)")
print(f"  A_required = {engine_result.get('A_required', 0):.2f} m2")
print(f"  A_design = {engine_result.get('A_design', 0):.2f} m2")
print(f"  N_tubes = {engine_result.get('N_tubes', 0)}")
print(f"  tube_passes = {engine_result.get('tube_passes', 0)}")
print(f"  D_shell = {engine_result.get('D_shell', 0):.3f} m = {engine_result.get('D_shell', 0)*1000:.0f} mm")
print(f"  L/D = {engine_result.get('length_diameter_ratio', 0):.1f}")
print(f"  baffle_spacing = {engine_result.get('baffle_spacing_m', 0):.3f} m")
print(f"  DP_tube = {engine_result.get('delta_P_tube', 0)/1000:.2f} kPa")
print(f"  DP_shell = {engine_result.get('delta_P_shell', 0)/1000:.2f} kPa")

h_tube_detail = engine_result.get('h_tube_detail', {})
h_shell_detail = engine_result.get('h_shell_detail', {})
print(f"  u_tube = {h_tube_detail.get('velocity', 0):.3f} m/s")
print(f"  u_shell = {h_shell_detail.get('velocity', 0):.3f} m/s")
print(f"  Re_tube = {h_tube_detail.get('Re', 0):.0f}")
print(f"  Re_shell = {h_shell_detail.get('Re', 0):.0f}")

# ============================================================
# 对比总结
# ============================================================
print("\n" + "=" * 80)
print("  [3] Device vs Engine 对比")
print("=" * 80)

comparisons = [
    ("Q (kW)", ds.get('Q_kW', 0), engine_result.get('Q', 0)/1000),
    ("LMTD (K)", ds.get('LMTD_K', 0), engine_result.get('LMTD', 0)),
    ("F factor", ds.get('F_factor', 0), engine_result.get('F', 0)),
    ("U (W/m2.K)", ds.get('U_W_m2K', 0), engine_result.get('U', 0)),
    ("A_req (m2)", ds.get('A_required_m2', 0), engine_result.get('A_required', 0)),
    ("A_design (m2)", ds.get('A_design_m2', 0), engine_result.get('A_design', 0)),
    ("N_tubes", ds.get('N_tubes', 0), engine_result.get('N_tubes', 0)),
    ("N_pass", ds.get('N_tube_passes', 0), engine_result.get('tube_passes', 0)),
    ("D_shell (m)", ds.get('D_shell_m', 0), engine_result.get('D_shell', 0)),
    ("L/D", ds.get('L_D_ratio', 0), engine_result.get('length_diameter_ratio', 0)),
    ("u_tube (m/s)", ht.get('u_tube_m_s', 0), h_tube_detail.get('velocity', 0)),
    ("u_shell (m/s)", ht.get('u_shell_m_s', 0), h_shell_detail.get('velocity', 0)),
    ("DP_tube (kPa)", pd.get('tube_Pa', 0)/1000, engine_result.get('delta_P_tube', 0)/1000),
    ("DP_shell (kPa)", pd.get('shell_Pa', 0)/1000, engine_result.get('delta_P_shell', 0)/1000),
]

print(f"\n{'参数':<20} {'Device层':>15} {'Engine层':>15} {'偏差':>10}")
print("-" * 65)
for name, dev_val, eng_val in comparisons:
    if eng_val != 0:
        diff_pct = abs(dev_val - eng_val) / abs(eng_val) * 100
    elif dev_val != 0:
        diff_pct = 100.0
    else:
        diff_pct = 0.0
    match = "OK" if diff_pct < 1.0 else "DIFF"
    print(f"  {name:<18} {dev_val:>15.3f} {eng_val:>15.3f} {diff_pct:>8.2f}% {match}")

print("\n" + "=" * 80)
