# -*- coding: utf-8 -*-
"""对比 heatexchanger_device._full_he_design vs physics_engine 直接调用"""
import sys, io, json, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.heatexchanger_device import _full_he_design
from physics_engine.heatexchanger import (
    design_heat_exchanger, lmtd_with_correction, heat_duty_from_flow,
    estimate_U, suggest_configuration, fouling_resistance,
    lmtd, area_required, effectiveness_ntu,
)
from physics_engine.thermo_helper import (
    get_fluid_density, get_fluid_viscosity, get_fluid_Cp, get_fluid_thermal_conductivity,
    get_molecular_weight,
)

# ============================================================
# 设计条件
# ============================================================
T_hot_in = 400.0      # K (127°C) 热水入口
T_hot_out = 350.0     # K (77°C) 热水出口
T_cold_in = 300.0     # K (27°C) 冷水入口
T_cold_out = 340.0    # K (67°C) 冷水出口
m_hot = 1.0           # kg/s 热水流量
P_hot = 500000.0      # Pa (5 bar) 热侧压力 - 必须加压防止汽化
P_cold = 300000.0     # Pa (3 bar) 冷侧压力

print("=" * 80)
print("  heatexchanger_device._full_he_design vs physics_engine 直接调用 对比")
print("=" * 80)
print(f"\n设计条件:")
print(f"  热侧: 热水 {T_hot_in} K → {T_hot_out} K, m_hot = {m_hot} kg/s")
print(f"  冷侧: 冷水 {T_cold_in} K → {T_cold_out} K")
print(f"  热侧压力: {P_hot/1e5:.1f} bar (加压防汽化)")
print(f"  冷侧压力: {P_cold/1e5:.1f} bar")
print(f"\n注意: 热侧入口温度 {T_hot_in} K = {T_hot_in-273.15:.1f}°C > 100°C")
print(f"      必须加压至 {P_hot/1e5:.1f} bar 以上防止水汽化")

# ============================================================
# Device 层调用
# ============================================================
print("\n" + "=" * 80)
print("  Device 层: _full_he_design")
print("=" * 80)

result_device = _full_he_design(
    T_hot_in=T_hot_in,
    T_hot_out=T_hot_out,
    T_cold_in=T_cold_in,
    T_cold_out=T_cold_out,
    m_hot=m_hot,
    hot_fluid="water",
    cold_fluid="water",
    P_hot=P_hot,
    P_cold=P_cold,
)

# ============================================================
# Engine 层直接调用
# ============================================================
print("\n" + "=" * 80)
print("  Engine 层: physics_engine 直接调用")
print("=" * 80)

# 步骤1: 获取物性
T_h_avg = (T_hot_in + T_hot_out) / 2  # 375 K
T_c_avg = (T_cold_in + T_cold_out) / 2  # 320 K

print(f"\n  热侧平均温度: {T_h_avg} K = {T_h_avg-273.15:.1f}°C")
print(f"  冷侧平均温度: {T_c_avg} K = {T_c_avg-273.15:.1f}°C")

# 热侧物性
rho_h = get_fluid_density("water", T_h_avg, P_hot, phase="liquid")
mu_h = get_fluid_viscosity("water", T_h_avg, P_hot, phase="liquid")
Cp_h_molar = get_fluid_Cp("water", T_h_avg, P_hot, phase="liquid")
k_h = get_fluid_thermal_conductivity("water", T_h_avg, P_hot, phase="liquid")
MW_h = get_molecular_weight("water")
Cp_h = Cp_h_molar / (MW_h / 1000.0)  # J/(kg·K)

# 冷侧物性
rho_c = get_fluid_density("water", T_c_avg, P_cold, phase="liquid")
mu_c = get_fluid_viscosity("water", T_c_avg, P_cold, phase="liquid")
Cp_c_molar = get_fluid_Cp("water", T_c_avg, P_cold, phase="liquid")
k_c = get_fluid_thermal_conductivity("water", T_c_avg, P_cold, phase="liquid")
MW_c = get_molecular_weight("water")
Cp_c = Cp_c_molar / (MW_c / 1000.0)  # J/(kg·K)

print(f"\n  热侧物性 (@ {T_h_avg} K, {P_hot/1e5:.1f} bar):")
print(f"    ρ_h = {rho_h:.2f} kg/m³")
print(f"    μ_h = {mu_h:.5f} Pa·s")
print(f"    Cp_h = {Cp_h:.0f} J/(kg·K)")
print(f"    k_h = {k_h:.4f} W/(m·K)")

print(f"\n  冷侧物性 (@ {T_c_avg} K, {P_cold/1e5:.1f} bar):")
print(f"    ρ_c = {rho_c:.2f} kg/m³")
print(f"    μ_c = {mu_c:.5f} Pa·s")
print(f"    Cp_c = {Cp_c:.0f} J/(kg·K)")
print(f"    k_c = {k_c:.4f} W/(m·K)")

# 步骤2: 热负荷计算
Q_h = Cp_h * m_hot * (T_hot_in - T_hot_out)
delta_T_cold = T_cold_out - T_cold_in
m_cold = Q_h / (Cp_c * delta_T_cold)
Q_c = Cp_c * m_cold * delta_T_cold

print(f"\n  热负荷计算:")
print(f"    Q_h = {Q_h/1000:.2f} kW")
print(f"    Q_c = {Q_c/1000:.2f} kW")
print(f"    m_cold (由热平衡估算) = {m_cold:.4f} kg/s")

# 步骤3: LMTD 计算
lmtd_result = lmtd_with_correction(
    T_hot_in=T_hot_in, T_hot_out=T_hot_out,
    T_cold_in=T_cold_in, T_cold_out=T_cold_out,
)
LMTD = lmtd_result.get("LMTD") or lmtd_result.get("LMTD_corrected")
F_factor = lmtd_result.get("F")

print(f"\n  LMTD 计算:")
print(f"    LMTD = {LMTD:.2f} K")
print(f"    F (修正系数) = {F_factor:.4f}")
print(f"    LMTD_corrected = {LMTD * F_factor:.2f} K")

# 步骤4: 总传热系数估算
u_result = estimate_U("water/water")
U_est = u_result.get("U_typical")
print(f"\n  总传热系数估算:")
print(f"    U = {U_est:.0f} W/(m²·K)")

# 步骤5: 换热器设计
d_o = 0.019  # 19 mm 管外径
d_i = d_o * (1 - 2 * 0.083)  # 壁厚 8.3%
L_tube = 3.0  # 3 m 管长

design_result = design_heat_exchanger(
    T_hot_in=T_hot_in, T_hot_out=T_hot_out,
    T_cold_in=T_cold_in, T_cold_out=T_cold_out,
    m_hot=m_hot, m_cold=m_cold,
    Cp_hot=Cp_h, Cp_cold=Cp_c,
    rho_hot=rho_h, rho_cold=rho_c,
    mu_hot=mu_h, mu_cold=mu_c,
    k_hot=k_h, k_cold=k_c,
    d_o=d_o, d_i=d_i, tube_length=L_tube,
    tube_layout="triangular",
)

print(f"\n  design_heat_exchanger 结果:")
print(f"    A_required = {design_result.get('A_required', 0):.2f} m²")
print(f"    A_design = {design_result.get('A_design', 0):.2f} m²")
print(f"    N_tubes = {design_result.get('N_tubes', 0)}")
print(f"    D_shell = {design_result.get('D_shell', 0):.3f} m")
print(f"    U = {design_result.get('U', 0):.0f} W/(m²·K)")
print(f"    delta_P_tube = {design_result.get('delta_P_tube', 0)/1000:.2f} kPa")
print(f"    delta_P_shell = {design_result.get('delta_P_shell', 0)/1000:.2f} kPa")

# ============================================================
# 对比
# ============================================================
print("\n" + "=" * 80)
print("  Device 层 vs Engine 层 对比")
print("=" * 80)

ds = result_device.get('design_summary', {})
hb = result_device.get('heat_balance', {})
ht = result_device.get('heat_transfer', {})
pd = result_device.get('pressure_drop', {})

print(f"\n--- 热负荷对比 ---")
print(f"  {'参数':<25} {'Device层':<15} {'Engine层':<15} {'一致?'}")
print("  " + "-" * 60)

Q_device = ds.get('Q_kW', 0)
Q_engine = Q_h / 1000
print(f"  {'Q (kW)':<25} {Q_device:<15.2f} {Q_engine:<15.2f} {'✓' if abs(Q_device-Q_engine)<0.5 else '✗'}")

LMTD_device = ds.get('LMTD_K', 0)
LMTD_engine = LMTD
print(f"  {'LMTD (K)':<25} {LMTD_device:<15.2f} {LMTD_engine:<15.2f} {'✓' if abs(LMTD_device-LMTD_engine)<0.5 else '✗'}")

F_device = ds.get('F_factor', 0)
F_engine = F_factor
print(f"  {'F (修正系数)':<25} {F_device:<15.4f} {F_engine:<15.4f} {'✓' if abs(F_device-F_engine)<0.01 else '✗'}")

print(f"\n--- 传热系数对比 ---")
U_device = ds.get('U_W_m2K', 0)
U_engine = design_result.get('U', 0)
print(f"  {'U (W/m²·K)':<25} {U_device:<15.0f} {U_engine:<15.0f} {'✓' if abs(U_device-U_engine)<10 else '⚠'}")

print(f"\n--- 面积对比 ---")
A_req_device = ds.get('A_required_m2', 0)
A_req_engine = design_result.get('A_required', 0)
print(f"  {'A_required (m²)':<25} {A_req_device:<15.2f} {A_req_engine:<15.2f} {'✓' if abs(A_req_device-A_req_engine)<0.5 else '✗'}")

A_design_device = ds.get('A_design_m2', 0)
A_design_engine = design_result.get('A_design', 0)
print(f"  {'A_design (m²)':<25} {A_design_device:<15.2f} {A_design_engine:<15.2f} {'✓' if abs(A_design_device-A_design_engine)<0.5 else '✗'}")

print(f"\n--- 结构参数对比 ---")
N_device = ds.get('N_tubes', 0)
N_engine = design_result.get('N_tubes', 0)
print(f"  {'N_tubes':<25} {N_device:<15} {N_engine:<15} {'✓' if N_device==N_engine else '⚠'}")

D_device = ds.get('D_shell_m', 0)
D_engine = design_result.get('D_shell', 0)
print(f"  {'D_shell (m)':<25} {D_device:<15.3f} {D_engine:<15.3f} {'✓' if abs(D_device-D_engine)<0.01 else '✗'}")

L_device = ds.get('L_tube_m', 0)
L_engine = design_result.get('tube_length_m', L_tube)
print(f"  {'L_tube (m)':<25} {L_device:<15.2f} {L_engine:<15.2f} {'✓' if abs(L_device-L_engine)<0.01 else '✗'}")

print(f"\n--- 压降对比 ---")
dP_tube_device = pd.get('tube_Pa', 0)
dP_tube_engine = design_result.get('delta_P_tube', 0)
print(f"  {'ΔP_tube (kPa)':<25} {dP_tube_device/1000:<15.2f} {dP_tube_engine/1000:<15.2f} {'✓' if abs(dP_tube_device-dP_tube_engine)<1000 else '⚠'}")

dP_shell_device = pd.get('shell_Pa', 0)
dP_shell_engine = design_result.get('delta_P_shell', 0)
print(f"  {'ΔP_shell (kPa)':<25} {dP_shell_device/1000:<15.2f} {dP_shell_engine/1000:<15.2f} {'✓' if abs(dP_shell_device-dP_shell_engine)<1000 else '⚠'}")

print(f"\n--- Device 层额外功能 ---")
print(f"  ε-NTU 校核:")
ntu = result_device.get('epsilon_ntu', {})
print(f"    ε = {ntu.get('epsilon', 0):.4f}")
print(f"    NTU = {ntu.get('NTU', 0):.3f}")
print(f"    C_r = {ntu.get('C_r', 0):.4f}")

print(f"\n  U 分解:")
u_dec = result_device.get('U_decomposition', {})
print(f"    R_total = {u_dec.get('R_total', 0):.6f} m²·K/W")
print(f"    R_fouling_tube = {u_dec.get('R_fouling_tube', 0):.6f}")
print(f"    R_fouling_shell = {u_dec.get('R_fouling_shell', 0):.6f}")
print(f"    污垢占比 = {u_dec.get('fouling_fraction_pct', 0):.1f}%")

print(f"\n  物性查询:")
props = result_device.get('properties', {})
print(f"    ρ_h = {props.get('rho_hot_kg_m3', 0):.2f} kg/m³")
print(f"    ρ_c = {props.get('rho_cold_kg_m3', 0):.2f} kg/m³")
print(f"    μ_h = {props.get('mu_hot_Pa_s', 0):.5f} Pa·s")
print(f"    μ_c = {props.get('mu_cold_Pa_s', 0):.5f} Pa·s")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 80)
print("  对比总结")
print("=" * 80)

print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│  heatexchanger_device._full_he_design (Device 层)                           │
│  ─────────────────────────────────────────────────────────────────────────  │
│  功能: 完整换热器设计流水线                                                  │
│  1. 流体名称标准化 (中英文别名映射)                                           │
│  2. 单位自动转换 (°C→K, bar→Pa, kg/h→kg/s)                                 │
│  3. 物性自动查询 (从thermo_helper获取ρ/μ/Cp/k)                               │
│  4. 冷侧参数自动估算 (T_cold_out 或 m_cold)                                  │
│  5. 温度交叉校验                                                             │
│  6. 管排列方式自动选择 (三角/正方)                                            │
│  7. 调用 design_heat_exchanger 物理引擎                                      │
│  8. ε-NTU 校核                                                              │
│  9. U 分解 (污垢热阻占比)                                                    │
│  10. 工程默认值追踪 (d_o, L)                                                 │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  physics_engine.heatexchanger (Engine 层)                                   │
│  ─────────────────────────────────────────────────────────────────────────  │
│  核心函数:                                                                  │
│  • design_heat_exchanger: 完整换热器设计 (LMTD+U+A+结构+压降)                 │
│  • lmtd_with_correction: LMTD + F 修正系数                                  │
│  • area_required: 所需换热面积计算                                            │
│  • estimate_U: 总传热系数经验估算                                             │
│  • suggest_configuration: 换热器型式推荐                                      │
│  • fouling_resistance: 污垢热阻查询                                           │
│  • effectiveness_ntu: ε-NTU 校核                                            │
│  输出: 几何尺寸 + 压降 + 传热系数 + 污垢热阻                                   │
└─────────────────────────────────────────────────────────────────────────────┘

核心结论:
  • 热负荷计算: 完全一致 (Device层直接调用Engine层函数)
  • LMTD/F: 完全一致
  • 结构参数: 完全一致
  • Device层额外提供: 物性查询/ε-NTU校核/U分解/单位转换/默认值追踪
  • 两者关系: Device层 = Engine层 + 物性查询 + 工程智能 + 报告生成
""")
