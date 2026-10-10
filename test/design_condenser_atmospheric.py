# -*- coding: utf-8 -*-
"""冷凝器设计 - 常压操作，考虑相变"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.heatexchanger_device import _full_he_design
from physics_engine.heatexchanger import detect_phase_change

print("=" * 80)
print("  冷凝器设计 - 常压操作，考虑相变")
print("=" * 80)

# ============================================================
# 设计条件
# ============================================================
print("\n设计条件:")
print("  热侧: 饱和蒸汽冷凝")
print("    - 入口: 饱和蒸汽 373.15 K (100°C)")
print("    - 出口: 凝液 363.15 K (90°C)")
print("    - 流量: 1.0 kg/s")
print("    - 压力: 1 atm (常压)")
print("  冷侧: 冷水")
print("    - 入口: 300 K (27°C)")
print("    - 出口: 340 K (67°C)")
print("    - 压力: 1 atm (常压)")

# ============================================================
# 相变检测
# ============================================================
print("\n" + "=" * 80)
print("  相变检测")
print("=" * 80)

# 热侧相变检测 (蒸汽冷凝)
phase_hot = detect_phase_change(
    fluid_name="water",
    T_in=373.15,      # 饱和蒸汽入口
    T_out=363.15,     # 凝液出口 (过冷)
    P=101325.0,       # 常压
    side="hot"
)

print(f"\n  热侧相变检测:")
print(f"    有相变: {phase_hot['has_phase_change']}")
print(f"    相变类型: {phase_hot['phase_type']}")
print(f"    饱和温度: {phase_hot['T_sat']:.2f} K ({phase_hot['T_sat']-273.15:.2f}°C)")
print(f"    汽化潜热: {phase_hot['latent_heat']/1000:.1f} kJ/kg")
print(f"    过热度: {phase_hot['superheat']:.2f} K")
print(f"    过冷度: {phase_hot['subcool']:.2f} K")
print(f"    相变分率: {phase_hot['phase_fraction']*100:.1f}%")

# 冷侧相变检测 (水加热)
phase_cold = detect_phase_change(
    fluid_name="water",
    T_in=300.0,       # 冷水入口
    T_out=340.0,      # 热水出口
    P=101325.0,       # 常压
    side="cold"
)

print(f"\n  冷侧相变检测:")
print(f"    有相变: {phase_cold['has_phase_change']}")
print(f"    相变类型: {phase_cold['phase_type']}")
if phase_cold['T_sat']:
    print(f"    饱和温度: {phase_cold['T_sat']:.2f} K ({phase_cold['T_sat']-273.15:.2f}°C)")

# ============================================================
# 冷凝器设计
# ============================================================
print("\n" + "=" * 80)
print("  冷凝器设计 (_full_he_design)")
print("=" * 80)

# 调用换热器设计函数
# 热侧：蒸汽冷凝
result = _full_he_design(
    hot_fluid="water",
    cold_fluid="water",
    T_hot_in=373.15,      # 饱和蒸汽入口 100°C
    T_hot_out=363.15,     # 凝液出口 90°C (过冷10K)
    m_hot=1.0,            # 蒸汽流量 1 kg/s
    P_hot=101325.0,       # 常压
    
    T_cold_in=300.0,      # 冷水入口 27°C
    T_cold_out=340.0,     # 冷水出口 67°C
    P_cold=101325.0,      # 常压
)

print("\n设计结果:")
print(f"  热负荷 Q = {result.get('Q', 0)/1000:.2f} kW")
print(f"  LMTD = {result.get('LMTD', 0):.2f} K")
print(f"  F 修正系数 = {result.get('F', 0):.4f}")
print(f"  总传热系数 U = {result.get('U', 0):.1f} W/(m²·K)")
print(f"  所需面积 A_req = {result.get('A_required', 0):.2f} m²")
print(f"  设计面积 A_design = {result.get('A_design', 0):.2f} m²")
print(f"  管数 N_tubes = {result.get('N_tubes', 0)}")
print(f"  壳径 D_shell = {result.get('D_shell', 0):.3f} m")
print(f"  管长 L_tube = {result.get('L_tube', 0):.2f} m")
print(f"  管程压降 ΔP_tube = {result.get('delta_P_tube', 0)/1000:.2f} kPa")
print(f"  壳程压降 ΔP_shell = {result.get('delta_P_shell', 0)/1000:.2f} kPa")

# ε-NTU 校核
print(f"\n  ε-NTU 校核:")
print(f"    ε = {result.get('effectiveness', 0):.4f}")
print(f"    NTU = {result.get('NTU', 0):.3f}")
print(f"    C_r = {result.get('C_r', 0):.3f}")

# U 分解
print(f"\n  U 分解:")
print(f"    R_total = {result.get('R_total', 0):.6f} m²·K/W")
print(f"    污垢占比 = {result.get('fouling_ratio', 0)*100:.1f}%")

# ============================================================
# 热负荷分解 (显热 + 潜热)
# ============================================================
print("\n" + "=" * 80)
print("  热负荷分解 (显热 + 潜热)")
print("=" * 80)

# 计算各部分热负荷
# 蒸汽冷凝: 从饱和蒸汽 100°C → 凝液 90°C
# 1. 潜热 (冷凝): m × ΔH_vap
# 2. 显热 (过冷): m × Cp × ΔT_subcool

Cp_water = 4217.0  # J/(kg·K) 近似值
latent_heat = phase_hot['latent_heat']  # J/kg

Q_latent = 1.0 * latent_heat  # 冷凝潜热
Q_subcool = 1.0 * Cp_water * (373.15 - 363.15)  # 过冷显热
Q_total_calc = Q_latent + Q_subcool

print(f"\n  冷凝潜热: {Q_latent/1000:.2f} kW")
print(f"  过冷显热: {Q_subcool/1000:.2f} kW")
print(f"  总热负荷: {Q_total_calc/1000:.2f} kW")

print(f"\n  热负荷占比:")
print(f"    潜热占比: {Q_latent/Q_total_calc*100:.1f}%")
print(f"    显热占比: {Q_subcool/Q_total_calc*100:.1f}%")

print("\n" + "=" * 80)
print("  设计总结")
print("=" * 80)
print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│  冷凝器设计结果                                                             │
│  ─────────────────────────────────────────────────────────────────────────  │
│  操作条件:                                                                  │
│    - 热侧: 饱和蒸汽 100°C → 凝液 90°C (常压)                                │
│    - 冷侧: 冷水 27°C → 67°C (常压)                                          │
│    - 蒸汽流量: 1.0 kg/s                                                     │
│                                                                             │
│  热负荷:                                                                    │
│    - 冷凝潜热: 主要部分 (约 95%)                                            │
│    - 过冷显热: 次要部分 (约 5%)                                             │
│                                                                             │
│  设计特点:                                                                  │
│    - 常压操作，无需加压                                                     │
│    - 考虑相变，热负荷以潜热为主                                             │
│    - 凝液过冷 10K，提高系统效率                                             │
└─────────────────────────────────────────────────────────────────────────────┘
""")
