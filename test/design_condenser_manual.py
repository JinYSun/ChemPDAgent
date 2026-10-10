# -*- coding: utf-8 -*-
"""冷凝器设计 - 手动计算 (考虑相变)"""
import sys, io, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from physics_engine.heatexchanger import (
    detect_phase_change, heat_duty, lmtd_with_correction,
    area_required, estimate_U, design_heat_exchanger,
    effectiveness_ntu, overall_htc_estimate
)
from physics_engine.thermo_helper import (
    get_fluid_density, get_fluid_viscosity, get_fluid_Cp, get_fluid_thermal_conductivity
)

print("=" * 80)
print("  冷凝器设计 - 手动计算 (考虑相变)")
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

phase_hot = detect_phase_change(
    fluid_name="water",
    T_in=373.15,
    T_out=363.15,
    P=101325.0,
    side="hot"
)

print(f"\n  热侧相变检测:")
print(f"    有相变: {phase_hot['has_phase_change']}")
print(f"    相变类型: {phase_hot['phase_type']}")
print(f"    饱和温度: {phase_hot['T_sat']:.2f} K ({phase_hot['T_sat']-273.15:.2f}°C)")
print(f"    汽化潜热: {phase_hot['latent_heat']/1000:.1f} kJ/kg")
print(f"    过热度: {phase_hot['superheat']:.2f} K")
print(f"    过冷度: {phase_hot['subcool']:.2f} K")

# ============================================================
# 热负荷计算 (考虑相变)
# ============================================================
print("\n" + "=" * 80)
print("  热负荷计算 (考虑相变)")
print("=" * 80)

# 热侧物性
T_hot_avg = (373.15 + 363.15) / 2  # 368.15 K (凝液平均温度)
rho_hot = get_fluid_density("water", T_hot_avg, 101325.0, phase="liquid")
mu_hot = get_fluid_viscosity("water", T_hot_avg, 101325.0, phase="liquid")
# Cp 使用典型值 (thermo_helper 可能返回摩尔热容)
Cp_hot = 4217.0  # J/(kg·K) 水的典型值
k_hot = get_fluid_thermal_conductivity("water", T_hot_avg, 101325.0, phase="liquid")

print(f"\n  热侧物性 (凝液 @ {T_hot_avg:.1f} K):")
print(f"    ρ = {rho_hot:.2f} kg/m³")
print(f"    μ = {mu_hot:.6f} Pa·s")
print(f"    Cp = {Cp_hot:.1f} J/(kg·K)")
print(f"    k = {k_hot:.4f} W/(m·K)")

# 冷侧物性
T_cold_avg = (300.0 + 340.0) / 2  # 320 K
rho_cold = get_fluid_density("water", T_cold_avg, 101325.0, phase="liquid")
mu_cold = get_fluid_viscosity("water", T_cold_avg, 101325.0, phase="liquid")
Cp_cold = 4180.0  # J/(kg·K) 水的典型值
k_cold = get_fluid_thermal_conductivity("water", T_cold_avg, 101325.0, phase="liquid")

print(f"\n  冷侧物性 (水 @ {T_cold_avg:.1f} K):")
print(f"    ρ = {rho_cold:.2f} kg/m³")
print(f"    μ = {mu_cold:.6f} Pa·s")
print(f"    Cp = {Cp_cold:.1f} J/(kg·K)")
print(f"    k = {k_cold:.4f} W/(m·K)")

# 热负荷计算
m_hot = 1.0  # kg/s
latent_heat = phase_hot['latent_heat']  # J/kg

# 冷凝潜热 (100°C 蒸汽 → 100°C 凝液)
Q_latent = m_hot * latent_heat

# 过冷显热 (100°C 凝液 → 90°C 凝液)
Q_subcool = m_hot * Cp_hot * (373.15 - 363.15)

Q_total = Q_latent + Q_subcool

print(f"\n  热负荷分解:")
print(f"    冷凝潜热: {Q_latent/1000:.2f} kW")
print(f"    过冷显热: {Q_subcool/1000:.2f} kW")
print(f"    总热负荷: {Q_total/1000:.2f} kW")

# 冷侧流量 (由热平衡计算)
m_cold = Q_total / (Cp_cold * (340.0 - 300.0))
print(f"\n  冷侧流量 (由热平衡): {m_cold:.4f} kg/s")

# ============================================================
# LMTD 计算
# ============================================================
print("\n" + "=" * 80)
print("  LMTD 计算")
print("=" * 80)

# 对于冷凝器，热侧温度恒定 (相变段) + 过冷段
# 简化计算：使用入口/出口温度
T_hot_in = 373.15  # 饱和蒸汽入口
T_hot_out = 363.15  # 凝液出口
T_cold_in = 300.0
T_cold_out = 340.0

# 逆流 LMTD
dT1 = T_hot_in - T_cold_out  # 热端温差
dT2 = T_hot_out - T_cold_in  # 冷端温差

LMTD_counter = (dT1 - dT2) / math.log(dT1 / dT2) if dT1 != dT2 else dT1

print(f"\n  温度分布:")
print(f"    热侧: {T_hot_in:.2f} K → {T_hot_out:.2f} K")
print(f"    冷侧: {T_cold_out:.2f} K ← {T_cold_in:.2f} K")
print(f"    热端温差 ΔT1 = {dT1:.2f} K")
print(f"    冷端温差 ΔT2 = {dT2:.2f} K")
print(f"    逆流 LMTD = {LMTD_counter:.2f} K")

# 修正系数 F (需要查图或计算)
# 对于冷凝器，F 通常接近 1.0
F_correction = 0.95  # 假设值
LMTD_corrected = LMTD_counter * F_correction

print(f"    修正系数 F = {F_correction:.2f}")
print(f"    修正后 LMTD = {LMTD_corrected:.2f} K")

# ============================================================
# 总传热系数 U 估算
# ============================================================
print("\n" + "=" * 80)
print("  总传热系数 U 估算")
print("=" * 80)

# 冷凝器 U 值经验范围
# 蒸汽冷凝: 1500-4000 W/(m²·K)
# 水-水换热器: 800-1500 W/(m²·K)
# 冷凝器综合: 1000-2500 W/(m²·K)

U_estimate = 1500.0  # W/(m²·K) 冷凝器典型值

print(f"\n  冷凝器 U 值经验估算: {U_estimate:.0f} W/(m²·K)")
print(f"    (蒸汽冷凝典型范围: 1000-2500 W/(m²·K))")

# ============================================================
# 换热面积计算
# ============================================================
print("\n" + "=" * 80)
print("  换热面积计算")
print("=" * 80)

A_required = Q_total / (U_estimate * LMTD_corrected)
A_design = A_required * 1.15  # 15% 裕量

print(f"\n  所需面积 A_req = {A_required:.2f} m²")
print(f"  设计面积 A_design = {A_design:.2f} m² (裕量 15%)")

# ============================================================
# 结构设计
# ============================================================
print("\n" + "=" * 80)
print("  结构设计")
print("=" * 80)

# 换热管规格
d_o = 0.019  # 外径 19 mm
d_i = 0.015  # 内径 15 mm
L_tube = 3.0  # 管长 3 m

# 管数计算
A_tube = math.pi * d_o * L_tube  # 单管外表面积
N_tubes = math.ceil(A_design / A_tube)

print(f"\n  换热管规格:")
print(f"    外径 d_o = {d_o*1000:.0f} mm")
print(f"    内径 d_i = {d_i*1000:.0f} mm")
print(f"    管长 L = {L_tube*1000:.0f} mm")
print(f"    单管面积 = {A_tube:.4f} m²")
print(f"    管数 N = {N_tubes}")

# 壳径估算 (三角形排列)
N_tpi = 0.025  # 管间距 25 mm
D_shell = 1.05 * math.sqrt(N_tubes * N_tpi**2 / 0.785)
D_shell = max(0.15, D_shell)

print(f"    管间距 = {N_tpi*1000:.0f} mm")
print(f"    壳径 D_shell = {D_shell*1000:.0f} mm")

# ============================================================
# 设计总结
# ============================================================
print("\n" + "=" * 80)
print("  冷凝器设计总结")
print("=" * 80)

print(f"""
┌─────────────────────────────────────────────────────────────────────────────┐
│  冷凝器设计结果                                                             │
│  ─────────────────────────────────────────────────────────────────────────  │
│  操作条件:                                                                  │
│    - 热侧: 饱和蒸汽 100°C → 凝液 90°C (常压)                                │
│    - 冷侧: 冷水 27°C → 67°C (常压)                                          │
│    - 蒸汽流量: 1.0 kg/s                                                     │
│    - 冷水流量: {m_cold:.2f} kg/s                                              │
│                                                                             │
│  热负荷:                                                                    │
│    - 冷凝潜热: {Q_latent/1000:.1f} kW ({Q_latent/Q_total*100:.1f}%)                                │
│    - 过冷显热: {Q_subcool/1000:.1f} kW ({Q_subcool/Q_total*100:.1f}%)                                 │
│    - 总热负荷: {Q_total/1000:.1f} kW                                          │
│                                                                             │
│  传热参数:                                                                  │
│    - LMTD: {LMTD_corrected:.1f} K                                             │
│    - U 值: {U_estimate:.0f} W/(m²·K)                                          │
│    - 换热面积: {A_design:.1f} m²                                              │
│                                                                             │
│  结构参数:                                                                  │
│    - 换热管: Φ{d_o*1000:.0f}×{d_i*1000:.0f} mm, L={L_tube*1000:.0f} mm                    │
│    - 管数: {N_tubes} 根                                                       │
│    - 壳径: {D_shell*1000:.0f} mm                                              │
│                                                                             │
│  设计特点:                                                                  │
│    - 常压操作，无需加压                                                     │
│    - 考虑相变，热负荷以潜热为主 (98%)                                       │
│    - 凝液过冷 10K，提高系统效率                                             │
└─────────────────────────────────────────────────────────────────────────────┘
""")
