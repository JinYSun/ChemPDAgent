# -*- coding: utf-8 -*-
"""对比 heatexchanger_device._full_he_design vs physics_engine 直接调用 - 含相变检测"""
import sys, io, json, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.heatexchanger_device import _full_he_design
from physics_engine.heatexchanger import (
    design_heat_exchanger, lmtd_with_correction, heat_duty_from_flow,
    estimate_U, suggest_configuration, fouling_resistance,
    lmtd, area_required, effectiveness_ntu, detect_phase_change,
)
from physics_engine.thermo_helper import (
    get_fluid_density, get_fluid_viscosity, get_fluid_Cp, get_fluid_thermal_conductivity,
    get_molecular_weight,
)

# ============================================================
# 案例1: 无相变 - 热水/冷水换热
# ============================================================
print("=" * 80)
print("  案例1: 无相变 - 热水/冷水换热")
print("=" * 80)

T_hot_in = 400.0      # K (127°C) 热水入口
T_hot_out = 350.0     # K (77°C) 热水出口
T_cold_in = 300.0     # K (27°C) 冷水入口
T_cold_out = 340.0    # K (67°C) 冷水出口
m_hot = 1.0           # kg/s
P_hot = 500000.0      # Pa (5 bar) - 加压防汽化
P_cold = 300000.0     # Pa (3 bar)

print(f"\n设计条件:")
print(f"  热侧: 热水 {T_hot_in} K → {T_hot_out} K, m = {m_hot} kg/s, P = {P_hot/1e5:.1f} bar")
print(f"  冷侧: 冷水 {T_cold_in} K → {T_cold_out} K, P = {P_cold/1e5:.1f} bar")

# 相变检测
print(f"\n--- 相变检测 (Engine层: detect_phase_change) ---")
phase_hot = detect_phase_change("water", T_hot_in, T_hot_out, P_hot, side="hot")
phase_cold = detect_phase_change("water", T_cold_in, T_cold_out, P_cold, side="cold")

print(f"\n  热侧:")
print(f"    T_sat = {phase_hot.get('T_sat', 'N/A')} K")
print(f"    有相变: {phase_hot.get('has_phase_change', False)}")
print(f"    相变类型: {phase_hot.get('phase_type', 'single')}")
if phase_hot.get('T_sat'):
    print(f"    过热度: {phase_hot.get('superheat', 0):.2f} K")
    print(f"    过冷度: {phase_hot.get('subcool', 0):.2f} K")

print(f"\n  冷侧:")
print(f"    T_sat = {phase_cold.get('T_sat', 'N/A')} K")
print(f"    有相变: {phase_cold.get('has_phase_change', False)}")
print(f"    相变类型: {phase_cold.get('phase_type', 'single')}")

# ============================================================
# 案例2: 热侧有相变 - 蒸汽冷凝
# ============================================================
print("\n" + "=" * 80)
print("  案例2: 热侧有相变 - 蒸汽冷凝")
print("=" * 80)

# 蒸汽在 5 bar 下饱和温度约 425 K (152°C)
T_steam_in = 440.0    # K (167°C) 过热蒸汽入口
T_steam_out = 400.0   # K (127°C) 冷凝液出口 (过冷)
m_steam = 0.5         # kg/s
P_steam = 500000.0    # Pa (5 bar)

T_cw_in = 300.0       # K (27°C) 冷却水入口
T_cw_out = 320.0      # K (47°C) 冷却水出口
P_cw = 300000.0       # Pa (3 bar)

print(f"\n设计条件:")
print(f"  热侧: 蒸汽 {T_steam_in} K → {T_steam_out} K, m = {m_steam} kg/s, P = {P_steam/1e5:.1f} bar")
print(f"  冷侧: 冷却水 {T_cw_in} K → {T_cw_out} K, P = {P_cw/1e5:.1f} bar")

# 相变检测
print(f"\n--- 相变检测 (Engine层: detect_phase_change) ---")
phase_steam = detect_phase_change("water", T_steam_in, T_steam_out, P_steam, side="hot")
phase_cw = detect_phase_change("water", T_cw_in, T_cw_out, P_cw, side="cold")

print(f"\n  热侧 (蒸汽):")
print(f"    T_sat = {phase_steam.get('T_sat', 'N/A')} K ({phase_steam.get('T_sat', 0) - 273.15:.1f}°C)")
print(f"    有相变: {phase_steam.get('has_phase_change', False)}")
print(f"    相变类型: {phase_steam.get('phase_type', 'single')}")
print(f"    潜热: {phase_steam.get('latent_heat', 0)/1000:.1f} kJ/kg")
print(f"    过热度: {phase_steam.get('superheat', 0):.2f} K")
print(f"    过冷度: {phase_steam.get('subcool', 0):.2f} K")
print(f"    相变分率: {phase_steam.get('phase_fraction', 0)*100:.1f}%")

print(f"\n  冷侧 (冷却水):")
print(f"    T_sat = {phase_cw.get('T_sat', 'N/A')} K")
print(f"    有相变: {phase_cw.get('has_phase_change', False)}")
print(f"    相变类型: {phase_cw.get('phase_type', 'single')}")

# 计算热负荷 (含相变)
if phase_steam.get('has_phase_change'):
    # 获取物性
    T_steam_avg = (T_steam_in + T_steam_out) / 2
    T_cw_avg = (T_cw_in + T_cw_out) / 2
    
    Cp_steam = get_fluid_Cp("water", T_steam_avg, P_steam, phase="liquid")
    MW_steam = get_molecular_weight("water")
    Cp_steam = Cp_steam / (MW_steam / 1000.0)  # J/(kg·K)
    
    Cp_cw = get_fluid_Cp("water", T_cw_avg, P_cw, phase="liquid")
    MW_cw = get_molecular_weight("water")
    Cp_cw = Cp_cw / (MW_cw / 1000.0)  # J/(kg·K)
    
    # 热负荷计算
    Q_sensible = m_steam * Cp_steam * abs(T_steam_in - T_steam_out)
    Q_latent = m_steam * phase_steam.get('latent_heat', 0) * phase_steam.get('phase_fraction', 0)
    Q_total_hot = Q_sensible + Q_latent
    
    # 冷侧流量估算
    delta_T_cw = T_cw_out - T_cw_in
    m_cw = Q_total_hot / (Cp_cw * delta_T_cw)
    Q_total_cw = m_cw * Cp_cw * delta_T_cw
    
    print(f"\n--- 热负荷计算 (含相变) ---")
    print(f"  热侧显热: {Q_sensible/1000:.2f} kW")
    print(f"  热侧潜热: {Q_latent/1000:.2f} kW")
    print(f"  热侧总热负荷: {Q_total_hot/1000:.2f} kW")
    print(f"  冷侧流量 (由热平衡): {m_cw:.4f} kg/s")
    print(f"  冷侧热负荷: {Q_total_cw/1000:.2f} kW")
    
    # LMTD 计算
    lmtd_result = lmtd_with_correction(
        T_hot_in=T_steam_in, T_hot_out=T_steam_out,
        T_cold_in=T_cw_in, T_cold_out=T_cw_out,
    )
    LMTD = lmtd_result.get("LMTD")
    F = lmtd_result.get("F")
    
    print(f"\n--- LMTD 计算 ---")
    print(f"  LMTD = {LMTD:.2f} K")
    print(f"  F (修正系数) = {F:.4f}")
    print(f"  LMTD_corrected = {LMTD * F:.2f} K")

# ============================================================
# 案例3: 冷侧有相变 - 液体沸腾
# ============================================================
print("\n" + "=" * 80)
print("  案例3: 冷侧有相变 - 液体沸腾")
print("=" * 80)

# 热水加热甲醇使其部分汽化
T_hot_in2 = 400.0    # K (127°C) 热水入口
T_hot_out2 = 370.0   # K (97°C) 热水出口
m_hot2 = 2.0         # kg/s
P_hot2 = 500000.0    # Pa (5 bar)

T_meth_in = 300.0    # K (27°C) 甲醇入口 (液相)
T_meth_out = 350.0   # K (77°C) 甲醇出口 (部分汽化)
P_meth = 300000.0    # Pa (3 bar)

print(f"\n设计条件:")
print(f"  热侧: 热水 {T_hot_in2} K → {T_hot_out2} K, m = {m_hot2} kg/s, P = {P_hot2/1e5:.1f} bar")
print(f"  冷侧: 甲醇 {T_meth_in} K → {T_meth_out} K, P = {P_meth/1e5:.1f} bar")

# 相变检测
print(f"\n--- 相变检测 (Engine层: detect_phase_change) ---")
phase_hot2 = detect_phase_change("water", T_hot_in2, T_hot_out2, P_hot2, side="hot")
phase_meth = detect_phase_change("methanol", T_meth_in, T_meth_out, P_meth, side="cold")

print(f"\n  热侧 (热水):")
print(f"    T_sat = {phase_hot2.get('T_sat', 'N/A')} K")
print(f"    有相变: {phase_hot2.get('has_phase_change', False)}")

print(f"\n  冷侧 (甲醇):")
print(f"    T_sat = {phase_meth.get('T_sat', 'N/A')} K ({phase_meth.get('T_sat', 0) - 273.15:.1f}°C)")
print(f"    有相变: {phase_meth.get('has_phase_change', False)}")
print(f"    相变类型: {phase_meth.get('phase_type', 'single')}")
print(f"    潜热: {phase_meth.get('latent_heat', 0)/1000:.1f} kJ/kg")
print(f"    过冷度: {phase_meth.get('subcool', 0):.2f} K")
print(f"    过热度: {phase_meth.get('superheat', 0):.2f} K")
print(f"    相变分率: {phase_meth.get('phase_fraction', 0)*100:.1f}%")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 80)
print("  相变检测功能总结")
print("=" * 80)

print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│  detect_phase_change (Engine 层新增函数)                                     │
│  ─────────────────────────────────────────────────────────────────────────  │
│  功能: 检测流体是否发生相变，并计算相关参数                                    │
│                                                                             │
│  输入:                                                                      │
│  • fluid_name: 流体名称 (water, methanol, ethanol...)                       │
│  • T_in, T_out: 入口/出口温度 (K)                                           │
│  • P: 操作压力 (Pa)                                                         │
│  • side: "hot" 或 "cold"                                                    │
│                                                                             │
│  输出:                                                                      │
│  • has_phase_change: 是否发生相变                                            │
│  • phase_type: "condensing" (冷凝) / "evaporating" (蒸发) / "single"        │
│  • T_sat: 饱和温度 (K)                                                      │
│  • latent_heat: 潜热 (J/kg)                                                 │
│  • superheat: 过热度 (K)                                                    │
│  • subcool: 过冷度 (K)                                                      │
│  • phase_fraction: 相变分率 (0-1)                                            │
│                                                                             │
│  计算原理:                                                                  │
│  1. 二分法求解饱和温度 (蒸气压 = 操作压力)                                     │
│  2. 比较 T_in/T_out 与 T_sat 判断相变类型                                     │
│  3. 调用 calc_enthalpy_vaporization 获取潜热                                 │
└─────────────────────────────────────────────────────────────────────────────┘

使用场景:
  • 冷凝器设计: 蒸汽 → 冷凝液 (热侧相变)
  • 再沸器设计: 液体 → 蒸汽 (冷侧相变)
  • 部分冷凝/部分汽化: 计算相变分率和潜热贡献
""")
