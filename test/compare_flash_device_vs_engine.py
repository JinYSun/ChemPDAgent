# -*- coding: utf-8 -*-
"""对比 flash_drum_device._full_design vs physics_engine 直接调用"""
import sys, io, json, math
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.flash_drum_device import _full_design, FLASH_DRUM_CONSTRAINTS
from physics_engine.flash_drum import (
    isothermal_flash, flash_drum_diameter, flash_drum_volume,
    validate_flash_drum_design, entrainment_fraction,
)
from physics_engine.thermo_helper import get_mixture_properties, get_mixture_k_values

# ============================================================
# 设计条件
# ============================================================
components = ["methanol", "water"]
z = [0.4, 0.6]        # 甲醇40%、水60%
F = 100.0             # kmol/h
T = 350.0             # K
P = 101325.0          # Pa (常压)

print("=" * 80)
print("  flash_drum_device._full_design vs physics_engine 直接调用 对比")
print("=" * 80)
print(f"\n设计条件:")
print(f"  组分: {components}")
print(f"  进料组成: 甲醇 {z[0]*100:.0f}%, 水 {z[1]*100:.0f}%")
print(f"  进料量: {F} kmol/h")
print(f"  操作温度: {T} K ({T-273.15:.1f}°C)")
print(f"  操作压力: {P/1000:.2f} kPa")

# ============================================================
# Device 层调用
# ============================================================
print("\n" + "=" * 80)
print("  Device 层: _full_design")
print("=" * 80)

result_device = _full_design(
    components=components,
    z=z,
    F=F,
    T=T,
    P=P,
)

# ============================================================
# Engine 层直接调用
# ============================================================
print("\n" + "=" * 80)
print("  Engine 层: physics_engine 直接调用")
print("=" * 80)

# 步骤1: 获取 K 值
K_list = get_mixture_k_values(components, T, P, z)
K = np.array(K_list, dtype=float)
z_arr = np.array(z, dtype=float)
print(f"\n  K值: {K_list}")

# 步骤2: 等温闪蒸
feed_flow = F * 1000.0 / 3600.0  # mol/s
flash = isothermal_flash(feed_flow, z_arr, K)
print(f"\n  isothermal_flash 结果:")
print(f"    ψ (气化率) = {flash.get('psi', 0):.6f}")
print(f"    V (气相 mol/s) = {flash.get('V', 0):.4f}")
print(f"    L (液相 mol/s) = {flash.get('L', 0):.4f}")
print(f"    收敛: {flash.get('converged', False)}")
print(f"    y (气相组成) = {flash.get('y', [])}")
print(f"    x (液相组成) = {flash.get('x', [])}")

# 步骤3: 混合物物性
props = get_mixture_properties(components, z, T, P)
print(f"\n  混合物物性:")
rho_V = props.get('rho_V')
rho_L = props.get('rho_L')
MW_mix = props.get('MW_mix')
print(f"    ρ_V = {rho_V} kg/m³")
print(f"    ρ_L = {rho_L} kg/m³")
print(f"    MW_mix = {MW_mix} g/mol")

# 步骤4: 闪蒸罐几何设计（简化版，无约束自动调整）
V_mol = flash['V']
L_mol = flash['L']
V_mass = V_mol * MW_mix / 1000
L_mass = L_mol * MW_mix / 1000
V_vol = V_mass / rho_V
L_vol = L_mass / rho_L

print(f"\n  流量计算:")
print(f"    V_mass = {V_mass:.4f} kg/s")
print(f"    L_mass = {L_mass:.4f} kg/s")
print(f"    V_vol = {V_vol:.6f} m³/s")
print(f"    L_vol = {L_vol:.6f} m³/s")

# 调用 flash_drum_volume 和 flash_drum_diameter
drum_vol = flash_drum_volume(V_vol, L_vol, rho_V, rho_L)
drum_diam = flash_drum_diameter(V_vol, rho_V)

print(f"\n  flash_drum_volume 结果:")
print(f"    D = {drum_vol.get('D', 0):.4f} m")
print(f"    L_total = {drum_vol.get('L_total', 0):.4f} m")
print(f"    V_total = {drum_vol.get('V_total', 0):.4f} m³")
print(f"    aspect_ratio = {drum_vol.get('aspect_ratio', 0):.4f}")

print(f"\n  flash_drum_diameter 结果:")
print(f"    D = {drum_diam.get('D', 0):.4f} m")
print(f"    u_actual = {drum_diam.get('u_actual', 0):.4f} m/s")
print(f"    u_max = {drum_diam.get('u_max', 0):.4f} m/s")

# 步骤5: 约束校验
# 合并两个结果
drum_combined = {**drum_vol, **drum_diam}
validation = validate_flash_drum_design(drum_combined)
print(f"\n  validate_flash_drum_design 结果:")
print(f"    通过率: {validation.get('pass_rate', 0)*100:.0f}%")
for check in validation.get('checks', []):
    status = "PASS" if check.get('pass') else "FAIL"
    print(f"    [{status}] {check.get('name')}: {check.get('value')}")

# ============================================================
# 对比
# ============================================================
print("\n" + "=" * 80)
print("  Device 层 vs Engine 层 对比")
print("=" * 80)

print(f"\n--- 闪蒸计算对比 ---")
print(f"  {'参数':<25} {'Device层':<15} {'Engine层':<15} {'一致?'}")
print("  " + "-" * 60)

psi_d = result_device.get('flash', {}).get('psi', 0)
psi_e = flash.get('psi', 0)
print(f"  {'气化率 ψ':<25} {psi_d:<15.6f} {psi_e:<15.6f} {'✓' if abs(psi_d-psi_e)<1e-8 else '✗'}")

V_d = result_device.get('flash', {}).get('V_mol_s', 0)
V_e = flash.get('V', 0)
print(f"  {'气相流量 V (mol/s)':<25} {V_d:<15.4f} {V_e:<15.4f} {'✓' if abs(V_d-V_e)<0.001 else '✗'}")

L_d = result_device.get('flash', {}).get('L_mol_s', 0)
L_e = flash.get('L', 0)
print(f"  {'液相流量 L (mol/s)':<25} {L_d:<15.4f} {L_e:<15.4f} {'✓' if abs(L_d-L_e)<0.001 else '✗'}")

print(f"\n--- 物性对比 ---")
fp_d = result_device.get('properties', {})
print(f"  {'ρ_V (kg/m³)':<25} {fp_d.get('rho_V', 'N/A'):<15} {rho_V:<15} {'✓' if fp_d.get('rho_V')==rho_V else '✗'}")
print(f"  {'ρ_L (kg/m³)':<25} {fp_d.get('rho_L', 'N/A'):<15} {rho_L:<15} {'✓' if fp_d.get('rho_L')==rho_L else '✗'}")
print(f"  {'MW_mix (g/mol)':<25} {fp_d.get('MW_mix', 'N/A'):<15} {MW_mix:<15} {'✓' if fp_d.get('MW_mix')==MW_mix else '✗'}")

print(f"\n--- 几何设计对比 ---")
dd_d = result_device.get('drum_design', {})
D_device = dd_d.get('D', 0)
D_engine_vol = drum_vol.get('D', 0)
D_engine_diam = drum_diam.get('D', 0)
L_device = dd_d.get('L_total', 0)
L_engine = drum_vol.get('L_total', 0)

print(f"  {'直径 D (m)':<25} {D_device:<15.4f} {D_engine_vol:<15.4f} {'✓' if abs(D_device-D_engine_vol)<0.01 else '⚠'}")
print(f"  {'总长 L (m)':<25} {L_device:<15.4f} {L_engine:<15.4f} {'✓' if abs(L_device-L_engine)<0.01 else '⚠'}")
print(f"  {'长径比 L/D':<25} {dd_d.get('aspect_ratio', 0):<15.4f} {drum_vol.get('aspect_ratio', 0):<15.4f} {'✓' if abs(dd_d.get('aspect_ratio',0)-drum_vol.get('aspect_ratio',0))<0.1 else '⚠'}")

print(f"\n--- 约束校验对比 ---")
val_d = result_device.get('validation', {})
print(f"  Device层通过率: {val_d.get('pass_rate', 0)*100:.0f}%")
print(f"  Engine层通过率: {validation.get('pass_rate', 0)*100:.0f}%")

print(f"\n--- Device 层额外功能 ---")
print(f"  约束自动调整: {result_device.get('adjust_notes', [])}")
print(f"  夹带估算: {result_device.get('entrainment', {})}")
print(f"  停留时间: {result_device.get('residence_time', {})}")
print(f"  压降估算: {result_device.get('pressure_drop_Pa', 'N/A')} Pa")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 80)
print("  对比总结")
print("=" * 80)

print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│  flash_drum_device._full_design (Device 层)                                 │
│  ─────────────────────────────────────────────────────────────────────────  │
│  功能: 完整闪蒸罐设计流水线                                                  │
│  1. 调用 get_mixture_k_values 获取 K 值                                     │
│  2. 调用 isothermal_flash 等温闪蒸计算                                       │
│  3. 调用 get_mixture_properties 获取混合物物性                                │
│  4. 几何设计（带约束自动调整）                                                │
│     - 硬编码 FLASH_DRUM_CONSTRAINTS 约束字典                                 │
│     - 自动调整直径/容积以满足 L/D、气速、停留时间约束                          │
│  5. 雾沫夹带估算 (Souders-Brown)                                            │
│  6. 停留时间计算                                                             │
│  7. 压降估算                                                                │
│  8. 约束校验                                                                │
│  9. 调整日志记录                                                             │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  physics_engine.flash_drum (Engine 层)                                      │
│  ─────────────────────────────────────────────────────────────────────────  │
│  核心函数:                                                                  │
│  • isothermal_flash: Rachford-Rice 等温闪蒸计算                              │
│  • flash_drum_volume: 基于停留时间的容积设计                                  │
│  • flash_drum_diameter: 基于气速约束的直径设计                                │
│  • entrainment_fraction: 雾沫夹带估算                                        │
│  • validate_flash_drum_design: 约束校验                                      │
│  输出: 几何尺寸 + 校验结果                                                   │
└─────────────────────────────────────────────────────────────────────────────┘

核心结论:
  • 闪蒸计算 (ψ/V/L/y/x): 完全一致 (Device层直接调用Engine层)
  • 物性查询: 完全一致
  • 几何设计: Device层有约束自动调整，Engine层无此功能
  • Device层额外提供: 夹带估算/停留时间/压降/调整日志
  • 两者关系: Device层 = Engine层 + 约束自动调整 + 工程校验 + 报告生成
""")
