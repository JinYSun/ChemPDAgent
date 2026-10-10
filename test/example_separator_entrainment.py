# -*- coding: utf-8 -*-
"""
立式气液分离器液滴夹带校核示例
"""
import sys, os, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from physics_engine.flash_drum import entrainment_fraction
import math

# ============================================================
# 已知条件
# ============================================================
D = 1.4              # m, 分离器内径
rho_G = 15.0         # kg/m³, 气相密度
rho_L = 650.0        # kg/m³, 液相密度
m_gas = 3200.0       # kg/h, 气相质量流量

print("=" * 70)
print("  立式气液分离器 — 液滴夹带校核")
print("=" * 70)

# ============================================================
# [步骤1] 空塔气速计算
# ============================================================
A_cross = math.pi * (D / 2) ** 2
Q_gas = m_gas / rho_G / 3600   # m³/s
u_actual = Q_gas / A_cross      # m/s

print(f"\n[步骤1] 空塔气速计算")
print(f"  分离器内径 D:         {D} m")
print(f"  截面积 A:             π/4 × {D}² = {A_cross:.4f} m²")
print(f"  气相质量流量:         {m_gas:.0f} kg/h")
print(f"  气相密度 rho_G:       {rho_G} kg/m³")
print(f"  液相密度 rho_L:       {rho_L} kg/m³")
print(f"  气相体积流量 Q_gas:   {m_gas}/{rho_G}/3600 = {Q_gas:.4f} m³/s")
print(f"  空塔气速 u:           {Q_gas:.4f}/{A_cross:.4f} = {u_actual:.2f} m/s")

# ============================================================
# [步骤2] Souders-Brown 气速校核 (标准方法)
# ============================================================
ent_result = entrainment_fraction(
    vapor_velocity=u_actual,
    gas_density=rho_G,
    liquid_density=rho_L,
    separator_diameter=D,
)

print(f"\n[步骤2] Souders-Brown 气速校核")
print(f"  函数: entrainment_fraction()")
print(f"  方法: {ent_result.get('method')}")
print(f"  最大允许气速 u_max:   {ent_result.get('u_max_m_s', 'N/A')} m/s")
print(f"  速度比 u/u_max:       {ent_result.get('velocity_ratio', 'N/A')}")
print(f"  夹带风险等级:         {ent_result.get('entrainment_risk', 'N/A')}")
print(f"  夹带分率估算:         {ent_result.get('entrainment_fraction', 'N/A')}")
print(f"  备注: {ent_result.get('note', '')}")

# ============================================================
# [步骤3] 带表面张力的补充校核（临界液滴直径）
# ============================================================
sigma = 0.025  # N/m, 假设典型烃类体系
ent_result_sigma = entrainment_fraction(
    vapor_velocity=u_actual,
    gas_density=rho_G,
    liquid_density=rho_L,
    surface_tension=sigma,
    separator_diameter=D,
)

print(f"\n[步骤3] 临界液滴直径估算（补充信息）")
print(f"  假设表面张力 sigma:   {sigma} N/m")
d_crit = ent_result_sigma.get('critical_droplet_diameter_um')
if d_crit is not None:
    print(f"  临界液滴直径 d_crit:  {d_crit} μm")
    print(f"  说明: 小于此尺寸的液滴可能被气流夹带")
else:
    print(f"  无法计算临界液滴直径")

# ============================================================
# [结论]
# ============================================================
print(f"\n{'='*70}")
print(f"  校核结论")
print(f"{'='*70}")
print(f"  分离器内径:           {D} m")
print(f"  空塔气速:             {u_actual:.2f} m/s")
print(f"  Souders-Brown 允许气速: {ent_result.get('u_max_m_s', 'N/A')} m/s")
print(f"  速度比:               {ent_result.get('velocity_ratio', 'N/A')}")
risk = ent_result.get('entrainment_risk', 'unknown')
if risk == 'low':
    print(f"  结论: 气速远低于夹带限，该分离器尺寸不会导致液滴夹带")
elif risk == 'medium':
    print(f"  结论: 夹带风险低，该分离器尺寸基本安全")
elif risk == 'high':
    print(f"  结论: 接近夹带限，建议增大罐径或加装除沫器")
else:
    print(f"  结论: 超允许气速，存在严重夹带风险！必须增大罐径")
print(f"{'='*70}")
