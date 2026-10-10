# -*- coding: utf-8 -*-
"""150m³ 乙醇储罐设计方案 - 卧式/立式/球形三种方案对比"""
import sys, io, json, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.storage_tank_device import _full_design

# ============================================================
# 设计条件
# ============================================================
V = 150.0           # 设计容量 150 m³
fluid = "ethanol"   # 储存介质：乙醇
T = 298.15          # 操作温度 25°C (常温)
P = 101325.0        # 常压操作

print("=" * 70)
print("  150 m³ 乙醇储罐设计方案")
print("=" * 70)
print(f"\n设计条件:")
print(f"  储存介质: 乙醇 (ethanol)")
print(f"  设计容量: {V} m³")
print(f"  操作温度: {T} K ({T - 273.15:.1f}°C)")
print(f"  操作压力: {P/1000:.2f} kPa (常压)")

results = {}

# ============================================================
# 方案一：立式圆柱储罐
# ============================================================
print("\n" + "=" * 70)
print("  方案一：立式圆柱储罐 (Vertical Cylindrical Tank)")
print("=" * 70)

result_v = _full_design(
    V=V,
    tank_type="vertical",
    fluid=fluid,
    T=T,
    P=P,
)
results["vertical"] = result_v

print(f"\n--- 几何尺寸 ---")
print(f"  罐体直径 D: {result_v.get('D', 0):.2f} m")
print(f"  罐体高度 H: {result_v.get('H', 0):.2f} m")
print(f"  高径比 H/D: {result_v.get('H_D_ratio', result_v.get('aspect_ratio', 0)):.2f}")
print(f"  总容积: {result_v.get('V_total', 0):.2f} m³")
print(f"  设计容积(含15%气相裕量): {result_v.get('design_volume_m3', 0):.2f} m³")

print(f"\n--- 流体信息 ---")
fp = result_v.get('fluid_properties', {})
print(f"  液相密度: {fp.get('rho_liquid_kg_m3', 'N/A')} kg/m³")
print(f"  分子量: {fp.get('MW_g_mol', 'N/A')} g/mol")
print(f"  饱和蒸气压: {fp.get('vapor_pressure_Pa', 'N/A')} Pa")
print(f"  液体质量: {result_v.get('liquid_mass_tonnes', 0):.1f} 吨")

print(f"\n--- 基础荷载 ---")
fnd = result_v.get('foundation', {})
print(f"  总重量: {fnd.get('total_weight_kN', 0):.1f} kN")
print(f"  基底压力: {fnd.get('base_pressure_kPa', 0):.1f} kPa")

print(f"\n--- 蒸发损耗 ---")
bl = result_v.get('breathing_losses', {})
if bl:
    print(f"  每次呼吸损耗: {bl.get('loss_per_breath_kg', 0):.4f} kg")
    print(f"  年损耗量: {bl.get('annual_loss_kg', 0):.1f} kg")

print(f"\n--- 约束校核 ---")
val = result_v.get('validation', {})
print(f"  通过率: {val.get('pass_rate', 0)*100:.0f}% ({val.get('n_pass', 0)}/{val.get('n_total', 0)})")
for check in val.get('checks', []):
    status = "PASS" if check.get('pass') else "FAIL"
    print(f"  [{status}] {check.get('name')}: {check.get('value')}")
for w in result_v.get('warnings', []):
    if '默认值' not in str(w):
        print(f"  [警告] {w}")

# ============================================================
# 方案二：卧式圆柱储罐
# ============================================================
print("\n" + "=" * 70)
print("  方案二：卧式圆柱储罐 (Horizontal Cylindrical Tank)")
print("=" * 70)

result_h = _full_design(
    V=V,
    tank_type="horizontal",
    fluid=fluid,
    T=T,
    P=P,
)
results["horizontal"] = result_h

print(f"\n--- 几何尺寸 ---")
print(f"  罐体直径 D: {result_h.get('D', 0):.2f} m")
print(f"  罐体长度 L: {result_h.get('L', 0):.2f} m")
print(f"  长径比 L/D: {result_h.get('L_D_ratio', result_h.get('aspect_ratio', 0)):.2f}")
print(f"  总容积: {result_h.get('V_total', 0):.2f} m³")
print(f"  设计容积(含15%气相裕量): {result_h.get('design_volume_m3', 0):.2f} m³")

print(f"\n--- 流体信息 ---")
fp = result_h.get('fluid_properties', {})
print(f"  液相密度: {fp.get('rho_liquid_kg_m3', 'N/A')} kg/m³")
print(f"  液体质量: {result_h.get('liquid_mass_tonnes', 0):.1f} 吨")

print(f"\n--- 基础荷载 ---")
fnd = result_h.get('foundation', {})
print(f"  总重量: {fnd.get('total_weight_kN', 0):.1f} kN")
print(f"  基底压力: {fnd.get('base_pressure_kPa', 0):.1f} kPa")

print(f"\n--- 蒸发损耗 ---")
bl = result_h.get('breathing_losses', {})
if bl:
    print(f"  每次呼吸损耗: {bl.get('loss_per_breath_kg', 0):.4f} kg")
    print(f"  年损耗量: {bl.get('annual_loss_kg', 0):.1f} kg")

print(f"\n--- 约束校核 ---")
val = result_h.get('validation', {})
print(f"  通过率: {val.get('pass_rate', 0)*100:.0f}% ({val.get('n_pass', 0)}/{val.get('n_total', 0)})")
for check in val.get('checks', []):
    status = "PASS" if check.get('pass') else "FAIL"
    print(f"  [{status}] {check.get('name')}: {check.get('value')}")
for w in result_h.get('warnings', []):
    if '默认值' not in str(w):
        print(f"  [警告] {w}")

# ============================================================
# 方案三：球形储罐
# ============================================================
print("\n" + "=" * 70)
print("  方案三：球形储罐 (Spherical Tank)")
print("=" * 70)

result_s = _full_design(
    V=V,
    tank_type="spherical",
    fluid=fluid,
    T=T,
    P=P,
)
results["spherical"] = result_s

print(f"\n--- 几何尺寸 ---")
print(f"  球体直径 D: {result_s.get('D', 0):.2f} m")
print(f"  总容积: {result_s.get('V_total', 0):.2f} m³")
print(f"  设计容积(含15%气相裕量): {result_s.get('design_volume_m3', 0):.2f} m³")

print(f"\n--- 流体信息 ---")
fp = result_s.get('fluid_properties', {})
print(f"  液相密度: {fp.get('rho_liquid_kg_m3', 'N/A')} kg/m³")
print(f"  液体质量: {result_s.get('liquid_mass_tonnes', 0):.1f} 吨")

print(f"\n--- 基础荷载 ---")
fnd = result_s.get('foundation', {})
print(f"  总重量: {fnd.get('total_weight_kN', 0):.1f} kN")
print(f"  基底压力: {fnd.get('base_pressure_kPa', 0):.1f} kPa")

print(f"\n--- 蒸发损耗 ---")
bl = result_s.get('breathing_losses', {})
if bl:
    print(f"  每次呼吸损耗: {bl.get('loss_per_breath_kg', 0):.4f} kg")
    print(f"  年损耗量: {bl.get('annual_loss_kg', 0):.1f} kg")

print(f"\n--- 约束校核 ---")
val = result_s.get('validation', {})
print(f"  通过率: {val.get('pass_rate', 0)*100:.0f}% ({val.get('n_pass', 0)}/{val.get('n_total', 0)})")
for check in val.get('checks', []):
    status = "PASS" if check.get('pass') else "FAIL"
    print(f"  [{status}] {check.get('name')}: {check.get('value')}")
for w in result_s.get('warnings', []):
    if '默认值' not in str(w):
        print(f"  [警告] {w}")

# ============================================================
# 三种方案对比
# ============================================================
print("\n" + "=" * 70)
print("  三种方案综合对比")
print("=" * 70)

# 提取关键参数
D_v = result_v.get('D', 0)
H_v = result_v.get('H', 0)
D_h = result_h.get('D', 0)
L_h = result_h.get('L', 0)
D_s = result_s.get('D', 0)

# 计算占地面积
area_v = math.pi * (D_v/2)**2 if D_v > 0 else 0
area_h = D_h * L_h if D_h > 0 and L_h > 0 else 0  # 卧式罐投影面积（近似）
area_s = math.pi * (D_s/2)**2 if D_s > 0 else 0

# 计算表面积（估算）
SA_v = result_v.get('surface_area_m2', math.pi * D_v * H_v + 2 * math.pi * (D_v/2)**2 if D_v > 0 else 0)
SA_h = result_h.get('surface_area_m2', math.pi * D_h * L_h + 2 * math.pi * (D_h/2)**2 if D_h > 0 else 0)
SA_s = result_s.get('surface_area_m2', math.pi * D_s**2 if D_s > 0 else 0)

print(f"\n{'参数':<25} {'立式罐':<15} {'卧式罐':<15} {'球罐':<15}")
print("-" * 70)
print(f"{'罐体直径 D (m)':<25} {D_v:<15.2f} {D_h:<15.2f} {D_s:<15.2f}")
print(f"{'高度/长度 (m)':<25} {H_v:<15.2f} {L_h:<15.2f} {'-':<15}")
print(f"{'长径比/高径比':<25} {result_v.get('H_D_ratio', result_v.get('aspect_ratio', 0)):<15.2f} {result_h.get('L_D_ratio', result_h.get('aspect_ratio', 0)):<15.2f} {'1.0':<15}")
print(f"{'总容积 (m³)':<25} {result_v.get('V_total', 0):<15.2f} {result_h.get('V_total', 0):<15.2f} {result_s.get('V_total', 0):<15.2f}")
print(f"{'液体质量 (吨)':<25} {result_v.get('liquid_mass_tonnes', 0):<15.1f} {result_h.get('liquid_mass_tonnes', 0):<15.1f} {result_s.get('liquid_mass_tonnes', 0):<15.1f}")
print(f"{'占地面积 (m²)':<25} {area_v:<15.1f} {area_h:<15.1f} {area_s:<15.1f}")
print(f"{'表面积 (m²)':<25} {SA_v:<15.1f} {SA_h:<15.1f} {SA_s:<15.1f}")

fnd_v = result_v.get('foundation', {})
fnd_h = result_h.get('foundation', {})
fnd_s = result_s.get('foundation', {})
print(f"{'总重量 (kN)':<25} {fnd_v.get('total_weight_kN', 0):<15.1f} {fnd_h.get('total_weight_kN', 0):<15.1f} {fnd_s.get('total_weight_kN', 0):<15.1f}")
print(f"{'基底压力 (kPa)':<25} {fnd_v.get('base_pressure_kPa', 0):<15.1f} {fnd_h.get('base_pressure_kPa', 0):<15.1f} {fnd_s.get('base_pressure_kPa', 0):<15.1f}")

bl_v = result_v.get('breathing_losses', {})
bl_h = result_h.get('breathing_losses', {})
bl_s = result_s.get('breathing_losses', {})
print(f"{'年蒸发损耗 (kg)':<25} {bl_v.get('annual_loss_kg', 0):<15.1f} {bl_h.get('annual_loss_kg', 0):<15.1f} {bl_s.get('annual_loss_kg', 0):<15.1f}")

val_v = result_v.get('validation', {})
val_h = result_h.get('validation', {})
val_s = result_s.get('validation', {})
print(f"{'校核通过率':<25} {val_v.get('pass_rate', 0)*100:<14.0f}% {val_h.get('pass_rate', 0)*100:<14.0f}% {val_s.get('pass_rate', 0)*100:<14.0f}%")

print("\n" + "-" * 70)
print("方案建议:")
print("  • 立式罐: 占地面积小，适合场地有限的情况，制造简单，成本较低")
print("  • 卧式罐: 适合低压/常压储存，便于运输和安装，可埋地")
print("  • 球罐: 表面积最小，蒸发损耗最低，但制造成本高，适合中高压储存")
print("-" * 70)
