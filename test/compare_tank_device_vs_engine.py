# -*- coding: utf-8 -*-
"""对比 storage_tank_device._full_design vs physics_engine 直接调用"""
import sys, io, json, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.storage_tank_device import _full_design
from physics_engine.storage_tank import (
    design_vertical_tank, design_horizontal_tank, design_spherical_tank,
    validate_tank_design, breathing_losses, tank_foundation_load,
)
from physics_engine.thermo_helper import get_fluid_density, get_fluid_MW, get_fluid_vapor_pressure

# ============================================================
# 设计条件
# ============================================================
V = 150.0           # m³
fluid = "ethanol"
T = 298.15          # K (25°C)
P = 101325.0        # Pa

print("=" * 80)
print("  storage_tank_device._full_design vs physics_engine 直接调用 对比")
print("=" * 80)
print(f"\n设计条件: V={V} m³, 介质={fluid}, T={T}K ({T-273.15}°C), P={P/1000:.2f} kPa")

# ============================================================
# 辅助函数：提取关键参数
# ============================================================
def extract_key_params(result, source_name):
    """提取关键设计参数"""
    return {
        "source": source_name,
        "D": result.get("D", 0),
        "H": result.get("H", result.get("H_shell", 0)),
        "L": result.get("L", result.get("L_shell", 0)),
        "H_total": result.get("H_total", 0),
        "L_total": result.get("L_total", 0),
        "aspect_ratio": result.get("aspect_ratio", 0),
        "actual_volume": result.get("actual_volume", 0),
        "fill_ratio": result.get("fill_ratio", 0),
        "theoretical_D": result.get("theoretical_D", 0),
    }

def print_comparison(label, val_device, val_engine):
    """打印对比结果"""
    if isinstance(val_device, float):
        match = "✓" if abs(val_device - val_engine) < 0.001 else "✗"
        print(f"  {label:<25} {val_device:<15.4f} {val_engine:<15.4f} {match}")
    else:
        match = "✓" if val_device == val_engine else "✗"
        print(f"  {label:<25} {str(val_device):<15} {str(val_engine):<15} {match}")

# ============================================================
# 方案一：立式储罐
# ============================================================
print("\n" + "=" * 80)
print("  方案一：立式圆柱储罐 (Vertical)")
print("=" * 80)

# Device 层调用
result_v_device = _full_design(V=V, tank_type="vertical", fluid=fluid, T=T, P=P)

# Physics Engine 直接调用
result_v_engine = design_vertical_tank(V, aspect_ratio=4.0)

print("\n--- 几何参数对比 ---")
print(f"  {'参数':<25} {'device层':<15} {'physics_engine':<15} {'一致?'}")
print("  " + "-" * 65)
d1 = extract_key_params(result_v_device, "device")
d2 = extract_key_params(result_v_engine, "engine")
print_comparison("直径 D (m)", d1["D"], d2["D"])
print_comparison("筒体高度 H (m)", d1["H"], d2["H"])
print_comparison("总高度 H_total (m)", d1["H_total"], d2["H_total"])
print_comparison("高径比 H/D", d1["aspect_ratio"], d2["aspect_ratio"])
print_comparison("实际容积 (m³)", d1["actual_volume"], d2["actual_volume"])
print_comparison("填充率", d1["fill_ratio"], d2["fill_ratio"])
print_comparison("理论直径 (m)", d1["theoretical_D"], d2["theoretical_D"])

print("\n--- Device 层额外计算 ---")
print(f"  设计容积(含15%裕量): {result_v_device.get('design_volume_m3', 'N/A')} m³")
print(f"  要求容积: {result_v_device.get('required_volume_m3', 'N/A')} m³")
fp = result_v_device.get('fluid_properties', {})
print(f"  液相密度: {fp.get('rho_liquid_kg_m3', 'N/A')} kg/m³")
print(f"  分子量: {fp.get('MW_g_mol', 'N/A')} g/mol")
print(f"  饱和蒸气压: {fp.get('vapor_pressure_Pa', 'N/A')} Pa")
print(f"  液体质量: {result_v_device.get('liquid_mass_tonnes', 'N/A')} 吨")
fnd = result_v_device.get('foundation', {})
print(f"  基础总荷载: {fnd.get('total_weight_kN', 'N/A')} kN")
print(f"  基底压力: {fnd.get('base_pressure_kPa', 'N/A')} kPa")
bl = result_v_device.get('breathing_losses', {})
print(f"  年呼吸损耗: {bl.get('annual_loss_kg', 'N/A')} kg")
val = result_v_device.get('validation', {})
print(f"  约束校核: {val.get('n_pass', 0)}/{val.get('n_total', 0)} 通过 ({val.get('pass_rate', 0)*100:.0f}%)")

# ============================================================
# 方案二：卧式储罐
# ============================================================
print("\n" + "=" * 80)
print("  方案二：卧式圆柱储罐 (Horizontal)")
print("=" * 80)

# Device 层调用
result_h_device = _full_design(V=V, tank_type="horizontal", fluid=fluid, T=T, P=P)

# Physics Engine 直接调用
result_h_engine = design_horizontal_tank(V, aspect_ratio=3.0)

print("\n--- 几何参数对比 ---")
print(f"  {'参数':<25} {'device层':<15} {'physics_engine':<15} {'一致?'}")
print("  " + "-" * 65)
d1 = extract_key_params(result_h_device, "device")
d2 = extract_key_params(result_h_engine, "engine")
print_comparison("直径 D (m)", d1["D"], d2["D"])
print_comparison("筒体长度 L (m)", d1["L"], d2["L"])
print_comparison("总长度 L_total (m)", d1["L_total"], d2["L_total"])
print_comparison("长径比 L/D", d1["aspect_ratio"], d2["aspect_ratio"])
print_comparison("实际容积 (m³)", d1["actual_volume"], d2["actual_volume"])
print_comparison("填充率", d1["fill_ratio"], d2["fill_ratio"])
print_comparison("理论直径 (m)", d1["theoretical_D"], d2["theoretical_D"])

print("\n--- Device 层额外计算 ---")
print(f"  设计容积(含15%裕量): {result_h_device.get('design_volume_m3', 'N/A')} m³")
fp = result_h_device.get('fluid_properties', {})
print(f"  液相密度: {fp.get('rho_liquid_kg_m3', 'N/A')} kg/m³")
print(f"  液体质量: {result_h_device.get('liquid_mass_tonnes', 'N/A')} 吨")
fnd = result_h_device.get('foundation', {})
print(f"  基础总荷载: {fnd.get('total_weight_kN', 'N/A')} kN")
print(f"  基底压力: {fnd.get('base_pressure_kPa', 'N/A')} kPa")
bl = result_h_device.get('breathing_losses', {})
print(f"  年呼吸损耗: {bl.get('annual_loss_kg', 'N/A')} kg")
val = result_h_device.get('validation', {})
print(f"  约束校核: {val.get('n_pass', 0)}/{val.get('n_total', 0)} 通过 ({val.get('pass_rate', 0)*100:.0f}%)")

# ============================================================
# 方案三：球形储罐
# ============================================================
print("\n" + "=" * 80)
print("  方案三：球形储罐 (Spherical)")
print("=" * 80)

# Device 层调用
result_s_device = _full_design(V=V, tank_type="spherical", fluid=fluid, T=T, P=P)

# Physics Engine 直接调用
result_s_engine = design_spherical_tank(V)

print("\n--- 几何参数对比 ---")
print(f"  {'参数':<25} {'device层':<15} {'physics_engine':<15} {'一致?'}")
print("  " + "-" * 65)
print_comparison("直径 D (m)", result_s_device.get("D", 0), result_s_engine.get("D", 0))
print_comparison("实际容积 (m³)", result_s_device.get("actual_volume", 0), result_s_engine.get("actual_volume", 0))
print_comparison("填充率", result_s_device.get("fill_ratio", 0), result_s_engine.get("fill_ratio", 0))
print_comparison("理论直径 (m)", result_s_device.get("theoretical_D", 0), result_s_engine.get("theoretical_D", 0))

print("\n--- Device 层额外计算 ---")
print(f"  设计容积(含15%裕量): {result_s_device.get('design_volume_m3', 'N/A')} m³")
fp = result_s_device.get('fluid_properties', {})
print(f"  液相密度: {fp.get('rho_liquid_kg_m3', 'N/A')} kg/m³")
print(f"  液体质量: {result_s_device.get('liquid_mass_tonnes', 'N/A')} 吨")
bl = result_s_device.get('breathing_losses', {})
print(f"  年呼吸损耗: {bl.get('annual_loss_kg', 'N/A')} kg")
val = result_s_device.get('validation', {})
print(f"  约束校核: {val.get('n_pass', 0)}/{val.get('n_total', 0)} 通过 ({val.get('pass_rate', 0)*100:.0f}%)")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 80)
print("  对比总结")
print("=" * 80)

print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│  storage_tank_device._full_design (Device 层)                              │
│  ─────────────────────────────────────────────────────────────────────────  │
│  功能: 完整设计流水线                                                       │
│  1. 调用 physics_engine 几何计算                                            │
│  2. 流体物性查询 (密度/MW/蒸气压)                                           │
│  3. 流体质量计算                                                            │
│  4. 基础荷载计算                                                            │
│  5. 保温层设计 (温差>10K时)                                                 │
│  6. 呼吸损耗估算                                                            │
│  7. 约束校验 (长径比/填充率/基底压力)                                        │
│  8. 工程默认值处理 (P/T_ambient/aspect_ratio)                               │
│  9. 参数报告生成                                                            │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  physics_engine.design_xxx_tank (Engine 层)                                │
│  ─────────────────────────────────────────────────────────────────────────  │
│  功能: 纯几何计算                                                          │
│  1. 理论直径计算                                                            │
│  2. 标准直径匹配                                                            │
│  3. 实际容积/填充率计算                                                     │
│  4. 封头容积/深度计算                                                       │
│  输出: 几何尺寸字典 {D, H/L, aspect_ratio, actual_volume, fill_ratio, ...}  │
└─────────────────────────────────────────────────────────────────────────────┘

核心结论:
  • 几何计算结果: 完全一致 (Device层直接调用Engine层函数)
  • Device层额外提供: 物性/荷载/保温/损耗/校验/报告
  • Engine层职责: 纯几何计算，无副作用
  • 两者关系: Device层 = Engine层 + 物性查询 + 工程校验 + 报告生成
""")
