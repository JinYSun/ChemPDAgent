# -*- coding: utf-8 -*-
"""
卧式储罐设计示例 — 150 m³ 乙醇常温储存
完整展示：物性方法 → 中间计算过程 → 设计方案结论
"""
import sys, os, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from physics_engine.thermo_helper import get_chemical_properties
from physics_engine.storage_tank import (
    effective_capacity_horizontal,
    horizontal_cylindrical_volume,
    head_volume,
    head_depth,
    tank_foundation_load,
    _STANDARD_HORIZONTAL_DIAMETERS,
)
import math

# ============================================================
# 物性计算 — 调用 thermo_helper.get_chemical_properties
# ============================================================
T_ambient = 298.15   # K (25°C 常温)
P_ambient = 101325.0 # Pa (常压)
component = "ethanol"

props = get_chemical_properties(component, T=T_ambient, P=P_ambient)

# ============================================================
# 设计条件
# ============================================================
V_liquid_required = 150.0     # m³, 有效液相储存容积
fill_fraction = 0.85          # 最高操作液位 85%D
vapor_fraction = 0.15         # 预留 15% 气相空间
head_type = "ellipsoidal"     # 标准椭圆封头 (2:1)

print("=" * 70)
print("  卧式储罐设计 — 150 m³ 乙醇常温储存")
print("=" * 70)

# ============================================================
# [物性] 乙醇物性数据 (thermo + CoolProp 双引擎)
# ============================================================
print(f"\n[物性数据] {component} @ T={T_ambient}K ({T_ambient-273.15:.1f}°C), P={P_ambient/1000:.1f}kPa")
print(f"  数据来源: thermo_helper.get_chemical_properties (thermo + CoolProp 双引擎)")
print(f"  分子量 MW:              {props['MW']:.2f} g/mol")
print(f"  常压沸点 Tb:            {props['Tb']:.2f} K ({props['Tb']-273.15:.1f}°C)")
print(f"  临界温度 Tc:            {props['Tc']:.2f} K")
print(f"  临界压力 Pc:            {props['Pc']/1e5:.2f} bar")
print(f"  液相密度 rho_L:         {props['rho_L']:.1f} kg/m³  [TH]")
print(f"  气相密度 rho_G:         {props['rho_G']:.4f} kg/m³  [CP]")
print(f"  液相黏度 mu_L:          {props['mu_L']*1e6:.1f} μPa·s  [TH]")
print(f"  液相导热系数 k_L:       {props['k_L']*1e3:.2f} mW/m/K  [TH]")
print(f"  饱和蒸气压 Psat:        {props['Psat']:.0f} Pa  [TH]")
print(f"  汽化潜热 Hvap:          {props['Hvap']:.0f} J/mol ({props['Hvap']/props['MW']*1000:.0f} J/g)  [TH]")
print(f"  表面张力 sigma:         {props['sigma']*1e3:.2f} mN/m  [CP]")
print(f"  液相摩尔体积 Vm_L:      {props['Vm_L']*1e6:.2f} cm³/mol  [TH]")

# 乙醇设计密度（取物性计算值）
ethanol_density = props['rho_L']  # kg/m³

# ============================================================
# [步骤1] 设计条件汇总
# ============================================================
V_total = V_liquid_required / (1 - vapor_fraction)
print(f"\n[步骤1] 设计条件")
print(f"  有效液相容积:         {V_liquid_required:.1f} m³")
print(f"  储存介质:             乙醇 ({component})")
print(f"  操作温度:             {T_ambient}K ({T_ambient-273.15:.1f}°C)")
print(f"  操作压力:             常压 ({P_ambient/1000:.1f} kPa)")
print(f"  介质密度 (物性计算):  {ethanol_density:.1f} kg/m³")
print(f"  气相预留比例:         {vapor_fraction*100:.0f}%")
print(f"  设计总几何容积:       V_total = {V_liquid_required:.1f}/{1-vapor_fraction:.2f} = {V_total:.2f} m³")
print(f"  最高操作液位:         {fill_fraction*100:.0f}% D")
print(f"  封头类型:             标准椭圆封头 (2:1)")

# ============================================================
# [步骤2] 标准直径系列遍历，筛选满足 L/D总 ∈ [3,6] 的方案
# ============================================================
print(f"\n[步骤2] 标准直径系列遍历")
print(f"  标准直径: {_STANDARD_HORIZONTAL_DIAMETERS}")
print(f"  约束条件: 3 <= L/D总 <= 6, V_effective >= {V_liquid_required} m³")
print()
print(f"  {'序号':>4} | {'D/m':>5} | {'L筒体/m':>7} | {'L总/m':>6} | "
      f"{'L/D总':>5} | {'V有效/m³':>8} | {'判定':>4}")
print(f"  {'-'*60}")

candidates = []
for idx, D_std in enumerate(_STANDARD_HORIZONTAL_DIAMETERS):
    V_one_head = head_volume(D_std, head_type)
    V_heads = 2 * V_one_head
    V_shell_needed = V_total - V_heads
    L_calc = max(2.0, V_shell_needed / (math.pi * (D_std / 2) ** 2)) if V_shell_needed > 0 else 2.0

    hd = head_depth(D_std, head_type)
    L_total_calc = L_calc + 2 * hd
    LD_total = L_total_calc / D_std if D_std > 0 else 0

    vol_result = horizontal_cylindrical_volume(D_std, L_calc, head_type=head_type)
    eff = effective_capacity_horizontal(D_std, L_calc, fill_fraction=fill_fraction, head_type=head_type)
    V_eff = eff["V_effective"]

    ok = 3 <= LD_total <= 6 and V_eff >= V_liquid_required
    mark = "OK" if ok else "--"
    print(f"  {idx+1:>4} | {D_std:>5.1f} | {L_calc:>7.2f} | {L_total_calc:>6.2f} | "
          f"{LD_total:>5.2f} | {V_eff:>8.2f} |  {mark}")

    if ok:
        candidates.append({
            "D": D_std, "L_shell": L_calc, "L_total": L_total_calc,
            "LD_total": LD_total, "V_total": vol_result["total_volume"],
            "V_effective": V_eff,
        })

print(f"  {'-'*60}")

if not candidates:
    print("\n  未找到满足条件的方案！")
    sys.exit(1)

# 选取 L/D总 最接近 4.0 的方案
best = min(candidates, key=lambda c: abs(c["LD_total"] - 4.0))
D = best["D"]
L_shell = best["L_shell"]
print(f"\n  选取方案: D = {D} m (L/D总 = {best['LD_total']:.2f}, 最接近 4.0)")

# ============================================================
# [步骤3] 封头容积计算
# ============================================================
V_head = head_volume(D, head_type)
print(f"\n[步骤3] 封头容积计算")
print(f"  函数: head_volume(D={D}, head_type='{head_type}')")
print(f"  公式: V_head = π/24 × D³")
print(f"  V_head = π/24 × {D}³ = {V_head:.4f} m³")

# ============================================================
# [步骤4] 筒体长度计算
# ============================================================
V_heads_total = 2 * V_head
V_shell_needed = V_total - V_heads_total
A_cross = math.pi * (D / 2) ** 2
print(f"\n[步骤4] 筒体长度计算")
print(f"  函数: horizontal_cylindrical_volume 反推")
print(f"  公式: L_shell = (V_total - 2×V_head) / (π/4 × D²)")
print(f"  L_shell = ({V_total:.2f} - 2×{V_head:.4f}) / (π/4 × {D}²)")
print(f"  L_shell = {V_shell_needed:.4f} / {A_cross:.4f} = {L_shell:.2f} m")

# ============================================================
# [步骤5] 总长计算
# ============================================================
hd = head_depth(D, head_type)
L_total = L_shell + 2 * hd
LD_total = L_total / D
print(f"\n[步骤5] 总长计算")
print(f"  函数: head_depth(D={D}, head_type='{head_type}')")
print(f"  head_depth = D/4 = {hd:.3f} m (椭圆封头)")
print(f"  L_total = L_shell + 2×head_depth = {L_shell:.2f} + 2×{hd:.3f} = {L_total:.2f} m")
print(f"  总长径比 L/D = {L_total:.2f}/{D} = {LD_total:.2f}")

# ============================================================
# [步骤6] 总容积计算
# ============================================================
vol_result = horizontal_cylindrical_volume(D, L_shell, head_type=head_type)
V_shell_vol = vol_result["shell_volume"]
V_head_total_vol = vol_result["head_volume_total"]
V_total_actual = vol_result["total_volume"]
print(f"\n[步骤6] 总容积计算")
print(f"  函数: horizontal_cylindrical_volume(D={D}, L={L_shell:.2f}, head_type='{head_type}')")
print(f"  V_shell = π/4 × {D}² × {L_shell:.2f} = {V_shell_vol:.2f} m³")
print(f"  V_head_total = 2×{V_head:.4f} = {V_head_total_vol:.4f} m³")
print(f"  V_total = {V_shell_vol:.2f} + {V_head_total_vol:.2f} = {V_total_actual:.2f} m³")

# ============================================================
# [步骤7] 有效容量计算 (基于 85% 液位)
# ============================================================
eff = effective_capacity_horizontal(
    D, L_shell,
    fill_fraction=fill_fraction, head_type=head_type,
    n_level_points=20,
)
V_effective = eff["V_effective"]
vapor_space = eff["vapor_space"]
effective_frac = eff["effective_fraction"]
print(f"\n[步骤7] 有效容量计算")
print(f"  函数: effective_capacity_horizontal(D={D}, L={L_shell:.2f}, fill_fraction={fill_fraction})")
print(f"  最高操作液位 h_max = {fill_fraction}×{D} = {eff['h_max']:.3f} m")
print(f"  V_effective = {V_effective:.2f} m³ (85%液位下液相容积)")
print(f"  有效占比 = {effective_frac*100:.1f}% (卧式罐弓形截面非线性: 85%液位→{effective_frac*100:.1f}%容积)")
print(f"  vapor_space = V_total - V_effective = {V_total_actual:.2f} - {V_effective:.2f} = {vapor_space:.2f} m³")

# ============================================================
# [步骤8] 基础荷载估算
# ============================================================
found = tank_foundation_load(D=D, H=L_shell, liquid_density=ethanol_density)
print(f"\n[步骤8] 基础荷载估算")
print(f"  函数: tank_foundation_load(D={D}, H={L_shell:.2f}, rho={ethanol_density:.1f})")
print(f"  液体重量:  {found['liquid_weight']/1000:.1f} kN")
print(f"  罐壁重量:  {found['shell_weight']/1000:.1f} kN")
print(f"  罐底重量:  {found['bottom_weight']/1000:.1f} kN")
print(f"  总重量:    {found['total_weight']/1000:.1f} kN")
print(f"  基底压力:  {found['base_pressure']/1000:.1f} kPa (地基承载力一般要求 < 200 kPa)")

# ============================================================
# [步骤9] 液位-容积对照表
# ============================================================
print(f"\n[步骤9] 液位-容积对照表")
print(f"  {'液位/D':>8} | {'液位高度/m':>10} | {'液相容积/m³':>12} | {'容积占比/%':>10}")
print(f"  {'-'*52}")
for row in eff["level_table"]:
    print(f"  {row['fill_fraction']:>8.2%} | {row['liquid_height_m']:>10.3f} | "
          f"{row['liquid_volume_m3']:>12.2f} | {row['volume_pct']:>10.1f}")

# ============================================================
# [结论] 设计方案汇总
# ============================================================
print(f"\n{'='*70}")
print(f"  设计方案结论")
print(f"{'='*70}")
print(f"  储存介质:             乙醇 (常温常压)")
print(f"  介质密度:             {ethanol_density:.1f} kg/m³ (thermo_helper 计算)")
print(f"  罐内径 D:             {D} m")
print(f"  筒体长度 L:           {L_shell:.2f} m")
print(f"  封头:                 标准椭圆封头 (2:1), 深度 {hd:.3f} m ×2")
print(f"  总长 (含封头):        {L_total:.2f} m")
print(f"  总长径比 L/D:         {LD_total:.2f} (工业常规 3~6)")
print(f"  总几何容积:           {V_total_actual:.2f} m³")
print(f"  有效液相容积:         {V_effective:.2f} m³ (> {V_liquid_required} m³ 需求)")
print(f"  气相空间:             {vapor_space:.2f} m³")
print(f"  最高操作液位:         {eff['h_max']:.3f} m ({fill_fraction*100:.0f}% D)")
print(f"  总重量 (满液):        {found['total_weight']/1000:.1f} kN")
print(f"  基底压力:             {found['base_pressure']/1000:.1f} kPa")
print(f"{'='*70}")
