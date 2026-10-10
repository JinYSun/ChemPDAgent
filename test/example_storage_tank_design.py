# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
"""
常压储罐设计计算
=================
已知条件：
  - 液体密度 ρ = 1200 kg/m³
  - 设计有效容积 V_eff = 500 m³（含 10% 气相空间裕量）
  - 高径比 H/D = 1.0 ~ 1.2
  - 材料许用应力 [σ] = 140 MPa
  - 腐蚀裕量 C = 2 mm

求解：罐体直径 D、高度 H、壁厚
"""
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from physics_engine.storage_tank import (
    vertical_cylindrical_volume,
    design_vertical_tank,
    tank_shell_thickness,
    tank_foundation_load,
    validate_tank_design,
)

print("=" * 70)
print("          常压立式圆柱储罐设计计算")
print("=" * 70)

# ============================================================
# 1. 已知参数
# ============================================================
rho = 1200.0          # 液体密度 kg/m³
V_liquid_req = 500.0  # 液相储存需求 m³ (用户要求的储存量)
vapor_fraction = 0.10 # 气相空间裕量 10%
H_D_min = 1.0         # 高径比下限
H_D_max = 1.2         # 高径比上限
allowable_stress = 140e6   # 许用应力 Pa (140 MPa)
corrosion_allowance = 0.002  # 腐蚀裕量 m (2 mm)
joint_efficiency = 0.85      # 焊缝系数（Q235-B 双面焊）
g = 9.81

print(f"\n【已知条件】")
print(f"  液体密度 rho          = {rho} kg/m3")
print(f"  液相储存需求 V_liq    = {V_liquid_req} m3")
print(f"  气相空间裕量          = {vapor_fraction*100:.0f}%")
print(f"  高径比 H/D           = {H_D_min} ~ {H_D_max}")
print(f"  材料许用应力 [σ]     = {allowable_stress/1e6:.0f} MPa")
print(f"  腐蚀裕量 C           = {corrosion_allowance*1000:.0f} mm")
print(f"  焊缝系数 φ           = {joint_efficiency}")

# ============================================================
# 2. 液相容积
# ============================================================
V_total = V_liquid_req / (1 - vapor_fraction)
print(f"\n{'='*70}")
print(f"【步骤1】总容积反推")
print(f"  V_total = V_liquid / (1 - 气相分率)")
print(f"  V_total = {V_liquid_req} / (1 - {vapor_fraction}) = {V_total:.2f} m3")

# ============================================================
# 3. 理论直径与高度（手算验证）
# ============================================================
print(f"\n{'='*70}")
print(f"【步骤2】理论直径与高度计算（H/D = 1.0 ~ 1.2）")

# 取 H/D 中间值 1.1 进行理论计算
H_D_design = (H_D_min + H_D_max) / 2  # 1.1
# V = π/4 * D² * H = π/4 * D² * (H/D * D) = π/4 * (H/D) * D³
# D = (4V / (π * H/D))^(1/3)
D_theo = (4 * V_total / (math.pi * H_D_design)) ** (1/3)
H_theo = H_D_design * D_theo

print(f"  取 H/D = {H_D_design}")
print(f"  由 V = pi/4 * D^2 * H = pi/4 * (H/D) * D^3")
print(f"  → D = (4V / (π × H/D))^(1/3)")
print(f"  理论直径 D_theo = {D_theo:.4f} m")
print(f"  理论高度 H_theo = {H_theo:.4f} m")

# 分别计算 H/D=1.0 和 H/D=1.2 的范围
for ratio_name, ratio_val in [("H/D=1.0", H_D_min), ("H/D=1.2", H_D_max)]:
    D_range = (4 * V_total / (math.pi * ratio_val)) ** (1/3)
    H_range = ratio_val * D_range
    print(f"  {ratio_name}: D = {D_range:.4f} m, H = {H_range:.4f} m")

# ============================================================
# 4. 调用 design_vertical_tank 匹配标准直径
# ============================================================
print(f"\n{'='*70}")
print(f"【步骤3】匹配标准直径系列（调用 design_vertical_tank）")

# 500 m3 是液相储存需求，需反推总容积
# V_total = 500 / 0.9 = 555.56 m3 (预留 10% 气相空间)
V_total_required = V_total  # 555.56 m3

# 以总容积 500 m3 匹配标准直径，确保 H/D 在 1.0~1.2
STANDARD_DS = [3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30]
best_D_std = None
best_H_std = None
best_aspect = None
best_vol = None

print(f"  标准直径系列: {STANDARD_DS} ... (m)")
D_min_range = (4*V_total_required/(math.pi*H_D_max))**(1/3)
D_max_range = (4*V_total_required/(math.pi*H_D_min))**(1/3)
print(f"  理论直径范围: D(H/D=1.2)={D_min_range:.4f} m ~ D(H/D=1.0)={D_max_range:.4f} m")
print(f"  (基于总容积 V_total = {V_total_required} m3 计算)")
print(f"\n  遍历标准直径，寻找 H/D 在 1.0~1.2 且总容积 >= {V_total_required} m3 的方案：")

for D_try in STANDARD_DS:
    # 以总容积反推高度
    H_try = V_total_required / (math.pi * (D_try/2)**2)
    aspect_try = H_try / D_try
    V_try = math.pi * (D_try/2)**2 * H_try
    # 检查：总容积 >= 500 且 H/D 在范围内
    ok_vol = V_try >= V_total_required - 0.01
    ok_aspect = H_D_min <= aspect_try <= H_D_max
    marker = ""
    if ok_aspect and ok_vol:
        marker = " <-- 满足"
    elif ok_vol:
        marker = " (容积OK, H/D超)"
    elif ok_aspect:
        marker = " (H/D OK, 容积不足)"
    print(f"    D={D_try:>2} m -> H={H_try:.4f} m, H/D={aspect_try:.4f}, V={V_try:.2f} m3{marker}")
    if ok_aspect and ok_vol and best_D_std is None:
        best_D_std = D_try
        best_H_std = H_try
        best_aspect = aspect_try
        best_vol = V_try

if best_D_std is None:
    print("\n  未找到精确满足条件的标准直径，取最接近方案")
    best_D_std = 8
    best_H_std = V_total_required / (math.pi * (best_D_std/2)**2)
    best_aspect = best_H_std / best_D_std
    best_vol = math.pi * (best_D_std/2)**2 * best_H_std

D_std = best_D_std
H_std = best_H_std
aspect_actual = best_aspect
V_actual = best_vol
fill_ratio = V_liquid_req / V_actual

print(f"\n  >>> 选定标准直径 D = {D_std} m")
print(f"  >>> 对应高度 H = {H_std:.4f} m")
print(f"  >>> 实际高径比 H/D = {aspect_actual:.4f}")
print(f"  >>> 实际总容积 V = {V_actual:.2f} m3")
print(f"  >>> 填充率 = {fill_ratio:.4f}")

# 验证气相空间
vapor_vol = V_actual - V_liquid_req
vapor_pct = vapor_vol / V_actual * 100
print(f"\n  气相空间体积 = {V_actual:.2f} - {V_liquid_req:.1f} = {vapor_vol:.2f} m3")
print(f"  气相空间占比 = {vapor_pct:.1f}% (要求 >= 10%)")
if vapor_pct >= 10:
    print(f"  [OK] 气相空间满足要求")
else:
    print(f"  [NG] 气相空间不足，需增大罐体")

# ============================================================
# 5. 罐壁厚度计算（GB 50341 / API 650）
# ============================================================
print(f"\n{'='*70}")
print(f"【步骤4】罐壁厚度计算（GB 50341 / API 650 变点法）")

result_thickness = tank_shell_thickness(
    D=D_std,
    H=H_std,
    liquid_density=rho,
    allowable_stress=allowable_stress,
    joint_efficiency=joint_efficiency,
    corrosion_allowance=corrosion_allowance,
)

print(f"  公式：t_calc = P × D / (2[σ]φ - P)")
print(f"  t_design = t_calc + C_corr")
print(f"  其中 P = rho * g * h_liquid（静水压力）")
print(f"  罐体分 {result_thickness['n_courses']} 圈板（每圈 1m 高度）\n")

print(f"  {'圈板':>4} | {'液柱高度(m)':>10} | {'静压P(Pa)':>12} | {'t_calc(mm)':>10} | {'t_design(mm)':>12} | {'选定(mm)':>8}")
print(f"  {'-'*4}-+-{'-'*10}-+-{'-'*12}-+-{'-'*10}-+-{'-'*12}-+-{'-'*8}")
for c in result_thickness["courses"]:
    print(f"  {c['course']:>4} | {c['h_liquid']:>10.2f} | {c['pressure_Pa']:>12.0f} | {c['t_calc_m']*1000:>10.3f} | {c['t_design_m']*1000:>12.3f} | {c['t_selected_mm']:>8.1f}")

print(f"\n  底层壁厚: {result_thickness['t_bottom_shell_mm']:.1f} mm")
print(f"  顶层壁厚: {result_thickness['t_top_shell_mm']:.1f} mm")
print(f"  罐底板厚: {result_thickness['t_bottom_plate_mm']:.1f} mm")

# ============================================================
# 6. 基础荷载校核
# ============================================================
print(f"\n{'='*70}")
print(f"【步骤5】基础荷载校核")

result_foundation = tank_foundation_load(
    D=D_std, H=H_std, liquid_density=rho
)

print(f"  液体重量 = {result_foundation['liquid_weight']/1000:.1f} kN")
print(f"  罐壁重量 = {result_foundation['shell_weight']/1000:.1f} kN")
print(f"  罐底重量 = {result_foundation['bottom_weight']/1000:.1f} kN")
print(f"  总重量   = {result_foundation['total_weight']/1000:.1f} kN")
print(f"  基底压力 = {result_foundation['base_pressure']/1000:.1f} kPa")
if result_foundation['base_pressure']/1000 < 200:
    print(f"  [OK] 基底压力 {result_foundation['base_pressure']/1000:.1f} kPa < 200 kPa（地基承载力满足）")
else:
    print(f"  [NG] 基底压力 {result_foundation['base_pressure']/1000:.1f} kPa 超限")

# ============================================================
# 7. 设计校验
# ============================================================
print(f"\n{'='*70}")
print(f"【步骤6】设计校验汇总")

design_summary = {
    "D": D_std,
    "H": H_std,
    "fill_ratio": fill_ratio,
    "base_pressure": result_foundation["base_pressure"],
}
result_validate = validate_tank_design(design_summary, tank_type="vertical_cylindrical")

for chk in result_validate["checks"]:
    status = "[OK]" if chk["pass"] else "[NG]"
    print(f"  {status} | {chk['name']}: {chk['value']:.4f} (限值: {chk['limit']})")
for w in result_validate["warnings"]:
    print(f"  ⚠ {w}")
print(f"  校验通过率: {result_validate['n_pass']}/{result_validate['n_total']}")

# ============================================================
# 8. 最终结论
# ============================================================
print(f"\n{'='*70}")
print(f"【计算结论】")
print(f"  罐体直径 D = {D_std} m（标准系列）")
print(f"  罐体高度 H = {H_std:.2f} m")
print(f"  实际高径比 H/D = {aspect_actual:.2f}（要求 1.0~1.2）")
print(f"  实际总容积 V = {V_actual:.2f} m³")
print(f"  液相容积 = {V_liquid_req:.1f} m³，气相空间 = {vapor_vol:.2f} m³（{vapor_pct:.1f}%）")
print(f"  底层壁厚 = {result_thickness['t_bottom_shell_mm']:.0f} mm")
print(f"  顶层壁厚 = {result_thickness['t_top_shell_mm']:.0f} mm")
print(f"  罐底板厚 = {result_thickness['t_bottom_plate_mm']:.0f} mm")
print(f"  设计标准 = {result_thickness['standard']}")
print(f"{'='*70}")
