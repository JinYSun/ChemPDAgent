# -*- coding: utf-8 -*-
"""
立式储罐设计示例 — 500 m³ 常压储罐
完整展示：设计条件 → 直径/高度计算 → 壁厚校核 → 基础荷载 → 结论
"""
import sys, os, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from physics_engine.storage_tank import (
    design_vertical_tank,
    vertical_cylindrical_volume,
    head_volume,
    head_depth,
    tank_shell_thickness,
    tank_foundation_load,
    validate_tank_design,
    _STANDARD_VERTICAL_DIAMETERS,
)
import math

# ============================================================
# 设计条件
# ============================================================
V_effective = 500.0         # m³, 有效容积
vapor_fraction = 0.10       # 10% 气相空间裕量
aspect_ratio_target = 1.1   # 高径比 H/D 目标值 (范围 1.0~1.2)
liquid_density = 1200.0     # kg/m³
allowable_stress = 140e6    # Pa (140 MPa)
corrosion_allowance = 0.002 # m (2 mm)
head_type_top = "conical"   # 锥顶 (常压储罐标准配置)
head_type_bottom = "flat"   # 平底

print("=" * 70)
print("  立式储罐设计 — 500 m³ 常压储罐")
print("=" * 70)

# ============================================================
# [步骤1] 设计总容积
# ============================================================
V_total_design = V_effective / (1 - vapor_fraction)
print(f"\n[步骤1] 设计条件")
print(f"  有效容积:             {V_effective:.1f} m³")
print(f"  气相裕量:             {vapor_fraction*100:.0f}%")
print(f"  设计总几何容积:       V_total = {V_effective:.1f}/{1-vapor_fraction:.2f} = {V_total_design:.2f} m³")
print(f"  目标高径比 H/D:       {aspect_ratio_target} (范围 1.0~1.2)")
print(f"  介质密度:             {liquid_density} kg/m³")
print(f"  材料许用应力:         {allowable_stress/1e6:.0f} MPa")
print(f"  腐蚀裕量:             {corrosion_allowance*1000:.0f} mm")
print(f"  顶封头:               锥顶 (conical)")
print(f"  底封头:               平底 (flat)")

# ============================================================
# [步骤2] 理论直径估算
# ============================================================
# 按纯圆柱体 + 高径比目标值估算
D_theo = (4 * V_total_design / (math.pi * aspect_ratio_target)) ** (1/3)
H_theo = aspect_ratio_target * D_theo
print(f"\n[步骤2] 理论直径估算")
print(f"  公式: D = (4×V_total/(π×H/D))^(1/3)")
print(f"  D_theo = (4×{V_total_design:.2f}/(π×{aspect_ratio_target}))^(1/3) = {D_theo:.2f} m")
print(f"  H_theo = {aspect_ratio_target}×{D_theo:.2f} = {H_theo:.2f} m")

# ============================================================
# [步骤3] 标准直径系列遍历，筛选 H/D ∈ [1.0, 1.2]
# ============================================================
print(f"\n[步骤3] 标准直径系列遍历")
print(f"  标准直径: {_STANDARD_VERTICAL_DIAMETERS}")
print(f"  约束条件: 1.0 <= H_total/D <= 1.2, V_total >= {V_total_design:.2f} m³")
print()
print(f"  {'序号':>4} | {'D/m':>5} | {'H筒体/m':>7} | {'H总/m':>6} | "
      f"{'H/D总':>5} | {'V总/m³':>8} | {'判定':>4}")
print(f"  {'-'*60}")

candidates = []
for idx, D_std in enumerate(_STANDARD_VERTICAL_DIAMETERS):
    # 用 design_vertical_tank 计算
    result = design_vertical_tank(
        volume_required=V_effective,
        aspect_ratio=aspect_ratio_target,
        use_standard=True,
        head_type=head_type_top,
    )
    # 对每个标准直径手动计算
    V_one_head = head_volume(D_std, head_type_top)
    V_heads = 2 * V_one_head
    V_shell_needed = V_total_design - V_heads
    if V_shell_needed <= 0:
        H_calc = 1.0
    else:
        H_calc = V_shell_needed / (math.pi * (D_std / 2) ** 2)

    hd_top = head_depth(D_std, head_type_top)
    hd_bottom = head_depth(D_std, head_type_bottom)  # flat = 0
    H_total = H_calc + hd_top + hd_bottom
    HD_total = H_total / D_std if D_std > 0 else 0

    vol_result = vertical_cylindrical_volume(D_std, H_calc, head_type=head_type_top)
    V_total_actual = vol_result["total_volume"]

    ok = 1.0 <= HD_total <= 1.2 and V_total_actual >= V_total_design
    mark = "OK" if ok else "--"
    print(f"  {idx+1:>4} | {D_std:>5.1f} | {H_calc:>7.2f} | {H_total:>6.2f} | "
          f"{HD_total:>5.2f} | {V_total_actual:>8.2f} |  {mark}")

    if ok:
        candidates.append({
            "D": D_std, "H_shell": H_calc, "H_total": H_total,
            "HD_total": HD_total, "V_total": V_total_actual,
        })

print(f"  {'-'*60}")

if not candidates:
    print("\n  未找到满足 H/D ∈ [1.0, 1.2] 的方案，放宽至 [0.8, 1.5] 重新筛选...")
    for idx, D_std in enumerate(_STANDARD_VERTICAL_DIAMETERS):
        V_one_head = head_volume(D_std, head_type_top)
        V_heads = 2 * V_one_head
        V_shell_needed = V_total_design - V_heads
        if V_shell_needed <= 0:
            H_calc = 1.0
        else:
            H_calc = V_shell_needed / (math.pi * (D_std / 2) ** 2)
        hd_top = head_depth(D_std, head_type_top)
        H_total = H_calc + hd_top
        HD_total = H_total / D_std if D_std > 0 else 0
        vol_result = vertical_cylindrical_volume(D_std, H_calc, head_type=head_type_top)
        V_total_actual = vol_result["total_volume"]
        ok = 0.8 <= HD_total <= 1.5 and V_total_actual >= V_total_design
        if ok:
            candidates.append({
                "D": D_std, "H_shell": H_calc, "H_total": H_total,
                "HD_total": HD_total, "V_total": V_total_actual,
            })
            print(f"  候选: D={D_std}m, H/D={HD_total:.2f}")

if not candidates:
    print("  仍未找到满足条件的方案！")
    sys.exit(1)

# 选取 H/D总 最接近 1.1 的方案
best = min(candidates, key=lambda c: abs(c["HD_total"] - aspect_ratio_target))
D = best["D"]
H_shell = best["H_shell"]
print(f"\n  选取方案: D = {D} m (H/D总 = {best['HD_total']:.2f}, 最接近 {aspect_ratio_target})")

# ============================================================
# [步骤4] 封头容积计算
# ============================================================
V_head_top = head_volume(D, head_type_top)
V_head_bottom = head_volume(D, head_type_bottom)
hd_top = head_depth(D, head_type_top)
hd_bottom = head_depth(D, head_type_bottom)
print(f"\n[步骤4] 封头计算")
print(f"  顶封头 (conical): V_head = π/24 × {D}³ = {V_head_top:.4f} m³, 深度 = {hd_top:.3f} m")
print(f"  底封头 (flat):    V_head = 0 m³, 深度 = {hd_bottom:.3f} m")

# ============================================================
# [步骤5] 筒体高度计算
# ============================================================
V_heads_total = V_head_top + V_head_bottom
V_shell_needed = V_total_design - V_heads_total
A_cross = math.pi * (D / 2) ** 2
print(f"\n[步骤5] 筒体高度计算")
print(f"  公式: H_shell = (V_total - V_heads) / (π/4 × D²)")
print(f"  H_shell = ({V_total_design:.2f} - {V_heads_total:.4f}) / (π/4 × {D}²)")
print(f"  H_shell = {V_shell_needed:.4f} / {A_cross:.4f} = {H_shell:.2f} m")

# ============================================================
# [步骤6] 总高计算
# ============================================================
H_total = H_shell + hd_top + hd_bottom
HD_total = H_total / D
print(f"\n[步骤6] 总高计算")
print(f"  H_total = H_shell + hd_top + hd_bottom = {H_shell:.2f} + {hd_top:.3f} + {hd_bottom:.3f} = {H_total:.2f} m")
print(f"  总高径比 H/D = {H_total:.2f}/{D} = {HD_total:.2f}")

# ============================================================
# [步骤7] 总容积计算
# ============================================================
vol_result = vertical_cylindrical_volume(D, H_shell, head_type=head_type_top)
V_shell_vol = vol_result["shell_volume"]
V_head_total_vol = vol_result["head_volume_total"]
V_total_actual = vol_result["total_volume"]
fill_ratio = V_effective / V_total_actual
print(f"\n[步骤7] 总容积计算")
print(f"  V_shell = π/4 × {D}² × {H_shell:.2f} = {V_shell_vol:.2f} m³")
print(f"  V_heads = {V_head_total_vol:.4f} m³")
print(f"  V_total = {V_shell_vol:.2f} + {V_head_total_vol:.4f} = {V_total_actual:.2f} m³")
print(f"  填充率 = {V_effective:.1f}/{V_total_actual:.2f} = {fill_ratio:.2%}")

# ============================================================
# [步骤8] 罐壁厚度设计 (GB 50341 / API 650)
# ============================================================
thickness = tank_shell_thickness(
    D=D, H=H_total,
    liquid_density=liquid_density,
    allowable_stress=allowable_stress,
    joint_efficiency=0.85,
    corrosion_allowance=corrosion_allowance,
)
print(f"\n[步骤8] 罐壁厚度设计 (GB 50341 / API 650)")
print(f"  函数: tank_shell_thickness(D={D}, H={H_total:.2f}, rho={liquid_density})")
print(f"  圈板数: {thickness['n_courses']}")
print(f"  {'圈板':>4} | {'液柱高/m':>8} | {'计算厚度/mm':>10} | {'设计厚度/mm':>10} | {'选定厚度/mm':>10}")
print(f"  {'-'*52}")
for c in thickness["courses"]:
    print(f"  {c['course']:>4} | {c['h_liquid']:>8.2f} | {c['t_calc_m']*1000:>10.2f} | "
          f"{c['t_design_m']*1000:>10.2f} | {c['t_selected_mm']:>10.1f}")
print(f"  底层壁厚: {thickness['t_bottom_shell_mm']:.1f} mm")
print(f"  顶层壁厚: {thickness['t_top_shell_mm']:.1f} mm")
print(f"  罐底板厚: {thickness['t_bottom_plate_mm']:.1f} mm")
print(f"  罐顶板厚: {thickness['t_roof_mm']:.1f} mm")

# ============================================================
# [步骤9] 基础荷载估算
# ============================================================
found = tank_foundation_load(D=D, H=H_shell, liquid_density=liquid_density)
print(f"\n[步骤9] 基础荷载估算")
print(f"  函数: tank_foundation_load(D={D}, H={H_shell:.2f}, rho={liquid_density})")
print(f"  液体重量:  {found['liquid_weight']/1000:.1f} kN")
print(f"  罐壁重量:  {found['shell_weight']/1000:.1f} kN")
print(f"  罐底重量:  {found['bottom_weight']/1000:.1f} kN")
print(f"  总重量:    {found['total_weight']/1000:.1f} kN")
print(f"  基底压力:  {found['base_pressure']/1000:.1f} kPa (地基承载力一般要求 < 200 kPa)")

# ============================================================
# [步骤10] 设计校验
# ============================================================
design_dict = {
    "D": D, "H": H_total,
    "fill_ratio": fill_ratio,
    "base_pressure": found["base_pressure"],
}
validation = validate_tank_design(design_dict, tank_type="vertical_cylindrical")
print(f"\n[步骤10] 设计校验")
print(f"  校验通过率: {validation['pass_rate']*100:.0f}%")
for chk in validation["checks"]:
    status = "PASS" if chk["pass"] else "FAIL"
    print(f"    {chk['name']}: {chk['value']:.2f} (限值 {chk['limit']}) [{status}]")
for w in validation["warnings"]:
    print(f"  [!] {w}")

# ============================================================
# [结论] 设计方案汇总
# ============================================================
print(f"\n{'='*70}")
print(f"  设计方案结论")
print(f"{'='*70}")
print(f"  储存介质密度:         {liquid_density} kg/m³")
print(f"  有效容积:             {V_effective:.1f} m³")
print(f"  罐内径 D:             {D} m")
print(f"  筒体高度 H:           {H_shell:.2f} m")
print(f"  总高 (含封头):        {H_total:.2f} m")
print(f"  总高径比 H/D:         {HD_total:.2f}")
print(f"  顶封头:               锥顶, 深度 {hd_top:.3f} m")
print(f"  底封头:               平底")
print(f"  总几何容积:           {V_total_actual:.2f} m³")
print(f"  填充率:               {fill_ratio:.2%}")
print(f"  底层壁厚:             {thickness['t_bottom_shell_mm']:.1f} mm")
print(f"  顶层壁厚:             {thickness['t_top_shell_mm']:.1f} mm")
print(f"  总重量 (满液):        {found['total_weight']/1000:.1f} kN")
print(f"  基底压力:             {found['base_pressure']/1000:.1f} kPa")
print(f"{'='*70}")
