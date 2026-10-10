# -*- coding: utf-8 -*-
"""甲醇-乙酸酯化反应反应器设计示例 - 固定床(PFR)和全混流(CSTR)"""
import sys, io, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.reactor_device import _full_reactor_design

# ============================================================
# 反应参数
# ============================================================
# 甲醇(A) + 乙酸(B) → 乙酸甲酯 + 水
# 二级不可逆双分子反应: r = k * C_A * C_B
# 
# 进料条件:
# - 总进料流量 F_total = 0.05 mol/s
# - 甲醇入口浓度 C_A0 = 20 mol/m³
# - 乙酸入口浓度 C_B0 = 15 mol/m³
# - 反应速率常数 k = 0.1 m³/(mol·s)
# - 乙酸目标转化率 X_B = 90%
# - 温度 T = 200°C = 473.15 K
# - 压力 P = 1 atm = 101325 Pa
#
# 化学计量关系: A + B → 产物 (1:1)
# 乙酸(B)为限量组分，甲醇(A)过量
# 乙酸转化率 X_B = 0.90 对应甲醇转化率 X_A = (C_B0/C_A0) * X_B = (15/20) * 0.90 = 0.675

F_total = 0.05        # mol/s 混合进料
C_A0 = 20.0           # mol/m³ 甲醇
C_B0 = 15.0           # mol/m³ 乙酸
X_B_target = 0.90     # 乙酸目标转化率
X_A_target = (C_B0 / C_A0) * X_B_target  # 甲醇转化率 = 0.675
k = 0.1               # m³/(mol·s) 二级反应速率常数
T = 200               # °C (会自动转为K)
P = 101325            # Pa (常压)

# 计算 F_A0
# F_total = F_A + F_B + 其他惰性组分
# 假设只有甲醇和乙酸，则 F_A + F_B = F_total
# C_A0 / C_B0 = F_A0 / F_B0 (同体积流量)
# F_A0 = F_total * C_A0 / (C_A0 + C_B0) = 0.05 * 20/35 = 0.02857 mol/s
F_A0 = F_total * C_A0 / (C_A0 + C_B0)

print("=" * 70)
print("  甲醇-乙酸酯化反应 - 反应器设计")
print("=" * 70)
print(f"\n反应: 甲醇(A) + 乙酸(B) → 乙酸甲酯 + 水")
print(f"动力学: r = k × C_A × C_B (二级不可逆双分子反应)")
print(f"\n进料条件:")
print(f"  总进料流量: {F_total} mol/s")
print(f"  甲醇(A)入口浓度: {C_A0} mol/m³")
print(f"  乙酸(B)入口浓度: {C_B0} mol/m³")
print(f"  甲醇进料量 F_A0: {F_A0:.4f} mol/s")
print(f"\n反应条件:")
print(f"  速率常数 k: {k} m³/(mol·s)")
print(f"  温度: {T}°C = {T+273.15} K")
print(f"  压力: {P/1000:.2f} kPa")
print(f"\n转化率:")
print(f"  乙酸目标转化率 X_B: {X_B_target*100:.1f}%")
print(f"  甲醇对应转化率 X_A: {X_A_target*100:.1f}%")

# ============================================================
# 1. 固定床反应器 (PFR) 设计
# ============================================================
print("\n" + "=" * 70)
print("  1. 固定床反应器 (PFR) 设计")
print("=" * 70)

result_pfr = _full_reactor_design(
    reactor_type="fixed_bed",
    X_target=X_A_target,        # 甲醇转化率
    k=k,
    delta_H_rxn=0.0,            # 反应焓未知，设为0
    T_in=T,
    F_A0=F_A0,
    C_A0=C_A0,
    C_B0=C_B0,                  # 乙酸浓度，启用双分子动力学
    n=2,                        # 二级反应
    P_in=P,
    rho_gas=800.0,              # 液相密度估算 (kg/m³)
    mu=5e-4,                    # 液相粘度 (Pa·s)
    Cp_mix=150.0,               # 混合热容 J/(mol·K)
)

print("\n--- 固定床反应器结果 ---")
print(f"反应器体积 V: {result_pfr.get('V', 0)*1000:.2f} L = {result_pfr.get('V', 0):.4f} m³")
print(f"空时 τ: {result_pfr.get('tau', 0):.1f} s")
print(f"床层直径 D: {result_pfr.get('D_bed', 0)*1000:.1f} mm")
print(f"床层长度 L: {result_pfr.get('L_bed', 0)*1000:.1f} mm")
print(f"长径比 L/D: {result_pfr.get('H_D_ratio', 0):.2f}")

if 'nth_order_analysis' in result_pfr:
    nth = result_pfr['nth_order_analysis']
    print(f"\n解析解分析:")
    print(f"  反应级数: {nth.get('order_label', 'N/A')}")
    print(f"  计算体积: {nth.get('V', 0)*1000:.2f} L")

if 'flow_and_transport' in result_pfr:
    ft = result_pfr['flow_and_transport']
    print(f"\n流动参数:")
    print(f"  体积流量 v0: {ft.get('v0_m3_s', 0)*1e6:.2f} cm³/s")
    print(f"  空塔气速: {ft.get('u_superficial_m_s', 0):.4f} m/s")
    print(f"  颗粒Re数: {ft.get('Re_particle', 0):.2f}")

if 'validation' in result_pfr:
    val = result_pfr['validation']
    print(f"\n约束校核:")
    print(f"  通过率: {val.get('pass_rate', 0)*100:.0f}%")
    for check in val.get('checks', []):
        status = "✓" if check.get('pass') else "✗"
        print(f"  {status} {check.get('name')}: {check.get('value')} (限值: {check.get('limit')})")

# ============================================================
# 2. 全混流反应器 (CSTR) 设计
# ============================================================
print("\n" + "=" * 70)
print("  2. 全混流反应器 (CSTR) 设计")
print("=" * 70)

result_cstr = _full_reactor_design(
    reactor_type="cstr",
    X_target=X_A_target,        # 甲醇转化率
    k=k,
    delta_H_rxn=0.0,            # 反应焓未知，设为0
    T_in=T,
    F_A0=F_A0,
    C_A0=C_A0,
    C_B0=C_B0,                  # 乙酸浓度，启用双分子动力学
    n=2,                        # 二级反应
    P_in=P,
    rho_gas=800.0,              # 液相密度估算 (kg/m³)
    mu=5e-4,                    # 液相粘度 (Pa·s)
    Cp_mix=150.0,               # 混合热容 J/(mol·K)
    num_tanks=1,                # 单釜
)

print("\n--- 全混流反应器结果 ---")
print(f"反应器体积 V: {result_cstr.get('V', result_cstr.get('V_total', 0))*1000:.2f} L = {result_cstr.get('V', result_cstr.get('V_total', 0)):.4f} m³")
print(f"空时 τ: {result_cstr.get('tau', 0):.1f} s")
print(f"釜体直径 D: {result_cstr.get('D_tank', result_cstr.get('D_bed', 0))*1000:.1f} mm")
print(f"釜体高度 H: {result_cstr.get('H_tank', result_cstr.get('L_bed', 0))*1000:.1f} mm")
print(f"高径比 H/D: {result_cstr.get('H_D_ratio', 0):.2f}")

if 'nth_order_analysis' in result_cstr:
    nth = result_cstr['nth_order_analysis']
    print(f"\n解析解分析:")
    print(f"  反应级数: {nth.get('order_label', 'N/A')}")
    print(f"  计算体积: {nth.get('V', 0)*1000:.2f} L")

if 'flow_and_transport' in result_cstr:
    ft = result_cstr['flow_and_transport']
    print(f"\n流动参数:")
    print(f"  体积流量 v0: {ft.get('v0_m3_s', 0)*1e6:.2f} cm³/s")

if 'validation' in result_cstr:
    val = result_cstr['validation']
    print(f"\n约束校核:")
    print(f"  通过率: {val.get('pass_rate', 0)*100:.0f}%")
    for check in val.get('checks', []):
        status = "✓" if check.get('pass') else "✗"
        print(f"  {status} {check.get('name')}: {check.get('value')} (限值: {check.get('limit')})")

# ============================================================
# 3. 对比分析
# ============================================================
print("\n" + "=" * 70)
print("  3. PFR vs CSTR 对比")
print("=" * 70)

V_pfr = result_pfr.get('V', 0)
V_cstr = result_cstr.get('V', result_cstr.get('V_total', 0))
ratio = V_cstr / V_pfr if V_pfr > 0 else 0

print(f"\n{'参数':<20} {'固定床(PFR)':<15} {'全混流(CSTR)':<15} {'CSTR/PFR':<10}")
print("-" * 60)
print(f"{'反应器体积 (L)':<20} {V_pfr*1000:<15.2f} {V_cstr*1000:<15.2f} {ratio:<10.2f}")
print(f"{'空时 τ (s)':<20} {result_pfr.get('tau', 0):<15.1f} {result_cstr.get('tau', 0):<15.1f}")

print("\n结论:")
if ratio > 1:
    print(f"  • CSTR 体积是 PFR 的 {ratio:.2f} 倍")
    print(f"  • 对于二级反应，PFR 效率更高")
else:
    print(f"  • PFR 和 CSTR 体积相近")
