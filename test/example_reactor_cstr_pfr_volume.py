# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
"""
液相一级不可逆反应 A -> B 反应器体积计算
==========================================
已知:
  k = 0.2 min^-1
  进料流量 v0 = 100 L/min
  进口浓度 C_A0 = 2 mol/L
  目标转化率 X = 0.85

求解: CSTR 体积 和 PFR 体积
"""
sys.path.insert(0, r"c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents")

from physics_engine.reactor import cstr_design, pfr_design

print("=" * 70)
print("       液相一级不可逆反应 A -> B -- 反应器体积计算")
print("=" * 70)

# ============================================================
# 1. 已知参数
# ============================================================
k_min = 0.2           # min^-1
v0_L_min = 100.0      # L/min
C_A0_mol_L = 2.0      # mol/L
X_target = 0.85

# 转换为 SI 一致单位 (mol, m3, s)
k = k_min / 60.0                # s^-1
v0 = v0_L_min / 60.0 / 1000.0   # m3/s
C_A0 = C_A0_mol_L * 1000.0      # mol/m3
F_A0 = v0 * C_A0                 # mol/s

print(f"\n[已知条件]")
print(f"  反应: A -> B, 液相一级不可逆")
print(f"  k = {k_min} min^-1 = {k:.6f} s^-1")
print(f"  v0 = {v0_L_min} L/min = {v0:.6f} m3/s")
print(f"  C_A0 = {C_A0_mol_L} mol/L = {C_A0:.0f} mol/m3")
print(f"  F_A0 = v0 * C_A0 = {F_A0:.4f} mol/s")
print(f"  X_target = {X_target}")

# ============================================================
# 2. CSTR 体积计算
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤1] CSTR 体积计算")

cstr_result = cstr_design(
    F_A0=F_A0,
    X_target=X_target,
    k=k,
    C_A0=C_A0,
    n=1.0,
    num_tanks=1,
)

V_CSTR_m3 = cstr_result["V_total"]
V_CSTR_L = V_CSTR_m3 * 1000.0
tau_CSTR = cstr_result["tau"]

# 解析公式验证
V_CSTR_analytical = v0 * X_target / (k * (1 - X_target))
tau_analytical = X_target / (k * (1 - X_target))

print(f"\n  CSTR 设计方程:")
print(f"    V = F_A0 * X / r_A = F_A0 * X / (k * C_A0 * (1-X))")
print(f"    V = v0 * X / (k * (1-X))")
print(f"    V = {v0:.6f} * {X_target} / ({k:.6f} * {1-X_target})")
print(f"    V = {V_CSTR_analytical:.4f} m3 = {V_CSTR_analytical*1000:.1f} L")
print(f"\n  函数计算结果:")
print(f"    V_CSTR = {V_CSTR_m3:.4f} m3 = {V_CSTR_L:.1f} L")
print(f"    空时 tau = {tau_CSTR:.1f} s = {tau_CSTR/60:.2f} min")
print(f"    出口浓度 C_A = {cstr_result['C_A_out']:.1f} mol/m3 = {cstr_result['C_A_out']/1000:.4f} mol/L")

# ============================================================
# 3. PFR 体积计算
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤2] PFR 体积计算")

pfr_result = pfr_design(
    F_A0=F_A0,
    X_target=X_target,
    k=k,
    C_A0=C_A0,
    n=1.0,
    epsilon=0.0,  # 液相
)

V_PFR_m3 = pfr_result["V"]
V_PFR_L = V_PFR_m3 * 1000.0
tau_PFR = pfr_result["tau"]

# 解析公式验证
V_PFR_analytical = v0 / k * (-1) * __import__('math').log(1 - X_target)
tau_PFR_analytical = -__import__('math').log(1 - X_target) / k

print(f"\n  PFR 设计方程:")
print(f"    V = (v0/k) * (-ln(1-X))")
print(f"    V = ({v0:.6f}/{k:.6f}) * (-ln(1-{X_target}))")
print(f"    V = {V_PFR_analytical:.4f} m3 = {V_PFR_analytical*1000:.1f} L")
print(f"\n  函数计算结果:")
print(f"    V_PFR = {V_PFR_m3:.4f} m3 = {V_PFR_L:.1f} L")
print(f"    空时 tau = {tau_PFR:.1f} s = {tau_PFR/60:.2f} min")

# ============================================================
# 4. 对比分析
# ============================================================
print(f"\n{'='*70}")
print(f"[对比分析]")

ratio = V_CSTR_m3 / V_PFR_m3 if V_PFR_m3 > 0 else float('inf')

print(f"  {'项目':<20} {'CSTR':>12} {'PFR':>12} {'比值':>10}")
print(f"  {'-'*56}")
print(f"  {'反应器体积 (m3)':<20} {V_CSTR_m3:>12.4f} {V_PFR_m3:>12.4f} {ratio:>10.2f}")
print(f"  {'反应器体积 (L)':<20} {V_CSTR_L:>12.1f} {V_PFR_L:>12.1f}")
print(f"  {'空时 tau (min)':<20} {tau_CSTR/60:>12.2f} {tau_PFR/60:>12.2f}")
print(f"  {'出口浓度 (mol/L)':<20} {cstr_result['C_A_out']/1000:>12.4f} {C_A0_mol_L*(1-X_target):>12.4f}")

print(f"\n  V_CSTR / V_PFR = {ratio:.2f}")
print(f"  物理解释: 对一级反应, V_CSTR/V_PFR = X/(1-X) / (-ln(1-X))")
print(f"           = {X_target}/{1-X_target} / {-__import__('math').log(1-X_target):.4f} = {ratio:.2f}")

# ============================================================
# 5. 最终结论
# ============================================================
print(f"\n{'='*70}")
print(f"[设计结论]")
print(f"  反应: A -> B, k = {k_min} min^-1, 一级不可逆")
print(f"  进料: v0 = {v0_L_min} L/min, C_A0 = {C_A0_mol_L} mol/L")
print(f"  目标转化率: X = {X_target*100:.0f}%")
print(f"")
print(f"  CSTR 所需体积: V = {V_CSTR_L:.1f} L ({V_CSTR_m3:.4f} m3)")
print(f"  PFR 所需体积:  V = {V_PFR_L:.1f} L ({V_PFR_m3:.4f} m3)")
print(f"  体积比: V_CSTR/V_PFR = {ratio:.2f}")
print(f"")
print(f"  PFR 体积更小, 对一级反应效率更高 (无返混, 推动力大)")
print(f"{'='*70}")
