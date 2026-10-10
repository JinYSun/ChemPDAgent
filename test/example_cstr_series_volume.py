# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
"""
液相一级不可逆反应 A -> 产物
多釜串联 CSTR 体积计算 (N = 2, 3, 4, 5)
==========================================
已知:
  F_A0 = 2.0 mol/s
  C_A0 = 1000 mol/m3
  k = 0.021 s^-1 (一级)
  X = 0.90
  dH_rxn = -50000 J/mol (放热)
"""
sys.path.insert(0, r"c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents")

from physics_engine.reactor import cstr_design, reaction_enthalpy

print("=" * 70)
print("       液相一级反应 A -> 产物 -- 多釜串联 CSTR 设计")
print("=" * 70)

# ============================================================
# 1. 已知参数
# ============================================================
F_A0 = 2.0          # mol/s
C_A0 = 1000.0       # mol/m3
k = 0.021           # s^-1
X_target = 0.90
dH_rxn = -50000.0   # J/mol

v0 = F_A0 / C_A0    # m3/s

print(f"\n[已知条件]")
print(f"  F_A0 = {F_A0} mol/s")
print(f"  C_A0 = {C_A0} mol/m3")
print(f"  v0 = F_A0/C_A0 = {v0:.6f} m3/s")
print(f"  k = {k} s^-1 (一级)")
print(f"  X_target = {X_target}")
print(f"  dH_rxn = {dH_rxn} J/mol (放热)")

# ============================================================
# 2. 单釜基准 (N=1)
# ============================================================
print(f"\n{'='*70}")
print(f"[基准] 单釜 CSTR (N=1)")

r1 = cstr_design(F_A0=F_A0, X_target=X_target, k=k, C_A0=C_A0, n=1.0, num_tanks=1)
V_single = r1["V_total"]
tau_single = r1["tau"]

# 解析验证: tau = X / (k*(1-X))
tau_analytical = X_target / (k * (1 - X_target))
V_analytical = v0 * tau_analytical

print(f"  解析: tau = X/(k*(1-X)) = {X_target}/({k}*{1-X_target}) = {tau_analytical:.2f} s")
print(f"  解析: V = v0 * tau = {v0:.4f} * {tau_analytical:.2f} = {V_analytical:.2f} m3")
print(f"  函数: V = {V_single:.4f} m3, tau = {tau_single:.2f} s")

# ============================================================
# 3. 多釜串联计算 (N = 2, 3, 4, 5)
# ============================================================
print(f"\n{'='*70}")
print(f"[计算] 多釜串联 CSTR (N = 2, 3, 4, 5)")

results = []
for N in [2, 3, 4, 5]:
    r = cstr_design(F_A0=F_A0, X_target=X_target, k=k, C_A0=C_A0, n=1.0, num_tanks=N)
    V_total = r["V_total"]
    V_per = r["V_per_tank"]
    tau_total = r["tau"]
    results.append({"N": N, "V_total": V_total, "V_per": V_per, "tau": tau_total})
    print(f"\n  N={N}:")
    print(f"    V_total = {V_total:.4f} m3 = {V_total*1000:.2f} L")
    print(f"    V_per   = {V_per:.4f} m3 = {V_per*1000:.2f} L")
    print(f"    tau     = {tau_total:.2f} s = {tau_total/60:.2f} min")

# ============================================================
# 4. 反应热
# ============================================================
print(f"\n{'='*70}")
print(f"[反应热]")

Q = reaction_enthalpy(delta_H_rxn=dH_rxn, conversion=X_target, molar_flow=F_A0)
print(f"  Q = -dH_rxn * F_A0 * X = -({dH_rxn}) * {F_A0} * {X_target}")
print(f"  Q = {Q:.0f} W = {Q/1000:.2f} kW (需移除热量)")

# ============================================================
# 5. 汇总对比
# ============================================================
print(f"\n{'='*70}")
print(f"[汇总对比]")

print(f"\n  {'N釜':<8} {'V_total(m3)':<14} {'V_per(m3)':<12} {'tau(s)':<10} {'V/V_PFR':<10}")
print(f"  {'-'*56}")

# PFR 基准: V_PFR = (v0/k)*(-ln(1-X))
import math
V_PFR = (v0 / k) * (-math.log(1 - X_target))
print(f"  {'PFR':<8} {V_PFR:<14.4f} {'-':<12} {'-':<10} {'1.00':<10}")
print(f"  {'1(CSTR)':<8} {V_single:<14.4f} {V_single:<12.4f} {tau_single:<10.2f} {V_single/V_PFR:<10.2f}")
for r in results:
    ratio = r["V_total"] / V_PFR
    print(f"  {r['N']:<8} {r['V_total']:<14.4f} {r['V_per']:<12.4f} {r['tau']:<10.2f} {ratio:<10.2f}")

# ============================================================
# 6. 结论
# ============================================================
print(f"\n{'='*70}")
print(f"[设计结论]")
print(f"  反应: A -> 产物, 一级不可逆, k = {k} s^-1")
print(f"  进料: F_A0 = {F_A0} mol/s, C_A0 = {C_A0} mol/m3, v0 = {v0:.4f} m3/s")
print(f"  目标转化率: X = {X_target*100:.0f}%")
print(f"  反应热: Q = {Q/1000:.2f} kW (放热)")
print(f"")
print(f"  PFR 基准体积: {V_PFR:.4f} m3")
print(f"")
print(f"  {'N':<4} {'V_total(m3)':<14} {'V_per(m3)':<12} {'V/V_PFR':<10}")
print(f"  {'-'*42}")
print(f"  {'1':<4} {V_single:<14.4f} {V_single:<12.4f} {V_single/V_PFR:<10.2f}")
for r in results:
    print(f"  {r['N']:<4} {r['V_total']:<14.4f} {r['V_per']:<12.4f} {r['V_total']/V_PFR:<10.2f}")
print(f"\n  随釜数增加, 总容积趋近 PFR, 但单釜体积递减")
print(f"{'='*70}")
