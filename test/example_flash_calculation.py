# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
"""
等温闪蒸计算
=============
已知:
  进料 F = 100 kmol/h
  组成 z = [甲烷 0.5, 乙烷 0.2, 丙烷 0.3]
  闪蒸温度 T = -20 degC = 253.15 K
  闪蒸压力 P = 1.2 MPa = 1200000 Pa

求解:
  气相分率 psi = V/F
  液相组成 x_i
  气相组成 y_i
"""
import math
sys.path.insert(0, r"c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents")

from physics_engine.thermo_helper import (
    get_mixture_properties,
    get_mixture_k_values,
)

print("=" * 70)
print("          等温闪蒸计算 (Rachford-Rice)")
print("=" * 70)

# ============================================================
# 1. 已知条件
# ============================================================
components = ['methane', 'ethane', 'propane']
z = [0.5, 0.2, 0.3]          # 进料摩尔分数
F = 100.0                     # 进料流量 kmol/h
T_C = -20.0                   # 闪蒸温度 degC
T = T_C + 273.15              # 闪蒸温度 K
P_MPa = 1.2                   # 闪蒸压力 MPa
P = P_MPa * 1e6               # 闪蒸压力 Pa

print(f"\n[已知条件]")
print(f"  进料流量 F   = {F} kmol/h")
print(f"  组分         = 甲烷/乙烷/丙烷")
print(f"  进料组成 z   = {z}")
print(f"  闪蒸温度 T   = {T_C} degC = {T:.2f} K")
print(f"  闪蒸压力 P   = {P_MPa} MPa = {P:.0f} Pa")

# ============================================================
# 2. 获取混合物物性与 K 值
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤1] 获取混合物物性与 K 值")

mix_props = get_mixture_properties(components, z, T, P)
print(f"  混合物平均分子量 MW_mix = {mix_props['MW_mix']:.4f} g/mol")
print(f"  K值计算方法: {mix_props['K_method']}")

Ks = get_mixture_k_values(components, T, P, z, method="auto")
print(f"  K值 = {Ks}")
for i, c in enumerate(components):
    print(f"    K_{c} = {Ks[i]:.6f}")

# 检查是否处于两相区
f_0 = sum(z[i] * (Ks[i] - 1) / (1 + 0 * (Ks[i] - 1)) for i in range(len(components)))
f_1 = sum(z[i] * (Ks[i] - 1) / (1 + 1.0 * (Ks[i] - 1)) for i in range(len(components)))
print(f"\n  两相区检验:")
print(f"    f(psi=0) = {f_0:.6f}  ({'>0 有液相' if f_0 > 0 else '<=0 全气相'})")
print(f"    f(psi=1) = {f_1:.6f}  ({'<0 有气相' if f_1 < 0 else '>=0 全液相'})")

if f_0 <= 0:
    print(f"  ==> 全气相 (所有 Ki > 1 或 f(0)<=0), 闪蒸无液相")
elif f_1 >= 0:
    print(f"  ==> 全液相 (所有 Ki < 1 或 f(1)>=0), 闪蒸无气相")
else:
    print(f"  ==> 处于两相区, 可进行闪蒸计算")

# ============================================================
# 3. Rachford-Rice 迭代求解
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤2] Rachford-Rice 迭代求解气相分率 psi = V/F")
print(f"  方程: f(psi) = Sum[ z_i*(K_i-1) / (1 + psi*(K_i-1)) ] = 0")
print(f"  求解区间: psi in [0, 1]")

def rachford_rice(psi, z, K):
    """Rachford-Rice 函数"""
    return sum(z[i] * (K[i] - 1) / (1 + psi * (K[i] - 1)) for i in range(len(z)))

def rachford_rice_deriv(psi, z, K):
    """Rachford-Rice 导数 (Newton-Raphson)"""
    return sum(-z[i] * (K[i] - 1)**2 / (1 + psi * (K[i] - 1))**2 for i in range(len(z)))

# Newton-Raphson 迭代
psi = 0.5  # 初始猜测
max_iter = 100
tol = 1e-10

print(f"\n  初始猜测 psi = {psi}")
print(f"  {'迭代':>4} | {'psi':>12} | {'f(psi)':>14} | {'df/dpsi':>14}")
print(f"  {'-'*4}-+-{'-'*12}-+-{'-'*14}-+-{'-'*14}")

for it in range(max_iter):
    f_val = rachford_rice(psi, z, Ks)
    df_val = rachford_rice_deriv(psi, z, Ks)

    if it < 5 or it % 10 == 0 or abs(f_val) < tol:
        print(f"  {it:>4} | {psi:>12.8f} | {f_val:>14.6e} | {df_val:>14.6e}")

    if abs(f_val) < tol:
        print(f"  ==> 收敛! 迭代 {it+1} 次")
        break

    # Newton 步
    delta = -f_val / df_val
    # 限制步长, 保持在 [0, 1] 区间
    psi_new = psi + delta
    if psi_new < 0:
        psi_new = psi / 2
    elif psi_new > 1:
        psi_new = (psi + 1) / 2
    psi = psi_new

psi_final = psi
print(f"\n  气相分率 psi = V/F = {psi_final:.6f}")

# ============================================================
# 4. 计算气液相组成
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤3] 计算气液相组成")
print(f"  x_i = z_i / [1 + psi*(K_i - 1)]")
print(f"  y_i = K_i * x_i")

x = []
y = []
for i in range(len(components)):
    xi = z[i] / (1 + psi_final * (Ks[i] - 1))
    yi = Ks[i] * xi
    x.append(xi)
    y.append(yi)

# 归一化验证
sum_x = sum(x)
sum_y = sum(y)
x_norm = [xi / sum_x for xi in x]
y_norm = [yi / sum_y for yi in y]

print(f"\n  {'组分':>8} | {'z_i':>8} | {'K_i':>10} | {'x_i(液相)':>10} | {'y_i(气相)':>10} | {'回收率':>8}")
print(f"  {'-'*8}-+-{'-'*8}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*8}")
for i, c in enumerate(components):
    recovery = psi_final * y_norm[i] / z[i] * 100  # 气相回收率
    print(f"  {c:>8} | {z[i]:>8.4f} | {Ks[i]:>10.6f} | {x_norm[i]:>10.6f} | {y_norm[i]:>10.6f} | {recovery:>7.2f}%")

print(f"  {'Sum':>8} | {sum(z):>8.4f} | {'':>10} | {sum(x_norm):>10.6f} | {sum(y_norm):>10.6f} |")

# ============================================================
# 5. 物料平衡
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤4] 物料平衡")

V = psi_final * F          # 气相流量 kmol/h
L = (1 - psi_final) * F    # 液相流量 kmol/h

print(f"  进料 F = {F} kmol/h")
print(f"  气相 V = psi * F = {psi_final:.6f} * {F} = {V:.4f} kmol/h")
print(f"  液相 L = (1-psi) * F = {1-psi_final:.6f} * {F} = {L:.4f} kmol/h")
print(f"  校验: V + L = {V + L:.4f} kmol/h (应等于 {F})")

print(f"\n  气相各组分流量 (kmol/h):")
for i, c in enumerate(components):
    print(f"    {c}: V*y_i = {V:.4f} * {y_norm[i]:.6f} = {V * y_norm[i]:.4f}")

print(f"\n  液相各组分流量 (kmol/h):")
for i, c in enumerate(components):
    print(f"    {c}: L*x_i = {L:.4f} * {x_norm[i]:.6f} = {L * x_norm[i]:.4f}")

# ============================================================
# 6. 最终结论
# ============================================================
print(f"\n{'='*70}")
print(f"[计算结论]")
print(f"  闪蒸条件: T = {T_C} degC, P = {P_MPa} MPa")
print(f"  K值方法: {mix_props['K_method']}")
print(f"  气相分率 psi = V/F = {psi_final:.6f} ({psi_final*100:.2f}%)")
print(f"  气相流量 V = {V:.4f} kmol/h")
print(f"  液相流量 L = {L:.4f} kmol/h")
print(f"")
print(f"  液相组成 x_i:")
for i, c in enumerate(components):
    print(f"    {c}: x = {x_norm[i]:.6f}")
print(f"  气相组成 y_i:")
for i, c in enumerate(components):
    print(f"    {c}: y = {y_norm[i]:.6f}")
print(f"{'='*70}")
