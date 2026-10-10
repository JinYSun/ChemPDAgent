# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
"""
CSTR 双分子反应设计 -- 甲醇 + 乙酸 (酯化反应)
================================================
已知:
  混合进料 F_total = 0.05 mol/s
  甲醇入口浓度 C_A0 = 20 mol/m3
  乙酸入口浓度 C_B0 = 15 mol/m3
  二级不可逆反应 r = k * C_A * C_B
  k = 0.1 m3/(mol*s)
  乙酸目标转化率 X_B = 0.90
  T = 200 degC, P = 1 atm

注意: 需先验证化学计量约束是否允许 X_B = 0.90
"""
sys.path.insert(0, r"c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents")

from physics_engine.reactor import cstr_design

print("=" * 70)
print("       CSTR 双分子反应设计 -- 甲醇 + 乙酸")
print("=" * 70)

# ============================================================
# 1. 已知参数
# ============================================================
F_total = 0.05        # 混合进料总摩尔流量 mol/s
C_A0 = 20.0           # 甲醇入口浓度 mol/m3
C_B0 = 15.0           # 乙酸入口浓度 mol/m3
k = 0.1               # 反应常数 m3/(mol*s)
X_target_B = 0.90     # 乙酸目标转化率
T = 200               # degC
P_atm = 1             # atm

print(f"\n[已知条件]")
print(f"  反应: CH3OH + CH3COOH -> CH3COOCH3 + H2O")
print(f"  二级不可逆, r = k * C_A * C_B")
print(f"  F_total = {F_total} mol/s (混合进料)")
print(f"  C_A0 (甲醇) = {C_A0} mol/m3")
print(f"  C_B0 (乙酸) = {C_B0} mol/m3")
print(f"  k = {k} m3/(mol*s)")
print(f"  X_target (乙酸) = {X_target_B}")
print(f"  T = {T} degC, P = {P_atm} atm")

# ============================================================
# 2. 物料平衡与化学计量约束
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤1] 物料平衡与化学计量约束")

# 总体积流量
C_total = C_A0 + C_B0  # mol/m3
v0 = F_total / C_total  # m3/s

# 各组分摩尔流量
F_A0 = v0 * C_A0  # 甲醇 mol/s
F_B0 = v0 * C_B0  # 乙酸 mol/s

print(f"  总浓度 C_total = {C_A0} + {C_B0} = {C_total} mol/m3")
print(f"  体积流量 v0 = F_total / C_total = {F_total} / {C_total} = {v0:.6f} m3/s = {v0*1e6:.1f} mL/s")
print(f"  F_A0 (甲醇) = v0 * C_A0 = {v0:.6f} * {C_A0} = {F_A0:.6f} mol/s")
print(f"  F_B0 (乙酸) = v0 * C_B0 = {v0:.6f} * {C_B0} = {F_B0:.6f} mol/s")

# 化学计量约束 (1:1 反应)
M_ratio = C_B0 / C_A0  # B/A 摩尔比
X_max_A = M_ratio       # A(甲醇)的最大转化率 (受限于B不足)
X_max_B = C_A0 / C_B0   # B(乙酸)的最大转化率... 不对

# 正确分析:
# 反应 A + B -> P, 1:1
# 消耗量: C_A0 * X_A = C_B0 * X_B
# 当 X_B = 1 时, 需要的 A 消耗 = C_B0, 但 A 最多提供 C_A0
# 所以 X_B 的最大值受限于: C_A0 * X_A_max = C_B0 => X_A_max = C_B0/C_A0
# 但 X_A 最大只能是 1, 所以如果 C_B0 < C_A0, 则 B 是限量组分
# B 全部消耗时: C_A0 * X_A = C_B0, X_A = C_B0/C_A0 = 0.75
# 此时 X_B = 1.0 (乙酸100%转化)
# 等等... 这里需要重新理解转化率的定义

# 在 cstr_design 中, X 是相对于 A 的转化率
# X = (F_A0 - F_A) / F_A0
# C_A_out = C_A0 * (1 - X)
# C_B_out = C_B0 - C_A0 * X  (因为 1:1 消耗)
# C_B_out >= 0 要求 X <= C_B0/C_A0 = 0.75

# 所以 A(甲醇)的转化率上限 X_A_max = C_B0/C_A0 = 0.75
# 当 X_A = 0.75 时, B(乙酸)全部消耗, X_B = 1.0

# 用户要求的是乙酸转化率 X_B = 0.90
# 乙酸消耗 = C_B0 * X_B = 15 * 0.90 = 13.5 mol/m3
# 甲醇消耗 = C_A0 * X_A = 13.5 => X_A = 13.5/20 = 0.675
# 甲醇剩余 = 20*(1-0.675) = 6.5 mol/m3 > 0, OK
# 乙酸剩余 = 15*(1-0.90) = 1.5 mol/m3 > 0, OK

# 所以实际上, 如果以 A(甲醇)为基准, X_A = 0.675 < X_max = 0.75, 是可行的!

print(f"\n  化学计量分析 (A + B -> P, 1:1):")
print(f"  C_B0/C_A0 = {C_B0}/{C_A0} = {M_ratio:.4f}")
print(f"  A(甲醇)转化率上限 = C_B0/C_A0 = {M_ratio:.4f} (B耗尽时)")
print(f"  B(乙酸)转化率上限 = 1.0 (A过量时)")

# 关键: 用户要求 X_B = 0.90, 换算为 X_A
# C_A0 * X_A = C_B0 * X_B
# X_A = C_B0 * X_B / C_A0 = 15 * 0.90 / 20 = 0.675
X_A_equiv = C_B0 * X_target_B / C_A0
print(f"\n  用户要求: X_B (乙酸) = {X_target_B}")
print(f"  等价于 X_A (甲醇) = C_B0 * X_B / C_A0 = {C_B0} * {X_target_B} / {C_A0} = {X_A_equiv:.4f}")
print(f"  X_A = {X_A_equiv:.4f} < X_A_max = {M_ratio:.4f} => 化学计量可行!")

# 验证出口浓度
CA_out = C_A0 * (1 - X_A_equiv)
CB_out = C_B0 * (1 - X_target_B)
print(f"\n  出口浓度验证:")
print(f"    C_A_out = {C_A0} * (1 - {X_A_equiv:.4f}) = {CA_out:.2f} mol/m3")
print(f"    C_B_out = {C_B0} * (1 - {X_target_B}) = {CB_out:.2f} mol/m3")
print(f"    两者均 > 0, 可行")

# ============================================================
# 3. CSTR 设计计算 (以甲醇 A 为基准)
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤2] CSTR 设计计算 (以甲醇 A 为基准)")

# 注意: cstr_design 中 X 是相对于 A 的转化率
# 我们需要传入 X_A_equiv
cstr_result = cstr_design(
    F_A0=F_A0,
    X_target=X_A_equiv,
    k=k,
    C_A0=C_A0,
    n=1.0,
    num_tanks=1,
    C_B0=C_B0,  # 关键! 启用双分子动力学
)

print(f"\n  调用 cstr_design(F_A0={F_A0:.6f}, X={X_A_equiv:.4f}, k={k}, C_A0={C_A0}, C_B0={C_B0})")
print(f"  动力学模式: {cstr_result.get('kinetics', 'N/A')}")

if 'error' in cstr_result:
    print(f"  [ERROR] {cstr_result['error']}")
else:
    V_CSTR = cstr_result["V_total"]
    tau_CSTR = cstr_result["tau"]
    CA_out_calc = cstr_result.get("C_A_out", 0)
    CB_out_calc = cstr_result.get("C_B_out", 0)

    # 手动验证
    r_out = k * CA_out_calc * CB_out_calc
    V_verify = F_A0 * X_A_equiv / r_out if r_out > 0 else float('inf')

    print(f"\n  出口浓度:")
    print(f"    C_A_out (甲醇) = {CA_out_calc:.4f} mol/m3")
    print(f"    C_B_out (乙酸) = {CB_out_calc:.4f} mol/m3")
    print(f"    反应速率 r_out = k * CA * CB = {k} * {CA_out_calc:.4f} * {CB_out_calc:.4f} = {r_out:.4f} mol/(m3*s)")

    print(f"\n  CSTR 设计方程:")
    print(f"    V = F_A0 * X_A / (k * C_A_out * C_B_out)")
    print(f"    V = {F_A0:.6f} * {X_A_equiv:.4f} / ({k} * {CA_out_calc:.4f} * {CB_out_calc:.4f})")
    print(f"    V = {F_A0 * X_A_equiv:.6f} / {r_out:.4f}")
    print(f"    V = {V_verify:.4f} m3 = {V_verify*1000:.2f} L")

    print(f"\n  计算结果:")
    print(f"    V_CSTR = {V_CSTR:.4f} m3 = {V_CSTR*1000:.2f} L")
    print(f"    空时 tau = {tau_CSTR:.2f} s = {tau_CSTR/60:.2f} min")

    # 转化率换算
    print(f"\n  转化率换算:")
    print(f"    甲醇转化率 X_A = {X_A_equiv:.4f} ({X_A_equiv*100:.1f}%)")
    print(f"    乙酸转化率 X_B = {X_target_B:.4f} ({X_target_B*100:.1f}%)")

# ============================================================
# 4. 对比: 如果直接以 X_A = 0.90 传入 (错误做法)
# ============================================================
print(f"\n{'='*70}")
print(f"[步骤3] 对比: 若误将 X=0.90 作为甲醇转化率传入")

cstr_wrong = cstr_design(
    F_A0=F_A0,
    X_target=0.90,
    k=k,
    C_A0=C_A0,
    n=1.0,
    num_tanks=1,
    C_B0=C_B0,
)
X_A_capped = cstr_wrong.get("kinetics", "")
print(f"  传入 X_target = 0.90")
print(f"  X_max = C_B0/C_A0 = {M_ratio:.4f}")
print(f"  函数自动将 X 限制为 X_max * 0.999 = {M_ratio*0.999:.4f}")
if 'error' in cstr_wrong:
    print(f"  [ERROR] {cstr_wrong['error']}")
else:
    print(f"  实际计算 X_A = {M_ratio*0.999:.4f} (乙酸几乎耗尽)")
    print(f"  V_CSTR = {cstr_wrong['V_total']:.4f} m3 = {cstr_wrong['V_total']*1000:.2f} L")
    print(f"  C_B_out = {cstr_wrong.get('C_B_out', 0):.4f} mol/m3 (接近零)")

# ============================================================
# 5. 最终结论
# ============================================================
print(f"\n{'='*70}")
print(f"[设计结论]")
print(f"  反应: CH3OH + CH3COOH -> CH3COOCH3 + H2O")
print(f"  二级不可逆, r = k * C_A * C_B, k = {k} m3/(mol*s)")
print(f"  进料: F_total = {F_total} mol/s, C_A0 = {C_A0}, C_B0 = {C_B0} mol/m3")
print(f"  T = {T} degC, P = {P_atm} atm")
print(f"")
print(f"  化学计量约束:")
print(f"    甲醇过量 (C_A0/C_B0 = {C_A0/C_B0:.2f})")
print(f"    乙酸为限量组分")
print(f"    X_A_max = {M_ratio:.4f} (甲醇转化率上限)")
print(f"")

if 'error' not in cstr_result:
    print(f"  目标: 乙酸转化率 X_B = {X_target_B*100:.0f}%")
    print(f"  等价甲醇转化率 X_A = {X_A_equiv:.4f} ({X_A_equiv*100:.1f}%)")
    print(f"  CSTR 体积 V = {V_CSTR:.4f} m3 = {V_CSTR*1000:.2f} L")
    print(f"  空时 tau = {tau_CSTR:.2f} s = {tau_CSTR/60:.2f} min")
    print(f"  出口: C_A = {CA_out_calc:.2f} mol/m3, C_B = {CB_out_calc:.2f} mol/m3")
else:
    print(f"  [ERROR] {cstr_result.get('error', '未知错误')}")

print(f"{'='*70}")
