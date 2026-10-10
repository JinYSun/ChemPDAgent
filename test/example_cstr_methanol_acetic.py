"""CSTR 全混流反应器设计 — 甲醇 + 乙酸 二级不可逆双分子反应
=========================================================

反应: A(甲醇) + B(乙酸) -> 产物
速率: r = k * C_A * C_B  (二级不可逆双分子)

设计条件:
  混合进料 F_total = 0.05 mol/s
  C_A0 = 20 mol/m3 (甲醇)
  C_B0 = 15 mol/m3 (乙酸)
  k = 0.1 m3/(mol.s)
  乙酸目标转化率 X_B = 90%
  T = 200 degC, P = 1 atm
"""

import sys
import os
import io
import math

if sys.stdout.encoding and sys.stdout.encoding.lower() in ('gbk', 'cp936'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

_project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from physics_engine.reactor import (
    cstr_design,
    pfr_design,
    reaction_rate,
    adiabatic_temperature_rise,
    design_reactor,
)


def sep(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def cstr_bimolecular(F_A0, C_A0, C_B0, k, X, v0=None):
    """双分子二级反应 CSTR 精确设计 (A + B -> P, r = k*CA*CB)

    CSTR 设计方程: V = F_A0 * X / (k * C_A * C_B)
    其中:
      C_A = C_A0 * (1 - X)
      C_B = C_B0 - C_A0 * X   (1:1 化学计量)

    X 为甲醇(A)转化率。若需以乙酸(B)转化率为基准，
    则 X_A = X_B * C_B0 / C_A0
    """
    if v0 is None:
        v0 = 1.0  # 归一化

    C_A = C_A0 * (1 - X)
    C_B = C_B0 - C_A0 * X

    if C_A <= 0 or C_B <= 0:
        return {"error": "出口浓度为零或负, 转化率超出极限"}

    r = k * C_A * C_B
    V = F_A0 * X / r
    tau = V / v0

    return {
        "V": V,
        "tau": tau,
        "C_A_out": C_A,
        "C_B_out": C_B,
        "r_out": r,
    }


def pfr_bimolecular(F_A0, C_A0, C_B0, k, X, v0=None):
    """双分子二级反应 PFR 精确设计 (解析解)

    PFR 设计方程 (A+B->P, r=k*CA*CB, CA0!=CB0):
      tau = 1/(k*(CA0-CB0)) * ln(CA_out*CB0 / (CA0*CB_out))
    其中:
      CA_out = CA0*(1-X)
      CB_out = CB0 - CA0*X
    """
    if v0 is None:
        v0 = 1.0

    CA_out = C_A0 * (1 - X)
    CB_out = C_B0 - C_A0 * X

    if CA_out <= 0 or CB_out <= 0:
        return {"V": float('inf'), "tau": float('inf')}

    if abs(C_B0 - C_A0) < 1e-10:
        # 等浓度特殊情况: tau = X / (k*CA0*(1-X))
        tau = X / (k * C_A0 * (1 - X))
    else:
        # 正确解析解: tau = 1/(k*(CA0-CB0)) * ln(CA_out*CB0 / (CA0*CB_out))
        tau = 1.0 / (k * (C_A0 - C_B0)) * math.log(CA_out * C_B0 / (C_A0 * CB_out))

    V = tau * v0
    return {"V": V, "tau": tau}


def main():
    # ================================================================
    # 1. 设计条件
    # ================================================================
    sep("1. 设计条件")

    F_total = 0.05          # 混合进料总摩尔流量 (mol/s)
    C_A0 = 20.0             # 甲醇入口浓度 (mol/m3)
    C_B0 = 15.0             # 乙酸入口浓度 (mol/m3)
    k = 0.1                 # 反应速率常数 m3/(mol.s)
    n = 2.0                 # 总反应级数
    X_B_target = 0.90       # 乙酸目标转化率 90%
    T_C = 200.0
    T_K = T_C + 273.15
    P_atm = 101325.0

    # 由总进料和浓度计算体积流量和各组分摩尔流量
    # F_total = F_A + F_B (混合进料中仅含甲醇和乙酸)
    # C_A0/C_B0 = F_A/F_B (浓度比=摩尔流量比)
    F_A0 = F_total * C_A0 / (C_A0 + C_B0)   # 甲醇摩尔流量 (mol/s)
    F_B0 = F_total * C_B0 / (C_A0 + C_B0)   # 乙酸摩尔流量 (mol/s)
    v0 = F_total / (C_A0 + C_B0)             # 体积流量 (m3/s)

    # 将乙酸转化率转换为甲醇转化率 (1:1 化学计量)
    # X_B = F_B0_reacted / F_B0 = F_A0_reacted / F_A0 * (F_A0/F_B0)
    # X_A = X_B * C_B0 / C_A0
    X_A_target = X_B_target * C_B0 / C_A0

    print(f"  反应: 甲醇(A) + 乙酸(B) -> 产物")
    print(f"  速率方程: r = k * C_A * C_B")
    print(f"  混合进料总流量 F_total = {F_total} mol/s")
    print(f"  甲醇摩尔流量 F_A0 = {F_A0:.5f} mol/s")
    print(f"  乙酸摩尔流量 F_B0 = {F_B0:.5f} mol/s")
    print(f"  体积流量 v0 = {v0:.6f} m3/s")
    print(f"  甲醇入口浓度 C_A0 = {C_A0} mol/m3")
    print(f"  乙酸入口浓度 C_B0 = {C_B0} mol/m3")
    print(f"  反应速率常数 k    = {k} m3/(mol.s)")
    print(f"  反应级数     n    = {n} (双分子二级)")
    print(f"  乙酸目标转化率 X_B = {X_B_target*100:.0f}%")
    print(f"  对应甲醇转化率 X_A = {X_A_target*100:.2f}%")
    print(f"  温度         T    = {T_C} degC ({T_K:.2f} K)")
    print(f"  压力         P    = 1 atm ({P_atm} Pa)")

    # ================================================================
    # 2. 反应动力学分析
    # ================================================================
    sep("2. 反应动力学分析")

    # 限制组分分析
    # 甲醇(A)过量，乙酸(B)为限制组分
    # 甲醇最大理论转化率 X_A_max = C_B0/C_A0 = 0.75
    # 对应乙酸最大转化率 X_B_max = 1.0 (乙酸可完全转化)
    X_A_max = C_B0 / C_A0  # 甲醇最大转化率 = 75%
    X_B_max = 1.0           # 乙酸理论上可100%转化(甲醇过量)

    X_A_actual = X_A_target
    if X_A_target >= X_A_max:
        X_A_actual = X_A_max * 0.95
        print(f"  [!] 甲醇转化率 {X_A_target*100:.1f}% >= 极限 {X_A_max*100:.1f}%, 不可达!")
        print(f"  [!] 自动调整为 X_A = {X_A_actual*100:.1f}%")

    # 对应的乙酸转化率
    X_B_actual = X_A_actual * C_A0 / C_B0

    r_in = k * C_A0 * C_B0
    C_A_out = C_A0 * (1 - X_A_actual)
    C_B_out = C_B0 - C_A0 * X_A_actual
    r_out = k * C_A_out * C_B_out

    print(f"  入口反应速率 r_in  = k*CA0*CB0 = {k}*{C_A0}*{C_B0} = {r_in:.2f} mol/(m3.s)")
    print(f"  出口甲醇浓度 C_A   = {C_A0}*(1-{X_A_actual:.4f}) = {C_A_out:.2f} mol/m3")
    print(f"  出口乙酸浓度 C_B   = {C_B0} - {C_A0}*{X_A_actual:.4f} = {C_B_out:.2f} mol/m3")
    print(f"  出口反应速率 r_out = {k}*{C_A_out:.2f}*{C_B_out:.2f} = {r_out:.4f} mol/(m3.s)")
    print(f"  速率比 r_in/r_out  = {r_in/r_out:.1f}")
    print(f"\n  限制组分: B(乙酸)")
    print(f"  甲醇最大理论转化率: {X_A_max*100:.1f}%")
    print(f"  乙酸最大理论转化率: {X_B_max*100:.1f}% (甲醇过量)")
    print(f"\n  实际计算: 乙酸转化率 X_B = {X_B_actual*100:.2f}%")
    print(f"  对应甲醇转化率 X_A = {X_A_actual*100:.2f}%")

    # ================================================================
    # 3. CSTR 精确设计 (双分子动力学)
    # ================================================================
    sep("3. CSTR 精确设计 (双分子 r=k*CA*CB)")

    # 以实际进料流量计算
    res = cstr_bimolecular(F_A0, C_A0, C_B0, k, X_A_actual, v0)

    V_cstr = res["V"]
    tau_cstr = res["tau"]

    print(f"\n  CSTR 设计方程: V = F_A0*X_A / (k * C_A_out * C_B_out)")
    print(f"  代入: V = {F_A0:.5f}*{X_A_actual:.4f} / ({k}*{C_A_out:.2f}*{C_B_out:.2f})")
    print(f"\n  反应器体积 V   = {V_cstr:.6f} m3 = {V_cstr*1000:.4f} L")
    print(f"  空间时间 tau   = {tau_cstr:.4f} s")
    print(f"  体积流量 v0    = {v0:.6f} m3/s")

    # 校验: F_B_reacted = F_B0 * X_B
    F_B_reacted = F_B0 * X_B_actual
    print(f"\n  乙酸反应量 = F_B0 * X_B = {F_B0:.5f} * {X_B_actual:.4f} = {F_B_reacted:.5f} mol/s")
    print(f"  甲醇反应量 = F_A0 * X_A = {F_A0:.5f} * {X_A_actual:.4f} = {F_A0*X_A_actual:.5f} mol/s")
    print(f"  (两者应相等, 1:1 化学计量验证)")

    # ================================================================
    # 4. PFR 精确设计 (解析解对比)
    # ================================================================
    sep("4. PFR 精确设计 (解析解对比)")

    pfr_res = pfr_bimolecular(F_A0, C_A0, C_B0, k, X_A_actual, v0)
    tau_pfr = pfr_res["tau"]
    V_pfr = pfr_res["V"]

    print(f"\n  PFR 解析解 (双分子 A+B, CA0!=CB0):")
    M = C_B0 / C_A0
    print(f"  M = CB0/CA0 = {M:.3f}")
    print(f"  tau = {tau_pfr:.4f} s")
    print(f"  V_PFR = {V_pfr:.6f} m3 = {V_pfr*1000:.4f} L")

    print(f"\n  对比:")
    print(f"    CSTR tau = {tau_cstr:.4f} s,  V_CSTR = {V_cstr*1000:.4f} L")
    print(f"    PFR  tau = {tau_pfr:.4f} s,  V_PFR  = {V_pfr*1000:.4f} L")
    if tau_pfr > 0:
        print(f"    V_CSTR/V_PFR = {tau_cstr/tau_pfr:.2f}")

    # ================================================================
    # 5. 不同乙酸转化率下的 CSTR 与 PFR 对比
    # ================================================================
    sep("5. 不同乙酸转化率下的 CSTR vs PFR 空间时间")

    # 乙酸转化率范围: 0 ~ 100% (对应甲醇转化率 0 ~ X_A_max=75%)
    X_B_list = [0.10, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99]
    print(f"\n  乙酸最大转化率 = 100% (甲醇过量)")
    print(f"  对应甲醇最大转化率 = {X_A_max*100:.1f}%\n")
    print(f"  {'X_B':>6s} | {'X_A':>6s} | {'tau_CSTR(s)':>12s} | {'tau_PFR(s)':>11s} | {'CSTR/PFR':>9s} | {'C_A_out':>8s} | {'C_B_out':>8s}")
    print(f"  {'-'*6}-+-{'-'*6}-+-{'-'*12}-+-{'-'*11}-+-{'-'*9}-+-{'-'*8}-+-{'-'*8}")

    for X_B in X_B_list:
        X_A = X_B * C_B0 / C_A0
        if X_A >= X_A_max:
            continue
        cstr_r = cstr_bimolecular(F_A0, C_A0, C_B0, k, X_A, v0)
        pfr_r = pfr_bimolecular(F_A0, C_A0, C_B0, k, X_A, v0)
        if "error" in cstr_r:
            print(f"  {X_B:5.2f} | {X_A:5.3f} | {'N/A':>12s} | {pfr_r['tau']:11.4f} | {'N/A':>9s} | {'N/A':>8s} | {'N/A':>8s}")
            continue
        ratio = cstr_r["tau"] / pfr_r["tau"] if pfr_r["tau"] > 0 else float('inf')
        print(f"  {X_B:6.3f} | {X_A:5.3f} | {cstr_r['tau']:12.4f} | {pfr_r['tau']:11.4f} | {ratio:9.2f} | "
              f"{cstr_r['C_A_out']:8.2f} | {cstr_r['C_B_out']:8.2f}")

    # ================================================================
    # 6. 使用 design_reactor 完整设计
    # ================================================================
    sep("6. design_reactor 完整 CSTR 设计")

    full_result = design_reactor(
        reactor_type="cstr",
        F_A0=F_A0,
        X_target=X_A_target,
        delta_H_rxn=-50000.0,
        k=k,
        C_A0=C_A0,
        n=n,
        T_in=T_K,
        P_in=P_atm,
        num_tanks=1,
        rho_cat=2000.0,
        particle_diameter=0.005,
        epsilon=0.4,
        d_tube_inner=0.038,
        L_tube=6.0,
        U=350.0,
        T_coolant=T_K,
        rho_gas=1.2,
        mu=1e-3,
        Cp_mix=150.0,
        C_B0=C_B0,
    )

    V_full = full_result.get('V_total', 0)
    tau_full = full_result.get('tau', 0)
    D_tank = full_result.get('D_tank', 0)
    H_tank = full_result.get('H_tank', 0)

    print(f"\n  进料: F_A0 = {F_A0:.5f} mol/s, v0 = {v0:.6f} m3/s")
    print(f"  反应器类型   : {full_result.get('reactor_type')}")
    print(f"  总体积       : {V_full:.6f} m3  ({V_full*1000:.4f} L)")
    print(f"  停留时间     : {tau_full:.2f} s")
    print(f"  釜体直径     : {D_tank:.4f} m  ({D_tank*1000:.1f} mm)")
    print(f"  釜体高度     : {H_tank:.4f} m  ({H_tank*1000:.1f} mm)")
    print(f"  出口温度     : {full_result.get('T_out', 0):.2f} K")

    # 校验
    val = full_result.get("validation", {})
    print(f"\n  设计校验:")
    for chk in val.get("checks", []):
        tag = "PASS" if chk["status"] == "pass" else ("WARN" if chk["status"] == "warn" else "FAIL")
        print(f"    [{tag}] {chk['name']}: {chk['message']}")

    # ================================================================
    # 7. 多釜串联体积对比
    # ================================================================
    sep("7. 多釜串联体积对比 (cstr_design, 双分子动力学)")

    print(f"\n  {'N':>3s} | {'V_total(L)':>11s} | {'V_per(L)':>9s} | {'tau(s)':>8s} | {'X_intermediates':>20s}")
    print(f"  {'-'*3}-+-{'-'*11}-+-{'-'*9}-+-{'-'*8}-+-{'-'*20}")

    for N in range(1, 7):
        res_n = cstr_design(F_A0, X_A_target, k, C_A0, n, num_tanks=N, C_B0=C_B0)
        X_int = res_n.get("X_intermediates", [])
        X_str = ", ".join([f"{x:.3f}" for x in X_int]) if X_int else "-"
        print(f"  {N:3d} | {res_n['V_total']*1000:11.4f} | {res_n['V_per_tank']*1000:9.4f} | {res_n['tau']:8.2f} | {X_str:>20s}")

    # ================================================================
    # 8. 设计结论
    # ================================================================
    sep("8. 设计结论")

    print(f"""
  反应: 甲醇(A) + 乙酸(B) -> 产物, r = k*CA*CB

  进料条件:
    混合进料总流量 F_total = {F_total} mol/s
    甲醇摩尔流量 F_A0 = {F_A0:.5f} mol/s
    乙酸摩尔流量 F_B0 = {F_B0:.5f} mol/s
    体积流量 v0 = {v0:.6f} m3/s

  进料浓度:
    C_A0 = {C_A0} mol/m3 (甲醇)
    C_B0 = {C_B0} mol/m3 (乙酸)
    CB0/CA0 = {C_B0/C_A0:.2f} (乙酸为限制组分, 甲醇过量)

  转化率:
    乙酸目标转化率 X_B = {X_B_target*100:.0f}%
    对应甲醇转化率 X_A = {X_A_target*100:.2f}%
    甲醇最大理论转化率 = {X_A_max*100:.1f}% (由化学计量比限制)

  CSTR 精确设计 (双分子动力学):
    反应器体积 V = {V_cstr:.6f} m3 = {V_cstr*1000:.4f} L
    空间时间 tau = {tau_cstr:.4f} s

  PFR 精确设计 (解析解):
    空间时间 tau = {tau_pfr:.4f} s
    反应器体积 V = {V_pfr:.6f} m3 = {V_pfr*1000:.4f} L

  CSTR/PFR 体积比 = {tau_cstr/tau_pfr:.2f}

  出口浓度:
    甲醇 C_A_out = {C_A_out:.2f} mol/m3
    乙酸 C_B_out = {C_B_out:.2f} mol/m3
""")

    print(f"{'='*70}")
    print(f"  计算完成!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
