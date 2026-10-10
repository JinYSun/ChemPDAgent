"""VCM 正确回流设计：recycle to extinction（总循环，无弛放）。
关键：HCl 略过量 = 新鲜乙炔 + 产品带走的 HCl 杂质，使每种原子都有出路。
反应器 ext = 0.9*min(限制反应物)（物理正确的"90%限制反应物"）。
外层自适应调 fresh_HCl 匹配塔底 HCl 杂质；内层收敛循环。全程项目工具(mixer+tray)。
"""
import warnings
warnings.filterwarnings("ignore")
from mixer.funcs import calc_outlet_molar_flows as mixer_outlet
from tray_distillation.funcs import calc_mass_balance
from pump.funcs import calc_outlet_pressure as pump_outlet
from heat_exchangers.funcs import (
    calc_heat_duty, calc_log_mean_temp_difference, calc_heat_transfer_area)

C = ["Acetylene", "Hydrogen Chloride", "Vinyl Chloride"]
FRESH_A = 1500.0            # 新鲜乙炔固定
CONV = 0.90
LK, HK = "Hydrogen Chloride", "Vinyl Chloride"
LIGHT = ["Acetylene", "Hydrogen Chloride"]; HEAVY = ["Vinyl Chloride"]
D_PUR, W_PUR = 0.990, 0.995
def full(d): return {c: d.get(c, 0.0) for c in C}

fresh_h = 1500.0           # 初值，外层自适应
for outer in range(200):
    FRESH = {"Acetylene": FRESH_A, "Hydrogen Chloride": fresh_h, "Vinyl Chloride": 0.0}
    recycle = {c: 0.0 for c in C}
    for it in range(20000):
        mx = mixer_outlet([full(FRESH), full(recycle)]); s = full(mx["outlet_molar_flows_mol_per_s"])
        # 反应器：90% × 限制反应物（物理正确，不会消耗超过存量）
        ext = CONV * min(s["Acetylene"], s["Hydrogen Chloride"])
        s_rx = {"Acetylene": s["Acetylene"]-ext, "Hydrogen Chloride": s["Hydrogen Chloride"]-ext,
                "Vinyl Chloride": s["Vinyl Chloride"]+ext}
        col_feed = {c: max(v, 1e-9) for c, v in s_rx.items()}
        mb = calc_mass_balance(col_feed, D_PUR, W_PUR, LK, HK, LIGHT, HEAVY)
        dist = full(mb["distillate_flows_mol_per_s"])
        diff = sum(abs(dist[c] - recycle[c]) for c in C)
        recycle = {c: 0.5*recycle[c] + 0.5*dist[c] for c in C}
        if diff < 1e-7: break
    bot = full(mb["bottoms_flows_mol_per_s"])
    # 外层：新鲜HCl 应 = 消耗(=ext,稳态等于新鲜乙炔) + 塔底HCl杂质
    #       但乙炔可能未完全转化 -> 用乙炔净累积调HCl
    a_accum = FRESH_A - (ext + bot["Acetylene"])   # 乙炔进 - (消耗+塔底)  >0 表示乙炔在系统里没出路
    fresh_h_new = fresh_h + 0.8 * a_accum          # 乙炔憋着->说明ext不够->需更多HCl支持反应
    if abs(fresh_h_new - fresh_h) < 1e-4: fresh_h = fresh_h_new; break
    fresh_h = fresh_h_new

sh = {"Acetylene":"C2H2","Hydrogen Chloride":"HCl","Vinyl Chloride":"C2H3Cl"}
def row(n,d): return f"{n:12s} "+" ".join(f"{sh[c]}={d[c]:9.2f}" for c in C)+f"  Σ={sum(d.values()):9.2f}"
print(f"外层收敛 outer={outer}, 反算 新鲜HCl={fresh_h:.2f} mol/s (乙炔={FRESH_A})")
print(f"内层循环 iters={it}, resid={diff:.1e}\n")
FRESH = {"Acetylene": FRESH_A, "Hydrogen Chloride": fresh_h, "Vinyl Chloride": 0.0}
print("物料平衡 (mol/s):")
for n,d in [("新鲜进料",full(FRESH)),("循环回流",recycle),("混合器出口",s),
            ("反应器出口",full(s_rx)),("塔顶D(全回流)",dist),("塔底W(产品)",bot)]:
    print(row(n,d))
print(f"  反应进度 ext={ext:.2f}")
print(f"\n>>> 产品VCM = {bot['Vinyl Chloride']:.2f} mol/s")
print(f"[校验-无弛放] C2H2:{FRESH_A} = 消耗{ext:.2f}+塔底{bot['Acetylene']:.3f} = {ext+bot['Acetylene']:.2f}")
print(f"              HCl :{fresh_h:.2f} = 消耗{ext:.2f}+塔底{bot['Hydrogen Chloride']:.2f} = {ext+bot['Hydrogen Chloride']:.2f}")
overall = ext / FRESH_A * 100
print(f"\n乙炔总转化率 = {overall:.2f}% (单程90%, 循环到耗尽)")
print(f"HCl 过量投料 = {fresh_h-FRESH_A:.2f} mol/s = 产品带走的HCl杂质(非浪费)")

print("\n=== 各设备参数 ===")
p1=pump_outlet(inlet_pressure_Pa=101325.0,target_pressure_Pa=0.15e6)
print(f"P1泵 常压→0.15MPa 压升{p1['pressure_rise_Pa']:.0f}Pa")
d1=calc_heat_duty(s,30+273.15,130+273.15,0.15e6); l1=calc_log_mean_temp_difference(30+273.15,130+273.15,250+273.15,220+273.15)
a1=calc_heat_transfer_area(s,30+273.15,130+273.15,0.15e6,250+273.15,220+273.15,utility_type="steam")
print(f"E1加热30→130℃ Q={d1['heat_duty_W']/1e6:.2f}MW LMTD={l1['log_mean_temp_difference_K']:.1f}K A={a1['heat_transfer_area_m2']:.0f}m²")
d2=calc_heat_duty(full(s_rx),130+273.15,120+273.15,0.15e6); l2=calc_log_mean_temp_difference(130+273.15,120+273.15,20+273.15,40+273.15)
a2=calc_heat_transfer_area(full(s_rx),130+273.15,120+273.15,0.15e6,20+273.15,40+273.15,utility_type="cooling_water")
print(f"E2冷却130→120℃ Q={d2['heat_duty_W']/1e6:.2f}MW LMTD={l2['log_mean_temp_difference_K']:.1f}K A={a2['heat_transfer_area_m2']:.0f}m²")
p2=pump_outlet(inlet_pressure_Pa=0.15e6,target_pressure_Pa=4.43e6)
print(f"P2泵 0.15→4.43MPa 压升{p2['pressure_rise_Pa']/1e6:.2f}MPa")
print(f"T1塔 D_total={mb['distillate_total_mol_per_s']:.1f} W_total={mb['bottoms_total_mol_per_s']:.1f}mol/s")
