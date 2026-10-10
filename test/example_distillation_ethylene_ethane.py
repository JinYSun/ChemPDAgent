# -*- coding: utf-8 -*-
"""乙烯-乙烷分类精馏塔设计 (C2 Splitter)
进料: 100 kmol/h, 乙烯40%/乙烷60%, 250K, 20bar
塔顶乙烯 > 98%, 塔底乙烯 < 0.5%
"""
import sys, io, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.distillation_device import _full_design

result = _full_design(
    components=["ethylene", "ethane"],  # 乙烯(LK)-乙烷(HK)体系
    z=[0.4, 0.6],                       # 进料摩尔组成: 乙烯40%, 乙烷60%
    F=100.0,                            # 进料流量 100 kmol/h
    x_D_spec=0.98,                      # 塔顶乙烯纯度 >= 98%
    x_B_spec=0.005,                     # 塔底乙烯 <= 0.5% (乙烷纯度 >= 99.5%)
    P=2000000.0,                        # 操作压力 20 bar = 2 MPa
    T_feed=250.0,                       # 进料温度 250 K (自动计算q值)
    RR_factor=2.0,                      # 回流比 = 2.0 x R_min
    tray_type='sieve',                  # 筛板塔
    flooding_fraction=0.80,             # 液泛分率 80%
    tray_spacing=0.60,                  # 塔板间距 0.60 m
)

print("=" * 70)
print("  乙烯-乙烷分类精馏塔设计结果 (C2 Splitter)")
print("=" * 70)

# === FUG 计算 ===
s = result.get("design_summary", {})
fug = result.get("fenske_underwood_gilliland", {})
print(f"\n--- FUG 设计 ---")
print(f"  相对挥发度 alpha   = {s.get('alpha_LK_HK', '?')}")
print(f"  最小理论板数 N_min = {s.get('N_min', '?')}")
print(f"  最小回流比 R_min   = {s.get('R_min', '?')}")
print(f"  实际回流比 R       = {s.get('R', '?')}")
print(f"  R/R_min            = {s.get('R_Rmin_ratio', '?')}")
print(f"  理论板数 N_theo    = {s.get('N_theoretical', '?')}")
print(f"  实际板数 N_actual  = {s.get('N_actual', '?')}")
print(f"  塔板效率 E0        = {s.get('tray_efficiency', '?')} {s.get('efficiency_note', '')}")
print(f"  进料位置 N_feed    = {s.get('feed_stage', '?')} (从塔顶数)")
print(f"  精馏段板数         = {s.get('N_rectifying', '?')}")
print(f"  提馏段板数         = {s.get('N_stripping', '?')}")
print(f"  进料热状态 q       = {s.get('q_state', '?')}")

# === 塔径塔高 ===
print(f"\n--- 塔体尺寸 ---")
print(f"  计算塔径 D_calc    = {s.get('D_column_calculated', '?')} m")
print(f"  标准塔径 D         = {s.get('D_column', '?')} m")
print(f"  塔高 H             = {s.get('H_column', '?')} m")
print(f"  H/D                = {s.get('H_D_ratio', '?')}")

# === 物料衡算 ===
mb = result.get("material_balance", {})
print(f"\n--- 物料衡算 ---")
print(f"  进料 F             = {mb.get('F_kmolh', '?')} kmol/h")
print(f"  塔顶 D             = {mb.get('D_kmolh', '?')} kmol/h")
print(f"  塔底 B             = {mb.get('B_kmolh', '?')} kmol/h")
print(f"  塔顶乙烯纯度       = {mb.get('x_D_LK_actual', '?')}")
print(f"  塔底乙烯含量       = {mb.get('x_B_LK_actual', '?')}")

# 组分分布
comp_dist = result.get("component_distribution", mb.get("component_distribution", []))
if comp_dist:
    print(f"\n--- 组分分布 ---")
    for c in comp_dist:
        print(f"  {c.get('component','?')}: 进料={c.get('feed_kmolh','?')} "
              f"塔顶={c.get('distillate_kmolh','?')} 塔底={c.get('bottoms_kmolh','?')} kmol/h")

# === 温度 ===
temp = result.get("temperatures", {})
print(f"\n--- 温度 ---")
print(f"  塔顶温度           = {temp.get('T_top_K', '?')} K ({temp.get('T_top_C', '?')} C)")
print(f"  塔底温度           = {temp.get('T_bottom_K', '?')} K ({temp.get('T_bottom_C', '?')} C)")
print(f"  进料泡点           = {temp.get('T_feed_bubble_K', '?')} K")
print(f"  进料露点           = {temp.get('T_feed_dew_K', '?')} K")

# === 热负荷 ===
hd = result.get("heat_duty", {})
print(f"\n--- 热负荷 ---")
print(f"  冷凝器 Qc          = {hd.get('Q_condenser_kW', '?')} kW")
print(f"  再沸器 Qr          = {hd.get('Q_reboiler_kW', '?')} kW")

# === 水力学 ===
hyd = result.get("hydraulics", {})
print(f"\n--- 水力学 ---")
print(f"  堰长 l_w           = {hyd.get('weir_length_m', '?')} m")
print(f"  堰高 h_w           = {hyd.get('weir_height_mm', '?')} mm")
print(f"  清液层高度 h_L     = {hyd.get('clear_liquid_height_mm', '?')} mm")
print(f"  单板压降 dP        = {hyd.get('pressure_drop_per_tray_Pa', '?')} Pa")
print(f"  降液管流速 u_dn    = {hyd.get('downcomer_velocity_m_s', '?')} m/s")
print(f"  空塔气速           = {hyd.get('vapor_velocity_m_s', '?')} m/s")

# === 校核 ===
val = result.get("validation", {})
print(f"\n--- 校核结果 ---")
for k, v in val.items():
    print(f"  {k}: {v}")

# === 物性 ===
props = result.get("properties", {})
print(f"\n--- 物性 ---")
print(f"  气相密度 rho_V     = {props.get('rho_vapor_kg_m3', '?')} kg/m3")
print(f"  液相密度 rho_L     = {props.get('rho_liquid_kg_m3', '?')} kg/m3")
print(f"  混合物粘度 mu      = {props.get('mu_feed_cP', props.get('mu_l_cP', '?'))} cP")

# === q 值 ===
q_info = result.get("q_info", result.get("q_auto_notes", []))
if q_info:
    print(f"\n--- q 值信息 ---")
    if isinstance(q_info, list):
        for note in q_info:
            print(f"  {note}")
    else:
        print(f"  {q_info}")

# === 警告 ===
warnings = result.get("warnings", [])
if warnings:
    print(f"\n--- 警告 ---")
    for w in warnings:
        print(f"  ! {w}")

# === 错误 ===
if "error" in result:
    print(f"\n--- 错误 ---")
    print(f"  {result['error']}")

# === 收敛 ===
print(f"\n--- 收敛状态 ---")
print(f"  converged = {result.get('converged', '?')}")
