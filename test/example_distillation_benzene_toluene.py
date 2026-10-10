# -*- coding: utf-8 -*-
"""苯-甲苯常压精馏塔设计示例"""
import sys, io, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.distillation_device import _full_design

result = _full_design(
    components=["benzene", "toluene"],  # 苯-甲苯体系
    z=[0.4, 0.6],                       # 进料摩尔组成：苯40%、甲苯60%
    F=100.0,                            # 进料流量 100 kmol/h
    x_D_spec=0.99,                      # 塔顶苯纯度 ≥ 0.99
    x_B_spec=0.01,                      # 塔底苯摩尔分数 ≤ 0.01 (甲苯纯度 ≥ 0.99)
    P=101325.0,                         # 常压操作 101.325 kPa
    q=1.0,                              # 饱和液体进料
    RR_factor=2.0,                      # 回流比 = 2.0 × R_min
    tray_type='sieve',                  # 筛板塔
    flooding_fraction=0.80,             # 液泛分率 80%
    tray_spacing=0.60,                  # 塔板间距 0.60 m
)

print("=" * 70)
print("  苯-甲苯常压精馏塔设计结果")
print("=" * 70)

# 设计摘要
s = result.get("design_summary", {})
print(f"\n--- 设计摘要 ---")
print(f"  相对挥发度 α     = {s.get('alpha_LK_HK', '?')}")
print(f"  最小回流比 R_min = {s.get('R_min', '?')}")
print(f"  实际回流比 R     = {s.get('R_actual', '?')}")
print(f"  最小理论板数 N_min = {s.get('N_min', '?')}")
print(f"  理论板数 N_theo  = {s.get('N_theoretical', '?')}")
print(f"  实际板数 N_actual = {s.get('N_actual', '?')}")
print(f"  塔板效率 E0      = {s.get('tray_efficiency', '?')}")
print(f"  进料位置 N_feed  = {s.get('N_feed', '?')}")

# 塔径
print(f"\n--- 塔径 ---")
print(f"  塔径 D           = {s.get('D_column_m', '?')} m")
print(f"  塔高 H           = {s.get('H_column_m', '?')} m")
print(f"  空塔气速 u       = {s.get('u_actual_m_s', '?')} m/s")
print(f"  液泛气速 u_flood = {s.get('u_flood_m_s', '?')} m/s")
print(f"  液泛分率         = {s.get('flooding_fraction', '?')}")

# 物料衡算
print(f"\n--- 物料衡算 ---")
print(f"  进料 F           = {s.get('F_kmolh', '?')} kmol/h")
print(f"  塔顶 D           = {s.get('D_kmolh', '?')} kmol/h")
print(f"  塔底 B           = {s.get('B_kmolh', '?')} kmol/h")

# 热负荷
print(f"\n--- 热负荷 ---")
print(f"  冷凝器 Qc        = {s.get('Q_condenser_kW', '?')} kW")
print(f"  再沸器 Qr        = {s.get('Q_reboiler_kW', '?')} kW")

# 水力学
print(f"\n--- 水力学 ---")
print(f"  堰长 l_w         = {s.get('weir_length_m', '?')} m")
print(f"  堰高 h_w         = {s.get('weir_height_mm', '?')} mm")
print(f"  清液层高度 h_L   = {s.get('clear_liquid_height_mm', '?')} mm")
print(f"  单板压降 ΔP      = {s.get('tray_pressure_drop_Pa', '?')} Pa")
print(f"  降液管流速 u_dn  = {s.get('downcomer_velocity_m_s', '?')} m/s")

# 校核
print(f"\n--- 校核结果 ---")
val = result.get("validation", {})
print(f"  液泛校核         = {val.get('flooding_check', '?')}")
print(f"  夹带校核         = {val.get('entrainment_check', '?')}")
print(f"  压降校核         = {val.get('pressure_drop_check', '?')}")
print(f"  漏液校核         = {val.get('weeping_check', '?')}")

# q 值说明
q_notes = result.get("q_auto_notes", [])
if q_notes:
    print(f"\n--- q 值计算 ---")
    for note in q_notes:
        print(f"  {note}")

# 温度
print(f"\n--- 温度 ---")
print(f"  塔顶温度 T_top   = {s.get('T_top_K', '?')} K")
print(f"  塔底温度 T_bot   = {s.get('T_bottom_K', '?')} K")
print(f"  进料泡点 T_bub   = {s.get('T_feed_bubble_K', '?')} K")

# 警告
warnings = result.get("warnings", [])
if warnings:
    print(f"\n--- 警告 ---")
    for w in warnings:
        print(f"  ⚠ {w}")

# 错误
if "error" in result:
    print(f"\n--- 错误 ---")
    print(f"  {result['error']}")
