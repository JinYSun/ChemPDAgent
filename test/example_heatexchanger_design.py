# -*- coding: utf-8 -*-
"""管壳式换热器设计示例"""
import sys, io, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.heatexchanger_device import _full_he_design

result = _full_he_design(
    T_hot_in=360.0,      # K  热侧入口温度
    T_hot_out=350.0,     # K  热侧出口温度
    m_hot=1.0,           # kg/s  热侧质量流量
    hot_fluid="water",   # 热流体
    T_cold_in=300.0,     # K  冷侧入口温度
    T_cold_out=320.0,    # K  冷侧出口温度
    cold_fluid="water",  # 冷流体
    P_hot=300000.0,      # Pa  热侧压力 (300 kPa，360K水远低于沸点，液态安全)
    P_cold=300000.0,     # Pa  冷侧压力
    d_o=0.015,           # m  管外径 15mm（减小管径→更多管数→更大壳径→提升壳程流速）
    L=1.5,               # m  管长 1.5m（短管降低管程压降，改善L/D）
)

print("=" * 70)
print("  管壳式换热器设计结果")
print("=" * 70)

# 设计摘要
s = result.get("design_summary", {})
print(f"\n--- 设计摘要 ---")
print(f"  热负荷 Q         = {s.get('Q_kW', '?')} kW")
print(f"  LMTD             = {s.get('LMTD_K', '?')} K")
print(f"  F 修正因子       = {s.get('F_factor', '?')}")
print(f"  总传热系数 U     = {s.get('U_W_m2K', '?')} W/m²·K")
print(f"  所需传热面积     = {s.get('A_required_m2', '?')} m²")
print(f"  设计传热面积     = {s.get('A_design_m2', '?')} m²")
print(f"  换热管数         = {s.get('N_tubes', '?')}")
print(f"  管程数           = {s.get('N_tube_passes', '?')}")
print(f"  壳径 D_shell     = {s.get('D_shell_m', '?')} m")
print(f"  管长 L           = {s.get('L_tube_m', '?')} m")
print(f"  管外径 d_o        = {s.get('d_o_m', '?')} m")
print(f"  管内径 d_i        = {s.get('d_i_m', '?')} m")
print(f"  管排列           = {s.get('tube_layout', '?')}")
print(f"  折流板数         = {s.get('N_baffles', '?')}")
print(f"  折流板间距       = {s.get('baffle_spacing_m', '?')} m")
print(f"  长径比 L/D       = {s.get('L_D_ratio', '?')}")

# 热平衡
hb = result.get("heat_balance", {})
print(f"\n--- 热平衡 ---")
print(f"  热侧放热 Q_hot   = {hb.get('Q_hot_kW', '?')} kW")
print(f"  冷侧吸热 Q_cold  = {hb.get('Q_cold_kW', '?')} kW")
print(f"  热平衡偏差       = {hb.get('imbalance_pct', '?')}%")

# ε-NTU
ntu = result.get("epsilon_ntu", {})
print(f"\n--- ε-NTU 校核 ---")
print(f"  效能 ε           = {ntu.get('epsilon', '?')}")
print(f"  NTU              = {ntu.get('NTU', '?')}")
print(f"  热容比 C_r       = {ntu.get('C_r', '?')}")

# 传热
ht = result.get("heat_transfer", {})
print(f"\n--- 传热详情 ---")
print(f"  管侧 h_tube      = {ht.get('h_tube_W_m2K', '?')} W/m²·K")
print(f"  壳侧 h_shell     = {ht.get('h_shell_W_m2K', '?')} W/m²·K")
print(f"  管侧流速         = {ht.get('u_tube_m_s', '?')} m/s")
print(f"  壳侧流速         = {ht.get('u_shell_m_s', '?')} m/s")
print(f"  管侧 Re          = {ht.get('Re_tube', '?')}")
print(f"  壳侧 Re          = {ht.get('Re_shell', '?')}")
print(f"  管侧流态         = {ht.get('tube_flow_regime', '?')}")
print(f"  壳侧流态         = {ht.get('shell_flow_regime', '?')}")

# 压降
dp = result.get("pressure_drop", {})
print(f"\n--- 压降 ---")
print(f"  管程压降         = {dp.get('tube_Pa', '?')} Pa")
print(f"  壳程压降         = {dp.get('shell_Pa', '?')} Pa")

# U 分解
ud = result.get("U_decomposition", {})
print(f"\n--- U 分解 ---")
print(f"  总热阻 R_total   = {ud.get('R_total', '?')} m²·K/W")
print(f"  管侧对流热阻    = {ud.get('R_conv_tube', '?')}")
print(f"  壳侧对流热阻    = {ud.get('R_conv_shell', '?')}")
print(f"  壁面热阻         = {ud.get('R_wall', '?')}")
print(f"  管侧污垢热阻    = {ud.get('R_fouling_tube', '?')}")
print(f"  壳侧污垢热阻    = {ud.get('R_fouling_shell', '?')}")
print(f"  污垢占比         = {ud.get('fouling_fraction_pct', '?')}%")

# 物性
props = result.get("properties", {})
print(f"\n--- 流体物性 ---")
print(f"  热侧 Cp          = {props.get('Cp_hot_J_kgK', '?')} J/(kg·K)")
print(f"  冷侧 Cp          = {props.get('Cp_cold_J_kgK', '?')} J/(kg·K)")
print(f"  热侧密度         = {props.get('rho_hot_kg_m3', '?')} kg/m³")
print(f"  冷侧密度         = {props.get('rho_cold_kg_m3', '?')} kg/m³")

# 校验
val = result.get("validation", {})
print(f"\n--- 校验 ---")
print(f"  校验结果         = {val}")
print(f"  收敛             = {result.get('converged', '?')}")

# 警告
warnings = result.get("warnings", [])
if warnings:
    print(f"\n--- 警告 ---")
    for w in warnings:
        print(f"  ⚠ {w}")

# 默认值参数
if "user_defaulted_params_report" in result:
    print(f"\n{result['user_defaulted_params_report']}")

# 转换备注
notes = result.get("conversion_notes", [])
if notes:
    print(f"\n--- 单位换算备注 ---")
    for n in notes:
        print(f"  {n}")
