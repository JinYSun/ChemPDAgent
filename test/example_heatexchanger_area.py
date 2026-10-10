"""
换热器传热面积计算示例
========================

利用 physics_engine.heatexchanger 物理引擎，完成：
  1. 热负荷计算 (heat_duty)
  2. 对数平均温差 LMTD 计算 (逆流)
  3. 传热面积计算 (area_required)
  4. 换热器结构参数估算

已知条件:
  热流体: 150°C → 90°C, 流量 5000 kg/h, Cp = 2.5 kJ/(kg·K)
  冷流体: 30°C → 80°C
  逆流操作
  总传热系数 U = 500 W/(m²·K)
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

from physics_engine.heatexchanger import (
    lmtd,
    heat_duty,
    area_required,
    suggest_configuration,
)


def sep(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def main():
    # ================================================================
    # 1. 已知条件
    # ================================================================
    sep("1. 已知条件")

    T_hot_in = 150.0       # 热流体入口温度 (°C)
    T_hot_out = 90.0       # 热流体出口温度 (°C)
    T_cold_in = 30.0       # 冷流体入口温度 (°C)
    T_cold_out = 80.0      # 冷流体出口温度 (°C)
    m_hot_kgh = 5000.0     # 热流体流量 (kg/h)
    m_hot = m_hot_kgh / 3600.0  # 热流体流量 (kg/s)
    Cp_hot = 2.5           # 热流体比热 (kJ/(kg·K))
    Cp_hot_J = Cp_hot * 1000.0  # 热流体比热 (J/(kg·K))
    U = 500.0              # 总传热系数 (W/(m²·K))

    print(f"  温度条件:")
    print(f"    热流体入口 T_h,in  = {T_hot_in} °C")
    print(f"    热流体出口 T_h,out = {T_hot_out} °C")
    print(f"    冷流体入口 T_c,in  = {T_cold_in} °C")
    print(f"    冷流体出口 T_c,out = {T_cold_out} °C")
    print(f"  热流体参数:")
    print(f"    流量 m_h           = {m_hot_kgh} kg/h = {m_hot:.4f} kg/s")
    print(f"    比热 Cp_h          = {Cp_hot} kJ/(kg·K) = {Cp_hot_J:.0f} J/(kg·K)")
    print(f"  操作条件:")
    print(f"    流动方式           = 逆流 (counter-current)")
    print(f"    总传热系数 U       = {U} W/(m²·K)")

    # ================================================================
    # 2. 热负荷计算 (引擎函数)
    # ================================================================
    sep("2. 热负荷计算 (heat_duty)")

    # 手动计算
    Q_hot_manual = m_hot * Cp_hot_J * (T_hot_in - T_hot_out)  # W

    print(f"\n  手动计算:")
    print(f"    Q = m_h × Cp_h × (T_h,in - T_h,out)")
    print(f"      = {m_hot:.4f} × {Cp_hot_J:.0f} × ({T_hot_in} - {T_hot_out})")
    print(f"      = {m_hot:.4f} × {Cp_hot_J:.0f} × {T_hot_in - T_hot_out}")
    print(f"      = {Q_hot_manual:.2f} W = {Q_hot_manual/1000:.4f} kW")

    # 引擎计算
    # heat_duty 的 m 单位为 kg/s, Cp 单位 J/(kg·K)
    result_hd = heat_duty(
        m_hot=m_hot, Cp_hot=Cp_hot_J,
        T_hot_in=T_hot_in, T_hot_out=T_hot_out,
    )

    print(f"\n  引擎 heat_duty() 计算:")
    print(f"    Q_hot = {result_hd['Q_hot']:.2f} W = {result_hd['Q_hot']/1000:.4f} kW")
    print(f"    热平衡偏差 = {result_hd['balance_error']:.2f}%")

    Q = result_hd['Q_hot']  # 热负荷 (W)

    # ================================================================
    # 3. 对数平均温差 LMTD (引擎函数)
    # ================================================================
    sep("3. 对数平均温差 LMTD 计算 (lmtd)")

    # 逆流端差
    dT1 = T_hot_in - T_cold_out   # 热端 (150 - 80 = 70°C)
    dT2 = T_hot_out - T_cold_in   # 冷端 (90 - 30 = 60°C)

    print(f"\n  逆流端差:")
    print(f"    ΔT1 = T_h,in - T_c,out = {T_hot_in} - {T_cold_out} = {dT1} °C")
    print(f"    ΔT2 = T_h,out - T_c,in = {T_hot_out} - {T_cold_in} = {dT2} °C")

    # 手动 LMTD
    if abs(dT1 - dT2) < 0.01:
        LMTD_manual = (dT1 + dT2) / 2.0
    else:
        LMTD_manual = (dT1 - dT2) / math.log(dT1 / dT2)

    print(f"\n  手动计算 LMTD:")
    print(f"    LMTD = (ΔT1 - ΔT2) / ln(ΔT1/ΔT2)")
    print(f"         = ({dT1} - {dT2}) / ln({dT1}/{dT2})")
    print(f"         = {dT1 - dT2} / {math.log(dT1/dT2):.6f}")
    print(f"         = {LMTD_manual:.4f} °C")

    # 引擎计算 (逆流, 不施加壳程修正 F=1)
    result_lmtd = lmtd(
        T_hot_in, T_hot_out, T_cold_in, T_cold_out,
        flow_arrangement="counter",
        apply_correction=False,  # 纯逆流, F=1
    )

    print(f"\n  引擎 lmtd() 计算:")
    print(f"    LMTD_raw = {result_lmtd['LMTD_raw']:.4f} °C")
    print(f"    ΔT1 = {result_lmtd['dT1']} °C, ΔT2 = {result_lmtd['dT2']} °C")
    print(f"    最小接近温度 = {result_lmtd['approach_temp_min']} °C")
    print(f"    有效: {result_lmtd['valid']}")
    if result_lmtd['warnings']:
        for w in result_lmtd['warnings']:
            print(f"    警告: {w}")

    LMTD_val = result_lmtd['LMTD_raw']

    # ================================================================
    # 4. 传热面积计算 (引擎函数)
    # ================================================================
    sep("4. 传热面积计算 (area_required)")

    # 手动计算
    A_manual = Q / (U * LMTD_val)

    print(f"\n  基本传热方程: Q = U × A × LMTD")
    print(f"    A = Q / (U × LMTD)")
    print(f"      = {Q:.2f} / ({U} × {LMTD_val:.4f})")
    print(f"      = {A_manual:.4f} m2")

    # 引擎计算 (含15%面积裕量)
    overspec = 15.0
    result_area = area_required(
        Q=Q, U=U, LMTD=LMTD_val,
        overspec_percent=overspec,
        d_o=0.025, tube_length=6.0,
        tube_passes=2, tube_layout="triangular",
    )

    print(f"\n  引擎 area_required() 计算:")
    print(f"    净传热面积 A_req = {result_area['A_required']:.4f} m2")
    print(f"    设计面积 A_design = {result_area['A_design']:.4f} m2 (含{overspec}%裕量)")
    print(f"    估算管数 (φ25×6000) = {result_area['N_tubes']} 根")
    print(f"    估算壳径 = {result_area['D_shell_standard']} mm")
    print(f"    长径比 L/D = {result_area['length_diameter_ratio']:.2f}")
    if result_area['warnings']:
        for w in result_area['warnings']:
            print(f"    警告: {w}")

    # ================================================================
    # 5. 换热器结构方案建议
    # ================================================================
    sep("5. 换热器结构方案建议 (suggest_configuration)")

    sug = suggest_configuration(T_hot_in, T_hot_out, T_cold_in, T_cold_out)
    print(f"\n  温度效率 P = (T_c,out - T_c,in) / (T_h,in - T_c,in)")
    P_val = (T_cold_out - T_cold_in) / (T_hot_in - T_cold_in)
    R_val = (T_hot_in - T_hot_out) / (T_cold_out - T_cold_in)
    print(f"    P = ({T_cold_out}-{T_cold_in}) / ({T_hot_in}-{T_cold_in}) = {P_val:.4f}")
    print(f"  热容比 R = (T_h,in - T_h,out) / (T_c,out - T_c,in)")
    print(f"    R = ({T_hot_in}-{T_hot_out}) / ({T_cold_out}-{T_cold_in}) = {R_val:.4f}")
    print(f"\n  推荐配置: {sug['recommended_config']}")
    print(f"    F 因子 = {sug['F_factor']:.4f}")
    print(f"    壳程数 = {sug['n_shells']}")
    print(f"    说明: {sug['reasoning']}")
    if sug.get('alternatives'):
        print(f"    可选方案: {sug['alternatives']}")

    # 若采用多壳程, 修正面积
    if sug['F_factor'] < 1.0 and sug['F_factor'] > 0:
        F_corr = sug['F_factor']
        LMTD_corr = LMTD_val * F_corr
        A_corr = Q / (U * LMTD_corr)
        print(f"\n  若采用 {sug['recommended_config']} (F={F_corr:.4f}):")
        print(f"    修正 LMTD = {LMTD_val:.4f} × {F_corr:.4f} = {LMTD_corr:.4f} °C")
        print(f"    修正面积 A = {A_corr:.4f} m2")
        print(f"    (比纯逆流面积增大 {A_corr/A_manual:.2f} 倍)")

    # ================================================================
    # 6. 不同 U 值下的面积敏感性分析
    # ================================================================
    sep("6. 总传热系数敏感性分析")

    U_values = [300, 400, 500, 600, 700, 800]
    print(f"\n  {'U (W/m2·K)':>12s} | {'A_req (m2)':>12s} | {'A_design (m2)':>14s}")
    print(f"  {'-'*12}-+-{'-'*12}-+-{'-'*14}")

    for u_val in U_values:
        a_req = Q / (u_val * LMTD_val)
        a_des = a_req * (1 + overspec / 100)
        marker = " <-- 本题" if u_val == U else ""
        print(f"  {u_val:>12.0f} | {a_req:>12.4f} | {a_des:>14.4f}{marker}")

    # ================================================================
    # 7. 设计结论
    # ================================================================
    sep("7. 设计结论")

    print(f"""
  温度条件:
    热流体: {T_hot_in}°C → {T_hot_out}°C (降温 {T_hot_in - T_hot_out}°C)
    冷流体: {T_cold_in}°C → {T_cold_out}°C (升温 {T_cold_out - T_cold_in}°C)
    逆流端差: ΔT1 = {dT1}°C, ΔT2 = {dT2}°C

  热负荷:
    Q = {Q:.2f} W = {Q/1000:.4f} kW

  传热推动力:
    LMTD (逆流) = {LMTD_val:.4f} °C

  传热面积:
    净面积 A_req   = {A_manual:.4f} m2
    设计面积 (含{overspec:.0f}%裕量) = {result_area['A_design']:.4f} m2

  结构估算 (φ25×6000 三角排列):
    管数 ≈ {result_area['N_tubes']} 根
    壳径 ≈ {result_area['D_shell_standard']} mm
""")

    print(f"{'='*70}")
    print(f"  计算完成!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
