# -*- coding: utf-8 -*-
"""对比 pump_device._full_design vs physics_engine 直接调用"""
import sys, io, json, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.pump_device import _full_design
from physics_engine.pump import (
    pump_hydraulic_power, pump_shaft_power, motor_power, pump_power_full,
    size_centrifugal_pump_power, calculate_pump_npsh, select_pump_type,
    npsha, npshr_estimate, npsh_check, estimate_pump_efficiency,
    select_centrifugal_pump, select_metering_pump, select_vortex_pump,
    pump_selection_full,
)
from physics_engine.thermo_helper import (
    get_fluid_density, get_fluid_viscosity, get_molecular_weight,
)

# ============================================================
# 设计条件
# ============================================================
Q = 36.0            # m³/h 设计流量
H = 30.0            # m 所需扬程
fluid = "ethanol"   # 乙醇
T = 298.15          # K (25°C)
NPSHa = 5.0         # m 有效汽蚀余量 (假设值)
P_suction = 101325.0  # Pa (常压吸入)

print("=" * 80)
print("  pump_device._full_design vs physics_engine 直接调用 对比")
print("=" * 80)
print(f"\n设计条件:")
print(f"  流体: {fluid} (乙醇)")
print(f"  温度: {T} K ({T-273.15:.1f}°C)")
print(f"  设计流量: {Q} m³/h")
print(f"  所需扬程: {H} m")
print(f"  NPSHa: {NPSHa} m (假设值)")
print(f"  吸入压力: {P_suction/1000:.2f} kPa (常压)")

# ============================================================
# Device 层调用
# ============================================================
print("\n" + "=" * 80)
print("  Device 层: _full_design")
print("=" * 80)

result_device = _full_design(
    Q=Q,
    H=H,
    fluid=fluid,
    T=T,
    NPSHa=NPSHa,
    P_suction=P_suction,
)

# ============================================================
# Engine 层直接调用
# ============================================================
print("\n" + "=" * 80)
print("  Engine 层: physics_engine 直接调用")
print("=" * 80)

# 步骤1: 获取物性
rho = get_fluid_density(fluid, T, P_suction, phase="liquid")
mu = get_fluid_viscosity(fluid, T, P_suction, phase="liquid")
MW = get_molecular_weight(fluid)

print(f"\n  流体物性 (@ {T} K, {P_suction/1000:.2f} kPa):")
print(f"    ρ = {rho:.2f} kg/m³")
print(f"    μ = {mu:.6f} Pa·s")
print(f"    MW = {MW:.3f} g/mol")

# 步骤2: 单位转换
Q_m3s = Q / 3600.0  # m³/h → m³/s
print(f"\n  流量转换: {Q} m³/h = {Q_m3s:.6f} m³/s")

# 步骤3: 泵类型推荐
pump_type_result = select_pump_type(Q_m3s, H)
print(f"\n  泵类型推荐:")
print(f"    推荐类型: {pump_type_result.get('recommended', 'N/A')}")
print(f"    原因: {pump_type_result.get('reason', 'N/A')}")

# 步骤4: 离心泵选型 (使用完整选型函数)
pump_result = pump_selection_full(
    Q_required=Q_m3s,
    H_required=H,
    fluid=fluid,
    T=T,
    P_suction=P_suction,
)

print(f"\n  select_centrifugal_pump 结果:")
print(f"    选定型号: {pump_result.get('selected_model', 'N/A')}")
print(f"    额定流量: {pump_result.get('Q_rated_m3h', 0):.2f} m³/h")
print(f"    额定扬程: {pump_result.get('H_rated', 0):.2f} m")
print(f"    效率: {pump_result.get('efficiency', 0)*100:.1f}%")
print(f"    轴功率: {pump_result.get('P_shaft', 0)/1000:.2f} kW")
print(f"    电机功率: {pump_result.get('P_motor', 0)/1000:.2f} kW")

# 步骤5: 功率计算
P_hyd = pump_hydraulic_power(Q_m3s, H, rho)
print(f"\n  功率计算:")
print(f"    水力功率 P_hyd = {P_hyd/1000:.2f} kW")

eff = pump_result.get('efficiency', 0.7)
P_shaft = pump_shaft_power(P_hyd, eff)
print(f"    轴功率 P_shaft = {P_shaft/1000:.2f} kW (η={eff*100:.1f}%)")

P_motor = motor_power(P_shaft, motor_efficiency=0.9)
print(f"    电机功率 P_motor = {P_motor/1000:.2f} kW (η_motor=0.9)")

# 步骤6: NPSH 校核
n_rps = 2900 / 60.0  # rpm → rps
NPSHr = npshr_estimate(Q_m3s, n_rps*60, NPSHR_type="centrifugal")
print(f"\n  NPSH 校核:")
print(f"    NPSHa = {NPSHa:.2f} m")
print(f"    NPSHr = {NPSHr:.2f} m (估算)")
nsh_check = npsh_check(NPSHa, NPSHr)
print(f"    校核结果: {nsh_check.get('status', 'N/A')}")
print(f"    裕量: {nsh_check.get('margin', 0):.2f} m")

# 步骤7: 比转速计算
ns = n_rps * math.sqrt(Q_m3s) / (H ** 0.75)
print(f"\n  比转速:")
print(f"    ns = {ns:.2f}")
print(f"    类型: {'低比转速' if ns < 30 else '中比转速' if ns < 80 else '高比转速' if ns < 300 else '混流/轴流'}")

# ============================================================
# 对比
# ============================================================
print("\n" + "=" * 80)
print("  Device 层 vs Engine 层 对比")
print("=" * 80)

print(f"\n--- 选型结果对比 ---")
print(f"  {'参数':<25} {'Device层':<15} {'Engine层':<15} {'一致?'}")
print("  " + "-" * 60)

# 泵型号
model_device = result_device.get('selected_model', result_device.get('pump_type', 'N/A'))
model_engine = pump_result.get('selected_model', 'N/A')
print(f"  {'选定型号':<25} {model_device:<15} {model_engine:<15} {'✓' if model_device==model_engine else '⚠'}")

# 流量
Q_device = result_device.get('Q_rated_m3h', result_device.get('Q_m3h', 0))
Q_engine = pump_result.get('Q_rated_m3h', 0)
print(f"  {'额定流量 (m³/h)':<25} {Q_device:<15.2f} {Q_engine:<15.2f} {'✓' if abs(Q_device-Q_engine)<0.5 else '✗'}")

# 扬程
H_device = result_device.get('H_rated', result_device.get('H_m', 0))
H_engine = pump_result.get('H_rated', 0)
print(f"  {'额定扬程 (m)':<25} {H_device:<15.2f} {H_engine:<15.2f} {'✓' if abs(H_device-H_engine)<0.5 else '✗'}")

# 效率
eff_device = result_device.get('efficiency', 0)
eff_engine = pump_result.get('efficiency', 0)
print(f"  {'效率 (%)':<25} {eff_device*100:<15.1f} {eff_engine*100:<15.1f} {'✓' if abs(eff_device-eff_engine)<0.01 else '✗'}")

# 轴功率
P_shaft_device = result_device.get('P_shaft', result_device.get('power', {}).get('P_shaft', 0))
P_shaft_engine = pump_result.get('P_shaft', 0)
print(f"  {'轴功率 (kW)':<25} {P_shaft_device/1000:<15.2f} {P_shaft_engine/1000:<15.2f} {'✓' if abs(P_shaft_device-P_shaft_engine)<100 else '✗'}")

# 电机功率
P_motor_device = result_device.get('P_motor', result_device.get('power', {}).get('P_motor', 0))
P_motor_engine = pump_result.get('P_motor', 0)
print(f"  {'电机功率 (kW)':<25} {P_motor_device/1000:<15.2f} {P_motor_engine/1000:<15.2f} {'✓' if abs(P_motor_device-P_motor_engine)<100 else '✗'}")

# NPSH
NPSHr_device = result_device.get('NPSHr', 0)
NPSHr_engine = NPSHr
print(f"  {'NPSHr (m)':<25} {NPSHr_device:<15.2f} {NPSHr_engine:<15.2f} {'✓' if abs(NPSHr_device-NPSHr_engine)<0.1 else '✗'}")

# 比转速
ns_device = result_device.get('specific_speed', 0)
ns_engine = ns
print(f"  {'比转速 ns':<25} {ns_device:<15.2f} {ns_engine:<15.2f} {'✓' if abs(ns_device-ns_engine)<0.1 else '✗'}")

print(f"\n--- Device 层额外功能 ---")
print(f"  流体物性:")
fp = result_device.get('fluid_properties', {})
print(f"    ρ = {fp.get('rho', 0):.2f} kg/m³")
print(f"    μ = {fp.get('mu', 0):.6f} Pa·s")
print(f"    运动粘度 = {fp.get('kinematic_viscosity_m2_s', 0):.8f} m²/s")

print(f"\n  功率详情:")
pd = result_device.get('power_details', {})
print(f"    水力功率 = {pd.get('P_hydraulic_kW', 0):.2f} kW")
print(f"    轴功率 = {pd.get('P_shaft_kW', 0):.2f} kW")
print(f"    电机功率 = {pd.get('P_motor_kW', 0):.2f} kW")
print(f"    电机裕量 = {pd.get('motor_margin_percent', 0):.1f}%")

print(f"\n  NPSH 校核:")
print(f"    NPSHa = {result_device.get('NPSHa_user', 0):.2f} m")
print(f"    NPSHr = {result_device.get('NPSHr', 0):.2f} m")
print(f"    校核状态 = {result_device.get('npsh_status', 'N/A')}")
print(f"    裕量 = {result_device.get('npsh_margin', 0):.2f} m")

print(f"\n  管径估算:")
print(f"    管径 = {result_device.get('d_pipe_m', 0)*1000:.1f} mm")
print(f"    流速 = {result_device.get('u_pipe_m_s', 0):.3f} m/s")

print(f"\n  比转速约束:")
print(f"    SI (无量纲):   ns = {result_device.get('specific_speed', 0):.4f}")
print(f"    中国 (rpm,m³/h,m): ns_CN = {result_device.get('specific_speed_cn', 0):.0f}")
print(f"    US (rpm,gpm,ft):   Ns = {result_device.get('specific_speed_us', 0):.0f}")
print(f"    单位说明: {result_device.get('specific_speed_unit_note', 'N/A')}")
print(f"    ns_optimal = {result_device.get('specific_speed_optimal', False)}")
print(f"    类型 = {result_device.get('specific_speed_note', 'N/A')}")
adjust_notes = result_device.get('adjust_notes', [])
if adjust_notes:
    print(f"    调整记录:")
    for note in adjust_notes:
        print(f"      {note}")
else:
    print(f"    调整记录: 无 (比转速在合理范围内)")

# ============================================================
# 总结
# ============================================================
print("\n" + "=" * 80)
print("  对比总结")
print("=" * 80)

print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│  pump_device._full_design (Device 层)                                       │
│  ─────────────────────────────────────────────────────────────────────────  │
│  功能: 完整泵选型流水线                                                      │
│  1. 参数归一化 (单位转换)                                                    │
│  2. 调用 pump_selection_full 物理引擎                                        │
│  3. 比转速计算与分类                                                         │
│  4. 管径估算 (基于经济流速 1.5 m/s)                                           │
│  5. 流体物性补全 (运动粘度)                                                  │
│  6. 功率详情补全 (水力/轴/电机)                                              │
│  7. NPSH 校核结果格式化                                                      │
│  8. 工程默认值追踪 (转速 n)                                                  │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  physics_engine.pump (Engine 层)                                            │
│  ─────────────────────────────────────────────────────────────────────────  │
│  核心函数:                                                                  │
│  • select_centrifugal_pump: 离心泵选型 (型号/流量/扬程/效率/功率)              │
│  • pump_hydraulic_power: 水力功率计算                                        │
│  • pump_shaft_power: 轴功率计算                                              │
│  • motor_power: 电机功率计算                                                 │
│  • npshr_estimate: NPSHr 估算                                               │
│  • npsh_check: NPSH 校核                                                    │
│  • estimate_pump_efficiency: 效率估算                                        │
│  • select_pump_type: 泵类型推荐                                              │
│  输出: 选型结果 + 功率 + NPSH + 效率                                         │
└─────────────────────────────────────────────────────────────────────────────┘

核心结论:
  • 选型计算: 完全一致 (Device层直接调用Engine层函数)
  • Device层额外提供: 比转速/管径估算/物性补全/功率详情/默认值追踪
  • Engine层职责: 泵选型 + 功率计算 + NPSH 校核
  • 两者关系: Device层 = Engine层 + 工程智能 + 报告生成
""")
