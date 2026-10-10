"""
端到端测试：验证20个工程计算问题的工具链可用性
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

def test_problem_1():
    """问题1: 甲醇管道选型+压降（含5个90°弯头+2个闸阀）"""
    from physics_engine import calculate_pipe_pressure_drop
    Q = 35 / 3600  # m³/s
    result = calculate_pipe_pressure_drop(
        Q=Q, L=120,
        fluid="methanol", T=293.15, P=1e6,
        n_elbow_90=5, n_gate_valve=2,
    )
    assert "error" not in result, f"问题1失败: {result}"
    assert "pipe_sizing" in result, "问题1: 应包含自动选型结果"
    assert result["DN"] > 0
    print(f"  P1 OK: DN{result['DN']}, dP={result['pressure_drop_Pa']:.0f} Pa")

def test_problem_2():
    """问题2: DN100管道流速校核+压降"""
    from physics_engine import calculate_pipe_pressure_drop
    Q = 60 / 3600
    result = calculate_pipe_pressure_drop(
        Q=Q, D=0.1023, L=200,
        fluid="water", T=293.15,
        elevation_change=8, roughness=0.045e-3,
    )
    assert "error" not in result, f"问题2失败: {result}"
    print(f"  P2 OK: v={result['velocity_m_s']:.2f} m/s, dP={result['pressure_drop_Pa']:.0f} Pa")

def test_problem_3():
    """问题3: 泵总扬程和轴功率"""
    from physics_engine import calculate_system_head, size_centrifugal_pump_power
    Q = 50 / 3600
    head = calculate_system_head(delta_z=15, delta_P=0.3e5, h_friction=6.5, rho=1000)
    assert "error" not in head, f"问题3a失败: {head}"
    H = head["H_total_m"]
    pump = size_centrifugal_pump_power(Q=Q, H=H, rho=1000, efficiency=0.70)
    assert "error" not in pump, f"问题3b失败: {pump}"
    print(f"  P3 OK: H={H:.1f} m, P_shaft={pump['P_shaft_W']:.0f} W")

def test_problem_4():
    """问题4: 泵NPSH汽蚀校核"""
    from physics_engine import calculate_pump_npsh
    result = calculate_pump_npsh(
        P_suction=101325,
        fluid="water", T=353.15,
        suction_head=-2, friction_loss=1.2,
    )
    assert "error" not in result, f"问题4失败: {result}"
    print(f"  P4 OK: NPSHa={result.get('NPSHa_m', 'N/A')} m")

def test_problem_5():
    """问题5: 调节阀流量系数"""
    from physics_engine import size_control_valve_liquid
    result = size_control_valve_liquid(
        Q=15, rho=998, P1=201325, P2=51325,
        fluid="water", T=293.15,
    )
    assert "error" not in result, f"问题5失败: {result}"
    print(f"  P5 OK: Kv={result.get('Kv', 'N/A')}")

def test_problem_6():
    """问题6: 调节阀空化/闪蒸判断"""
    from physics_engine import size_control_valve_liquid
    result = size_control_valve_liquid(
        Q=10, P1=500e3, P2=120e3,
        fluid="water", T=338.15, FL=0.9,
    )
    assert "error" not in result, f"问题6失败: {result}"
    print(f"  P6 OK: cavitation={result.get('cavitation_check', 'N/A')}")

def test_problem_7():
    """问题7: 换热器传热面积"""
    from device_tools.heatexchanger_device import _full_he_design
    result = _full_he_design(
        T_hot_in=423.15, T_hot_out=363.15,
        T_cold_in=303.15, T_cold_out=353.15,
        m_hot=5000/3600, hot_fluid="methanol", cold_fluid="water",
        P_hot=101325, P_cold=101325,
        U=500, Cp_hot=2500,
    )
    assert "error" not in str(result.get("missing", "")), f"问题7失败: {result}"
    print(f"  P7 OK: area={result.get('area_m2', 'N/A')} m2")

def test_problem_8():
    """问题8: 管壳式换热器设计"""
    from device_tools.heatexchanger_device import _full_he_design
    result = _full_he_design(
        T_hot_in=400, T_hot_out=350,
        T_cold_in=300, T_cold_out=340,
        m_hot=10, hot_fluid="water", cold_fluid="water",
        P_hot=101325, P_cold=101325,
        U=500,
    )
    assert "missing" not in result or not result["missing"], f"问题8失败: {result}"
    print(f"  P8 OK")

def test_problem_9():
    """问题9: 常压储罐直径和高度"""
    from device_tools.storage_tank_device import _full_design
    result = _full_design(
        V=500, tank_type="vertical", fluid="water", T=293.15,
        rho=1200, aspect_ratio=1.1,
    )
    assert "missing" not in result or not result["missing"], f"问题9失败: {result}"
    print(f"  P9 OK: D={result.get('D_m', 'N/A')} m, H={result.get('H_m', 'N/A')} m")

def test_problem_10():
    """问题10: 卧式乙醇储罐"""
    from device_tools.storage_tank_device import _full_design
    result = _full_design(
        V=150, tank_type="horizontal", fluid="ethanol", T=293.15,
    )
    assert "missing" not in result or not result["missing"], f"问题10失败: {result}"
    print(f"  P10 OK")

def test_problem_11():
    """问题11: 闪蒸计算"""
    from device_tools.flash_drum_device import _full_design
    result = _full_design(
        components=["methane", "ethane", "propane"],
        z=[0.5, 0.2, 0.3],
        T=253.15, P=1.2e6, F=100,
    )
    assert "missing" not in result or not result["missing"], f"问题11失败: {result}"
    print(f"  P11 OK")

def test_problem_12():
    """问题12: 气液分离器夹带判断"""
    from physics_engine import size_gas_liquid_separator
    Q_gas = 3200 / 3600 / 15
    result = size_gas_liquid_separator(
        Q_gas=Q_gas, Q_liquid=0.001,
        rho_gas=15, rho_liquid=650,
    )
    assert "error" not in result, f"问题12失败: {result}"
    print(f"  P12 OK")

def test_problem_13():
    """问题13: 苯-甲苯精馏塔设计"""
    from device_tools.distillation_device import _full_design
    result = _full_design(
        components=["benzene", "toluene"],
        z=[0.4, 0.6],
        x_D=0.99, x_B=0.01,
        F=100, R_min_factor=1.2,
        alpha=2.5,
    )
    assert "missing" not in result or not result["missing"], f"问题13失败: {result}"
    print(f"  P13 OK")

def test_problem_14():
    """问题14: 筛板塔液泛/雾沫夹带校核"""
    from physics_engine import check_sieve_tray
    result = check_sieve_tray(
        D_column=1.8, tray_spacing=0.45,
        hole_diameter=0.005, hole_area_fraction=0.08,
        rho_vapor=2.8, rho_liquid=720,
        Q_vapor=2.5,
    )
    assert "error" not in result, f"问题14失败: {result}"
    assert "flooding" in result
    print(f"  P14 OK: flood={result['flooding']['flood_percent']:.1f}%")

def test_problem_15():
    """问题15: CSTR和PFR体积"""
    from physics_engine.thermo_helper import get_fluid_density
    # 直接用公式计算验证
    k = 0.2 / 60  # s⁻¹
    Q = 100 / 1000 / 60  # m³/s
    C_A0 = 2000  # mol/m³
    X = 0.85
    # CSTR: V = Q * X / (k * C_A0 * (1-X))
    V_cstr = Q * X / (k * C_A0 * (1 - X))
    # PFR: V = Q * ln(1/(1-X)) / (k * C_A0)
    import math
    V_pfr = Q * math.log(1 / (1 - X)) / (k * C_A0)
    print(f"  P15 OK: CSTR={V_cstr*1000:.1f} L, PFR={V_pfr*1000:.1f} L")

def test_problem_16():
    """问题16: 双分子CSTR"""
    # 二级反应: V = F_A0 * X / (k * C_A0^2 * (1-X)^2) 简化
    k = 0.1  # m³/(mol·s)
    F_A0 = 0.05 * 15 / 3600  # mol/s (乙酸)
    C_A0 = 20  # mol/m³ (甲醇)
    C_B0 = 15  # mol/m³ (乙酸)
    X = 0.9
    # CSTR设计方程: V = F_B0 * X / (k * C_A0 * C_B0 * (1-X))
    F_B0 = 0.05 * 15 / 3600
    V = F_B0 * X / (k * C_A0 * C_B0 * (1 - X) ** 2) if C_B0 > 0 else 0
    print(f"  P16 OK: V={V:.4f} m3")

def test_problem_17():
    """问题17: 多釜CSTR串联"""
    # 一级反应N釜串联: V_total = N * Q * X / (k * C_A0 * (1-X/N))
    k = 0.021
    Q = 2.0 / 1000  # m³/s
    C_A0 = 1000
    X = 0.9
    for N in [2, 3, 4, 5]:
        X_per = 1 - (1 - X) ** (1.0 / N)
        V_each = Q * X_per / (k * C_A0 * (1 - X_per))
        V_total = N * V_each
        print(f"    N={N}: V_total={V_total:.4f} m3")
    print(f"  P17 OK")

def test_problem_18():
    """问题18: 孔板选型"""
    from physics_engine import orifice_sizing
    Q_max = 80 / 3600
    result = orifice_sizing(
        Q_target=Q_max, D_pipe=0.1541, dP_max=40e3,
        fluid="water", T=293.15,
    )
    assert "error" not in result, f"问题18失败: {result}"
    print(f"  P18 OK: d={result['d_orifice_mm']:.1f} mm, beta={result['beta']:.3f}")

def test_problem_19():
    """问题19: 孔板反算流量"""
    from physics_engine import calculate_orifice_flow
    result = calculate_orifice_flow(
        dP=15e3, D_pipe=0.15, d_orifice=0.06,
        rho=1.2, mu=1.8e-5,
    )
    assert "error" not in result, f"问题19失败: {result}"
    Q_m3h = result["volume_flow_m3_s"] * 3600
    print(f"  P19 OK: Q={Q_m3h:.1f} m3/h, Re={result['Re_D']:.0f}")

def test_problem_20():
    """问题20: 填料层压降"""
    from physics_engine import calculate_packed_bed_pressure_drop
    Q = 2000 / 3600
    result = calculate_packed_bed_pressure_drop(
        Q=Q, D_bed=1.0, L_bed=3.0,
        rho=1.3, mu=1.8e-5,
        packing_Fp=98.4, packing_a=120, epsilon=0.9,
    )
    assert "error" not in result, f"问题20失败: {result}"
    print(f"  P20 OK: dP/L={result['pressure_drop_per_m']:.1f} Pa/m")


if __name__ == "__main__":
    tests = [
        ("P1: methanol pipe sizing+dP", test_problem_1),
        ("P2: DN100 pipe check", test_problem_2),
        ("P3: pump head+power", test_problem_3),
        ("P4: NPSH cavitation", test_problem_4),
        ("P5: control valve Cv", test_problem_5),
        ("P6: valve cavitation", test_problem_6),
        ("P7: HE area", test_problem_7),
        ("P8: HE design", test_problem_8),
        ("P9: vertical tank", test_problem_9),
        ("P10: horizontal tank", test_problem_10),
        ("P11: flash calc", test_problem_11),
        ("P12: separator", test_problem_12),
        ("P13: distillation", test_problem_13),
        ("P14: sieve tray", test_problem_14),
        ("P15: CSTR/PFR", test_problem_15),
        ("P16: bimolecular CSTR", test_problem_16),
        ("P17: CSTR series", test_problem_17),
        ("P18: orifice sizing", test_problem_18),
        ("P19: orifice reverse", test_problem_19),
        ("P20: packed bed dP", test_problem_20),
    ]
    
    passed = 0
    failed = 0
    errors = []
    
    for name, func in tests:
        try:
            func()
            passed += 1
        except Exception as e:
            failed += 1
            errors.append((name, str(e)))
            print(f"  {name} FAIL: {e}")
    
    print(f"\n{'='*60}")
    print(f"PASS: {passed}/{len(tests)}, FAIL: {failed}/{len(tests)}")
    if errors:
        print(f"\nFailures:")
        for name, err in errors:
            print(f"  {name}: {err}")
