"""physics_engine 物理计算工具注册模块

将 physics_engine 中的关键函数注册为 UniversalToolRegistry 工具，
供 planner/executor 在校核、验算、物性查询场景下调用。
"""
from tools.registry import UniversalToolRegistry, ToolCapability


def register_physics_tools(registry: UniversalToolRegistry) -> UniversalToolRegistry:
    """注册所有 physics_engine 物理计算工具"""

    # ================================================================
    # 热力学物性查询 (category: thermo)
    # ================================================================
    from physics_engine import (
        get_fluid_density,
        get_fluid_viscosity,
        get_fluid_Cp,
        get_fluid_vapor_pressure,
        get_fluid_MW,
        get_fluid_thermal_conductivity,
        get_fluid_surface_tension,
    )

    registry.register(ToolCapability(
        name="get_fluid_density",
        description="查询纯流体或混合物的质量密度 (kg/m³)。支持液相/气相，自动路由 CoolProp HEOS → thermo。",
        category="thermo",
        parameters={
            "fluid": {"type": "str|list", "required": True, "description": "流体英文名（纯流体）或组分名列表（混合物），如 'water' 或 ['benzene', 'toluene']"},
            "T": {"type": "number", "required": True, "description": "温度 (K)"},
            "P": {"type": "number", "required": False, "description": "压力 (Pa)，默认 101325"},
            "phase": {"type": "string", "required": False, "description": "相态: 'liquid'(默认) 或 'gas'"},
            "z": {"type": "list", "required": False, "description": "混合物摩尔分率列表（混合物必填）"},
        },
        tags=["density", "thermo", "物性"],
        aliases=["fluid_density", "density_calc", "calc_density"],
    ), get_fluid_density)

    registry.register(ToolCapability(
        name="get_fluid_viscosity",
        description="查询流体动力粘度 (Pa·s)。液相用 Arrhenius 混合，气相用 Wilke 方法。",
        category="thermo",
        parameters={
            "fluid": {"type": "str|list", "required": True, "description": "流体英文名或组分名列表"},
            "T": {"type": "number", "required": True, "description": "温度 (K)"},
            "P": {"type": "number", "required": False, "description": "压力 (Pa)，默认 101325"},
            "phase": {"type": "string", "required": False, "description": "相态: 'liquid'(默认) 或 'gas'"},
            "z": {"type": "list", "required": False, "description": "混合物摩尔分率列表"},
        },
        tags=["viscosity", "thermo", "物性"],
        aliases=["fluid_viscosity", "viscosity_calc", "calc_viscosity"],
    ), get_fluid_viscosity)

    registry.register(ToolCapability(
        name="get_fluid_Cp",
        description="查询流体定压摩尔热容 Cp (J/mol/K)。",
        category="thermo",
        parameters={
            "fluid": {"type": "str|list", "required": True, "description": "流体英文名或组分名列表"},
            "T": {"type": "number", "required": True, "description": "温度 (K)"},
            "P": {"type": "number", "required": False, "description": "压力 (Pa)，默认 101325"},
            "phase": {"type": "string", "required": False, "description": "相态: 'liquid'(默认) 或 'gas'"},
            "z": {"type": "list", "required": False, "description": "混合物摩尔分率列表"},
        },
        tags=["heat_capacity", "Cp", "thermo", "物性"],
        aliases=["fluid_Cp", "heat_capacity", "calc_Cp"],
    ), get_fluid_Cp)

    registry.register(ToolCapability(
        name="get_fluid_vapor_pressure",
        description="查询纯流体饱和蒸气压或混合物泡点压力 (Pa)。CoolProp HEOS 优先，fallback 到 thermo Wagner/Antoine。",
        category="thermo",
        parameters={
            "fluid": {"type": "str|list", "required": True, "description": "流体英文名或组分名列表"},
            "T": {"type": "number", "required": True, "description": "温度 (K)"},
            "z": {"type": "list", "required": False, "description": "混合物摩尔分率列表（混合物泡点必填）"},
        },
        tags=["vapor_pressure", "Psat", "thermo", "物性"],
        aliases=["fluid_vapor_pressure", "vapor_pressure_calc", "calc_vapor_pressure", "Psat"],
    ), get_fluid_vapor_pressure)

    registry.register(ToolCapability(
        name="get_fluid_MW",
        description="查询流体分子量 (g/mol)。纯流体直接查库，混合物按摩尔分率加权。",
        category="thermo",
        parameters={
            "fluid": {"type": "str|list", "required": True, "description": "流体英文名或组分名列表"},
            "z": {"type": "list", "required": False, "description": "混合物摩尔分率列表"},
        },
        tags=["molecular_weight", "MW", "thermo", "物性"],
        aliases=["fluid_MW", "molecular_weight", "calc_MW"],
    ), get_fluid_MW)

    registry.register(ToolCapability(
        name="get_fluid_thermal_conductivity",
        description="查询流体导热系数 (W/m/K)。",
        category="thermo",
        parameters={
            "fluid": {"type": "str|list", "required": True, "description": "流体英文名或组分名列表"},
            "T": {"type": "number", "required": True, "description": "温度 (K)"},
            "P": {"type": "number", "required": False, "description": "压力 (Pa)，默认 101325"},
            "phase": {"type": "string", "required": False, "description": "相态: 'liquid'(默认) 或 'gas'"},
            "z": {"type": "list", "required": False, "description": "混合物摩尔分率列表"},
        },
        tags=["thermal_conductivity", "thermo", "物性"],
        aliases=["fluid_thermal_conductivity", "thermal_conductivity_calc"],
    ), get_fluid_thermal_conductivity)

    registry.register(ToolCapability(
        name="get_fluid_surface_tension",
        description="查询液体表面张力 (N/m)。纯流体直接查库，混合物用 Macleod-Sugden Parachor 法。",
        category="thermo",
        parameters={
            "fluid": {"type": "str|list", "required": True, "description": "流体英文名或组分名列表"},
            "T": {"type": "number", "required": True, "description": "温度 (K)"},
            "z": {"type": "list", "required": False, "description": "混合物摩尔分率列表"},
        },
        tags=["surface_tension", "thermo", "物性"],
        aliases=["fluid_surface_tension", "surface_tension_calc"],
    ), get_fluid_surface_tension)

    # ================================================================
    # 管道计算 (category: pipe)
    # ================================================================
    from physics_engine import (
        calculate_pipe_pressure_drop,
        select_pipe_diameter,
    )

    registry.register(ToolCapability(
        name="calculate_pipe_pressure_drop",
        description="单相直管段压降计算 (Darcy-Weisbach)。支持自动选径(D不填时自动选DN)、局部阻力(弯头/阀门个数)、物性自动查询。",
        category="pipe",
        parameters={
            "Q": {"type": "number", "required": True, "description": "体积流量 (m³/s)"},
            "L": {"type": "number", "required": True, "description": "管道长度 (m)"},
            "D": {"type": "number", "required": False, "description": "管道内径 (m)，不填时自动选型"},
            "rho": {"type": "number", "required": False, "description": "流体密度 (kg/m³)，给定后跳过物性查询"},
            "mu": {"type": "number", "required": False, "description": "动力粘度 (Pa·s)，给定后跳过物性查询"},
            "roughness": {"type": "number", "required": False, "description": "管壁粗糙度 (m)，默认 4.5e-5（碳钢）"},
            "elevation_change": {"type": "number", "required": False, "description": "高程变化 (m)，默认 0"},
            # ── 局部阻力 ──
            "K_fittings": {"type": "number", "required": False, "description": "局部阻力系数之和 (无量纲)"},
            "n_elbow_90": {"type": "int", "required": False, "description": "90°弯头个数 (K=0.3/个)"},
            "n_elbow_45": {"type": "int", "required": False, "description": "45°弯头个数 (K=0.15/个)"},
            "n_gate_valve": {"type": "int", "required": False, "description": "闸阀个数 (K=0.15/个)"},
            "n_globe_valve": {"type": "int", "required": False, "description": "截止阀个数 (K=6.0/个)"},
            "n_check_valve": {"type": "int", "required": False, "description": "止回阀个数 (K=2.0/个)"},
            # ── 可选流体信息（用于自动查询物性） ──
            "fluid": {"type": "string", "required": False, "description": "(选填) 流体名称，当 rho/mu 未提供时用于自动查询"},
            "T": {"type": "number", "required": False, "description": "(选填) 温度 K，当 rho/mu 未提供时用于自动查询"},
            "P": {"type": "number", "required": False, "description": "(选填) 压力 Pa，默认 101325"},
            "phase": {"type": "string", "required": False, "description": "(选填) 相态: 'liquid'(默认)/'gas'"},
            "z": {"type": "list", "required": False, "description": "(选填) 混合物摩尔分率列表"},
            # ── 自动选径参数 ──
            "fluid_type": {"type": "string", "required": False, "description": "自动选径时的流体类型: 'liquid'(默认)/'gas'等"},
        },
        tags=["pressure_drop", "pipe", "管道", "压降校核", "自动选径"],
        aliases=[
            "pipe_pressure_drop_calc",
            "pipe_pressure_drop_calculation",
            "pressure_drop_verification",
            "pressure_drop_calc",
            "calc_pipe_pressure_drop",
        ],
    ), calculate_pipe_pressure_drop)

    registry.register(ToolCapability(
        name="select_pipe_diameter",
        description="根据流量选择标准管径（符合工艺标准流速限制的公称口径匹配）。返回 DN、外径、壁厚、实际流速等。",
        category="pipe",
        parameters={
            "Q": {"type": "number", "required": True, "description": "体积流量 (m³/s)"},
            "fluid_type": {"type": "string", "required": False, "description": "流体类型: 'liquid'(默认)/'gas'/'steam'/'slurry' 等"},
            "max_velocity": {"type": "number", "required": False, "description": "最大允许流速 (m/s)，覆盖默认推荐值"},
            "min_velocity": {"type": "number", "required": False, "description": "最小允许流速 (m/s)，覆盖默认推荐值"},
        },
        tags=["pipe_sizing", "DN", "管道"],
        aliases=[
            "pipe_diameter_selection",
            "pipe_sizing",
            "select_pipe_dn",
            "calc_pipe_diameter",
        ],
    ), select_pipe_diameter)

    # ================================================================
    # 泵计算 (category: pump)
    # ================================================================
    from physics_engine import (
        size_centrifugal_pump_power,
        calculate_pump_npsh,
        select_pump_type,
        calculate_system_head,
    )

    registry.register(ToolCapability(
        name="calculate_system_head",
        description="系统总扬程计算。H_total = Δz + ΔP/(ρg) + h_friction。支持物性自动查询。",
        category="pump",
        parameters={
            "delta_z": {"type": "number", "required": False, "description": "排出液面与吸入液面高差 (m)，默认 0"},
            "delta_P": {"type": "number", "required": False, "description": "排出与吸入空间压差 (Pa)，默认 0"},
            "h_friction": {"type": "number", "required": False, "description": "管路总阻力损失 (m液柱)，默认 0"},
            "rho": {"type": "number", "required": False, "description": "流体密度 (kg/m³)，给定后跳过查询"},
            "fluid": {"type": "string", "required": False, "description": "(选填) 流体名称，用于自动查询密度"},
            "T": {"type": "number", "required": False, "description": "(选填) 温度 K"},
            "P": {"type": "number", "required": False, "description": "(选填) 压力 Pa，默认 101325"},
        },
        tags=["system_head", "扬程", "泵"],
        aliases=["calc_system_head", "total_head", "system_head_calc"],
    ), calculate_system_head)

    registry.register(ToolCapability(
        name="size_centrifugal_pump_power",
        description="离心泵功率估算。由流量、扬程计算水力功率、轴功率、电机功率。支持物性自动查询。",
        category="pump",
        parameters={
            "Q": {"type": "number", "required": True, "description": "体积流量 (m³/s)"},
            "H": {"type": "number", "required": True, "description": "扬程 (m)"},
            "rho": {"type": "number", "required": False, "description": "流体密度 (kg/m³)，默认 1000"},
            "delta_P": {"type": "number", "required": False, "description": "压差 (Pa)，可选，若提供则自动换算 H"},
            "efficiency": {"type": "number", "required": False, "description": "泵效率 (0~1)，不指定则自动估算"},
            "fluid": {"type": "string", "required": False, "description": "(选填) 流体名称，用于自动查询密度"},
            "T": {"type": "number", "required": False, "description": "(选填) 温度 K"},
        },
        tags=["pump_power", "功率", "泵"],
        aliases=["pump_power_calc", "calc_pump_power", "centrifugal_pump_power"],
    ), size_centrifugal_pump_power)

    registry.register(ToolCapability(
        name="calculate_pump_npsh",
        description="泵 NPSH 汽蚀校核。计算可用 NPSHa，对比必需 NPSHr，判断汽蚀风险。支持物性自动查询。",
        category="pump",
        parameters={
            "P_suction": {"type": "number", "required": True, "description": "泵入口绝对压力 (Pa)"},
            "P_vapor": {"type": "number", "required": False, "description": "液体饱和蒸气压 (Pa)，给定后跳过查询"},
            "rho": {"type": "number", "required": False, "description": "液体密度 (kg/m³)，给定后跳过查询"},
            "suction_head": {"type": "number", "required": False, "description": "吸入高度 (m)，正=倒灌，默认 0"},
            "friction_loss": {"type": "number", "required": False, "description": "吸入管路摩擦损失 (m)，默认 0"},
            "NPSH_required": {"type": "number", "required": False, "description": "必需汽蚀余量 (m)，由泵样本给出"},
            "Q": {"type": "number", "required": False, "description": "流量 (m³/s)，用于估算 NPSHr"},
            "n": {"type": "number", "required": False, "description": "转速 (rpm)，默认 2900"},
            "margin": {"type": "number", "required": False, "description": "安全余量 (m)，默认 0.5"},
            "fluid": {"type": "string", "required": False, "description": "(选填) 流体名称，用于自动查询"},
            "T": {"type": "number", "required": False, "description": "(选填) 温度 K"},
        },
        tags=["NPSH", "cavitation", "汽蚀", "泵"],
        aliases=["pump_npsh_calc", "npsh_check", "npsh_verification"],
    ), calculate_pump_npsh)

    registry.register(ToolCapability(
        name="select_pump_type",
        description="按比转速 n_s 推荐泵型（离心/混流/轴流等）。",
        category="pump",
        parameters={
            "Q": {"type": "number", "required": True, "description": "流量 (m³/s)"},
            "H": {"type": "number", "required": True, "description": "扬程 (m)"},
            "n": {"type": "number", "required": False, "description": "转速 (rpm)，默认 2900"},
        },
        tags=["pump_type", "比转速", "泵"],
        aliases=["pump_type_selection", "select_pump", "pump_selection"],
    ), select_pump_type)

    # ================================================================
    # 阀门计算 (category: valve)
    # ================================================================
    from physics_engine import size_control_valve_liquid, size_safety_valve_gas

    registry.register(ToolCapability(
        name="size_control_valve_liquid",
        description="液体调节阀 Cv 值计算与口径选择 (IEC 60534)。含空化/闪蒸判断。",
        category="valve",
        parameters={
            "Q": {"type": "number", "required": True, "description": "体积流量 (m³/h)"},
            "rho": {"type": "number", "required": False, "description": "液体密度 (kg/m³)，给定后跳过物性查询"},
            "P1": {"type": "number", "required": False, "description": "阀前压力 (Pa)"},
            "P2": {"type": "number", "required": False, "description": "阀后压力 (Pa)"},
            "Pv": {"type": "number", "required": False, "description": "液体饱和蒸气压 (Pa)，给定后跳过物性查询"},
            # ── 可选流体信息（用于自动查询物性） ──
            "fluid": {"type": "string", "required": False, "description": "(选填) 流体名称，当 rho/Pv 未提供时用于自动查询"},
            "T": {"type": "number", "required": False, "description": "(选填) 温度 K，当 rho/Pv 未提供时用于自动查询"},
            "P": {"type": "number", "required": False, "description": "(选填) 压力 Pa，默认 101325"},
        },
        tags=["Cv", "valve", "调节阀"],
        aliases=["control_valve_sizing", "valve_cv_calc", "size_control_valve"],
    ), size_control_valve_liquid)

    registry.register(ToolCapability(
        name="size_safety_valve_gas",
        description="气体安全阀泄放量计算与口径选择 (API 526)。",
        category="valve",
        parameters={
            "W": {"type": "number", "required": True, "description": "泄放质量流量 (kg/s)"},
            "P_set": {"type": "number", "required": True, "description": "安全阀设定压力 (Pa)"},
            "P_back": {"type": "number", "required": True, "description": "背压 (Pa)"},
            "T": {"type": "number", "required": True, "description": "泄放温度 (K)"},
            "MW": {"type": "number", "required": True, "description": "气体分子量 (g/mol)"},
        },
        tags=["safety_valve", "泄放", "阀"],
        aliases=["safety_valve_sizing", "psv_sizing", "relief_valve_calc"],
    ), size_safety_valve_gas)

    # ================================================================
    # 仪表计算 (category: instrument)
    # ================================================================
    from physics_engine import calculate_orifice_flow, orifice_sizing

    registry.register(ToolCapability(
        name="calculate_orifice_flow",
        description="孔板流量计流量计算（ISO 5167）。由差压反算流量，支持物性自动查询。",
        category="instrument",
        parameters={
            "dP": {"type": "number", "required": True, "description": "孔板差压 (Pa)"},
            "D_pipe": {"type": "number", "required": True, "description": "管道内径 (m)"},
            "d_orifice": {"type": "number", "required": False, "description": "孔板开孔直径 (m)，可用 beta 替代"},
            "rho": {"type": "number", "required": False, "description": "流体密度 (kg/m³)，给定后跳过查询"},
            "mu": {"type": "number", "required": False, "description": "动力粘度 (Pa·s)，给定后跳过查询"},
            "beta": {"type": "number", "required": False, "description": "直径比 d/D，可替代 d_orifice"},
            "fluid": {"type": "string", "required": False, "description": "(选填) 流体名称，用于自动查询"},
            "T": {"type": "number", "required": False, "description": "(选填) 温度 K"},
            "P": {"type": "number", "required": False, "description": "(选填) 压力 Pa，默认 101325"},
            "phase": {"type": "string", "required": False, "description": "(选填) 相态: 'liquid'(默认)/'gas'"},
        },
        tags=["orifice", "流量计", "仪表"],
        aliases=["orifice_flow_calc", "calc_orifice_flow"],
    ), calculate_orifice_flow)

    registry.register(ToolCapability(
        name="orifice_sizing",
        description="孔板孔径设计（由目标流量反推孔径和直径比β，基于ISO 5167二分法精确定位）。支持物性自动查询。",
        category="instrument",
        parameters={
            "Q_target": {"type": "number", "required": True, "description": "目标体积流量 (m³/s)"},
            "D_pipe": {"type": "number", "required": True, "description": "管道内径 (m)"},
            "dP_max": {"type": "number", "required": True, "description": "最大允许差压 (Pa)"},
            "rho": {"type": "number", "required": False, "description": "流体密度 (kg/m³)，给定后跳过查询"},
            "mu": {"type": "number", "required": False, "description": "动力粘度 (Pa·s)，给定后跳过查询"},
            "fluid": {"type": "string", "required": False, "description": "(选填) 流体名称，用于自动查询"},
            "T": {"type": "number", "required": False, "description": "(选填) 温度 K"},
            "P": {"type": "number", "required": False, "description": "(选填) 压力 Pa，默认 101325"},
            "phase": {"type": "string", "required": False, "description": "(选填) 相态: 'liquid'(默认)/'gas'"},
        },
        tags=["orifice", "孔板选型", "仪表"],
        aliases=["orifice_plate_sizing", "orifice_design", "calc_orifice_diameter"],
    ), orifice_sizing)

    # ================================================================
    # 填料床压降 (category: packed_bed)
    # ================================================================
    from physics_engine import calculate_packed_bed_pressure_drop, select_packed_tower_diameter

    registry.register(ToolCapability(
        name="calculate_packed_bed_pressure_drop",
        description=(
            "填料床/固定床压降计算（Ergun方程修正版），支持气液两相湿压降。\n"
            "【推荐用法】指定 packing_name 即可自动获取填料参数(Fp/ε/dp)，无需手动填 d_particle。\n"
            "【自动塔径】D_bed 未提供时自动基于泛点气速法估算（需提供 Q_liq 和 rho_liq）。\n"
            "【湿压降】提供 Q_liq 时自动计算持液量和 Leva 修正。\n"
            "【注意】Q 必须为数值(m³/s)，如 2000 Nm³/h 应转换为 2000/3600=0.556 m³/s。\n"
            "【填料名称】可用: step_cascade_ring_50, pall_ring_50, pall_ring_38, pall_ring_25, "
            "intalox_saddle_50, raschig_ring_50 等。"
        ),
        category="packed_bed",
        parameters={
            "Q": {"type": "number", "required": True, "description": "气相体积流量 (m³/s)，必须为数值。如2000Nm³/h需转换为2000/3600"},
            "D_bed": {"type": "number", "required": False, "description": "床层直径 (m)，不填时自动基于泛点气速估算（需提供Q_liq和rho_liq）"},
            "L_bed": {"type": "number", "required": False, "description": "床层/填料层高度 (m)，默认3.0"},
            "rho": {"type": "number", "required": False, "description": "气相密度 (kg/m³)，给定后跳过查询"},
            "mu": {"type": "number", "required": False, "description": "气相动力粘度 (Pa·s)，给定后跳过查询"},
            "packing_name": {"type": "string", "required": False, "description": "【推荐】填料名称，自动查表获取Fp/ε/dp。50mm阶梯环填'step_cascade_ring_50'，50mm鲍尔环填'pall_ring_50'"},
            "d_particle": {"type": "number", "required": False, "description": "填料公称直径 (m)，如50mm填料填0.05。已用packing_name时无需填此项"},
            "packing_Fp": {"type": "number", "required": False, "description": "填料因子 (m⁻¹)，如50mmPP阶梯环Fp=98.4。已用packing_name时无需填此项"},
            "epsilon": {"type": "number", "required": False, "description": "空隙率，散装填料默认0.90。已用packing_name时无需填此项"},
            "packing_a": {"type": "number", "required": False, "description": "填料比表面积 (m²/m³)。已用packing_name时无需填此项"},
            "Q_liq": {"type": "number", "required": False, "description": "液相体积流量 (m³/s)，提供时计算气液两相湿压降并自动估算塔径"},
            "rho_liq": {"type": "number", "required": False, "description": "液相密度 (kg/m³)"},
            "mu_liq": {"type": "number", "required": False, "description": "液相动力粘度 (Pa·s)，默认0.001"},
            "flooding_factor": {"type": "number", "required": False, "description": "泛点因子（操作气速/泛点气速），默认0.7"},
            "fluid": {"type": "string", "required": False, "description": "(选填) 气相流体名称，用于自动查询物性，如'air'"},
            "T": {"type": "number", "required": False, "description": "(选填) 温度 K，与fluid配合自动查询物性"},
            "P": {"type": "number", "required": False, "description": "(选填) 压力 Pa，默认 101325"},
            "phase": {"type": "string", "required": False, "description": "(选填) 相态: 'gas'(默认)/'liquid'"},
        },
        tags=["packed_bed", "pressure_drop", "填料", "床层", "吸收塔"],
        aliases=["packed_bed_dP", "packed_bed_calc", "bed_pressure_drop"],
    ), calculate_packed_bed_pressure_drop)

    registry.register(ToolCapability(
        name="select_packed_tower_diameter",
        description=(
            "填料塔直径选择（基于GPDC泛点气速法）。"
            "使用Sherwood-Fair关联式估算泛点气速，按设计因子确定操作气速，"
            "计算所需塔径并匹配标准系列。"
            "通常先调用此工具确定塔径，再调用calculate_packed_bed_pressure_drop计算压降。"
        ),
        category="packed_bed",
        parameters={
            "Q_gas": {"type": "number", "required": True, "description": "气相体积流量 (m³/s)"},
            "Q_liquid": {"type": "number", "required": True, "description": "液相体积流量 (m³/s)"},
            "rho_gas": {"type": "number", "required": True, "description": "气相密度 (kg/m³)"},
            "rho_liquid": {"type": "number", "required": True, "description": "液相密度 (kg/m³)"},
            "mu_liquid": {"type": "number", "required": False, "description": "液相黏度 (Pa·s)，默认0.001"},
            "packing_name": {"type": "string", "required": False, "description": "填料名称，默认'pall_ring_50'"},
            "flooding_factor": {"type": "number", "required": False, "description": "泛点因子，设计操作气速/泛点气速，默认0.7"},
            "F_packing": {"type": "number", "required": False, "description": "填料因子 Fp (m⁻¹)，若提供则覆盖 packing_name"},
        },
        tags=["packed_bed", "tower_diameter", "flooding", "填料塔"],
        aliases=["packed_tower_diameter", "tower_sizing"],
    ), select_packed_tower_diameter)

    # ================================================================
    # 储罐辅助计算 (category: storage_tank)
    # ================================================================
    from physics_engine import optimize_tank_aspect_ratio, calculate_liquid_level_residence_time

    registry.register(ToolCapability(
        name="optimize_tank_aspect_ratio",
        description="储罐径高比优化（最小表面积设计）。返回最优 D/H 及对应尺寸。",
        category="storage_tank",
        parameters={
            "volume_required": {"type": "number", "required": True, "description": "所需有效容积 (m³)"},
            "tank_type": {"type": "string", "required": False, "description": "储罐类型: 'vertical'(默认)/'horizontal'"},
            "head_type": {"type": "string", "required": False, "description": "封头类型: 'ellipsoidal'(默认)/'torispherical'/'hemispherical'"},
        },
        tags=["aspect_ratio", "optimization", "储罐"],
        aliases=["tank_aspect_ratio", "optimize_tank", "tank_H_D_opt"],
    ), optimize_tank_aspect_ratio)

    registry.register(ToolCapability(
        name="calculate_liquid_level_residence_time",
        description="储罐液位停留时间计算。校验给定流量下液体在罐内的停留时间是否满足工艺要求。",
        category="storage_tank",
        parameters={
            "volume": {"type": "number", "required": True, "description": "罐体有效容积 (m³)"},
            "Q_in": {"type": "number", "required": True, "description": "入口流量 (m³/s)"},
            "Q_out": {"type": "number", "required": True, "description": "出口流量 (m³/s)"},
            "tank_type": {"type": "string", "required": False, "description": "储罐类型: 'vertical'(默认)/'horizontal'"},
            "D": {"type": "number", "required": False, "description": "罐径 (m)，不指定则由优化得出"},
        },
        tags=["residence_time", "停留时间", "储罐"],
        aliases=["residence_time_calc", "tank_residence_time"],
    ), calculate_liquid_level_residence_time)

    # ================================================================
    # 分离器计算 (category: separator)
    # ================================================================
    from physics_engine import size_gas_liquid_separator

    registry.register(ToolCapability(
        name="size_gas_liquid_separator",
        description="气液分离器尺寸设计/校核（Souders-Brown / Watkins 关联式）。可计算最小直径，也可校核已有直径的夹带风险。",
        category="separator",
        parameters={
            "Q_gas": {"type": "number", "required": True, "description": "气相体积流量 (m³/s)"},
            "Q_liquid": {"type": "number", "required": True, "description": "液相体积流量 (m³/s)，用户未提供时可假设极小值如 0.001"},
            "rho_gas": {"type": "number", "required": True, "description": "气相密度 (kg/m³)"},
            "rho_liquid": {"type": "number", "required": True, "description": "液相密度 (kg/m³)"},
            "mu_gas": {"type": "number", "required": False, "description": "气相粘度 (Pa·s)，默认 2e-5"},
            "D_existing": {"type": "number", "required": False, "description": "已有分离器直径 (m)，提供时自动校核夹带风险并输出速度比"},
            "method": {"type": "string", "required": False, "description": "计算方法: 'watkins'(默认，即Souders-Brown法)/'residence_time'(停留时间法)。**不建议传此参数，使用默认值即可**"},
            "residence_time_liquid": {"type": "number", "required": False, "description": "液相停留时间 (s)，默认 300"},
        },
        tags=["separator", "气液分离", "Souders-Brown"],
        aliases=["gas_liquid_separator_sizing", "separator_sizing"],
    ), size_gas_liquid_separator)

    # ================================================================
    # 精馏塔板校核 (category: distillation)
    # ================================================================
    from physics_engine import check_sieve_tray

    registry.register(ToolCapability(
        name="check_sieve_tray",
        description="筛板精馏塔综合校核：液泛、雾沫夹带、压降、漏液四项校核。",
        category="distillation",
        parameters={
            "D_column": {"type": "number", "required": True, "description": "塔径 (m)"},
            "tray_spacing": {"type": "number", "required": True, "description": "板间距 (m)"},
            "hole_diameter": {"type": "number", "required": True, "description": "孔径 (m)"},
            "hole_area_fraction": {"type": "number", "required": True, "description": "开孔率 (-)"},
            "rho_vapor": {"type": "number", "required": True, "description": "气相密度 (kg/m³)"},
            "rho_liquid": {"type": "number", "required": True, "description": "液相密度 (kg/m³)"},
            "Q_vapor": {"type": "number", "required": True, "description": "气相体积流量 (m³/s)"},
            "sigma": {"type": "number", "required": False, "description": "表面张力 (N/m)，可选"},
            "mu_liquid": {"type": "number", "required": False, "description": "液相黏度 (Pa·s)，可选"},
            "L_mass": {"type": "number", "required": False, "description": "液相质量流量 (kg/s)，可选"},
            "weir_height": {"type": "number", "required": False, "description": "堰高 (m)，默认 0.050"},
        },
        tags=["sieve_tray", "液泛", "雾沫夹带", "精馏"],
        aliases=["sieve_tray_check", "tray_flooding_check", "check_tray"],
    ), check_sieve_tray)

    # ================================================================
    # 换热器壳侧压降 (category: heat_exchanger)
    # ================================================================
    from physics_engine import calculate_shell_side_pressure_drop

    registry.register(ToolCapability(
        name="calculate_shell_side_pressure_drop",
        description="管壳式换热器壳侧压降计算（Bell-Delaware 法简化版）。",
        category="heat_exchanger",
        parameters={
            "m": {"type": "number", "required": True, "description": "壳程质量流量 (kg/s)"},
            "rho": {"type": "number", "required": True, "description": "壳程流体密度 (kg/m³)"},
            "mu": {"type": "number", "required": True, "description": "壳程流体粘度 (Pa·s)"},
            "D_shell": {"type": "number", "required": True, "description": "壳体内径 (m)"},
            "d_o": {"type": "number", "required": True, "description": "管外径 (m)"},
            "pitch": {"type": "number", "required": True, "description": "管间距 (m)"},
            "layout": {"type": "string", "required": False, "description": "管束排列: 'triangular'(默认)/'square'"},
            "baffle_spacing": {"type": "number", "required": False, "description": "折流板间距 (m)"},
        },
        tags=["shell_side", "pressure_drop", "换热器"],
        aliases=["shell_side_dP", "shell_pressure_drop", "hex_shell_dP"],
    ), calculate_shell_side_pressure_drop)

    return registry
