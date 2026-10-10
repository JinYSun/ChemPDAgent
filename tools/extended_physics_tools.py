"""扩展物理计算工具注册模块

将 physics_engine 中所有物性计算和设备参数计算函数注册为 UniversalToolRegistry 工具。
分为:
  - 物性计算 (calc_XXX): 单流体 (T, P) 查询, 返回标量 → 包装为 {"result": value}
  - 设备参数计算 (pipe, pump, HEX, distillation, reactor, ...): 返回 dict
  - 相平衡计算 (bubble/dew/flash/K-values): 返回 dict
"""
from tools.registry import UniversalToolRegistry, ToolCapability


# ── 通用参数模板 ──
_FL = {"type": "string", "required": True, "description": "流体名称 (英文名), 如 'water'"}
_T  = {"type": "number", "required": True, "description": "温度 (K)"}
_P  = {"type": "number", "required": False, "description": "压力 (Pa), 默认 101325"}
_PH = {"type": "string", "required": False, "description": "相态: 'liquid'(默认) / 'gas'"}
_MOL = {"type": "boolean", "required": False, "description": "true=摩尔基 (默认), false=质量基"}

_TP  = {"name": _FL, "T": _T, "P": _P}
_TPH = {"name": _FL, "T": _T, "P": _P, "phase": _PH}
_TOnly = {"name": _FL, "T": _T}


def _scalar_wrapper(func, result_key="result"):
    """将返回标量的函数包装为返回 dict 的工具函数"""
    def wrapper(**kwargs):
        r = func(**kwargs)
        if isinstance(r, dict):
            return r
        return {result_key: r}
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def _reg(registry, name, func, desc, cat, params=None, tags=None, aliases=None, wrap=True):
    """注册一个物性工具 (标量自动包装, params 可省略则使用 **kwargs 模式)"""
    if params is None:
        params = {}
    registry.register(ToolCapability(
        name=name, description=desc, category=cat,
        parameters=params, tags=tags or [], aliases=aliases or [],
    ), _scalar_wrapper(func) if wrap else func)


def register_extended_physics_tools(registry: UniversalToolRegistry):
    """注册所有扩展物理计算工具到 UniversalToolRegistry"""
    from physics_engine import thermo_helper as th
    from physics_engine import (
        pipe, pump, heatexchanger as hx, distillation as dist,
        flash_drum as fd, reactor as rx, packed_bed as pb,
        storage_tank as st, common,
    )

    # ================================================================
    # ① 热力学基础物性 (category: thermo_property)
    # ================================================================
    _S = "thermo_property"
    _reg(registry, "calc_vapor_pressure", th.calc_vapor_pressure,
         "纯流体饱和蒸气压 Psat (Pa), CoolProp HEOS → thermo Wagner/Antoine",
         _S, _TOnly, ["Psat", "蒸气压"], ["vapor_pressure", "Psat"])
    _reg(registry, "calc_sublimation_pressure", th.calc_sublimation_pressure,
         "升华压 Psub (Pa), T < 三相点, Clausius-Clapeyron", _S, _TOnly, ["升华压"])
    _reg(registry, "calc_liquid_density", th.calc_liquid_density,
         "液相质量密度 (kg/m³), CoolProp HEOS → thermo, 超临界降级",
         _S, _TP, ["density", "液相密度"], ["liquid_density", "rho_l"])
    _reg(registry, "calc_gas_density", th.calc_gas_density,
         "气相质量密度 (kg/m³), CoolProp HEOS → thermo → 理想气体, 超临界降级",
         _S, _TP, ["density", "气相密度"], ["gas_density", "rho_g"])
    _reg(registry, "calc_solid_density", th.calc_solid_density,
         "固相密度 (kg/m³)", _S, _TOnly, ["固相密度"])
    _reg(registry, "calc_molar_volume", th.calc_molar_volume,
         "摩尔体积 Vm (m³/mol), CoolProp → EOS", _S, _TPH, ["摩尔体积"])
    _reg(registry, "calc_compressibility_factor", th.calc_compressibility_factor,
         "压缩因子 Z = PV/(RT) (无量纲)", _S, _TPH, ["压缩因子", "Z"])
    _reg(registry, "calc_fugacity", th.calc_fugacity,
         "逸度 f (Pa), CoolProp → EOS", _S, _TPH, ["逸度", "fugacity"])
    _reg(registry, "calc_fugacity_coefficient", th.calc_fugacity_coefficient,
         "逸度系数 φ = f/P (无量纲)", _S, _TPH, ["逸度系数"])
    _reg(registry, "calc_second_virial_coefficient", th.calc_second_virial_coefficient,
         "第二维里系数 B (m³/mol), Tsonopoulos 关联式", _S, _TOnly, ["virial"])
    _reg(registry, "get_critical_constants", th.get_critical_constants,
         "临界常数 (Tc K, Pc Pa, Vc m³/mol, Zc, ω, MW g/mol)", _S, {"name": _FL},
         ["临界常数", "Tc", "Pc", "omega"], wrap=True)
    _reg(registry, "classify_regime", th.classify_regime,
         "判断热力学区域 (ideal_gas/liquid/gas/supercritical)", _S, _TP, ["regime"])

    # ================================================================
    # ② 热力学能量性质 (category: thermo_energy)
    # ================================================================
    _E = "thermo_energy"
    _MOL_P = {"name": _FL, "T": _T, "P": _P, "molar": _MOL}
    _reg(registry, "calc_heat_capacity", th.calc_heat_capacity,
         "定压热容 Cp (J/mol/K 或 J/kg/K), CoolProp → DIPPR", _E, _TPH,
         ["Cp", "热容"], ["Cp", "heat_capacity"])
    _reg(registry, "calc_Cv", th.calc_Cv,
         "定容热容 Cv (J/mol/K), CoolProp → Cv=Cp-R", _E, _TPH, ["Cv"])
    _reg(registry, "calc_isentropic_exponent", th.calc_isentropic_exponent,
         "绝热指数 γ=Cp/Cv (无量纲)", _E, _TPH, ["gamma", "绝热指数"])
    _reg(registry, "calc_enthalpy", th.calc_enthalpy,
         "相对焓 H (J/mol 或 J/kg), CoolProp Hmolar → thermo", _E, _MOL_P, ["enthalpy", "焓"])
    _reg(registry, "calc_entropy", th.calc_entropy,
         "相对熵 S (J/mol/K 或 J/kg/K)", _E, _MOL_P, ["entropy", "熵"])
    _reg(registry, "calc_internal_energy", th.calc_internal_energy,
         "内能 U (J/mol 或 J/kg), CoolProp → U=H-PV", _E, _MOL_P, ["内能"])
    _reg(registry, "calc_gibbs_energy", th.calc_gibbs_energy,
         "Gibbs 自由能 G=H-TS (J/mol 或 J/kg)", _E, _MOL_P, ["Gibbs", "自由能"])
    _reg(registry, "calc_helmholtz_energy", th.calc_helmholtz_energy,
         "Helmholtz 自由能 A=U-TS (J/mol 或 J/kg)", _E, _MOL_P, ["Helmholtz"])
    _reg(registry, "calc_enthalpy_vaporization", th.calc_enthalpy_vaporization,
         "汽化焓 ΔHvap (J/mol), CoolProp → Watson/DIPPR, T≥Tc→None",
         _E, _TOnly, ["汽化焓", "Hvap"], ["Hvap"])
    _reg(registry, "calc_enthalpy_fusion", th.calc_enthalpy_fusion,
         "融化焓 ΔHfus (J/mol)", _E, {"name": _FL}, ["融化焓"])
    _reg(registry, "calc_enthalpy_sublimation", th.calc_enthalpy_sublimation,
         "升华焓 ΔHsub (J/mol)", _E, {"name": _FL}, ["升华焓"])
    _reg(registry, "get_formation_enthalpy", th.get_formation_enthalpy,
         "标准生成焓 ΔHf° (J/mol)", _E,
         {"fluid": _FL, "T": _T, "P": _P, "phase": _PH}, ["生成焓"])
    _reg(registry, "get_formation_gibbs", th.get_formation_gibbs,
         "标准 Gibbs 生成能 ΔGf° (J/mol)", _E,
         {"name": _FL, "phase": _PH}, ["Gibbs 生成能"])
    _reg(registry, "get_combustion_enthalpy", th.get_combustion_enthalpy,
         "燃烧焓 ΔHc (J/mol), 负值=放热", _E,
         {"name": _FL, "phase": _PH, "higher": {"type": "boolean", "required": False,
             "description": "true=高位热值 (默认), false=低位"}}, ["燃烧焓"])

    # ================================================================
    # ③ 输运物性 (category: thermo_transport)
    # ================================================================
    _TR = "thermo_transport"
    _reg(registry, "calc_viscosity_liquid", th.calc_viscosity_liquid,
         "液相动力粘度 μ (Pa·s), CoolProp → DIPPR, 超临界降级",
         _TR, _TP, ["viscosity", "液相粘度"], ["mu_l", "viscosity_liquid"])
    _reg(registry, "calc_viscosity_gas", th.calc_viscosity_gas,
         "气相动力粘度 μ (Pa·s), CoolProp → Chapman-Enskog",
         _TR, _TP, ["viscosity", "气相粘度"], ["mu_g", "viscosity_gas"])
    _reg(registry, "calc_kinematic_viscosity", th.calc_kinematic_viscosity,
         "运动粘度 ν=μ/ρ (m²/s)", _TR, _TPH, ["运动粘度", "nu"])
    _reg(registry, "calc_thermal_conductivity_liquid", th.calc_thermal_conductivity_liquid,
         "液相导热系数 λ (W/m/K)", _TR, _TP, ["导热系数", "液相"],
         ["k_l", "lambda_l", "thermal_conductivity_liquid"])
    _reg(registry, "calc_thermal_conductivity_gas", th.calc_thermal_conductivity_gas,
         "气相导热系数 λ (W/m/K)", _TR, _TP, ["导热系数", "气相"],
         ["k_g", "lambda_g", "thermal_conductivity_gas"])
    _reg(registry, "calc_surface_tension", th.calc_surface_tension,
         "表面张力 σ (N/m), CoolProp → DIPPR, T≥Tc→0",
         _TR, _TOnly, ["表面张力", "sigma"], ["surface_tension"])
    _DIFF_PARAMS = {"solute": {"type": "string", "required": True, "description": "溶质名称"},
                    "solvent": {"type": "string", "required": True, "description": "溶剂名称"},
                    "T": _T}
    _reg(registry, "calc_diffusion_liquid_wilke_chang", th.calc_diffusion_liquid_wilke_chang,
         "液相无限稀释扩散系数 D°_AB (m²/s), Wilke-Chang 方程, 超临界→None",
         _TR, _DIFF_PARAMS, ["扩散系数", "Wilke-Chang"], ["D_AB_wilke"])
    _reg(registry, "calc_diffusion_liquid_hayduk_minhas", th.calc_diffusion_liquid_hayduk_minhas,
         "液相扩散系数 (m²/s), Hayduk-Minhas (水/非水自动切换), 超临界→None",
         _TR, _DIFF_PARAMS, ["扩散系数", "Hayduk-Minhas"])
    _reg(registry, "calc_diffusion_gas", th.calc_diffusion_gas,
         "气相二元扩散系数 D_AB (m²/s), Chapman-Enskog",
         _TR, {"solute": {"type": "string", "required": True, "description": "组分A"},
               "solvent": {"type": "string", "required": True, "description": "组分B"},
               "T": _T, "P": _P}, ["气相扩散系数"])
    _reg(registry, "calc_diffusion_self", th.calc_diffusion_self,
         "自扩散系数 (m²/s)", _TR, _TP, ["自扩散"])

    # ================================================================
    # ④ 推导物性 (category: thermo_derived)
    # ================================================================
    _D = "thermo_derived"
    _reg(registry, "calc_prandtl_number", th.calc_prandtl_number,
         "Prandtl 数 Pr=μ·Cp_mass/λ (无量纲)", _D, _TPH, ["Pr", "Prandtl"])
    _reg(registry, "calc_thermal_diffusivity", th.calc_thermal_diffusivity,
         "热扩散率 α=λ/(ρ·Cp_mass) (m²/s)", _D, _TPH, ["热扩散率", "alpha"])
    _reg(registry, "calc_isothermal_compressibility", th.calc_isothermal_compressibility,
         "等温压缩率 κ_T=-(1/V)(∂V/∂P)_T (1/Pa)", _D, _TPH, ["压缩率"])
    _reg(registry, "calc_isentropic_compressibility", th.calc_isentropic_compressibility,
         "等熵压缩率 κ_s=κ_T/γ (1/Pa)", _D, _TPH, ["等熵压缩率"])
    _reg(registry, "calc_isobaric_expansion", th.calc_isobaric_expansion,
         "等压膨胀系数 β=(1/V)(∂V/∂T)_P (1/K)", _D, _TPH, ["膨胀系数", "beta"])
    _reg(registry, "calc_joule_thomson_coefficient", th.calc_joule_thomson_coefficient,
         "Joule-Thomson 系数 μ_JT (K/Pa)", _D, _TPH, ["JT系数", "节流"])
    _reg(registry, "calc_speed_of_sound", th.calc_speed_of_sound,
         "声速 c (m/s), CoolProp → c=√(γ/(ρ·κ_s))", _D, _TPH, ["声速", "sound"])
    _reg(registry, "calc_parachor", th.calc_parachor,
         "Parachor (m⁵·N^(1/4)/mol)", _D, _TOnly, ["Parachor"])

    # ================================================================
    # ⑤ 特征温度 / 安全 / 分子描述符
    # ================================================================
    _reg(registry, "get_triple_point", th.get_triple_point,
         "三相点 (Tt K, Pt Pa)", _S, {"name": _FL}, ["三相点"])
    _reg(registry, "get_normal_boiling_point", th.get_normal_boiling_point,
         "标准沸点 Tb (K, 1 atm)", _S, {"name": _FL}, ["沸点", "Tb"])
    _reg(registry, "get_normal_melting_point", th.get_normal_melting_point,
         "标准熔点 Tm (K)", _S, {"name": _FL}, ["熔点", "Tm"])
    _reg(registry, "get_acentric_factor", th.get_acentric_factor,
         "偏心因子 ω (无量纲)", _S, {"name": _FL}, ["omega", "偏心因子"])
    _reg(registry, "get_dipole_moment", th.get_dipole_moment,
         "偶极矩 (Debye)", _S, {"name": _FL}, ["偶极矩"])
    _reg(registry, "get_safety_properties", th.get_safety_properties,
         "安全性质 (闪点, 自燃温度, LEL/UEL, GHS)", _S, {"name": _FL},
         ["安全", "闪点", "LEL"], wrap=True)
    _reg(registry, "get_molecular_descriptors", th.get_molecular_descriptors,
         "分子描述符 (SMILES, InChI, HBA/HBD, PSA)", _S, {"name": _FL},
         ["SMILES", "分子描述符"], wrap=True)
    _reg(registry, "get_reduced_properties", th.get_reduced_properties,
         "对比性质 (Tr=T/Tc, Pr=P/Pc)", _S, _TP, ["对比温度", "对比压力"])

    # ================================================================
    # ⑥ 相平衡 / 混合物 (category: phase_equilibrium)
    # ================================================================
    _VLE = "phase_equilibrium"
    _COMP = {"type": "list", "required": True, "description": "组分名称列表"}
    _Z = {"type": "list", "required": True, "description": "摩尔分率列表"}
    _reg(registry, "calc_henry_constant", th.calc_henry_constant,
         "Henry 常数 (Pa), 气体在水中的溶解度", _VLE,
         {"name": _FL, "T": {"type": "number", "required": False, "description": "温度 K, 默认 298.15"}},
         ["Henry", "溶解度"])
    _reg(registry, "calc_activity_coefficients_unifac", th.calc_activity_coefficients_unifac,
         "UNIFAC 活度系数 γ_i", _VLE,
         {"components": _COMP, "T": _T, "xs": _Z}, ["UNIFAC", "活度系数"], wrap=True)
    _reg(registry, "calc_bubble_point_T", th.calc_bubble_point_T,
         "泡点温度 (K), 给定 P 和液相组成", _VLE,
         {"components": _COMP, "x": _Z, "P": _P}, ["泡点温度"], wrap=True)
    _reg(registry, "calc_dew_point_T", th.calc_dew_point_T,
         "露点温度 (K), 给定 P 和气相组成", _VLE,
         {"components": _COMP, "y": _Z, "P": _P}, ["露点温度"], wrap=True)
    _reg(registry, "calc_bubble_point_P", th.calc_bubble_point_P,
         "泡点压力 (Pa), 给定 T 和液相组成", _VLE,
         {"components": _COMP, "x": _Z, "T": _T}, ["泡点压力"], wrap=True)
    _reg(registry, "calc_dew_point_P", th.calc_dew_point_P,
         "露点压力 (Pa), 给定 T 和气相组成", _VLE,
         {"components": _COMP, "y": _Z, "T": _T}, ["露点压力"], wrap=True)
    _reg(registry, "calc_flash_rachford_rice", th.calc_flash_rachford_rice,
         "Rachford-Rice 闪蒸计算 (气液分率和组成)", _VLE,
         {"components": _COMP, "z": _Z, "K": {"type": "list", "required": True,
             "description": "K值列表"}}, ["闪蒸", "Rachford-Rice"], wrap=True)
    _reg(registry, "get_mixture_k_values", th.get_mixture_k_values,
         "K 值计算 (auto/henry/eos_pr/unifac/wilson/nrtl)", _VLE,
         {"components": _COMP, "T": _T, "P": _P, "zs": _Z,
          "method": {"type": "string", "required": False, "description": "方法: auto(默认)"}},
         ["K值", "相平衡常数"], wrap=True)
    _reg(registry, "get_chemical_properties", th.get_chemical_properties,
         "组分综合物性查询 (T,P 下全部可查物性)", _VLE,
         {"component_name": _FL, "T": {"type": "number", "required": False, "description": "温度 K"},
          "P": {"type": "number", "required": False, "description": "压力 Pa"}},
         ["综合物性"], wrap=True)
    _reg(registry, "get_mixture_properties", th.get_mixture_properties,
         "混合物综合物性 (给定组成和 T,P)", _VLE,
         {"components": _COMP, "z": _Z, "T": _T, "P": _P},
         ["混合物物性"], wrap=True)

    # ================================================================
    # ⑦ 管道计算扩展 (category: pipe)
    # ================================================================
    _reg(registry, "churchill_friction", pipe.churchill_friction,
         "Churchill 摩擦因子 f (全 Re 范围通用)", "pipe",
         {"Re": {"type": "number", "required": True, "description": "雷诺数"},
          "roughness": {"type": "number", "required": False, "description": "绝对粗糙度 (m)"},
          "diameter": {"type": "number", "required": False, "description": "管道内径 (m)"}},
         ["摩擦因子", "Churchill"])
    _reg(registry, "calculate_two_phase_pressure_drop", pipe.calculate_two_phase_pressure_drop,
         "两相流管道压降 (Lockhart-Martinelli)", "pipe", wrap=False)
    _reg(registry, "calculate_pipe_fitting_pressure_drop", pipe.calculate_pipe_fitting_pressure_drop,
         "管件局部阻力压降", "pipe", wrap=False)

    # ================================================================
    # ⑧ 泵计算扩展 (category: pump)
    # ================================================================
    _reg(registry, "estimate_pump_efficiency", pump.estimate_pump_efficiency,
         "泵效率估算 (η)", "pump",
         {"Q": {"type": "number", "required": True, "description": "流量 (m³/s)"},
          "H": {"type": "number", "required": True, "description": "扬程 (m)"},
          "n": {"type": "number", "required": False, "description": "转速 rpm, 默认 2900"}},
         ["泵效率"])
    _reg(registry, "suggest_pump_type", pump.suggest_pump_type,
         "推荐泵型 (按流量/扬程/粘度)", "pump", wrap=False)
    _reg(registry, "npshr_estimate", pump.npshr_estimate,
         "必需汽蚀余量 NPSHr 估算 (m)", "pump",
         {"Q": {"type": "number", "required": True, "description": "流量 (m³/s)"},
          "n": {"type": "number", "required": True, "description": "转速 (rpm)"},
          "NPSHR_type": {"type": "string", "required": False,
              "description": "泵型: centrifugal(默认)"}}, ["NPSHr"])

    # ================================================================
    # ⑨ 换热器参数计算 (category: heat_exchanger)
    # ================================================================
    _HX = "heat_exchanger"
    _reg(registry, "lmtd", hx.lmtd,
         "对数平均温差 LMTD (K), 含逆流/并流/交叉流", _HX, wrap=False)
    _reg(registry, "lmtd_with_correction", hx.lmtd_with_correction,
         "LMTD + 修正因子 F (多壳程/交叉流)", _HX, wrap=False)
    _reg(registry, "heat_duty_from_flow", hx.heat_duty_from_flow,
         "热负荷 Q (W) 由流量和温差计算", _HX, wrap=False)
    _reg(registry, "heat_exchanger_area", hx.heat_exchanger_area,
         "换热面积 A (m²) 由 Q/U/LMTD 计算", _HX, wrap=False)
    _reg(registry, "estimate_U", hx.estimate_U,
         "总传热系数 U 典型范围估算", _HX, wrap=False)
    _reg(registry, "effectiveness_ntu", hx.effectiveness_ntu,
         "ε-NTU 法换热器计算", _HX, wrap=False)
    _reg(registry, "fouling_resistance", hx.fouling_resistance,
         "污垢热阻查询 (TEMA 标准)", _HX, wrap=False)
    _reg(registry, "select_exchanger_type", hx.select_exchanger_type,
         "换热器类型推荐 (管壳式/板式/...)", _HX, wrap=False)
    _reg(registry, "condensate_heat_transfer", hx.condensate_heat_transfer,
         "冷凝传热系数 (Nusselt 膜状冷凝)", _HX, wrap=False)
    _reg(registry, "boiling_heat_transfer", hx.boiling_heat_transfer,
         "沸腾传热系数 (Rohsenow 核态沸腾)", _HX, wrap=False)
    _reg(registry, "detect_phase_change", hx.detect_phase_change,
         "检测流体是否发生相变", _HX, wrap=False)

    # ================================================================
    # ⑩ 精馏塔参数计算 (category: distillation)
    # ================================================================
    _DST = "distillation"
    _reg(registry, "fenske_minimum_stages", dist.fenske_minimum_stages,
         "Fenske 最小理论板数 N_min", _DST, wrap=False)
    _reg(registry, "underwood_minimum_reflux", dist.underwood_minimum_reflux,
         "Underwood 最小回流比 R_min", _DST, wrap=False)
    _reg(registry, "gilliland_correlation", dist.gilliland_correlation,
         "Gilliland 关联式求实际理论板数 N", _DST, wrap=False)
    _reg(registry, "calculate_column_diameter", dist.calculate_column_diameter,
         "精馏塔塔径计算 (Souders-Brown flooding)", _DST, wrap=False)
    _reg(registry, "calculate_q_factor", dist.calculate_q_factor,
         "进料热状态 q 因子计算", _DST, wrap=False)
    _reg(registry, "estimate_tray_efficiency", dist.estimate_tray_efficiency,
         "塔板效率估算 (O'Connell)", _DST, wrap=False)
    _reg(registry, "calculate_tray_hydraulics", dist.calculate_tray_hydraulics,
         "塔板水力学计算 (液泛/夹带/压降/漏液)", _DST, wrap=False)
    _reg(registry, "bubble_point_temperature", dist.bubble_point_temperature,
         "泡点温度计算 (Antoine + 迭代)", _DST, wrap=False)
    _reg(registry, "dew_point_temperature", dist.dew_point_temperature,
         "露点温度计算 (Antoine + 迭代)", _DST, wrap=False)
    _reg(registry, "relative_volatility_from_pressure", dist.relative_volatility_from_pressure,
         "相对挥发度 α 由饱和压计算", _DST, wrap=False)
    _reg(registry, "feed_stage_estimate", dist.feed_stage_estimate,
         "进料板位置估算 (Kirkbride)", _DST, wrap=False)
    _reg(registry, "calculate_condenser_duty", dist.calculate_condenser_duty,
         "冷凝器热负荷 (W)", _DST, wrap=False)
    _reg(registry, "calculate_reboiler_duty", dist.calculate_reboiler_duty,
         "再沸器热负荷 (W)", _DST, wrap=False)

    # ================================================================
    # ⑪ 闪蒸罐参数计算 (category: flash_drum)
    # ================================================================
    _FL2 = "flash_drum"
    _reg(registry, "solve_rachford_rice", fd.solve_rachford_rice,
         "Rachford-Rice 方程求解 (气相分率 ψ)", _FL2, wrap=False)
    _reg(registry, "isothermal_flash", fd.isothermal_flash,
         "等温闪蒸计算", _FL2, wrap=False)
    _reg(registry, "flash_drum_diameter", fd.flash_drum_diameter,
         "闪蒸罐直径计算 (Souders-Brown)", _FL2, wrap=False)
    _reg(registry, "flash_drum_volume", fd.flash_drum_volume,
         "闪蒸罐容积计算", _FL2, wrap=False)
    _reg(registry, "entrainment_fraction", fd.entrainment_fraction,
         "雾沫夹带分率计算", _FL2, wrap=False)

    # ================================================================
    # ⑫ 反应器参数计算 (category: reactor)
    # ================================================================
    _RX = "reactor"
    _reg(registry, "reaction_rate", rx.reaction_rate,
         "反应速率 r (mol/m³/s), n 级反应动力学", _RX, wrap=False)
    _reg(registry, "reaction_enthalpy", rx.reaction_enthalpy,
         "反应焓 ΔH_rxn (J/mol), 由生成焓计算", _RX, wrap=False)
    _reg(registry, "adiabatic_temperature_rise", rx.adiabatic_temperature_rise,
         "绝热温升 ΔT_ad (K)", _RX, wrap=False)
    _reg(registry, "selectivity_analysis", rx.selectivity_analysis,
         "选择性分析 (平行/连串反应)", _RX, wrap=False)
    _reg(registry, "pfr_design", rx.pfr_design,
         "PFR 反应器设计 (体积/管数/床层长度)", _RX, wrap=False)
    _reg(registry, "cstr_design", rx.cstr_design,
         "CSTR 反应器设计 (体积/停留时间)", _RX, wrap=False)
    _reg(registry, "multitubular_design", rx.multitubular_design,
         "列管式反应器设计", _RX, wrap=False)
    _reg(registry, "fixed_bed_dimensions", rx.fixed_bed_dimensions,
         "固定床尺寸 (床高/截面积)", _RX, wrap=False)
    _reg(registry, "fixed_bed_pressure_drop_ergun", rx.fixed_bed_pressure_drop_ergun,
         "固定床 Ergun 压降 (Pa)", _RX, wrap=False)
    _reg(registry, "fixed_bed_catalyst_inventory", rx.fixed_bed_catalyst_inventory,
         "固定床催化剂装填量 (kg)", _RX, wrap=False)
    _reg(registry, "minimum_fluidization_velocity", rx.minimum_fluidization_velocity,
         "最小流化速度 u_mf (m/s)", _RX, wrap=False)
    _reg(registry, "terminal_velocity", rx.terminal_velocity,
         "终端速度 u_t (m/s)", _RX, wrap=False)
    _reg(registry, "check_fluidization_regime", rx.check_fluidization_regime,
         "流化状态判断 (散式/聚式/节涌)", _RX, wrap=False)
    _reg(registry, "estimate_catalyst_properties", rx.estimate_catalyst_properties,
         "催化剂典型物性估算", _RX, wrap=False)
    _reg(registry, "calculate_thiele_modulus", rx.calculate_thiele_modulus,
         "Thiele 模数 φ (内扩散影响)", _RX, wrap=False)
    _reg(registry, "effectiveness_factor", rx.effectiveness_factor,
         "内扩散有效因子 η", _RX, wrap=False)
    _reg(registry, "pore_diffusion_correction", rx.pore_diffusion_correction,
         "孔扩散修正 (Weisz-Prater)", _RX, wrap=False)
    _reg(registry, "catalyst_deactivation", rx.catalyst_deactivation,
         "催化剂失活动力学", _RX, wrap=False)
    _reg(registry, "check_conversion_feasibility", rx.check_conversion_feasibility,
         "转化率可达性约束检查 (化学计量/体积/速率/停留时间)", _RX, wrap=False)
    _reg(registry, "reactor_heat_removal_check", rx.reactor_heat_removal_check,
         "反应器移热能力校核", _RX, wrap=False)

    # ================================================================
    # ⑬ 储罐参数计算 (category: storage_tank)
    # ================================================================
    _ST = "storage_tank"
    _reg(registry, "tank_shell_thickness", st.tank_shell_thickness,
         "储罐壁厚计算 (API 650)", _ST, wrap=False)
    _reg(registry, "breathing_losses", st.breathing_losses,
         "储罐呼吸损耗 (API 2000)", _ST, wrap=False)
    _reg(registry, "insulation_thickness", st.insulation_thickness,
         "保温层厚度计算", _ST, wrap=False)
    _reg(registry, "tank_foundation_load", st.tank_foundation_load,
         "储罐基础载荷计算", _ST, wrap=False)
    _reg(registry, "tank_relief_valve_sizing", st.tank_relief_valve_sizing,
         "储罐泄压阀选型 (API 2000)", _ST, wrap=False)
    _reg(registry, "fire_case_relief_check", st.fire_case_relief_check,
         "火灾工况泄放量校核 (API 521)", _ST, wrap=False)

    # ================================================================
    # ⑭ 填料床扩展 (category: packed_bed)
    # ================================================================
    _reg(registry, "packed_bed_flooding_check", pb.packed_bed_flooding_check,
         "填料塔泛点校核", "packed_bed", wrap=False)

    # ================================================================
    # ⑮ 通用衡算 (category: engineering)
    # ================================================================
    _EN = "engineering"
    _reg(registry, "material_balance", common.material_balance,
         "物料衡算 (质量/摩尔)", _EN, wrap=False)
    _reg(registry, "energy_balance", common.energy_balance,
         "能量衡算 (焓变/热负荷)", _EN, wrap=False)
    _reg(registry, "clausius_clapeyron", common.clausius_clapeyron,
         "Clausius-Clapeyron 蒸气压方程", _EN, wrap=False)
    _reg(registry, "antoine_equation", common.antoine_equation,
         "Antoine 蒸气压方程", _EN, wrap=False)
    _reg(registry, "raoults_law", common.raoults_law,
         "Raoult 定律 (理想溶液 VLE)", _EN, wrap=False)
