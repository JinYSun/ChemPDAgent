"""physics_engine — 化工物理计算引擎 (thermo + CoolProp 双引擎版)

内嵌了原 physics_toolbox 的全部计算逻辑，无外部路径依赖。
双引擎: thermo (Chemical, Mixture, PRMIX, SRKMIX, UNIFAC) + CoolProp (124 种纯流体 NIST 级精度)
"""

from .thermo_helper import (
    HAS_COOLPROP,
    clear_chem_cache,
    # --- 临界常数 & 相态判断 ---
    get_critical_constants,
    get_molecular_weight,
    classify_regime,
    get_reduced_properties,
    # --- 蒸气压 & 升华压 ---
    calc_vapor_pressure,
    calc_sublimation_pressure,
    # --- 密度 / 摩尔体积 ---
    calc_liquid_density,
    calc_gas_density,
    calc_solid_density,
    calc_molar_volume,
    # --- 压缩因子 / 逸度 ---
    calc_compressibility_factor,
    calc_fugacity,
    calc_fugacity_coefficient,
    # --- 热容 ---
    calc_heat_capacity,
    calc_Cv,
    calc_isentropic_exponent,
    # --- 焓 / 熵 / 内能 / Gibbs / Helmholtz ---
    calc_enthalpy,
    calc_entropy,
    calc_internal_energy,
    calc_gibbs_energy,
    calc_helmholtz_energy,
    get_formation_enthalpy,
    get_formation_gibbs,
    get_combustion_enthalpy,
    # --- 汽化/融化/升华焓 ---
    calc_enthalpy_vaporization,
    calc_enthalpy_fusion,
    calc_enthalpy_sublimation,
    # --- 粘度 ---
    calc_viscosity_liquid,
    calc_viscosity_gas,
    calc_kinematic_viscosity,
    # --- 导热系数 / 表面张力 ---
    calc_thermal_conductivity_liquid,
    calc_thermal_conductivity_gas,
    calc_surface_tension,
    # --- 推导物性 ---
    calc_thermal_diffusivity,
    calc_prandtl_number,
    calc_joule_thomson_coefficient,
    calc_isobaric_expansion,
    calc_isothermal_compressibility,
    calc_isentropic_compressibility,
    calc_speed_of_sound,
    calc_second_virial_coefficient,
    calc_parachor,
    calc_permittivity,
    # --- 扩散系数 ---
    calc_diffusion_gas,
    calc_diffusion_liquid_wilke_chang,
    calc_diffusion_liquid_hayduk_minhas,
    calc_diffusion_self,
    # --- 分子信息 ---
    get_triple_point,
    get_normal_boiling_point,
    get_normal_melting_point,
    get_dipole_moment,
    get_refractive_index,
    get_standard_entropy,
    get_critical_compressibility,
    get_acentric_factor,
    # --- 安全/环境 ---
    get_safety_properties,
    get_molecular_descriptors,
    # --- 相平衡 ---
    calc_henry_constant,
    calc_activity_coefficients_unifac,
    calc_bubble_point_T,
    calc_dew_point_T,
    calc_bubble_point_P,
    calc_dew_point_P,
    calc_flash_rachford_rice,
    get_mixture_k_values,
    # --- 综合属性 & 门面函数 ---
    get_chemical_properties,
    get_mixture_properties,
    get_fluid_density,
    get_fluid_viscosity,
    get_fluid_Cp,
    get_fluid_vapor_pressure,
    get_fluid_MW,
    get_fluid_thermal_conductivity,
    get_fluid_surface_tension,
)

from .common import (
    material_balance,
    energy_balance,
    clausius_clapeyron,
    antoine_equation,
    raoults_law,
    safe_float,
    safe_int,
)

# --- 管道计算 ---
from .pipe import (
    calculate_pipe_pressure_drop,
    calculate_two_phase_pressure_drop,
    calculate_pipe_fitting_pressure_drop,
    select_pipe_diameter,
)

# --- 填料床压降 ---
from .packed_bed import (
    calculate_packed_bed_pressure_drop,
    select_packed_tower_diameter,
    packed_bed_flooding_check,
    calculate_packed_bed_pressure_drop_bs,
)

# --- 泵计算 ---
from .pump import (
    size_centrifugal_pump_power,
    calculate_pump_npsh,
    calculate_pump_system_curve,
    match_pump_system_curve,
    select_pump_type,
    estimate_pump_efficiency,
    calculate_system_head,
)

# --- 阀门计算 ---
from .valve import (
    size_control_valve_liquid,
    size_safety_valve_gas,
)

# --- 仪表计算 ---
from .instrument import (
    calculate_orifice_flow,
    orifice_sizing,
)

# --- 储罐计算 ---
from .storage_tank import (
    optimize_tank_aspect_ratio,
    calculate_liquid_level_residence_time,
)

# --- 分离器计算 ---
from .flash_drum import (
    size_gas_liquid_separator,
    calculate_stokes_settling_efficiency,
)

# --- 精馏塔板校核 ---
from .distillation import (
    check_sieve_tray,
)

# --- 精馏塔板校核 ---
from .distillation import (
    check_sieve_tray,
)

# --- 换热器壳侧压降 ---
from .heatexchanger import (
    calculate_shell_side_pressure_drop,
)

__all__ = [
    "HAS_COOLPROP", "clear_chem_cache",
    # 临界常数 & 相态判断
    "get_critical_constants", "get_molecular_weight",
    "classify_regime", "get_reduced_properties",
    # 蒸气压 & 升华压
    "calc_vapor_pressure", "calc_sublimation_pressure",
    # 密度 / 摩尔体积
    "calc_liquid_density", "calc_gas_density", "calc_solid_density", "calc_molar_volume",
    # 压缩因子 / 逸度
    "calc_compressibility_factor", "calc_fugacity", "calc_fugacity_coefficient",
    # 热容
    "calc_heat_capacity", "calc_Cv", "calc_isentropic_exponent",
    # 焓 / 熵 / 内能 / Gibbs / Helmholtz
    "calc_enthalpy", "calc_entropy", "calc_internal_energy",
    "calc_gibbs_energy", "calc_helmholtz_energy",
    "get_formation_enthalpy", "get_formation_gibbs", "get_combustion_enthalpy",
    # 汽化/融化/升华焓
    "calc_enthalpy_vaporization", "calc_enthalpy_fusion", "calc_enthalpy_sublimation",
    # 粘度
    "calc_viscosity_liquid", "calc_viscosity_gas", "calc_kinematic_viscosity",
    # 导热系数 / 表面张力
    "calc_thermal_conductivity_liquid", "calc_thermal_conductivity_gas", "calc_surface_tension",
    # 推导物性
    "calc_thermal_diffusivity", "calc_prandtl_number",
    "calc_joule_thomson_coefficient", "calc_isobaric_expansion",
    "calc_isothermal_compressibility", "calc_isentropic_compressibility",
    "calc_speed_of_sound", "calc_second_virial_coefficient",
    "calc_parachor", "calc_permittivity",
    # 扩散系数
    "calc_diffusion_gas", "calc_diffusion_liquid_wilke_chang",
    "calc_diffusion_liquid_hayduk_minhas", "calc_diffusion_self",
    # 分子信息
    "get_triple_point", "get_normal_boiling_point", "get_normal_melting_point",
    "get_dipole_moment", "get_refractive_index", "get_standard_entropy",
    "get_critical_compressibility", "get_acentric_factor",
    # 安全/环境
    "get_safety_properties", "get_molecular_descriptors",
    # 相平衡
    "calc_henry_constant",
    "calc_activity_coefficients_unifac",
    "calc_bubble_point_T", "calc_dew_point_T",
    "calc_bubble_point_P", "calc_dew_point_P",
    "calc_flash_rachford_rice", "get_mixture_k_values",
    # 综合属性 & 门面函数
    "get_chemical_properties", "get_mixture_properties",
    "get_fluid_density", "get_fluid_viscosity", "get_fluid_Cp",
    "get_fluid_vapor_pressure", "get_fluid_MW",
    "get_fluid_thermal_conductivity", "get_fluid_surface_tension",
    # common
    "material_balance", "energy_balance",
    "clausius_clapeyron", "antoine_equation", "raoults_law",
    # 管道计算
    "calculate_pipe_pressure_drop",
    "calculate_two_phase_pressure_drop",
    "calculate_pipe_fitting_pressure_drop",
    "select_pipe_diameter",
    # 填料床压降
    "calculate_packed_bed_pressure_drop",
    "select_packed_tower_diameter",
    "packed_bed_flooding_check",
    "calculate_packed_bed_pressure_drop_bs",
    # 泵计算
    "size_centrifugal_pump_power",
    "calculate_pump_npsh",
    "calculate_pump_system_curve",
    "match_pump_system_curve",
    "select_pump_type",
    "estimate_pump_efficiency",
    "calculate_system_head",
    # 阀门计算
    "size_control_valve_liquid",
    "size_safety_valve_gas",
    # 仪表计算
    "calculate_orifice_flow",
    "orifice_sizing",
    # 储罐计算
    "optimize_tank_aspect_ratio",
    "calculate_liquid_level_residence_time",
    # 分离器计算
    "size_gas_liquid_separator",
    "calculate_stokes_settling_efficiency",
    # 精馏塔板校核
    "check_sieve_tray",
    # 换热器壳侧压降
    "calculate_shell_side_pressure_drop",
]
