"""设备级工具函数：每个设备作为一个独立工具，返回该设备的全部计算参数。

每个函数调用对应设备的函数模块中的所有 calc_* 函数，
将所有返回值合并为一个字典返回给 LLM。
"""

from typing import Dict, List, Any

# ─────────────────────────────────────────────────────────────────
# 导入各设备的函数模块
# ─────────────────────────────────────────────────────────────────
from pump.funcs import calc_outlet_pressure as pump_calc_outlet_pressure
from pressure_reducing_valve.funcs import calc_outlet_pressure as prv_calc_outlet_pressure
from flash_drum.funcs import calc_flash_drum
from heat_exchangers.funcs import (
    calc_heat_duty,
    calc_log_mean_temp_difference,
    calc_overall_heat_transfer_coefficient,
    calc_heat_transfer_area,
)
from stor.funcs import (
    calc_reaction_extent as stor_calc_reaction_extent,
    calc_outlet_molar_flows as stor_calc_outlet_molar_flows,
    calc_outlet_total_molar_flow as stor_calc_outlet_total_molar_flow,
    calc_outlet_volumetric_flow as stor_calc_outlet_volumetric_flow,
    calc_outlet_concentrations as stor_calc_outlet_concentrations,
)
from cstr.funcs import (
    calc_reaction_extent as cstr_calc_reaction_extent,
    calc_outlet_molar_flows as cstr_calc_outlet_molar_flows,
    calc_outlet_total_molar_flow as cstr_calc_outlet_total_molar_flow,
    calc_outlet_volumetric_flow as cstr_calc_outlet_volumetric_flow,
    calc_outlet_concentrations as cstr_calc_outlet_concentrations,
    calc_fluid_density as cstr_calc_fluid_density,
    calc_fluid_viscosity as cstr_calc_fluid_viscosity,
    calc_archimedes_number,
    calc_minimum_fluidization_velocity,
    calc_reaction_rate,
    calc_cstr_volume,
    calc_residence_time as cstr_calc_residence_time,
)
from pfr.funcs import (
    calc_reaction_extent as pfr_calc_reaction_extent,
    calc_outlet_molar_flows as pfr_calc_outlet_molar_flows,
    calc_outlet_total_molar_flow as pfr_calc_outlet_total_molar_flow,
    calc_outlet_volumetric_flow as pfr_calc_outlet_volumetric_flow,
    calc_outlet_concentrations as pfr_calc_outlet_concentrations,
    calc_inlet_total_molar_flow,
    calc_inlet_volumetric_flow_stp,
    calc_inlet_volumetric_flow_operating,
    calc_bed_volume,
    calc_catalyst_mass,
    calc_single_tube_diameter,
    calc_tube_number,
    calc_bed_length_single_tube,
    calc_bed_length_multitube,
)
from tray_distillation.funcs import (
    calc_mass_balance,
    calc_operating_conditions,
    calc_feed_thermal_condition,
    calc_min_reflux_ratio,
    calc_min_theoretical_stages,
    calc_theoretical_stages,
    calc_actual_stages,
    calc_feed_stage,
)
from mixer.funcs import (
    calc_outlet_molar_flows as mixer_calc_outlet_molar_flows,
    calc_outlet_total_molar_flow as mixer_calc_outlet_total_molar_flow,
    calc_outlet_composition as mixer_calc_outlet_composition,
    calc_outlet_pressure as mixer_calc_outlet_pressure,
    calc_outlet_temperature as mixer_calc_outlet_temperature,
    calc_outlet_volumetric_flow as mixer_calc_outlet_volumetric_flow,
    calc_mixer_state as mixer_calc_state,
)


# ═════════════════════════════════════════════════════════════════
# 设备 1：换热器
# ═════════════════════════════════════════════════════════════════
def run_heat_exchanger(
    process_fluid_molar_flows_mol_per_s: Dict[str, float],
    process_fluid_temp_in_K: float,
    process_fluid_temp_out_K: float,
    process_fluid_pressure_Pa: float,
    utility_fluid_temp_in_K: float,
    utility_fluid_temp_out_K: float,
    utility_type: str = "cooling_water",
) -> Dict[str, Any]:
    """运行换热器设备，返回热负荷、对数平均温差、总传热系数、换热面积等全部参数。"""
    result: Dict[str, Any] = {}

    duty = calc_heat_duty(
        process_fluid_molar_flows_mol_per_s=process_fluid_molar_flows_mol_per_s,
        process_fluid_temp_in_K=process_fluid_temp_in_K,
        process_fluid_temp_out_K=process_fluid_temp_out_K,
        process_fluid_pressure_Pa=process_fluid_pressure_Pa,
    )
    result.update(duty)

    lmtd = calc_log_mean_temp_difference(
        process_fluid_temp_in_K=process_fluid_temp_in_K,
        process_fluid_temp_out_K=process_fluid_temp_out_K,
        utility_fluid_temp_in_K=utility_fluid_temp_in_K,
        utility_fluid_temp_out_K=utility_fluid_temp_out_K,
    )
    result.update(lmtd)

    k_result = calc_overall_heat_transfer_coefficient(
        process_fluid_molar_flows_mol_per_s=process_fluid_molar_flows_mol_per_s,
        process_fluid_temp_in_K=process_fluid_temp_in_K,
        process_fluid_temp_out_K=process_fluid_temp_out_K,
        process_fluid_pressure_Pa=process_fluid_pressure_Pa,
        utility_type=utility_type,
    )
    result.update(k_result)

    area = calc_heat_transfer_area(
        process_fluid_molar_flows_mol_per_s=process_fluid_molar_flows_mol_per_s,
        process_fluid_temp_in_K=process_fluid_temp_in_K,
        process_fluid_temp_out_K=process_fluid_temp_out_K,
        process_fluid_pressure_Pa=process_fluid_pressure_Pa,
        utility_fluid_temp_in_K=utility_fluid_temp_in_K,
        utility_fluid_temp_out_K=utility_fluid_temp_out_K,
        utility_type=utility_type,
    )
    result.update(area)

    return result


# ═════════════════════════════════════════════════════════════════
# 设备 2：化学计量反应器
# ═════════════════════════════════════════════════════════════════
def run_stoichiometric_reactor(
    inlet_molar_flows_mol_per_s: Dict[str, float],
    stoichiometric_coefficients: Dict[str, float],
    key_component: str,
    key_component_conversion: float,
    reactor_temp_K: float,
    reactor_pressure_Pa: float,
) -> Dict[str, Any]:
    """运行化学计量反应器设备，返回反应进度、出口流量、出口总流量、出口体积流量、出口浓度。"""
    result: Dict[str, Any] = {}

    extent = stor_calc_reaction_extent(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    result.update(extent)

    outlet_flows = stor_calc_outlet_molar_flows(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    result["outlet_molar_flows_mol_per_s"] = outlet_flows

    total_flow = stor_calc_outlet_total_molar_flow(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    result.update(total_flow)

    vol_flow = stor_calc_outlet_volumetric_flow(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        reactor_temp_K=reactor_temp_K,
        reactor_pressure_Pa=reactor_pressure_Pa,
    )
    result.update(vol_flow)

    concentrations = stor_calc_outlet_concentrations(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        reactor_temp_K=reactor_temp_K,
        reactor_pressure_Pa=reactor_pressure_Pa,
    )
    result["outlet_concentrations_mol_per_m3"] = concentrations

    return result


# ═════════════════════════════════════════════════════════════════
# 设备 3：CSTR 反应器
# ═════════════════════════════════════════════════════════════════
def run_cstr_reactor(
    inlet_molar_flows_mol_per_s: Dict[str, float],
    stoichiometric_coefficients: Dict[str, float],
    key_component: str,
    key_component_conversion: float,
    reactor_temp_K: float,
    reactor_pressure_Pa: float,
    catalyst_particle_diameter_m: float,
    catalyst_particle_density_kg_per_m3: float,
    arrhenius_pre_exponential_factor: float,
    activation_energy_J_per_mol: float,
    reaction_order: int = 1,
) -> Dict[str, Any]:
    """运行 CSTR 反应器设备，返回 12 个底层工具结果合并后的全部参数。"""
    result: Dict[str, Any] = {}

    # 1. 反应进度
    extent = cstr_calc_reaction_extent(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    result.update(extent)

    # 2. 出口摩尔流量（dict-of-components）
    outlet_flows = cstr_calc_outlet_molar_flows(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    result["outlet_molar_flows_mol_per_s"] = outlet_flows

    # 3. 出口总摩尔流量
    total_flow = cstr_calc_outlet_total_molar_flow(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    result.update(total_flow)

    # 4. 出口体积流量
    vol_flow = cstr_calc_outlet_volumetric_flow(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
    )
    result.update(vol_flow)

    # 5. 出口浓度（dict-of-components）
    concentrations = cstr_calc_outlet_concentrations(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
    )
    result["outlet_concentrations_mol_per_m3"] = concentrations

    # 6. 流体密度
    density = cstr_calc_fluid_density(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
    )
    result.update(density)

    # 7. 流体粘度
    viscosity = cstr_calc_fluid_viscosity(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
    )
    result.update(viscosity)

    # 8. 阿基米德数
    archimedes = calc_archimedes_number(
        particle_diameter_m=catalyst_particle_diameter_m,
        particle_density_kg_per_m3=catalyst_particle_density_kg_per_m3,
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
    )
    result.update(archimedes)

    # 9. 最小流化速度
    umf = calc_minimum_fluidization_velocity(
        particle_diameter_m=catalyst_particle_diameter_m,
        particle_density_kg_per_m3=catalyst_particle_density_kg_per_m3,
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
    )
    result.update(umf)

    # 10. 反应速率
    rate = calc_reaction_rate(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
        arrhenius_pre_exponential_factor=arrhenius_pre_exponential_factor,
        activation_energy_J_mol=activation_energy_J_per_mol,
        reaction_order=reaction_order,
    )
    result.update(rate)

    # 11. CSTR 体积
    volume = calc_cstr_volume(
        method="kinetic",
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
        arrhenius_pre_exponential_factor=arrhenius_pre_exponential_factor,
        activation_energy_J_mol=activation_energy_J_per_mol,
        reaction_order=reaction_order,
    )
    result.update(volume)

    # 12. 停留时间
    residence = cstr_calc_residence_time(
        method="kinetic",
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
        arrhenius_pre_exponential_factor=arrhenius_pre_exponential_factor,
        activation_energy_J_mol=activation_energy_J_per_mol,
        reaction_order=reaction_order,
    )
    result.update(residence)

    return result


# ═════════════════════════════════════════════════════════════════
# 设备 4：PFR 反应器
# ═════════════════════════════════════════════════════════════════
def run_pfr_reactor(
    inlet_molar_flows_mol_per_s: Dict[str, float],
    stoichiometric_coefficients: Dict[str, float],
    key_component: str,
    key_component_conversion: float,
    reactor_temp_K: float,
    reactor_pressure_Pa: float,
    catalyst_particle_density_kg_per_m3: float,
    catalyst_bed_void_fraction: float,
    arrhenius_pre_exponential_factor: float,
    activation_energy_J_per_mol: float,
    reaction_order: int = 1,
    tube_inner_diameter_m: float = 0.05,
    superficial_velocity_m_per_s: float = 1.0,
) -> Dict[str, Any]:
    """运行 PFR 反应器设备，返回 14 个底层工具结果合并后的全部参数。"""
    result: Dict[str, Any] = {}

    # 1. 反应进度
    extent = pfr_calc_reaction_extent(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    result.update(extent)

    # 2. 出口摩尔流量
    outlet_flows = pfr_calc_outlet_molar_flows(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    result["outlet_molar_flows_mol_per_s"] = outlet_flows

    # 3. 出口总摩尔流量
    total_out = pfr_calc_outlet_total_molar_flow(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
    )
    result.update(total_out)

    # 4. 出口体积流量
    vol_out = pfr_calc_outlet_volumetric_flow(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
    )
    result.update(vol_out)

    # 5. 出口浓度
    concentrations = pfr_calc_outlet_concentrations(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
    )
    result["outlet_concentrations_mol_per_m3"] = concentrations

    # 6. 入口总摩尔流量
    inlet_total = calc_inlet_total_molar_flow(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
    )
    result.update(inlet_total)

    # 7. 入口体积流量（标准状态）
    inlet_stp = calc_inlet_volumetric_flow_stp(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
    )
    result.update(inlet_stp)

    # 8. 入口体积流量（操作工况）
    inlet_op = calc_inlet_volumetric_flow_operating(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
    )
    result.update(inlet_op)

    # 9. 床层体积（动力学法）
    bed_vol = calc_bed_volume(
        method="kinetic",
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
        arrhenius_pre_exponential_factor=arrhenius_pre_exponential_factor,
        activation_energy_J_mol=activation_energy_J_per_mol,
        reaction_order=reaction_order,
    )
    result.update(bed_vol)

    # 10. 催化剂质量
    cat_mass = calc_catalyst_mass(
        particle_density_kg_per_m3=catalyst_particle_density_kg_per_m3,
        bed_void_fraction=catalyst_bed_void_fraction,
        method="kinetic",
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
        arrhenius_pre_exponential_factor=arrhenius_pre_exponential_factor,
        activation_energy_J_mol=activation_energy_J_per_mol,
        reaction_order=reaction_order,
    )
    result.update(cat_mass)

    # 11. 单管直径
    tube_diam = calc_single_tube_diameter(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
        superficial_velocity_m_per_s=superficial_velocity_m_per_s,
    )
    result.update(tube_diam)

    # 12. 管数
    n_tubes = calc_tube_number(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
        tube_inner_diameter_m=tube_inner_diameter_m,
        superficial_velocity_m_per_s=superficial_velocity_m_per_s,
    )
    result.update(n_tubes)

    # 13. 单管床层长度
    bed_len_single = calc_bed_length_single_tube(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
        superficial_velocity_m_per_s=superficial_velocity_m_per_s,
        method="kinetic",
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        arrhenius_pre_exponential_factor=arrhenius_pre_exponential_factor,
        activation_energy_J_mol=activation_energy_J_per_mol,
        reaction_order=reaction_order,
    )
    result.update(bed_len_single)

    # 14. 多管床层长度
    bed_len_multi = calc_bed_length_multitube(
        inlet_molar_flows=inlet_molar_flows_mol_per_s,
        temperature_K=reactor_temp_K,
        pressure_Pa=reactor_pressure_Pa,
        tube_inner_diameter_m=tube_inner_diameter_m,
        superficial_velocity_m_per_s=superficial_velocity_m_per_s,
        method="kinetic",
        stoichiometric_coefficients=stoichiometric_coefficients,
        key_component=key_component,
        key_component_conversion=key_component_conversion,
        arrhenius_pre_exponential_factor=arrhenius_pre_exponential_factor,
        activation_energy_J_mol=activation_energy_J_per_mol,
        reaction_order=reaction_order,
    )
    result.update(bed_len_multi)

    return result


# ═════════════════════════════════════════════════════════════════
# 设备 5：板式精馏塔
# ═════════════════════════════════════════════════════════════════
def run_distillation_column(
    feed_molar_flows_mol_per_s: Dict[str, float],
    distillate_purity: float,
    bottoms_purity: float,
    light_key_component: str,
    heavy_key_component: str,
    light_components: List[str],
    heavy_components: List[str],
    feed_temp_K: float,
    feed_pressure_Pa: float,
    reflux_factor: float = 1.2,
) -> Dict[str, Any]:
    """运行精馏塔设备，返回扁平字典（所有 value 为 number/string/bool 或 {组分:数值}）。"""
    result: Dict[str, Any] = {}

    # 1. 物料衡算（原函数已返回扁平 dict）
    balance = calc_mass_balance(
        feed_molar_flows=feed_molar_flows_mol_per_s,
        distillate_purity=distillate_purity,
        bottoms_purity=bottoms_purity,
        light_key_component=light_key_component,
        heavy_key_component=heavy_key_component,
        light_components=light_components,
        heavy_components=heavy_components,
    )
    result.update(balance)

    # 2. 操作条件
    conditions = calc_operating_conditions(
        feed_molar_flows=feed_molar_flows_mol_per_s,
        distillate_purity=distillate_purity,
        bottoms_purity=bottoms_purity,
        light_key_component=light_key_component,
        heavy_key_component=heavy_key_component,
        light_components=light_components,
        heavy_components=heavy_components,
    )
    result.update(conditions)

    # 3. 进料热状态
    feed_q = calc_feed_thermal_condition(
        feed_molar_flows=feed_molar_flows_mol_per_s,
        distillate_purity=distillate_purity,
        bottoms_purity=bottoms_purity,
        light_key_component=light_key_component,
        heavy_key_component=heavy_key_component,
        light_components=light_components,
        heavy_components=heavy_components,
        feed_temp_K=feed_temp_K,
        feed_pressure_Pa=feed_pressure_Pa,
    )
    result.update(feed_q)

    # 4. 最小回流比
    min_reflux = calc_min_reflux_ratio(
        feed_molar_flows=feed_molar_flows_mol_per_s,
        distillate_purity=distillate_purity,
        bottoms_purity=bottoms_purity,
        light_key_component=light_key_component,
        heavy_key_component=heavy_key_component,
        light_components=light_components,
        heavy_components=heavy_components,
        feed_temp_K=feed_temp_K,
        feed_pressure_Pa=feed_pressure_Pa,
        reflux_factor=reflux_factor,
    )
    result.update(min_reflux)

    # 5. 最小理论板数
    min_stages = calc_min_theoretical_stages(
        feed_molar_flows=feed_molar_flows_mol_per_s,
        distillate_purity=distillate_purity,
        bottoms_purity=bottoms_purity,
        light_key_component=light_key_component,
        heavy_key_component=heavy_key_component,
        light_components=light_components,
        heavy_components=heavy_components,
    )
    result.update(min_stages)

    # 6. 理论板数
    theo_stages = calc_theoretical_stages(
        feed_molar_flows=feed_molar_flows_mol_per_s,
        distillate_purity=distillate_purity,
        bottoms_purity=bottoms_purity,
        light_key_component=light_key_component,
        heavy_key_component=heavy_key_component,
        light_components=light_components,
        heavy_components=heavy_components,
        feed_temp_K=feed_temp_K,
        feed_pressure_Pa=feed_pressure_Pa,
        reflux_factor=reflux_factor,
    )
    result.update(theo_stages)

    # 7. 实际塔板数
    actual = calc_actual_stages(
        feed_molar_flows=feed_molar_flows_mol_per_s,
        distillate_purity=distillate_purity,
        bottoms_purity=bottoms_purity,
        light_key_component=light_key_component,
        heavy_key_component=heavy_key_component,
        light_components=light_components,
        heavy_components=heavy_components,
        feed_temp_K=feed_temp_K,
        feed_pressure_Pa=feed_pressure_Pa,
        reflux_factor=reflux_factor,
    )
    result.update(actual)

    # 8. 进料板位置
    feed_stage = calc_feed_stage(
        feed_molar_flows=feed_molar_flows_mol_per_s,
        distillate_purity=distillate_purity,
        bottoms_purity=bottoms_purity,
        light_key_component=light_key_component,
        heavy_key_component=heavy_key_component,
        light_components=light_components,
        heavy_components=heavy_components,
        feed_temp_K=feed_temp_K,
        feed_pressure_Pa=feed_pressure_Pa,
        reflux_factor=reflux_factor,
    )
    result.update(feed_stage)

    # 9. 塔底压力修正：总压降按修正后的物理塔板数累加
    #    ΔP_total = N_actual × ΔP_tray；P_bottom = P_top + ΔP_total；
    #    T_bottom 在修正后压力下同步重算。
    from tray_distillation.funcs import (
        _DP_PER_TRAY_ATMOSPHERIC_Pa as _DP_ATM,
        _DP_PER_TRAY_PRESSURIZED_Pa as _DP_PRESS,
        _compute_mass_balance as _td_mb,
        _normalize_comp_ids,
        _solve_bubble_T_multicomp,
        _solve_temp_at_pressure,
    )
    top_P = result.get("column_top_pressure_Pa", 101325.0)
    N_act = result.get("N_actual", 0)
    dp_per_tray = _DP_ATM if top_P <= 200000.0 else _DP_PRESS
    dp_total = (N_act or 0) * dp_per_tray
    bot_P_corr = top_P + dp_total
    result["pressure_drop_per_tray_Pa"] = dp_per_tray
    result["column_pressure_drop_Pa"] = dp_total
    result["pressure_drop_tray_basis"] = N_act
    if N_act and bot_P_corr != result.get("column_bottom_pressure_Pa"):
        result["column_bottom_pressure_Pa"] = bot_P_corr
        mb2 = _td_mb(feed_molar_flows_mol_per_s, distillate_purity, bottoms_purity,
                     light_key_component, heavy_key_component,
                     light_components, heavy_components)
        if mb2["is_multicomponent"]:
            W_comps = list(mb2["bottoms_zs"].keys())
            W_fracs = list(mb2["bottoms_zs"].values())
            bot_T_corr = _solve_bubble_T_multicomp(W_comps, W_fracs, bot_P_corr)
        else:
            chem_l, chem_h = _normalize_comp_ids(mb2["component_ids"])
            bot_T_corr = _solve_temp_at_pressure(
                bot_P_corr, mb2["mole_fractions"]["W"], chem_l, chem_h, is_dew=False)
        result["column_bottom_temperature_K"] = bot_T_corr

    return result


# ═════════════════════════════════════════════════════════════
# 设备 5b：回收率 + 质量分数规格的精馏/分离塔（流程级简化衡算）
# ═════════════════════════════════════════════════════════════
def _mw_map(components, MW=None):
    """g/mol 分子量字典；未提供则从 thermo 查询。"""
    if MW:
        return {c: float(MW[c]) for c in components}
    from thermo import Chemical
    out = {}
    for c in components:
        m = getattr(Chemical(c), 'MW', None)
        if not m:
            raise ValueError(f"无法获取组分 {c} 的分子量，请通过 MW 参数显式提供")
        out[c] = float(m)
    return out


def run_column_recovery(
    feed_molar_flows_mol_per_s: Dict[str, float],
    recover_to_top: Dict[str, float],
    *,
    top_mass_fraction: Any = None,
    bottom_mass_fraction: Any = None,
    swing_component: str = None,
    MW: Dict[str, float] = None,
) -> Dict[str, Any]:
    """按「关键组分回收率 + 产品质量分数(wt%)」规格做塔顶/塔底逐组分物料衡算。

    与 calc_mass_balance（仅支持摩尔纯度 Fenske 型分裂）互补，适用于题目直接给出
    回收率与 wt% 纯度、而塔内逐组分去向可由质守恒唯一确定的分离塔。

    Parameters
    ----------
    feed_molar_flows_mol_per_s : {组分: mol/s}
    recover_to_top : {组分: 该组分进料去塔顶的分率}；未列组分默认 0.0（全部去塔底）。
        可先给出“清晰分割”组分的 0/1，分配组分可留空由 wt% 规格反算。
    top_mass_fraction : (组分名, 值) 可选——强制塔顶该组分质量分数，通过调节
        swing_component 的塔顶流量满足（需同时给出 swing_component）。
    bottom_mass_fraction : (组分名, 值) 可选——同理强制塔底质量分数。
    swing_component : 受 wt% 规格调节的“摆动组分”名。
    MW : 可选 g/mol 字典，缺省从 thermo 查询。

    Returns 扁平 dict，含 top/bottom 逐组分摩尔流量、总量、质量流量及回收率。
    """
    feed = {c: float(v) for c, v in feed_molar_flows_mol_per_s.items() if float(v) > 0}
    comps = list(feed)
    mw = _mw_map(comps, MW)
    # 初始按 recover_to_top 分配（未列组分→塔底）
    top = {c: feed[c] * float(recover_to_top.get(c, 0.0)) for c in comps}

    def _apply_wt(spec):
        target_comp, frac = spec[0], float(spec[1])
        if not (0 < frac < 1):
            raise ValueError('质量分数规格必须在 (0,1)')
        if swing_component not in comps:
            raise ValueError('需指定存在于进料中的 swing_component')
        if target_comp == swing_component:
            raise ValueError('wt% 目标组分不能就是 swing_component')
        fixed_mass_excl_swing = sum(top[c] * mw[c] for c in comps if c != swing_component)
        target_mass = top[target_comp] * mw[target_comp]
        if target_mass <= 0:
            raise ValueError('wt% 目标组分在塔顶流量为零，规格不可行')
        total_mass = target_mass / frac
        swing_mass = total_mass - fixed_mass_excl_swing
        swing_mol = swing_mass / mw[swing_component]
        if swing_mol < -1e-9:
            raise ValueError(f'塔顶规格不可行：需 swing 质量 {swing_mass:.4g} g/s (<0)')
        top[swing_component] = min(max(swing_mol, 0.0), feed[swing_component])

    if top_mass_fraction is not None:
        _apply_wt(top_mass_fraction)
    # 塔底：先由质量守恒得 bottom，再按需用 bottom_mass_fraction 反算 swing
    bottom = {c: feed[c] - top[c] for c in comps}
    if bottom_mass_fraction is not None:
        target_comp, frac = bottom_mass_fraction[0], float(bottom_mass_fraction[1])
        if swing_component not in comps:
            raise ValueError('需指定存在于进料中的 swing_component')
        # 塔底总质量 = 目标组分质量 / frac；塔底 swing 质量 = 总 - 其余固定
        fixed_bot_mass_excl_swing = sum(bottom[c] * mw[c] for c in comps
                                        if c not in (target_comp, swing_component))
        target_mass = bottom[target_comp] * mw[target_comp]
        if target_mass <= 0:
            raise ValueError('塔底 wt% 目标组分流量为零，规格不可行')
        bot_total_mass = target_mass / frac
        swing_bot_mass = bot_total_mass - target_mass - fixed_bot_mass_excl_swing
        swing_bot_mol = swing_bot_mass / mw[swing_component]
        if swing_bot_mol < -1e-9:
            raise ValueError(f'塔底规格不可行：需 swing 质量 {swing_bot_mass:.4g} g/s (<0)')
        bottom[swing_component] = min(max(swing_bot_mol, 0.0), feed[swing_component])
        top[swing_component] = feed[swing_component] - bottom[swing_component]

    top = {c: v for c, v in top.items() if v > 0}
    bottom = {c: v for c, v in bottom.items() if v > 0}
    top_mass = sum(v * mw[c] for c, v in top.items())
    bot_mass = sum(v * mw[c] for c, v in bottom.items())
    result: Dict[str, Any] = {
        "top_flows_mol_per_s": top,
        "bottom_flows_mol_per_s": bottom,
        "top_total_mol_per_s": sum(top.values()),
        "bottom_total_mol_per_s": sum(bottom.values()),
        "top_mass_flow_kg_per_s": top_mass / 1000.0,
        "bottom_mass_flow_kg_per_s": bot_mass / 1000.0,
    }
    # 回收率报告（每个进料关键组分去塔顶的分率）
    result["recoveries_to_top"] = {
        c: (top.get(c, 0.0) / feed[c]) for c in feed if feed[c] > 0}
    return result


# ═════════════════════════════════════════════════════════════════
# 设备 6：泵
# ═════════════════════════════════════════════════════════════════
def run_pump(
    inlet_pressure_Pa: float,
    target_pressure_Pa: float,
) -> Dict[str, Any]:
    """运行泵设备，返回出口压力和压升。"""
    return pump_calc_outlet_pressure(
        inlet_pressure_Pa=inlet_pressure_Pa,
        target_pressure_Pa=target_pressure_Pa,
    )


def run_pressure_reducing_valve(
    inlet_pressure_Pa: float,
    target_pressure_Pa: float,
) -> Dict[str, Any]:
    """运行减压阀设备，返回出口压力和压降。"""
    return prv_calc_outlet_pressure(
        inlet_pressure_Pa=inlet_pressure_Pa,
        target_pressure_Pa=target_pressure_Pa,
    )


# ═════════════════════════════════════════════════════════════════
# 设备 7：闪蒸罐
# ═════════════════════════════════════════════════════════════════
def run_flash_drum(
    flash_pressure_Pa: float,
    inlet_molar_flows_mol_per_s: Dict[str, float],
    key_component: str,
    key_component_recovery: float,
    recovery_phase: str = "vapor",
) -> Dict[str, Any]:
    """运行闪蒸罐设备，反算闪蒸温度并返回气液出口流量分配。"""
    return calc_flash_drum(
        flash_pressure_Pa=flash_pressure_Pa,
        inlet_molar_flows_mol_per_s=inlet_molar_flows_mol_per_s,
        key_component=key_component,
        key_component_recovery=key_component_recovery,
        recovery_phase=recovery_phase,
    )


# ═════════════════════════════════════════════════════════════════
# 设备 8：混合器
# ═════════════════════════════════════════════════════════════════
def run_mixer(
    inlet_molar_flows_list: List[Dict[str, float]],
    inlet_temperatures_K: List[float],
    inlet_pressures_Pa: List[float],
    **options,
) -> Dict[str, Any]:
    """运行混合器设备，返回出口流量、组成、压力、绝热温度与体积流量。

    优先调用 mixer.funcs.calc_mixer_state 一次性获得物料/能量互相一致的出口状态
    （mixer 模块自身推荐流程计算走此接口，避免逐组分独立 TP 闪蒸对两相纯物流
    气化率的歧义与重复闪蒸）。当一次求解不可用时，回退到分离的 calc_outlet_* 接口。
    options 透传给 calc_mixer_state（如 thermodynamic_model / outlet_pressure_Pa /
    pressure_drop_Pa / inlet_vapor_fractions / inlet_enthalpies_J_per_mol）。
    """
    try:
        return mixer_calc_state(
            inlet_molar_flows_list, inlet_temperatures_K, inlet_pressures_Pa,
            **options)
    except Exception:
        if options:
            raise

    result: Dict[str, Any] = {}

    result.update(mixer_calc_outlet_molar_flows(inlet_molar_flows_list))
    result.update(mixer_calc_outlet_total_molar_flow(inlet_molar_flows_list))
    result.update(mixer_calc_outlet_composition(inlet_molar_flows_list))

    pressure = mixer_calc_outlet_pressure(inlet_pressures_Pa)
    result.update(pressure)

    temperature = mixer_calc_outlet_temperature(
        inlet_molar_flows_list, inlet_temperatures_K, inlet_pressures_Pa)
    result.update(temperature)

    result.update(mixer_calc_outlet_volumetric_flow(
        inlet_molar_flows_list,
        outlet_temperature_K=temperature["outlet_temperature_K"],
        outlet_pressure_Pa=pressure["outlet_pressure_Pa"],
    ))

    return result
