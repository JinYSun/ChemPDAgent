"""闪蒸罐完整设计 MCP 工具"""
import math
import numpy as np
from physics_engine.flash_drum import (
    isothermal_flash, flash_drum_diameter, flash_drum_volume,
    validate_flash_drum_design, entrainment_fraction,
)
from physics_engine.thermo_helper import get_mixture_properties, get_mixture_k_values
from physics_engine.common import _REQUIRED, check_required_params
from tools.base import ToolRegistry
from tools.registry import UniversalToolRegistry, ToolCapability

# ============================================================
# 闪蒸罐设计约束（硬编码）
# ============================================================
FLASH_DRUM_CONSTRAINTS = {
    "L_D_min": 2.0,           # 长径比下限
    "L_D_max": 5.0,           # 长径比上限
    "L_D_target": 3.0,        # 目标长径比
    "u_max": 3.0,             # 最大气速 (m/s)
    "tau_vapor_min": 5.0,     # 最小气相停留时间 (s)
    "tau_liquid_min": 180.0,  # 最小液相停留时间 (s)
    "tau_vapor_default": 10.0,    # 默认气相停留时间 (s)
    "tau_liquid_default": 300.0,  # 默认液相停留时间 (s)
}


def _full_design(components: list = _REQUIRED, z: list = _REQUIRED,
                 T: float = _REQUIRED, P: float = _REQUIRED,
                 F: float = _REQUIRED,
                 # ── 可选：用户给定物性/设计参数（给定后跳过内部计算/查询） ──
                 K_values: list = None,    # 气液平衡常数列表
                 rho_V: float = None,      # 气相密度 kg/m³
                 rho_L: float = None,      # 液相密度 kg/m³
                 sigma: float = None,      # 表面张力 N/m
                 **kwargs) -> dict:
    """闪蒸罐完整设计流水线"""
    # 检查必填参数（无默认值，缺失时触发 human-in-the-loop）
    missing = check_required_params(locals(), {
        "components": "组分名称列表，如 ['methanol', 'water']",
        "z": "进料摩尔分数列表，如 [0.4, 0.6]",
        "T": "闪蒸温度 (K)",
        "P": "闪蒸压力 (Pa)",
        "F": "进料摩尔流量 (kmol/h)",
    })
    if missing:
        return missing

    # 1. 获取 K 值 (气液平衡常数) — 用户给定值优先
    conversion_notes = []
    if K_values is not None:
        K_list = [float(kv) for kv in K_values]
        conversion_notes.append(f"[用户给定] 气液平衡常数 K_values = {K_list}")
    else:
        K_list = get_mixture_k_values(components, T, P, z)
    K = np.array(K_list, dtype=float)
    z_arr = np.array(z, dtype=float)

    # 2. 等温闪蒸计算 (进料单位 kmol/h → mol/s)
    feed_flow = F * 1000.0 / 3600.0  # mol/s
    flash = isothermal_flash(feed_flow, z_arr, K)

    # 3. 混合物物性（用户给定值优先）
    props = get_mixture_properties(components, z, T, P)

    # 4. 闪蒸罐几何设计（带约束自动调整）
    drum_design = {}
    skip_notes = []
    adjust_notes = []
    # ── 密度：用户给定值优先 ──
    if rho_V is not None:
        rho_V_val = float(rho_V)
        conversion_notes.append(f"[用户给定] 气相密度 rho_V = {rho_V_val:.2f} kg/m³")
    else:
        rho_V_val = props.get('rho_V')
    if rho_L is not None:
        rho_L_val = float(rho_L)
        conversion_notes.append(f"[用户给定] 液相密度 rho_L = {rho_L_val:.2f} kg/m³")
    else:
        rho_L_val = props.get('rho_L')
    MW_mix = props.get('MW_mix')

    if rho_V_val is None or rho_L_val is None:
        skip_notes.append("[跳过] 物性数据不完整（缺气相/液相密度），闪蒸罐几何设计已跳过")
    elif rho_V_val > 0 and rho_L_val > 0 and flash.get('psi', 0) > 0:
        V_mol = flash['V']  # mol/s
        L_mol = flash['L']  # mol/s
        if MW_mix is None:
            MW_mix = 50.0
            skip_notes.append("[估算] MW_mix 未从物性数据库获取，使用估算值 50 g/mol")
        V_mass = V_mol * MW_mix / 1000  # kg/s
        L_mass = L_mol * MW_mix / 1000  # kg/s
        V_vol = V_mass / rho_V_val if rho_V_val > 0 else 0.001  # m³/s
        L_vol = L_mass / rho_L_val if rho_L_val > 0 else 0.001  # m³/s
        
        try:
            # 获取约束参数
            C = FLASH_DRUM_CONSTRAINTS
            u_max_allowed = C["u_max"]
            L_D_target = C["L_D_target"]
            L_D_min = C["L_D_min"]
            L_D_max = C["L_D_max"]
            tau_vapor = C["tau_vapor_default"]
            tau_liquid = C["tau_liquid_default"]
            
            # 步骤1: 基于停留时间计算容积和初始直径
            drum_vol = flash_drum_volume(V_vol, L_vol, rho_V_val, rho_L_val, 
                                         vapor_holdup_time=tau_vapor,
                                         liquid_holdup_time=tau_liquid)
            D_vol = drum_vol["D"]
            L_vol_total = drum_vol["L_total"]
            
            # 步骤2: 基于气速约束计算最小直径
            drum_diam = flash_drum_diameter(V_vol, rho_V_val, max_vapor_velocity=u_max_allowed)
            D_vel = drum_diam["D"]
            
            # 步骤3: 取两者较大值（满足所有约束）
            D_final = max(D_vol, D_vel)
            if D_vel > D_vol:
                adjust_notes.append(f"[调整] 气速约束主导，直径从 {D_vol:.3f}m 增大至 {D_vel:.3f}m")
            
            # 步骤4: 基于最终直径重新计算长度（保持容积需求）
            V_total_req = drum_vol["V_total"]
            A_cross = math.pi * D_final**2 / 4
            L_required = V_total_req / A_cross
            
            # 添加分离高度
            disengagement_height = 1.0
            L_total_final = L_required + disengagement_height
            
            # 步骤5: 重新计算气速和停留时间
            u_actual = V_vol / A_cross
            V_vapor_space = (L_total_final - disengagement_height) * A_cross * 0.3
            V_liquid_space = (L_total_final - disengagement_height) * A_cross * 0.6
            tau_vapor_actual = V_vapor_space / V_vol if V_vol > 0 else 0
            tau_liquid_actual = V_liquid_space / L_vol if L_vol > 0 else 0
            
            # 步骤6: 检查停留时间约束，必要时增加容积
            if tau_vapor_actual < C["tau_vapor_min"]:
                scale_factor = C["tau_vapor_min"] / tau_vapor_actual
                V_total_req *= scale_factor
                adjust_notes.append(f"[调整] 气相停留时间不足，容积放大 {scale_factor:.2f} 倍")
                # 重新计算尺寸
                A_cross = math.pi * D_final**2 / 4
                L_required = V_total_req / A_cross
                L_total_final = L_required + disengagement_height
                V_vapor_space = (L_total_final - disengagement_height) * A_cross * 0.3
                V_liquid_space = (L_total_final - disengagement_height) * A_cross * 0.6
                tau_vapor_actual = V_vapor_space / V_vol if V_vol > 0 else 0
                tau_liquid_actual = V_liquid_space / L_vol if L_vol > 0 else 0
            
            if tau_liquid_actual < C["tau_liquid_min"]:
                scale_factor = C["tau_liquid_min"] / tau_liquid_actual
                V_total_req *= scale_factor
                adjust_notes.append(f"[调整] 液相停留时间不足，容积放大 {scale_factor:.2f} 倍")
                # 重新计算尺寸
                A_cross = math.pi * D_final**2 / 4
                L_required = V_total_req / A_cross
                L_total_final = L_required + disengagement_height
                V_vapor_space = (L_total_final - disengagement_height) * A_cross * 0.3
                V_liquid_space = (L_total_final - disengagement_height) * A_cross * 0.6
                tau_vapor_actual = V_vapor_space / V_vol if V_vol > 0 else 0
                tau_liquid_actual = V_liquid_space / L_vol if L_vol > 0 else 0
            
            # 步骤7: 最后检查长径比约束，必要时调整直径（在所有容积调整完成后）
            L_D_actual = L_total_final / D_final
            if L_D_actual > L_D_max:
                # 长径比过高，需要增大直径使 L/D 降到目标值
                # 近似: D = (4 * V_total / (pi * L_D_target))^(1/3)
                D_final = (4 * V_total_req / (math.pi * L_D_target)) ** (1/3)
                A_cross = math.pi * D_final**2 / 4
                L_required = V_total_req / A_cross
                L_total_final = L_required + disengagement_height
                L_D_actual = L_total_final / D_final
                u_actual = V_vol / A_cross
                # 重新计算停留时间
                V_vapor_space = (L_total_final - disengagement_height) * A_cross * 0.3
                V_liquid_space = (L_total_final - disengagement_height) * A_cross * 0.6
                tau_vapor_actual = V_vapor_space / V_vol if V_vol > 0 else 0
                tau_liquid_actual = V_liquid_space / L_vol if L_vol > 0 else 0
                adjust_notes.append(f"[调整] 长径比超限，直径增大至 {D_final:.3f}m，L/D={L_D_actual:.2f}")
            elif L_D_actual < L_D_min and D_final > 0:
                # 长径比过低（太胖），可适当减小直径
                D_final = (4 * V_total_req / (math.pi * L_D_target)) ** (1/3)
                A_cross = math.pi * D_final**2 / 4
                L_required = V_total_req / A_cross
                L_total_final = L_required + disengagement_height
                L_D_actual = L_total_final / D_final
                u_actual = V_vol / A_cross
                V_vapor_space = (L_total_final - disengagement_height) * A_cross * 0.3
                V_liquid_space = (L_total_final - disengagement_height) * A_cross * 0.6
                tau_vapor_actual = V_vapor_space / V_vol if V_vol > 0 else 0
                tau_liquid_actual = V_liquid_space / L_vol if L_vol > 0 else 0
                adjust_notes.append(f"[调整] 长径比过低，直径调整为 {D_final:.3f}m，L/D={L_D_actual:.2f}")
            
            # 步骤8: 液相段和气相段高度分配
            L_liquid = V_liquid_space / A_cross if A_cross > 0 else 0
            L_vapor = L_total_final - L_liquid
            
            # 组装最终结果
            drum_design = {
                "V_total": V_total_req,
                "V_vapor": V_vapor_space,
                "V_liquid": V_liquid_space,
                "D": round(D_final, 4),
                "L_liquid": round(L_liquid, 4),
                "L_vapor": round(L_vapor, 4),
                "L_total": round(L_total_final, 4),
                "aspect_ratio": round(L_total_final / D_final, 4) if D_final > 0 else 0,
                "disengagement_height": disengagement_height,
                "tau_vapor": round(tau_vapor_actual, 1),
                "tau_liquid": round(tau_liquid_actual, 1),
                "u_actual": round(u_actual, 4),
                "u_max": u_max_allowed,
                "A_cross": round(A_cross, 6),
                "V_mass_kg_s": round(V_mass, 4),
                "L_mass_kg_s": round(L_mass, 4),
                "V_vol_m3_s": round(V_vol, 6),
                "L_vol_m3_s": round(L_vol, 6),
                "MW_mix_g_mol": round(MW_mix, 2),
            }
        except Exception as e:
            skip_notes.append(f"[错误] 闪蒸罐几何设计失败: {e}")

    # 5. 雾沫夹带估算
    entrainment = None
    u_actual = drum_design.get('u_actual', 0)
    if u_actual > 0:
        mu_L = props.get('mu_L')      # 不提供默认值
        sigma = props.get('sigma')     # 不提供默认值
        ent_result = entrainment_fraction(
            vapor_velocity=u_actual,
            liquid_viscosity=mu_L,
            surface_tension=sigma,
            liquid_density=rho_L,
            gas_density=rho_V,
        )
        ent_frac = ent_result.get('entrainment_fraction')
        entrainment = {
            'method': ent_result.get('method'),
            'entrainment_fraction': round(ent_frac, 5) if ent_frac is not None else None,
            'entrainment_note': ent_result.get('note', ''),
        }

    # 6. 停留时间
    residence_time = None
    if drum_design:
        V_total = drum_design.get('V_total', 0)
        V_vol = drum_design.get('V_vol_m3_s', 0)
        L_vol = drum_design.get('L_vol_m3_s', 0)
        V_vapor = drum_design.get('V_vapor', 0)
        V_liquid = drum_design.get('V_liquid', 0)
        residence_time = {
            'tau_vapor_s': round(V_vapor / V_vol, 1) if V_vol > 0 else None,
            'tau_liquid_s': round(V_liquid / L_vol, 1) if L_vol > 0 else None,
        }

    # 7. 压降估算（简化 K 值法）
    pressure_drop_Pa = None
    if rho_V and u_actual:
        # 进口速度头损失
        u_inlet = u_actual * 1.5  # 入口缩小效应
        dP = 0.5 * rho_V * u_inlet ** 2 * 1.5  # 当量阻力系数
        pressure_drop_Pa = round(dP, 1)

    # 8. 格式化输出
    psi = flash.get('psi', 0)
    result = {
        'flash': {
            'psi': round(psi, 6),
            'V_mol_s': round(flash.get('V', 0), 4),
            'L_mol_s': round(flash.get('L', 0), 4),
            'F_mol_s': round(feed_flow, 4),
            'F_kmol_h': F,
            'y': flash.get('y', []).tolist() if hasattr(flash.get('y', []), 'tolist') else flash.get('y', []),
            'x': flash.get('x', []).tolist() if hasattr(flash.get('x', []), 'tolist') else flash.get('x', []),
            'converged': flash.get('converged', False),
        },
        'temperatures': {
            'T_flash_K': T,
            'T_flash_C': round(T - 273.15, 1),
            'P_Pa': P,
            'P_bar': round(P / 1e5, 3),
        },
        'properties': props,
        'drum_design': drum_design,
        'entrainment': entrainment,
        'residence_time': residence_time,
        'pressure_drop_Pa': pressure_drop_Pa,
        'K_values': K_list,
        'skip_notes': skip_notes,
        'adjust_notes': adjust_notes,
    }

    # 9. 约束校验
    if drum_design:
        validation = validate_flash_drum_design(drum_design)
        result['validation'] = validation
        result['warnings'] = validation.get('warnings', [])
        result['converged'] = validation.get('valid', True)
    else:
        result['validation'] = None
        result['warnings'] = skip_notes
        result['converged'] = flash.get('converged', False)

    # 添加 conversion_notes
    if conversion_notes:
        result['conversion_notes'] = conversion_notes

    return result


def register_flash_drum_device(registry: ToolRegistry) -> ToolRegistry:
    registry.register_function(name="flash_drum_design",
        description="闪蒸罐完整设计。封装完整设计流水线: 物性查询→计算→校验→报告。",
        func=_full_design,
        param_descriptions={
            "components": "(必填) 组分名称列表，无默认值",
            "z": "(必填) 进料摩尔分数列表，无默认值",
            "T": "(必填) 闪蒸温度 (K)，无默认值",
            "P": "(必填) 闪蒸压力 (Pa)，无默认值",
            "F": "(必填) 进料摩尔流量 (kmol/h)，无默认值"
        },
        category="device", tags=["flash_drum", "flash", "闪蒸"])
    return registry


def register_to_universal(registry: UniversalToolRegistry) -> None:
    """注册到通用注册中心"""
    capability = ToolCapability(
        name="flash_drum_design",
        description="闪蒸罐完整设计。计算气化分率、罐径、罐高等参数。",
        category="device_design",
        parameters={
            "components": {"type": "list", "required": True, "description": "(必填) 组分英文名列表，如 ['methanol', 'water']，无默认值"},
            "z": {"type": "list", "required": True, "description": "(必填) 进料摩尔分数列表，如 [0.4, 0.6]，无默认值"},
            "T": {"type": "number", "required": True, "description": "(必填) 闪蒸温度 K，无默认值"},
            "P": {"type": "number", "required": True, "description": "(必填) 闪蒸压力 Pa，无默认值"},
            "F": {"type": "number", "required": True, "description": "(必填) 进料摩尔流量 kmol/h，无默认值"},
            # ── 可选物性参数 ──
            "K_values": {"type": "list", "required": False, "description": "(选填) 气液平衡常数列表，给定后跳过内部计算"},
            "rho_V": {"type": "number", "required": False, "description": "(选填) 气相密度 kg/m³，给定后跳过物性查询"},
            "rho_L": {"type": "number", "required": False, "description": "(选填) 液相密度 kg/m³，给定后跳过物性查询"},
            "sigma": {"type": "number", "required": False, "description": "(选填) 表面张力 N/m，给定后跳过物性查询"},
        },
        examples=[
            {
                "description": "甲醇-水闪蒸",
                "input": {
                    "components": ["methanol", "water"],
                    "z": [0.4, 0.6],
                    "T": 350,
                    "P": 101325,
                },
            },
        ],
        tags=["flash_drum", "闪蒸", "气液分离"],
    )
    registry.register(capability, _full_design)