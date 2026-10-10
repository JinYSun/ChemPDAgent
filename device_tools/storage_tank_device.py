"""储罐完整设计 MCP 工具"""
import math
from physics_engine.storage_tank import (
    design_vertical_tank, design_horizontal_tank, design_spherical_tank,
    validate_tank_design, breathing_losses, insulation_thickness, tank_foundation_load,
)
from physics_engine.thermo_helper import get_fluid_density, get_fluid_MW, get_fluid_vapor_pressure
from physics_engine.common import _REQUIRED, check_required_params
from tools.base import ToolRegistry
from tools.registry import UniversalToolRegistry, ToolCapability

# ============================================================
# 储罐标准参数集（工程常用值）
# ============================================================
# 操作压力 (Pa) — 常用规格
STANDARD_PRESSURES = {
    101325.0: "常压 (1 atm)",
    200000.0: "低压 (2 bar)",
    500000.0: "中压 (5 bar)",
    1000000.0: "较高压 (10 bar)",
}
DEFAULT_PRESSURE = 101325.0  # 常压

# 环境温度 (K) — 常用规格
STANDARD_T_AMBIENT = {
    293.15: "20°C (温带室内)",
    298.15: "25°C (标准温度)",
    303.15: "30°C (夏季/热带)",
    278.15: "5°C (寒冷地区)",
}
DEFAULT_T_AMBIENT = 298.15  # 25°C

# 长径比 (L/D) — 按罐型分类
STANDARD_ASPECT_RATIO = {
    "vertical":   {"default": 4.0, "range": [2.0, 3.0, 4.0, 5.0, 6.0], "note": "立式罐常用 3-5"},
    "horizontal": {"default": 3.0, "range": [2.0, 3.0, 4.0, 5.0], "note": "卧式罐常用 2-4"},
    "spherical":  {"default": 1.0, "range": [1.0], "note": "球罐固定为 1.0"},
}

def _full_design(V: float = _REQUIRED, tank_type: str = _REQUIRED,
                 fluid: str = _REQUIRED, T: float = _REQUIRED,
                 P: float = None,
                 T_ambient: float = None,
                 aspect_ratio: float = None,
                 # ── 可选：用户给定物性参数（给定后跳过内部查询） ──
                 rho: float = None,        # 液体密度 kg/m³
                 MW: float = None,         # 分子量 g/mol
                 Pv: float = None,         # 饱和蒸气压 Pa
                 **kwargs) -> dict:
    # 检查必填参数（核心工艺参数必须由用户提供，P/T_ambient/aspect_ratio 有工程默认值）
    missing = check_required_params(locals(), {
        "V": "设计体积 (m³)",
        "tank_type": "罐类型: 'vertical'/'horizontal'/'spherical'",
        "fluid": "储存介质名称，如 'water', 'methanol'",
        "T": "操作温度 (K)",
    })
    if missing:
        return missing

    # P / T_ambient / aspect_ratio 默认值处理
    _user_defaulted_params = {}
    conversion_notes = []

    if P is None or P is _REQUIRED:
        P = DEFAULT_PRESSURE
        _user_defaulted_params["P"] = f"{P:.0f} Pa ({STANDARD_PRESSURES.get(P, '自定义')})，可选: {[f'{v:.0f}Pa({d})' for v, d in STANDARD_PRESSURES.items()]}"
        conversion_notes.append(f"[默认值] P 使用常压 {P:.0f} Pa")
    P = float(P)

    if T_ambient is None or T_ambient is _REQUIRED:
        T_ambient = DEFAULT_T_AMBIENT
        _user_defaulted_params["T_ambient"] = f"{T_ambient:.2f} K ({STANDARD_T_AMBIENT.get(T_ambient, '自定义')})，可选: {[f'{v:.2f}K({d})' for v, d in STANDARD_T_AMBIENT.items()]}"
        conversion_notes.append(f"[默认值] T_ambient 使用标准温度 {T_ambient:.2f} K ({T_ambient - 273.15:.1f}°C)")
    T_ambient = float(T_ambient)

    if aspect_ratio is None or aspect_ratio is _REQUIRED:
        ar_info = STANDARD_ASPECT_RATIO.get(tank_type, STANDARD_ASPECT_RATIO["vertical"])
        aspect_ratio = ar_info["default"]
        _user_defaulted_params["aspect_ratio"] = f"{aspect_ratio}（{tank_type}型，可选: {ar_info['range']}，{ar_info['note']}）"
        conversion_notes.append(f"[默认值] aspect_ratio 使用 {aspect_ratio}（{tank_type}型标准值）")
    aspect_ratio = float(aspect_ratio)

    # 1. 储罐几何设计
    design_kwargs = {}
    if tank_type != 'spherical' and aspect_ratio is not None:
        design_kwargs['aspect_ratio'] = aspect_ratio
    if tank_type == 'vertical':
        result = design_vertical_tank(V, **design_kwargs)
    elif tank_type == 'horizontal':
        result = design_horizontal_tank(V, **design_kwargs)
    else:
        result = design_spherical_tank(V, **design_kwargs)

    # 2. 设计容积、理论尺寸
    design_vol = V / 0.85  # 预留 15% 气相空间
    result['design_volume_m3'] = round(design_vol, 4)
    result['required_volume_m3'] = V
    result['tank_type'] = tank_type

    # 3. 流体物性查询（用户给定值优先）
    fluid_props = {}
    rho_val = None
    MW_val = None
    # ── 密度：用户给定值优先 ──
    if rho is not None:
        rho_val = float(rho)
        fluid_props['rho_liquid_kg_m3'] = round(rho_val, 1)
        conversion_notes.append(f"[用户给定] 液体密度 rho = {rho_val} kg/m³")
    else:
        try:
            rho_val = get_fluid_density(fluid, T, P, phase="liquid")
            if rho_val:
                fluid_props['rho_liquid_kg_m3'] = round(rho_val, 1)
        except Exception:
            pass
    # ── 分子量：用户给定值优先 ──
    if MW is not None:
        MW_val = float(MW)
        fluid_props['MW_g_mol'] = round(MW_val, 2)
        conversion_notes.append(f"[用户给定] 分子量 MW = {MW_val} g/mol")
    else:
        try:
            MW_val = get_fluid_MW(fluid)
            if MW_val:
                fluid_props['MW_g_mol'] = round(MW_val, 2)
        except Exception:
            pass
    # ── 饱和蒸气压：用户给定值优先 ──
    if Pv is not None:
        fluid_props['vapor_pressure_Pa'] = round(float(Pv), 1)
        conversion_notes.append(f"[用户给定] 饱和蒸气压 Pv = {Pv} Pa")
    else:
        try:
            P_vap = get_fluid_vapor_pressure(fluid, T)
            if P_vap:
                fluid_props['vapor_pressure_Pa'] = round(P_vap, 1)
        except Exception:
            pass
    result['fluid_properties'] = fluid_props

    # 4. 流体质量
    if rho_val:
        result['liquid_mass_kg'] = round(V * rho_val, 1)
        result['liquid_mass_tonnes'] = round(V * rho_val / 1000, 3)

    # 5. 基础荷载计算
    D = result.get('D', 0)
    H = result.get('H', result.get('L', 0))
    if D > 0 and H > 0 and rho_val:
        foundation = tank_foundation_load(D, H, rho_val)
        result['foundation'] = {
            'total_weight_N': round(foundation['total_weight'], 0),
            'total_weight_kN': round(foundation['total_weight'] / 1000, 1),
            'liquid_weight_N': round(foundation['liquid_weight'], 0),
            'shell_weight_N': round(foundation['shell_weight'], 0),
            'bottom_weight_N': round(foundation['bottom_weight'], 0),
            'base_pressure_Pa': round(foundation['base_pressure'], 0),
            'base_pressure_kPa': round(foundation['base_pressure'] / 1000, 1),
        }

    # 6. 保温层设计（当 T ≠ T_ambient 时）
    if D > 0 and abs(T - T_ambient) > 10:
        ins = insulation_thickness(D, T, T_ambient)
        result['insulation'] = {
            'thickness_m': round(ins['thickness'], 4),
            'thickness_mm': round(ins['thickness'] * 1000, 1),
            'heat_loss_W_m2': round(ins['heat_loss'], 1),
            'surface_temperature_K': round(ins['surface_temperature'], 1),
            'insulation_k_W_mK': ins['insulation_k'],
        }

    # 7. 蒸发损耗（呼吸损耗）估算
    P_vap = fluid_props.get('vapor_pressure_Pa')
    if D > 0 and P_vap and MW:
        try:
            bl = breathing_losses(D, T + 10, T - 10, P_vap, MW)
            result['breathing_losses'] = {
                'loss_per_breath_kg': round(bl['loss_per_breath'], 4),
                'annual_loss_kg': round(bl['annual_loss'], 1),
                'breaths_per_day': bl['n_breaths_per_day'],
            }
        except Exception:
            pass

    # 8. 约束校验
    validate_type = {'vertical': 'vertical_cylindrical', 'horizontal': 'horizontal_cylindrical'}.get(tank_type, 'spherical')
    validation = validate_tank_design(result, validate_type)
    result['validation'] = validation
    result['warnings'] = validation.get('warnings', []) + conversion_notes
    result['converged'] = validation.get('valid', True)

    # 9. 用户默认参数报告
    if _user_defaulted_params:
        _desc_map = {"P": "操作压力", "T_ambient": "环境温度", "aspect_ratio": "长径比(L/D)"}
        lines = ["\n---\n⚙️ **储罐默认参数提醒**（以下参数用户未提供，已使用工程默认值）：\n"]
        for _k, _v in _user_defaulted_params.items():
            _label = _desc_map.get(_k, _k)
            lines.append(f"- **{_label}** = {_v}")
        result["user_defaulted_params_report"] = "\n".join(lines)

    return result

def register_storage_tank_device(registry: ToolRegistry) -> ToolRegistry:
    registry.register_function(name="storage_tank_design",
        description="储罐完整设计。封装完整设计流水线: 物性查询→计算→校验→报告。",
        func=_full_design,
        param_descriptions={
            "V": "(必填) 设计体积 (m³)，无默认值", 
            "tank_type": "(必填) 罐类型: vertical/horizontal/spherical，无默认值", 
            "fluid": "(必填) 储存介质，无默认值", 
            "T": "(必填) 温度 (K)，无默认值",
            "P": "(选填) 操作压力 Pa，默认 101325(常压)，可选: 101325/200000/500000/1000000",
            "T_ambient": "(选填) 环境温度 K，默认 298.15(25°C)，可选: 278.15/293.15/298.15/303.15",
            "aspect_ratio": "(选填) 储罐的高径比或长径比 (L/D)，默认按罐型自动选取。如遇长径比校验失败可调整此值修复，【严禁】传入 force_aspect_ratio 等未知参数"
        },
        category="device", tags=["storage_tank", "tank", "储罐"])
    return registry

def register_to_universal(registry: UniversalToolRegistry) -> None:
    capability = ToolCapability(
        name="storage_tank_design",
        description="储罐完整设计。计算直径、高度、实际体积等参数。",
        category="device_design",
        parameters={
            "V": {"type": "number", "required": True, "description": "(必填) 设计体积 m³，无默认值"},
            "tank_type": {"type": "string", "required": True, "description": "(必填) 罐类型: vertical/horizontal/spherical，无默认值"},
            "fluid": {"type": "string", "required": True, "description": "(必填) 储存介质，无默认值"},
            "T": {"type": "number", "required": True, "description": "(必填) 温度 K，无默认值"},
            "P": {"type": "number", "required": False, "description": "(选填) 操作压力 Pa，默认 101325(常压)，可选: 101325/200000/500000/1000000"},
            "T_ambient": {"type": "number", "required": False, "description": "(选填) 环境温度 K，默认 298.15(25°C)，可选: 278.15/293.15/298.15/303.15"},
            "aspect_ratio": {"type": "number", "required": False, "description": "(选填) 目标高径比或长径比 (L/D)，默认按罐型自动选取。如遇到形状比例警告可调大或调小此参数修复，切勿传 force_aspect_ratio"},
            # ── 可选物性参数 ──
            "rho": {"type": "number", "required": False, "description": "(选填) 液体密度 kg/m³，给定后跳过物性查询"},
            "MW": {"type": "number", "required": False, "description": "(选填) 分子量 g/mol，给定后跳过物性查询"},
            "Pv": {"type": "number", "required": False, "description": "(选填) 饱和蒸气压 Pa，给定后跳过物性查询"},
        },
        examples=[{"description": "立式储罐", "input": {"V": 100, "tank_type": "vertical", "fluid": "water"}}],
        tags=["storage_tank", "储罐", "球罐"],
    )
    registry.register(capability, _full_design)