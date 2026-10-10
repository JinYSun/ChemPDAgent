"""
thermo_helper.py — 统一物性计算引擎 (thermo + CoolProp 双引擎版)

■ 引擎分工策略
  ┌─────────────────────────────────────────────────────────────────┐
  │  物性                │  主引擎          │  备用引擎              │
  ├─────────────────────────────────────────────────────────────────┤
  │  密度/摩尔体积       │  CoolProp(高压优先)│  thermo              │
  │  焓 / 熵 / 内能      │  CoolProp         │  thermo              │
  │  Gibbs / Helmholtz   │  CoolProp         │  thermo 推导          │
  │  Cp / Cv / γ         │  CoolProp         │  thermo              │
  │  声速                │  CoolProp (精确)  │  thermo EOS (估算)   │
  │  压缩因子 Z          │  CoolProp         │  thermo EOS          │
  │  逸度 / 逸度系数     │  CoolProp AS      │  thermo EOS          │
  │  等温/等熵压缩率     │  CoolProp         │  thermo EOS          │
  │  膨胀系数 β          │  CoolProp         │  thermo EOS          │
  │  JT 系数             │  CoolProp         │  thermo EOS          │
  │  蒸气压              │  thermo (更广覆盖)│  CoolProp            │
  │  升华压              │  thermo           │  —                   │
  │  汽化焓              │  thermo (更广覆盖)│  CoolProp            │
  │  动力粘度 (液/气)    │  CoolProp+相态检查│  thermo              │
  │  导热系数 (液/气)    │  CoolProp+相态检查│  thermo              │
  │  表面张力            │  CoolProp (两相区)│  thermo              │
  │  Prandtl 数          │  CoolProp         │  thermo              │
  │  扩散系数 (气)       │  Chapman-Enskog   │  FSG (高压用Takahashi)│
  │  扩散系数 (液,水溶液)│  Hayduk-Laudie    │  Wilke-Chang         │
  │  扩散系数 (液,有机)  │  Hayduk-Minhas    │  Wilke-Chang         │
  │  安全/环境属性       │  thermo           │  —                   │
  │  分子描述符          │  thermo           │  —                   │
  └─────────────────────────────────────────────────────────────────┘

■ CoolProp 覆盖流体：~124 种纯流体 (NIST REFPROP 级别精度)
  不在 CoolProp 中的物质自动回退到 thermo

■ 扩散系数补充 (CoolProp 不提供):
  气相: Chapman-Enskog 理论 (精度约±5%)；高压时附加 Takahashi 修正
  液相: 水溶液优先 Hayduk-Laudie，有机溶剂优先 Hayduk-Minhas，通用兜底 Wilke-Chang

■ 液相声速: CoolProp 直接输出 (精度 <0.1%), 完全替代 PR EOS 估算路线

"""
from __future__ import annotations

import math
import warnings
from functools import lru_cache
from typing import Dict, List, Optional, Tuple, Union

# ============================================================
# 引擎导入
# ============================================================
try:
    from thermo import Chemical, Mixture, PRMIX, SRKMIX
    from thermo.interaction_parameters import IPDB as _IPDB
    from thermo import UNIFAC_gammas
    IPDB = _IPDB
    UNIFAC_gammas_fn = UNIFAC_gammas
except ImportError as _e:
    raise ImportError(f"thermo 库未安装或版本不兼容: {_e}")

try:
    import CoolProp.CoolProp as _CP
    from CoolProp.CoolProp import PropsSI as _PropsSI
    from CoolProp.CoolProp import FluidsList as _FluidsList
    _COOLPROP_FLUIDS: set = set(f.lower() for f in _FluidsList())
    HAS_COOLPROP = True
except ImportError:
    HAS_COOLPROP = False
    _COOLPROP_FLUIDS = set()

R_GAS   = 8.314462618    # J/(mol·K)
N_AV    = 6.02214076e23  # Avogadro
K_B     = 1.380649e-23   # Boltzmann

USE_THERMO = True

# ============================================================
# 单位换算
# ============================================================
class Units:
    @staticmethod
    def Pa_to_bar(P):  return P / 1e5
    @staticmethod
    def Pa_to_atm(P):  return P / 101325.0
    @staticmethod
    def Pa_to_MPa(P):  return P / 1e6
    @staticmethod
    def Pa_to_psi(P):  return P / 6894.757
    @staticmethod
    def bar_to_Pa(P):  return P * 1e5
    @staticmethod
    def atm_to_Pa(P):  return P * 101325.0
    @staticmethod
    def MPa_to_Pa(P):  return P * 1e6
    @staticmethod
    def psi_to_Pa(P):  return P * 6894.757
    @staticmethod
    def K_to_C(T):     return T - 273.15
    @staticmethod
    def C_to_K(T):     return T + 273.15
    @staticmethod
    def K_to_F(T):     return (T - 273.15) * 9/5 + 32
    @staticmethod
    def F_to_K(T):     return (T - 32) * 5/9 + 273.15

# ============================================================
# 工具函数
# ============================================================
def _safe_float(val, allow_none: bool = True) -> Optional[float]:
    if val is None:
        return None
    try:
        v = float(val)
        return None if (math.isnan(v) or math.isinf(v)) else v
    except (TypeError, ValueError):
        return None


@lru_cache(maxsize=512)
def _get_chem_cached(name: str, T: float, P: float) -> Optional[Chemical]:
    try:
        return Chemical(name, T=T, P=P)
    except Exception:
        return None


def _get_chem(name: str, T: float = 298.15, P: float = 101325.0) -> Optional[Chemical]:
    return _get_chem_cached(name, T, P)


def clear_chem_cache() -> None:
    _get_chem_cached.cache_clear()


@lru_cache(maxsize=256)
def _get_coolprop_name(name: str) -> Optional[str]:
    """将任意化学品名称/CAS 转换为 CoolProp 可识别的流体名称。"""
    if not HAS_COOLPROP:
        return None
    chem = _get_chem_cached(name, 298.15, 101325.0)
    if chem is not None:
        cas = getattr(chem, 'CAS', None)
        if cas:
            try:
                _PropsSI('molar_mass', 'T', 300, 'P', 101325, cas)
                return cas
            except Exception:
                pass
    for candidate in [name, name.capitalize(), name.lower()]:
        try:
            _PropsSI('molar_mass', 'T', 300, 'P', 101325, candidate)
            return candidate
        except Exception:
            pass
    return None


def _cp_prop(prop: str, T: float, P: float, fluid_cp: str,
             fallback=None):
    """CoolProp PropsSI 封装，带异常处理"""
    if not HAS_COOLPROP or fluid_cp is None:
        return fallback
    try:
        val = _PropsSI(prop, 'T', T, 'P', P, fluid_cp)
        v = _safe_float(val)
        return v if v is not None else fallback
    except Exception:
        return fallback


def _cp_phase_index(T: float, P: float, fluid_cp: str) -> Optional[int]:
    """
    获取 CoolProp 相态索引。
    CoolProp 7.x 相态编码:
      0=liquid, 1=supercritical, 2=supercritical_gas,
      3=supercritical_liquid, 5=gas, 6=twophase,
      7=unknown, 8=not_imposed
    """
    if not HAS_COOLPROP or fluid_cp is None:
        return None
    try:
        return int(_PropsSI('Phase', 'T', T, 'P', P, fluid_cp))
    except Exception:
        return None


_CP_LIQUID_PHASES = frozenset({0, 3, 6})   # liquid, supercritical_liquid, twophase
_CP_GAS_PHASES    = frozenset({2, 5, 8})   # supercritical_gas, gas, not_imposed
# phase 1 (supercritical) 气液不可区分


def _cp_fugacity(T: float, P: float, fluid_cp: str) -> Tuple[Optional[float], Optional[float]]:
    """用 CoolProp AbstractState 计算逸度 (Pa) 和逸度系数 φ。"""
    if not HAS_COOLPROP or fluid_cp is None:
        return None, None
    try:
        AS = _CP.AbstractState('HEOS', fluid_cp)
        AS.update(_CP.PT_INPUTS, P, T)
        f = _safe_float(AS.fugacity(0))
        phi = _safe_float(AS.fugacity_coefficient(0))
        return f, phi
    except Exception:
        return None, None


def _cp_surface_tension(T: float, fluid_cp: str) -> Optional[float]:
    """CoolProp 表面张力 (N/m)，需在两相区内用 Q=0 调用"""
    if not HAS_COOLPROP or fluid_cp is None:
        return None
    try:
        val = _PropsSI('surface_tension', 'T', T, 'Q', 0, fluid_cp)
        v = _safe_float(val)
        return v if v is not None and v >= 0 else None
    except Exception:
        return None


# ============================================================
# 临界常数 & 相态判断
# ============================================================
def get_critical_constants(name: str) -> Optional[Tuple[float, float, float, float, float, float]]:
    """(Tc K, Pc Pa, Vc m³/mol, Zc -, omega -, MW g/mol)"""
    chem = _get_chem(name)
    if chem is None:
        return None
    try:
        return (float(chem.Tc), float(chem.Pc), float(chem.Vc),
                float(chem.Zc), float(chem.omega), float(chem.MW))
    except Exception:
        return None


def get_molecular_weight(name: str) -> Optional[float]:
    chem = _get_chem(name)
    return _safe_float(chem.MW) if chem else None


class PropertyRegime:
    IDEAL_GAS     = "ideal_gas"
    LOW_P         = "low_pressure"
    MODERATE_P    = "moderate_pressure"
    HIGH_P        = "high_pressure"
    NEAR_CRITICAL = "near_critical"
    SUPERCRITICAL = "supercritical"
    LIQUID        = "liquid"
    SOLID         = "solid"


def classify_regime(name: str, T: float, P: float) -> str:
    crit = get_critical_constants(name)
    if crit is None:
        return PropertyRegime.LOW_P if P < 5e5 else PropertyRegime.HIGH_P
    Tc, Pc = crit[0], crit[1]
    Tr = T / Tc if Tc > 0 else 1.0
    Pr = P / Pc if Pc > 0 else 0.0
    if Tr > 1.0 and Pr > 1.0:
        return PropertyRegime.SUPERCRITICAL
    if 0.9 < Tr < 1.1 and 0.8 < Pr < 1.2:
        return PropertyRegime.NEAR_CRITICAL
    if Tr < 0.85:
        return PropertyRegime.LIQUID
    if Pr < 0.01:
        return PropertyRegime.IDEAL_GAS
    if Pr < 0.1:
        return PropertyRegime.LOW_P
    if Pr < 0.5:
        return PropertyRegime.MODERATE_P
    return PropertyRegime.HIGH_P


def get_reduced_properties(name: str, T: float, P: float) -> Tuple[float, float]:
    crit = get_critical_constants(name)
    if crit is None:
        return 1.0, 0.0
    return (T / crit[0] if crit[0] > 0 else 1.0), (P / crit[1] if crit[1] > 0 else 0.0)


# ============================================================
# 蒸气压 & 升华压
# ============================================================
def calc_vapor_pressure(name: str, T: float) -> Optional[float]:
    """
    蒸气压 Psat (Pa)

    路由:
      1. thermo VaporPressure callable (覆盖最广，Antoine/Wagner/DIPPR/Yaws 自动路由)
      2. chem.Psat 属性
      3. CoolProp (仅 T < Tc 两相区，作为交叉验证备选)

    注：thermo 内部已实现 Antoine → DIPPR → Wagner → Lee-Kesler 自动路由，
    无需在此层重复实现各关联式切换。
    """
    chem = _get_chem(name, T=T)
    if chem is not None:
        vp_obj = getattr(chem, 'VaporPressure', None)
        if vp_obj is not None and callable(vp_obj):
            try:
                val = _safe_float(vp_obj(T))
                if val and val > 0:
                    return val
            except Exception:
                pass
        val = _safe_float(getattr(chem, 'Psat', None))
        if val and val > 0:
            return val
    cp_name = _get_coolprop_name(name)
    if cp_name:
        try:
            val = _safe_float(_PropsSI('P', 'T', T, 'Q', 0, cp_name))
            if val and val > 0:
                return val
        except Exception:
            pass
    return None


def calc_sublimation_pressure(name: str, T: float) -> Optional[float]:
    """升华压 Psub (Pa)，T < 三相点温度"""
    chem = _get_chem(name, T=T)
    if chem is None:
        return None
    sp_obj = getattr(chem, 'SublimationPressure', None)
    if sp_obj is not None and callable(sp_obj):
        try:
            val = _safe_float(sp_obj(T))
            if val and val > 0:
                return val
        except Exception:
            pass
    Tt = _safe_float(getattr(chem, 'Tt', None))
    Pt = _safe_float(getattr(chem, 'Pt', None))
    Hsub = _safe_float(getattr(chem, 'Hsubm', None))
    if Tt and Pt and Hsub and T > 0:
        try:
            return Pt * math.exp(-Hsub / R_GAS * (1.0 / T - 1.0 / Tt))
        except Exception:
            pass
    return None


# ============================================================
# 密度 / 摩尔体积
# ============================================================
def calc_liquid_density(name: str, T: float, P: float = 101325.0) -> Optional[float]:
    """
    液相质量密度 (kg/m³)

    路由:
      1. CoolProp HEOS — 仅当 CoolProp 判定为液相/两相/超临界时采用
      2. thermo VolumeLiquid
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        if phase is None or phase in _CP_LIQUID_PHASES or phase == 1:
            val = _cp_prop('D', T, P, cp_name)
            if val and val > 0:
                return val
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'rhol', None))
        if val and val > 0:
            return val
    return None


def calc_gas_density(name: str, T: float, P: float = 101325.0) -> Optional[float]:
    """
    气相质量密度 (kg/m³)

    路由:
      1. CoolProp HEOS — 仅当 CoolProp 判定为气相/超临界气时采用
      2. thermo (含 EOS 路由)
      3. 理想气体 (仅 Pr < 0.01 时)
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        if phase is None or phase in _CP_GAS_PHASES or phase == 1:
            val = _cp_prop('D', T, P, cp_name)
            if val and val > 0:
                return val
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'rhog', None))
        if val and val > 0:
            return val
        if classify_regime(name, T, P) == PropertyRegime.IDEAL_GAS:
            MW = _safe_float(chem.MW)
            if MW and T > 0:
                return P * (MW / 1000.0) / (R_GAS * T)
    return None


def calc_solid_density(name: str, T: float) -> Optional[float]:
    """固相质量密度 (kg/m³)"""
    chem = _get_chem(name, T=T)
    if chem is None:
        return None
    val = _safe_float(getattr(chem, 'rhos', None))
    return val if val and val > 0 else None


def calc_molar_volume(name: str, T: float, P: float,
                      phase: str = "liquid") -> Optional[float]:
    """摩尔体积 (m³/mol)"""
    MW = get_molecular_weight(name)
    if MW is None or MW <= 0:
        return None
    if phase == "solid":
        rho = calc_solid_density(name, T)
    elif phase in ("gas", "vapor"):
        rho = calc_gas_density(name, T, P)
    else:
        rho = calc_liquid_density(name, T, P)
    if rho and rho > 0:
        return (MW / 1000.0) / rho
    return None


# ============================================================
# 压缩因子 Z
# ============================================================
def calc_compressibility_factor(name: str, T: float, P: float,
                                phase: str = "gas") -> Optional[float]:
    """
    压缩因子 Z = PV/(nRT)

    路由:
      1. CoolProp 直接输出 (基于 Helmholtz EOS)
      2. thermo EOS — 多根问题处理：
         [FIX-4] 立方 EOS 在两相区有三个实根，本函数按 phase 参数选取
         对应根（gas→最大Z根，liquid→最小Z根），并在两相区条件下
         通过 _vprint 输出警告（两相区单相Z值无物理意义）。
      3. 维里方程 (低压气体)
      4. 理想气体 Z=1

    注意：若 T/P 处于两相区，单相 Z 无直接物理意义，需通过闪蒸计算确定
    各相的 Z 值后再使用。
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_prop('Z', T, P, cp_name)
        if val is not None:
            return val
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        eos = getattr(chem, 'eos', None)
        if eos:
            # [FIX-4] 检查两相区并给出警告
            Psat = calc_vapor_pressure(name, T)
            crit = get_critical_constants(name)
            if (Psat is not None and crit is not None
                    and T < crit[0] and abs(P - Psat) / max(Psat, 1.0) < 0.05):
                # 接近饱和线，Z 多根问题显著
                import os
                if os.getenv("FUNCS_VERBOSE", "0") == "1":
                    print(
                        f"  [calc_compressibility_factor] 警告：T={T:.2f}K, P={P:.0f}Pa "
                        f"接近饱和线 (Psat={Psat:.0f}Pa)，立方 EOS 存在多根，"
                        f"已按 phase='{phase}' 选根，两相区内单相 Z 无直接物理意义。"
                    )
            # 按请求相态选根
            z_attr = 'Z_g' if phase in ("gas", "vapor") else 'Z_l'
            z = _safe_float(getattr(eos, z_attr, None))
            if z is not None and z > 0:
                return z
        regime = classify_regime(name, T, P)
        if regime == PropertyRegime.IDEAL_GAS:
            return 1.0
        if regime == PropertyRegime.LOW_P and phase in ("gas", "vapor"):
            B = _safe_float(getattr(chem, 'Bvirial', None))
            if B is not None:
                return 1.0 + B * P / (R_GAS * T)
    return None


# ============================================================
# 逸度 & 逸度系数
# ============================================================
def calc_fugacity(name: str, T: float, P: float,
                  phase: str = "gas") -> Optional[float]:
    """逸度 f (Pa)"""
    cp_name = _get_coolprop_name(name)
    f, _ = _cp_fugacity(T, P, cp_name)
    if f is not None:
        return f
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        eos = getattr(chem, 'eos', None)
        if eos:
            attr = 'fugacity_g' if phase in ("gas","vapor") else 'fugacity_l'
            val = _safe_float(getattr(eos, attr, None))
            if val is not None:
                return val
    if classify_regime(name, T, P) == PropertyRegime.IDEAL_GAS:
        return P
    return None


def calc_fugacity_coefficient(name: str, T: float, P: float,
                              phase: str = "gas") -> Optional[float]:
    """逸度系数 φ = f/P"""
    cp_name = _get_coolprop_name(name)
    _, phi = _cp_fugacity(T, P, cp_name)
    if phi is not None:
        return phi
    f = calc_fugacity(name, T, P, phase)
    if f is not None and P > 0:
        return f / P
    if classify_regime(name, T, P) == PropertyRegime.IDEAL_GAS:
        return 1.0
    return None


# ============================================================
# 热容 Cp / Cv / γ
# ============================================================
def calc_heat_capacity(name: str, T: float, P: float = 101325.0,
                       phase: str = "liquid") -> Optional[float]:
    """
    定压热容 Cp (J/mol/K)

    路由:
      1. CoolProp Cpmolar（含相态检查，避免液相请求得到气相 Cp）
      2. thermo DIPPR 关联式

    [FIX-12] 近临界区 (Tr∈[0.9,1.1]) Cp 急剧发散，返回值时输出 verbose 警告。
    CoolProp 在近临界区处理比 thermo 可靠，但仍可能出现极大值。
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        cp_phase = _cp_phase_index(T, P, cp_name)
        phase_ok = True
        if cp_phase is not None:
            if phase in ("liquid",) and cp_phase not in _CP_LIQUID_PHASES and cp_phase != 1:
                phase_ok = False
            elif phase in ("gas", "vapor") and cp_phase not in _CP_GAS_PHASES and cp_phase != 1:
                phase_ok = False
        if phase_ok:
            val = _cp_prop('Cpmolar', T, P, cp_name)
            if val and val > 0:
                # [FIX-12] 近临界区 Cp 发散警告
                crit = get_critical_constants(name)
                if crit and crit[0] > 0:
                    Tr = T / crit[0]
                    if 0.9 < Tr < 1.1:
                        import os
                        if os.getenv("FUNCS_VERBOSE", "0") == "1":
                            print(
                                f"  [calc_heat_capacity] 警告：{name} 处于近临界区"
                                f" (Tr={Tr:.3f})，Cp={val:.1f} J/mol/K 可能急剧发散，"
                                f"请谨慎使用。"
                            )
                return val
    chem = _get_chem(name, T=T, P=P)
    if chem is None:
        return None
    if phase in ("gas", "vapor"):
        val = _safe_float(getattr(chem, 'Cpgm', None))
    elif phase == "solid":
        val = _safe_float(getattr(chem, 'Cpsm', None))
    else:
        val = _safe_float(getattr(chem, 'Cplm', None))
    return val if val and val > 0 else None


def calc_Cv(name: str, T: float, P: float = 101325.0,
            phase: str = "gas") -> Optional[float]:
    """定容热容 Cv (J/mol/K)"""
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_prop('Cvmolar', T, P, cp_name)
        if val and val > 0:
            return val
    chem = _get_chem(name, T=T, P=P)
    if chem is None:
        return None
    val = _safe_float(getattr(chem, 'Cvgm', None))
    if val and val > 0:
        return val
    if phase in ("gas", "vapor"):
        cp = calc_heat_capacity(name, T, P, "gas")
        return (cp - R_GAS) if cp else None
    eos = getattr(chem, 'eos', None)
    cp_l = calc_heat_capacity(name, T, P, "liquid")
    if eos and cp_l:
        V_l   = _safe_float(getattr(eos, 'V_l', None))
        dPdV  = _safe_float(getattr(eos, 'dP_dV_l', None))
        beta  = _safe_float(getattr(chem, 'isobaric_expansion_l', None))
        if V_l and dPdV and dPdV != 0 and beta is not None:
            kT = -1.0 / (V_l * dPdV)
            if kT > 0:
                cv_l = cp_l - T * V_l * beta**2 / kT
                return cv_l if cv_l > 0 else None
    return None


def calc_isentropic_exponent(name: str, T: float, P: float = 101325.0,
                             phase: str = "gas") -> Optional[float]:
    """绝热指数 γ = Cp/Cv"""
    if phase in ("gas", "vapor"):
        cp_name = _get_coolprop_name(name)
        if cp_name:
            val = _cp_prop('isentropic_expansion_coefficient', T, P, cp_name)
            if val and 1.0 < val < 20.0:
                return val
        cp = calc_heat_capacity(name, T, P, "gas")
        cv = calc_Cv(name, T, P, "gas")
        if cp and cv and cv > 0:
            return cp / cv
        return None
    else:
        cp = calc_heat_capacity(name, T, P, "liquid")
        cv = calc_Cv(name, T, P, "liquid")
        if cp and cv and cv > 0:
            return cp / cv
        return None


# ============================================================
# 焓 / 熵 / 内能 / Gibbs / Helmholtz
# ============================================================
def calc_enthalpy(name: str, T: float, P: float = 101325.0,
                  molar: bool = True) -> Optional[float]:
    """相对焓 H (J/mol 或 J/kg)"""
    cp_name = _get_coolprop_name(name)
    if cp_name:
        key = 'Hmolar' if molar else 'H'
        val = _cp_prop(key, T, P, cp_name)
        if val is not None:
            return val
    chem = _get_chem(name, T=T, P=P)
    if chem is None:
        return None
    val = _safe_float(getattr(chem, 'Hm', None))
    if val is not None and not molar:
        MW = _safe_float(chem.MW)
        return val / (MW / 1000.0) if MW else None
    return val


def calc_entropy(name: str, T: float, P: float = 101325.0,
                 molar: bool = True) -> Optional[float]:
    """相对熵 S (J/mol/K 或 J/kg/K)"""
    cp_name = _get_coolprop_name(name)
    if cp_name:
        key = 'Smolar' if molar else 'S'
        val = _cp_prop(key, T, P, cp_name)
        if val is not None:
            return val
    chem = _get_chem(name, T=T, P=P)
    if chem is None:
        return None
    return _safe_float(getattr(chem, 'Sm', None))


def calc_internal_energy(name: str, T: float, P: float = 101325.0,
                         molar: bool = True) -> Optional[float]:
    """内能 U (J/mol 或 J/kg)"""
    cp_name = _get_coolprop_name(name)
    if cp_name:
        key = 'Umolar' if molar else 'U'
        val = _cp_prop(key, T, P, cp_name)
        if val is not None:
            return val
    H = calc_enthalpy(name, T, P, molar=True)
    MW = get_molecular_weight(name)
    rho = calc_gas_density(name, T, P) or calc_liquid_density(name, T, P)
    if H is not None and MW and rho and rho > 0:
        Vm = (MW / 1000.0) / rho
        U = H - P * Vm
        return U if molar else U / (MW / 1000.0)
    return None


def calc_gibbs_energy(name: str, T: float, P: float = 101325.0,
                      molar: bool = True) -> Optional[float]:
    """Gibbs 自由能 G = H - TS (J/mol 或 J/kg)"""
    cp_name = _get_coolprop_name(name)
    if cp_name:
        key = 'Gmolar' if molar else 'G'
        val = _cp_prop(key, T, P, cp_name)
        if val is not None:
            return val
    H = calc_enthalpy(name, T, P, molar)
    S = calc_entropy(name, T, P, molar)
    if H is not None and S is not None:
        return H - T * S
    return None


def calc_helmholtz_energy(name: str, T: float, P: float = 101325.0,
                          molar: bool = True) -> Optional[float]:
    """Helmholtz 自由能 A = U - TS (J/mol 或 J/kg)"""
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val_mass = _cp_prop('A', T, P, cp_name)
        if val_mass is not None:
            if molar:
                MW = get_molecular_weight(name)
                if MW:
                    return val_mass * (MW / 1000.0)
            else:
                return val_mass
    U = calc_internal_energy(name, T, P, molar)
    S = calc_entropy(name, T, P, molar)
    if U is not None and S is not None:
        return U - T * S
    return None


# ============================================================
# 相变焓
# ============================================================
def calc_enthalpy_vaporization(name: str, T: float) -> Optional[float]:
    """
    汽化焓 ΔHvap (J/mol)

    路由:
      1. thermo EnthalpyVaporization (覆盖最广，Watson/DIPPR 106)
      2. CoolProp: ΔHvap = H_gas(T,Q=1) - H_liq(T,Q=0)
      T ≥ Tc → 返回 None
    """
    crit = get_critical_constants(name)
    if crit and T >= crit[0]:
        return None
    chem = _get_chem(name, T=T)
    if chem is not None:
        val = _safe_float(getattr(chem, 'Hvapm', None))
        if val and val > 0:
            return val
    cp_name = _get_coolprop_name(name)
    if cp_name:
        try:
            Hg = _PropsSI('Hmolar', 'T', T, 'Q', 1, cp_name)
            Hl = _PropsSI('Hmolar', 'T', T, 'Q', 0, cp_name)
            dH = _safe_float(Hg - Hl) if Hg and Hl else None
            if dH and dH > 0:
                return dH
        except Exception:
            pass
    return None


def calc_enthalpy_fusion(name: str) -> Optional[float]:
    """融化焓 ΔHfus (J/mol)"""
    chem = _get_chem(name)
    return _safe_float(getattr(chem, 'Hfusm', None)) if chem else None


def calc_enthalpy_sublimation(name: str) -> Optional[float]:
    """升华焓 ΔHsub (J/mol)"""
    chem = _get_chem(name)
    return _safe_float(getattr(chem, 'Hsubm', None)) if chem else None


# ============================================================
# 生成焓 / Gibbs 生成能 / 燃烧焓
# ============================================================
def get_formation_enthalpy(
    fluid: Union[str, List[str]],
    T: float = 298.15, P: float = 101325.0,
    phase: str = "gas", z: Optional[List[float]] = None,
) -> Optional[float]:
    """标准摩尔生成焓 ΔHf° (J/mol)"""
    if isinstance(fluid, list):
        if z is None:
            z = [1.0 / len(fluid)] * len(fluid)
        vals = []
        for n in fluid:
            chem = _get_chem(n, T, P)
            if chem is None:
                return None
            hf = _safe_float(getattr(chem, 'Hfgm' if phase == "gas" else 'Hfm', None))
            if hf is None:
                hf = _safe_float(getattr(chem, 'Hfm' if phase == "gas" else 'Hfgm', None))
            if hf is None:
                return None
            vals.append(hf)
        return sum(z[i] * vals[i] for i in range(len(fluid)))
    else:
        chem = _get_chem(fluid, T, P)
        if chem is None:
            return None
        hf = _safe_float(getattr(chem, 'Hfgm' if phase == "gas" else 'Hfm', None))
        if hf is None:
            hf = _safe_float(getattr(chem, 'Hfm' if phase == "gas" else 'Hfgm', None))
        return hf


def get_formation_gibbs(name: str, phase: str = "gas") -> Optional[float]:
    """标准 Gibbs 生成自由能 ΔGf° (J/mol)"""
    chem = _get_chem(name)
    if chem is None:
        return None
    return _safe_float(getattr(chem, 'Gfgm' if phase == "gas" else 'Gfm', None))


def get_combustion_enthalpy(name: str, phase: str = "gas",
                            higher: bool = True) -> Optional[float]:
    """燃烧焓 ΔHc (J/mol)，负值=放热"""
    chem = _get_chem(name)
    if chem is None:
        return None
    attr = ('Hcgm' if higher else 'Hcgm_lower') if phase in ("gas","vapor") else \
           ('Hcm'  if higher else 'Hcm_lower')
    return _safe_float(getattr(chem, attr, None))


# ============================================================
# 粘度
# [FIX-2] 液气粘度均增加 CoolProp 相态一致性检查
# ============================================================
def calc_viscosity_liquid(name: str, T: float, P: float = 101325.0) -> Optional[float]:
    """
    液相动力粘度 μ (Pa·s)

    路由:
      1. CoolProp viscosity — 仅当 CoolProp 判定为液相/两相时采用
         [FIX-2] 增加相态检查：若 CoolProp 判定为气相，跳过并回退 thermo，
         避免气相区返回气相粘度作为液相粘度使用。
      2. thermo DIPPR 101 / Vogel
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        # [FIX-2] 仅在液相/两相/超临界区使用 CoolProp 粘度作为液相粘度
        if phase is None or phase in _CP_LIQUID_PHASES or phase == 1:
            val = _cp_prop('viscosity', T, P, cp_name)
            if val and val > 0:
                return val
        # 若 CoolProp 判定为气相，跳过，回退 thermo
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'mul', None))
        if val and val > 0:
            return val
    return None


def calc_viscosity_gas(name: str, T: float, P: float = 101325.0) -> Optional[float]:
    """
    气相动力粘度 μ (Pa·s)

    路由:
      1. CoolProp viscosity — 仅当 CoolProp 判定为气相/超临界气时采用
         [FIX-2] 增加相态检查：若 CoolProp 判定为液相，跳过并回退 thermo。
      2. thermo Chapman-Enskog / Reichenberg
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        # [FIX-2] 仅在气相/超临界区使用 CoolProp 粘度作为气相粘度
        if phase is None or phase in _CP_GAS_PHASES or phase == 1:
            val = _cp_prop('viscosity', T, P, cp_name)
            if val and val > 0:
                return val
        # 若 CoolProp 判定为液相，跳过，回退 thermo
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'mug', None))
        if val and val > 0:
            return val
    return None


def calc_kinematic_viscosity(name: str, T: float, P: float = 101325.0,
                             phase: str = "liquid") -> Optional[float]:
    """运动粘度 ν = μ/ρ (m²/s)"""
    if phase in ("gas", "vapor"):
        mu = calc_viscosity_gas(name, T, P)
        rho = calc_gas_density(name, T, P)
    else:
        mu = calc_viscosity_liquid(name, T, P)
        rho = calc_liquid_density(name, T, P)
    return (mu / rho) if mu and rho and rho > 0 else None


# ============================================================
# 导热系数
# [FIX-3] 液气导热系数均增加 CoolProp 相态一致性检查
# ============================================================
def calc_thermal_conductivity_liquid(name: str, T: float,
                                     P: float = 101325.0) -> Optional[float]:
    """
    液相导热系数 λ (W/m/K)

    路由:
      1. CoolProp conductivity — 仅当 CoolProp 判定为液相/两相时采用
         [FIX-3] 增加相态检查，气相区 CoolProp 返回气相导热系数，跳过回退 thermo。
      2. thermo DIPPR 100 / Sato-Riedel
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        # [FIX-3] 仅液相/两相/超临界使用
        if phase is None or phase in _CP_LIQUID_PHASES or phase == 1:
            val = _cp_prop('conductivity', T, P, cp_name)
            if val and val > 0:
                return val
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'kl', None))
        if val and val > 0:
            return val
    return None


def calc_thermal_conductivity_gas(name: str, T: float,
                                  P: float = 101325.0) -> Optional[float]:
    """
    气相导热系数 λ (W/m/K)

    路由:
      1. CoolProp conductivity — 仅当 CoolProp 判定为气相/超临界气时采用
         [FIX-3] 增加相态检查，液相区跳过回退 thermo。
      2. thermo Eucken 法
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        # [FIX-3] 仅气相/超临界气使用
        if phase is None or phase in _CP_GAS_PHASES or phase == 1:
            val = _cp_prop('conductivity', T, P, cp_name)
            if val and val > 0:
                return val
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'kg', None))
        if val and val > 0:
            return val
    return None


# ============================================================
# 表面张力
# [FIX-9] T >= Tc 时返回 None 而非 0.0
# ============================================================
def calc_surface_tension(name: str, T: float) -> Optional[float]:
    """
    表面张力 σ (N/m)

    路由:
      1. CoolProp surface_tension (两相区界面，Q=0 调用)
      2. thermo DIPPR 106 / Macleod-Sugden

    [FIX-9] T >= Tc 时返回 None（超临界流体无气液界面，表面张力物理上不存在）。
    原代码返回 0.0 会导致 Macleod-Sugden 等下游计算中得到错误的零值而非跳过。
    """
    crit = get_critical_constants(name)
    if crit and T >= crit[0]:
        # [FIX-9] 超临界无气液界面，返回 None 而非 0.0
        return None
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_surface_tension(T, cp_name)
        if val is not None and val >= 0:
            return val
    chem = _get_chem(name, T=T)
    if chem is not None:
        val = _safe_float(getattr(chem, 'sigma', None))
        if val is not None and val >= 0:
            return val
    return None


# ============================================================
# 推导物性
# ============================================================
def calc_prandtl_number(name: str, T: float, P: float = 101325.0,
                        phase: str = "liquid") -> Optional[float]:
    """Prandtl 数 Pr = μ·Cp_mass/λ"""
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_prop('Prandtl', T, P, cp_name)
        if val and val > 0:
            return val
    if phase in ("gas", "vapor"):
        mu = calc_viscosity_gas(name, T, P)
        lam = calc_thermal_conductivity_gas(name, T, P)
        cp = calc_heat_capacity(name, T, P, "gas")
    else:
        mu = calc_viscosity_liquid(name, T, P)
        lam = calc_thermal_conductivity_liquid(name, T, P)
        cp = calc_heat_capacity(name, T, P, "liquid")
    MW = get_molecular_weight(name)
    if mu and lam and cp and MW and lam > 0:
        return mu * (cp / (MW / 1000.0)) / lam
    return None


def calc_thermal_diffusivity(name: str, T: float, P: float = 101325.0,
                             phase: str = "liquid") -> Optional[float]:
    """热扩散率 α = λ/(ρ·Cp_mass) (m²/s)"""
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        attr = 'alphag' if phase in ("gas","vapor") else 'alphal'
        val = _safe_float(getattr(chem, attr, None))
        if val and val > 0:
            return val
    if phase in ("gas", "vapor"):
        lam = calc_thermal_conductivity_gas(name, T, P)
        rho = calc_gas_density(name, T, P)
        cp  = calc_heat_capacity(name, T, P, "gas")
    else:
        lam = calc_thermal_conductivity_liquid(name, T, P)
        rho = calc_liquid_density(name, T, P)
        cp  = calc_heat_capacity(name, T, P, "liquid")
    MW = get_molecular_weight(name)
    if lam and rho and cp and MW and MW > 0:
        return lam / (rho * (cp / (MW / 1000.0)))
    return None


def calc_isothermal_compressibility(name: str, T: float, P: float,
                                    phase: str = "liquid") -> Optional[float]:
    """等温压缩率 κ_T = -(1/V)(∂V/∂P)_T (1/Pa)"""
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_prop('isothermal_compressibility', T, P, cp_name)
        if val and val > 0:
            return val
    if classify_regime(name, T, P) == PropertyRegime.IDEAL_GAS:
        return 1.0 / P if P > 0 else None
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        eos = getattr(chem, 'eos', None)
        if eos:
            V = _safe_float(getattr(eos, 'V_g' if phase in ("gas","vapor") else 'V_l', None))
            dPdV = _safe_float(getattr(eos, 'dP_dV_g' if phase in ("gas","vapor") else 'dP_dV_l', None))
            if V and dPdV and dPdV != 0:
                return -1.0 / (V * dPdV)
    return None


def calc_isentropic_compressibility(name: str, T: float, P: float,
                                    phase: str = "liquid") -> Optional[float]:
    """等熵压缩率 κ_s = κ_T / γ (1/Pa)"""
    kT = calc_isothermal_compressibility(name, T, P, phase)
    gamma = calc_isentropic_exponent(name, T, P, phase)
    if kT and gamma and gamma > 0:
        return kT / gamma
    return None


def calc_isobaric_expansion(name: str, T: float, P: float,
                            phase: str = "liquid") -> Optional[float]:
    """等压热膨胀系数 β = (1/V)(∂V/∂T)_P (1/K)"""
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_prop('isobaric_expansion_coefficient', T, P, cp_name)
        if val is not None:
            return val
    if classify_regime(name, T, P) == PropertyRegime.IDEAL_GAS:
        return 1.0 / T if T > 0 else None
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        attr = 'isobaric_expansion_g' if phase in ("gas","vapor") else 'isobaric_expansion_l'
        return _safe_float(getattr(chem, attr, None))
    return None


def calc_joule_thomson_coefficient(name: str, T: float, P: float,
                                   phase: str = "gas") -> Optional[float]:
    """Joule-Thomson 系数 μ_JT = (∂T/∂P)_H (K/Pa)"""
    if classify_regime(name, T, P) == PropertyRegime.IDEAL_GAS:
        return 0.0
    chem = _get_chem(name, T=T, P=P)
    if chem is None:
        return None
    return _safe_float(getattr(chem, 'JTg' if phase in ("gas","vapor") else 'JTl', None))


def calc_speed_of_sound(name: str, T: float, P: float = 101325.0,
                        phase: str = "gas") -> Optional[float]:
    """声速 c (m/s)"""
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_prop('speed_of_sound', T, P, cp_name)
        if val and val > 0:
            return val
    MW = get_molecular_weight(name)
    if MW and phase in ("gas", "vapor"):
        gamma = calc_isentropic_exponent(name, T, P, "gas")
        Z = calc_compressibility_factor(name, T, P, "gas") or 1.0
        if gamma:
            return math.sqrt(gamma * Z * R_GAS * T / (MW / 1000.0))
    return None


def calc_second_virial_coefficient(name: str, T: float) -> Optional[float]:
    """第二维里系数 B (m³/mol)"""
    chem = _get_chem(name, T=T)
    return _safe_float(getattr(chem, 'Bvirial', None)) if chem else None


def calc_parachor(name: str, T: float) -> Optional[float]:
    """Parachor，用于表面张力混合规则"""
    chem = _get_chem(name, T=T)
    return _safe_float(getattr(chem, 'Parachor', None)) if chem else None


def calc_permittivity(name: str, T: float) -> Optional[float]:
    """相对介电常数（液相）"""
    chem = _get_chem(name, T=T)
    return _safe_float(getattr(chem, 'permittivity', None)) if chem else None


# ============================================================
# 扩散系数
# ============================================================

def calc_diffusion_gas(
    name_A: str, name_B: str, T: float, P: float = 101325.0
) -> Optional[float]:
    """
    气相二元扩散系数 D_AB (m²/s) — Chapman-Enskog 理论 + 高压 Takahashi 修正

    低压路由 (P ≤ 10 bar, 建议范围):
        Chapman-Enskog 公式 (Reid et al., 5th Ed. Eq.11-3.2)
        D_AB [cm²/s] = 1.8583e-3 * T^1.5 * √(1/MA + 1/MB)
                       / (P_bar * σ_AB² * Ω_D)
        参数: σ in Å, P in bar, M in g/mol, T in K

    高压修正 (P > 10 bar):
        [FIX-1] 新增 Takahashi (1974) 压力修正：
        D * P / (D * P)_low = f(Tr, Pr)  [基于对应态]
        近似式: (D*P)_HP / (D*P)_LP ≈ exp(-0.0152 * Pr^1.1)
        仅在 Lennard-Jones 参数可用时尝试 Chapman-Enskog；
        若参数缺失（大分子有机物），直接返回 None 并建议使用 FSG。

    精度:
        低压: ±5% (非极性), ±10% (极性)
        高压修正后: ±15%（工程估算）

    [FIX-1] 删除原代码中重复的废弃 SI 单位路径（D_AB 变量被丢弃），
    仅保留 bar/Å 单位路径，并补充高压修正。
    """
    chem_A = _get_chem(name_A)
    chem_B = _get_chem(name_B)
    if chem_A is None or chem_B is None:
        return None
    MA = _safe_float(chem_A.MW)
    MB = _safe_float(chem_B.MW)
    sigma_A = _safe_float(getattr(chem_A, 'molecular_diameter', None))  # Å
    sigma_B = _safe_float(getattr(chem_B, 'molecular_diameter', None))
    eps_A   = _safe_float(getattr(chem_A, 'Stockmayer', None))  # ε/k_B in K
    eps_B   = _safe_float(getattr(chem_B, 'Stockmayer', None))
    if not all([MA, MB, sigma_A, sigma_B, eps_A, eps_B]):
        return None

    sigma_AB = (sigma_A + sigma_B) / 2.0   # Å
    eps_AB   = math.sqrt(eps_A * eps_B)    # K
    T_star   = T / eps_AB

    # Neufeld-Janzen-Aziz 碰撞积分 Ω_D
    A, B, C, D, E, F = 1.06036, 0.15610, 0.19300, 0.47635, 1.03587, 1.52996
    G, H = 1.76474, 3.89411
    try:
        Omega_D = (A / T_star**B + C / math.exp(D * T_star)
                   + E / math.exp(F * T_star) + G / math.exp(H * T_star))
    except Exception:
        return None

    # [FIX-1] 只保留 bar/Å 单位路径（原代码有两条路径，前一条的 D_AB 被丢弃）
    # Chapman-Enskog: P in bar, σ in Å, M in g/mol, T in K → D in cm²/s
    P_bar = P / 1e5
    D_low_cm2_s = (1.8583e-3 * T**1.5 * math.sqrt(1.0/MA + 1.0/MB)
                   / (P_bar * sigma_AB**2 * Omega_D))
    D_low_m2_s = D_low_cm2_s * 1e-4  # cm²/s → m²/s

    # [FIX-1] 高压 Takahashi 修正（P > 10 bar）
    # 参考: Takahashi, S. (1974), J. Chem. Eng. Japan, 7(6), 417–420
    # 对应态修正：(D*P)_HP / (D*P)_LP ≈ exp(-α * Pr^β)
    # 对非极性双组分混合物 α≈0.0152, β≈1.1（工程近似）
    # 使用混合物伪临界参数（Kay 规则）
    if P_bar > 10.0:
        crit_A = get_critical_constants(name_A)
        crit_B = get_critical_constants(name_B)
        if crit_A and crit_B:
            # Kay 规则伪临界：简单算术平均（1:1 混合，工程近似）
            Tc_mix = 0.5 * (crit_A[0] + crit_B[0])
            Pc_mix = 0.5 * (crit_A[1] + crit_B[1])
            if Tc_mix > 0 and Pc_mix > 0:
                Tr_mix = T / Tc_mix
                Pr_mix = P / Pc_mix
                # Takahashi 近似修正因子
                try:
                    correction = math.exp(-0.0152 * Pr_mix**1.1)
                    # 同时考虑密度增加对 D*P 乘积的影响
                    # D_HP ≈ D_low * correction（因为 D_HP * P_HP = D_low * P_low * factor）
                    # 此处 correction 已含压力比，直接返回 D_HP
                    D_low_m2_s = D_low_m2_s * correction
                except Exception:
                    pass  # 修正失败，退回低压值

    return D_low_m2_s


def _get_Vb_cm3_mol(solute: str) -> Optional[float]:
    """
    溶质在正常沸点下的液相摩尔体积 V_b (cm³/mol)。

    [FIX-10] 单位注释补充：
      CoolProp 返回 rho [kg/m³]，MW [g/mol]
      Vb [cm³/mol] = (MW [g/mol] / 1000) [kg/mol] / rho [kg/m³] × 1e6 [cm³/m³]
                   = MW [g/mol] / rho [kg/m³] × 1000
      即 MW / rho * 1000（数值上与原代码相同，此处澄清推导）

    路由:
      1. CoolProp 两相边界液相密度（Q=0）
      2. thermo 在 Tb-δT 处的液相密度（避免边界判定为气相）
    """
    chem_A = _get_chem(solute)
    if chem_A is None:
        return None
    MW_A = _safe_float(getattr(chem_A, 'MW', None))  # g/mol
    Tb_A = _safe_float(getattr(chem_A, 'Tb', None))  # K
    if not MW_A or not Tb_A:
        return None

    # 方法1: CoolProp 两相边界液相密度
    cp_name = _get_coolprop_name(solute)
    if cp_name and HAS_COOLPROP:
        try:
            rho_tb = _safe_float(_PropsSI('D', 'T', Tb_A, 'Q', 0, cp_name))
            if rho_tb and rho_tb > 0:
                # [FIX-10] 单位：(g/mol) / (kg/m³) × 1000 = cm³/mol
                return MW_A / rho_tb * 1000.0
        except Exception:
            pass

    # 方法2: thermo 在 Tb-δT 处取液相密度
    for dT in [1.0, 2.0, 5.0]:
        chem_near = _get_chem(solute, T=Tb_A - dT)
        if chem_near:
            rho_tb = _safe_float(getattr(chem_near, 'rhol', None))
            if rho_tb and rho_tb > 0:
                return MW_A / rho_tb * 1000.0
    return None


def calc_diffusion_liquid_wilke_chang(
    solute: str, solvent: str, T: float
) -> Optional[float]:
    """
    液相无限稀释扩散系数 D°_AB (m²/s) — Wilke-Chang 方程（通用兜底）

    公式 (Wilke-Chang, AIChE J. 1955):
        D°_AB = 7.4e-8 * (φ·M_B)^0.5 * T / (μ_B [cP] · V_A^0.6 [cm³/mol])

    溶剂缔合因子 φ（[FIX-6] 重写，消除原代码 `or` 链的优先级 bug）：
      水=2.6, 甲醇=1.9, 乙醇=1.5, 其他=1.0
      原代码：`phi = dict1.get(cas, 0) or dict2.get(name, 1.0)`
      当 CAS 未命中时 dict1 返回 0，`0 or dict2.get(...)` 才走名称匹配——
      逻辑上碰巧有效，但当 CAS 命中却不在字典时返回 0，导致 phi=0 进而除零。
      [FIX-6] 改为三段显式逻辑：CAS优先 → 名称次之 → 默认1.0。

    适用范围：普通小分子液体，无限稀释。
    精度：±10%（水为溶剂），±25%（有机溶剂）。
    对于水溶液，优先使用 calc_diffusion_liquid_hayduk_laudie（精度更优）。
    """
    chem_B = _get_chem(solvent, T=T)
    if chem_B is None:
        return None
    M_B  = _safe_float(getattr(chem_B, 'MW', None))
    mu_B = calc_viscosity_liquid(solvent, T)
    if not M_B or not mu_B or mu_B <= 0:
        return None
    mu_B_cP = mu_B * 1000.0  # Pa·s → cP

    V_A = _get_Vb_cm3_mol(solute)
    if V_A is None or V_A <= 0:
        return None

    # [FIX-6] 缔合因子 φ：三段显式逻辑，消除 `or` 链 bug
    _CAS_ASSOC = {
        '7732-18-5': 2.6,  # 水
        '67-56-1':   1.9,  # 甲醇
        '64-17-5':   1.5,  # 乙醇
    }
    _NAME_ASSOC = {
        'water': 2.6, 'h2o': 2.6,
        'methanol': 1.9, 'methyl alcohol': 1.9,
        'ethanol': 1.5, 'ethyl alcohol': 1.5,
    }
    cas_B = getattr(chem_B, 'CAS', None) or ''
    if cas_B in _CAS_ASSOC:
        phi = _CAS_ASSOC[cas_B]
    else:
        phi = _NAME_ASSOC.get(solvent.lower().strip(), 1.0)

    D_cm2_s = 7.4e-8 * math.sqrt(phi * M_B) * T / (mu_B_cP * V_A**0.6)
    return D_cm2_s * 1e-4  # cm²/s → m²/s


def calc_diffusion_liquid_hayduk_laudie(
    solute: str, solvent_water: str, T: float
) -> Optional[float]:
    """
    液相无限稀释扩散系数 D°_AB (m²/s) — Hayduk-Laudie 方程（水溶液专用）

    [FIX-7] 新增函数。建议文档指出：对于溶质溶于水的体系，Hayduk-Laudie
    精度优于 Wilke-Chang，应优先使用本函数。

    公式 (Hayduk & Laudie, AIChE J. 1974):
        D°_AB = 13.26e-9 / (μ_water^1.14 [cP] · V_A^0.589 [cm³/mol])

    适用范围：溶质溶于水，无限稀释，T 在 278–323 K 范围内精度最优。
    精度：±10%（水溶液），优于 Wilke-Chang（±10-20%）。

    参数:
        solute       — 溶质名称（thermo 可识别）
        solvent_water — 溶剂名称（须为水，若非水则返回 None 并建议使用其他方法）
        T            — 温度 (K)
    """
    # 验证溶剂是否为水
    chem_solvent = _get_chem(solvent_water, T=T)
    if chem_solvent is None:
        return None
    _WATER_CAS = {'7732-18-5'}
    _WATER_NAMES = {'water', 'h2o', '水'}
    cas_solvent = getattr(chem_solvent, 'CAS', '') or ''
    is_water = (cas_solvent in _WATER_CAS
                or solvent_water.lower().strip() in _WATER_NAMES)
    if not is_water:
        # 非水溶剂，Hayduk-Laudie 不适用
        import os
        if os.getenv("FUNCS_VERBOSE", "0") == "1":
            print(
                f"  [calc_diffusion_liquid_hayduk_laudie] 溶剂 '{solvent_water}' 非水，"
                f"Hayduk-Laudie 仅适用于水溶液，返回 None。"
                f"请使用 calc_diffusion_liquid_hayduk_minhas 或 calc_diffusion_liquid_wilke_chang。"
            )
        return None

    mu_water = calc_viscosity_liquid(solvent_water, T)
    if not mu_water or mu_water <= 0:
        return None
    mu_water_cP = mu_water * 1000.0  # Pa·s → cP

    V_A = _get_Vb_cm3_mol(solute)
    if V_A is None or V_A <= 0:
        return None

    # Hayduk-Laudie 公式
    D_m2_s = 13.26e-9 / (mu_water_cP**1.14 * V_A**0.589)
    return D_m2_s  # 直接 m²/s


def calc_diffusion_liquid_hayduk_minhas(
    solute: str, solvent: str, T: float
) -> Optional[float]:
    """
    液相无限稀释扩散系数 D°_AB (m²/s) — Hayduk-Minhas 方程（非水有机溶剂优先）

    [FIX-7] 本函数适用于非水有机溶剂（烃类、酮、酯等）。
    若检测到溶剂为水，输出 verbose 警告并建议改用 calc_diffusion_liquid_hayduk_laudie。

    非水系统公式 (Hayduk & Minhas, Can. J. Chem. Eng. 1982):
        D°_AB = 13.3e-8 * T^1.47 * μ_B^(10.2/V_A - 0.791) / V_A^0.71

    精度：±10%（非极性有机溶剂），优于 Wilke-Chang。
    不适合：水溶液体系（请使用 calc_diffusion_liquid_hayduk_laudie）。
    """
    # [FIX-7] 检测是否为水溶液体系
    chem_solvent = _get_chem(solvent, T=T)
    if chem_solvent is not None:
        _WATER_CAS = {'7732-18-5'}
        _WATER_NAMES = {'water', 'h2o', '水'}
        cas_solvent = getattr(chem_solvent, 'CAS', '') or ''
        is_water = (cas_solvent in _WATER_CAS
                    or solvent.lower().strip() in _WATER_NAMES)
        if is_water:
            import os
            if os.getenv("FUNCS_VERBOSE", "0") == "1":
                print(
                    f"  [calc_diffusion_liquid_hayduk_minhas] 检测到溶剂为水，"
                    f"Hayduk-Minhas 非水公式对水溶液精度较差，"
                    f"建议改用 calc_diffusion_liquid_hayduk_laudie。"
                    f"仍继续计算，结果置信度偏低。"
                )

    mu_B = calc_viscosity_liquid(solvent, T)
    if not mu_B or mu_B <= 0:
        return None
    mu_B_cP = mu_B * 1000.0

    V_A = _get_Vb_cm3_mol(solute)
    if V_A is None or V_A <= 0:
        return None

    try:
        exp_mu = 10.2 / V_A - 0.791
        D_cm2_s = 13.3e-8 * T**1.47 * mu_B_cP**exp_mu / V_A**0.71
    except Exception:
        return None
    return D_cm2_s * 1e-4  # cm²/s → m²/s


def calc_diffusion_self(name: str, T: float, P: float = 101325.0) -> Optional[float]:
    """自扩散系数 D_self (m²/s) — 气相近似（A-A Chapman-Enskog）"""
    return calc_diffusion_gas(name, name, T, P)


# ============================================================
# 独立属性函数
# ============================================================
def get_triple_point(name: str) -> Optional[Tuple[float, float]]:
    """三相点 (Tt K, Pt Pa)"""
    chem = _get_chem(name)
    if chem is None:
        return None
    Tt = _safe_float(getattr(chem, 'Tt', None))
    Pt = _safe_float(getattr(chem, 'Pt', None))
    return (Tt, Pt) if Tt and Pt else None


def get_normal_boiling_point(name: str) -> Optional[float]:
    """正常沸点 Tb (K)"""
    chem = _get_chem(name)
    return _safe_float(getattr(chem, 'Tb', None)) if chem else None


def get_normal_melting_point(name: str) -> Optional[float]:
    """正常熔点 Tm (K)"""
    chem = _get_chem(name)
    return _safe_float(getattr(chem, 'Tm', None)) if chem else None


def get_dipole_moment(name: str) -> Optional[float]:
    """偶极矩 μ (Debye)"""
    chem = _get_chem(name)
    return _safe_float(getattr(chem, 'dipole', None)) if chem else None


def get_refractive_index(name: str) -> Optional[float]:
    """折射率 nD (589 nm，液相 20°C)"""
    chem = _get_chem(name)
    return _safe_float(getattr(chem, 'RI', None)) if chem else None


def get_standard_entropy(name: str, phase: str = "gas") -> Optional[float]:
    """标准摩尔熵 S° (J/mol/K)"""
    chem = _get_chem(name)
    if chem is None:
        return None
    val = _safe_float(getattr(chem, 'S0gm' if phase == "gas" else 'S0m', None))
    if val is None:
        val = _safe_float(getattr(chem, 'S0m' if phase == "gas" else 'S0gm', None))
    return val


def get_critical_compressibility(name: str) -> Optional[float]:
    """临界压缩因子 Zc"""
    chem = _get_chem(name)
    return _safe_float(getattr(chem, 'Zc', None)) if chem else None


def get_acentric_factor(name: str) -> Optional[float]:
    """Pitzer 偏心因子 ω"""
    chem = _get_chem(name)
    return _safe_float(getattr(chem, 'omega', None)) if chem else None


# ============================================================
# 安全 / 环境属性
# ============================================================
def get_safety_properties(name: str) -> Dict:
    """闪点、自燃点、爆炸极限、GWP、ODP、logP 等"""
    chem = _get_chem(name)
    if chem is None:
        raise ValueError(f"无法识别组分 '{name}'")
    return {
        "Tflash":        _safe_float(getattr(chem, 'Tflash', None)),
        "Tautoignition": _safe_float(getattr(chem, 'Tautoignition', None)),
        "LFL":           _safe_float(getattr(chem, 'LFL', None)),
        "UFL":           _safe_float(getattr(chem, 'UFL', None)),
        "GWP":           _safe_float(getattr(chem, 'GWP', None)),
        "ODP":           _safe_float(getattr(chem, 'ODP', None)),
        "logP":          _safe_float(getattr(chem, 'logP', None)),
        "TWA":           getattr(chem, 'TWA', None),
        "STEL":          getattr(chem, 'STEL', None),
        "Carcinogen":    getattr(chem, 'Carcinogen', None),
        "Skin":          getattr(chem, 'Skin', None),
    }


def get_molecular_descriptors(name: str) -> Dict:
    """偶极矩、折射率、分子直径、UNIFAC 基团等"""
    chem = _get_chem(name)
    if chem is None:
        raise ValueError(f"无法识别组分 '{name}'")
    return {
        "dipole":             _safe_float(getattr(chem, 'dipole', None)),
        "RI":                 _safe_float(getattr(chem, 'RI', None)),
        "molecular_diameter": _safe_float(getattr(chem, 'molecular_diameter', None)),
        "Stockmayer":         _safe_float(getattr(chem, 'Stockmayer', None)),
        "Zc":                 _safe_float(getattr(chem, 'Zc', None)),
        "omega":              _safe_float(getattr(chem, 'omega', None)),
        "StielPolar":         _safe_float(getattr(chem, 'StielPolar', None)),
        "UNIFAC_groups":      getattr(chem, 'UNIFAC_groups', None),
        "smiles":             getattr(chem, 'smiles', None),
        "formula":            getattr(chem, 'formula', None),
    }


# ============================================================
# K 值 & 相平衡
# ============================================================

# ============================================================
# Wilson / NRTL 二元参数数据库
# 格式: {(comp_i, comp_j): {model: (Lambda_ij [J/mol], Lambda_ji [J/mol], alpha_NRTL)}}
# Lambda_ij = a_ij + b_ij * T  (当前简化为常数)
# alpha_NRTL: 非随机性参数 (0.2~0.47)
# ============================================================
_WILSON_NRTL_DB: Dict[Tuple[str, str], Dict] = {
    ("ethanol", "water"): {
        "wilson": {"Lambda_12": 526.06, "Lambda_21": 1447.86},   # J/mol, Gmeling 1977
        "nrtl":   {"tau_12": 631.05, "tau_21": 1197.41, "alpha": 0.3},  # J/mol
    },
    ("methanol", "water"): {
        "wilson": {"Lambda_12": 1197.41, "Lambda_21": 415.27},
        "nrtl":   {"tau_12": 351.26, "tau_21": 1367.17, "alpha": 0.3},
    },
    ("acetone", "water"): {
        "wilson": {"Lambda_12": 2627.73, "Lambda_21": -225.17},
        "nrtl":   {"tau_12": 2167.35, "tau_21": 273.14, "alpha": 0.3},
    },
    ("benzene", "toluene"): {
        "wilson": {"Lambda_12": 0.0, "Lambda_21": 0.0},  # 近理想体系
        "nrtl":   {"tau_12": 0.0, "tau_21": 0.0, "alpha": 0.3},
    },
}


def _resolve_wilson_nrtl_key(components: List[str]) -> Optional[Tuple[str, str]]:
    """在数据库中查找组分对（忽略顺序）"""
    if len(components) != 2:
        return None
    c1, c2 = components[0].lower().strip(), components[1].lower().strip()
    for (a, b) in _WILSON_NRTL_DB:
        if (c1 == a and c2 == b) or (c1 == b and c2 == a):
            return (a, b)
    return None


def _calc_gamma_wilson(components: List[str], T: float, xs: List[float]) -> Optional[List[float]]:
    """Wilson 活度系数计算（仅二元体系）"""
    key = _resolve_wilson_nrtl_key(components)
    if key is None or "wilson" not in _WILSON_NRTL_DB.get(key, {}):
        return None
    params = _WILSON_NRTL_DB[key]["wilson"]
    L12, L21 = params["Lambda_12"], params["Lambda_21"]
    # 体积比近似 (用临界体积)
    V1, V2 = 1.0, 1.0
    for i, name in enumerate(components):
        crit = get_critical_constants(name)
        if crit and crit[3] > 0:
            if i == 0: V1 = crit[3]
            else: V2 = crit[3]
    Lambda_12 = (V2 / V1) * math.exp(-L12 / (R_GAS * T)) if V1 > 0 else 1.0
    Lambda_21 = (V1 / V2) * math.exp(-L21 / (R_GAS * T)) if V2 > 0 else 1.0
    x1, x2 = xs[0], xs[1]
    denom1 = x1 + x2 * Lambda_12
    denom2 = x1 * Lambda_21 + x2
    if abs(denom1) < 1e-15 or abs(denom2) < 1e-15:
        return [1.0, 1.0]
    ln_gamma1 = -math.log(denom1) + x2 * (Lambda_12 / denom1 - Lambda_21 / denom2)
    ln_gamma2 = -math.log(denom2) - x1 * (Lambda_12 / denom1 - Lambda_21 / denom2)
    return [math.exp(ln_gamma1), math.exp(ln_gamma2)]


def _calc_gamma_nrtl(components: List[str], T: float, xs: List[float]) -> Optional[List[float]]:
    """NRTL 活度系数计算（仅二元体系）"""
    key = _resolve_wilson_nrtl_key(components)
    if key is None or "nrtl" not in _WILSON_NRTL_DB.get(key, {}):
        return None
    params = _WILSON_NRTL_DB[key]["nrtl"]
    tau_12 = params["tau_12"] / (R_GAS * T)
    tau_21 = params["tau_21"] / (R_GAS * T)
    alpha = params["alpha"]
    G12 = math.exp(-alpha * tau_12)
    G21 = math.exp(-alpha * tau_21)
    x1, x2 = xs[0], xs[1]
    denom1 = x1 + x2 * G21
    denom2 = x2 + x1 * G12
    if abs(denom1) < 1e-15 or abs(denom2) < 1e-15:
        return [1.0, 1.0]
    ln_gamma1 = x2**2 * (tau_21 * (G21 / denom1)**2 + tau_12 * G12 / denom2**2)
    ln_gamma2 = x1**2 * (tau_12 * (G12 / denom2)**2 + tau_21 * G21 / denom1**2)
    return [math.exp(ln_gamma1), math.exp(ln_gamma2)]


def _can_use_wilson(components: List[str]) -> bool:
    return _resolve_wilson_nrtl_key(components) is not None and len(components) == 2


def _can_use_nrtl(components: List[str]) -> bool:
    key = _resolve_wilson_nrtl_key(components)
    return key is not None and "nrtl" in _WILSON_NRTL_DB.get(key, {}) and len(components) == 2


def _can_use_unifac(components: List[str]) -> bool:
    if UNIFAC_gammas_fn is None:
        return False
    for name in components:
        chem = _get_chem(name)
        if not chem or not getattr(chem, 'UNIFAC_groups', None):
            return False
    return True


# ============================================================
# Henry 常数数据库
# ============================================================
_T_HENRY_REF = 298.15

_HENRY_WATER_DB: Dict[str, Tuple[float, float, str]] = {
    'CO2':  (1.648e8, 2400.0, 'carbon dioxide'),
    'H2S':  (1.024e8, 2100.0, 'hydrogen sulfide'),
    'SO2':  (1.250e5, 3120.0, 'sulfur dioxide'),
    'NH3':  (5.900e0, 4110.0, 'ammonia'),
    'N2':   (1.608e9, 1300.0, 'nitrogen'),
    'O2':   (4.349e9, 1700.0, 'oxygen'),
    'H2':   (7.835e9,  500.0, 'hydrogen'),
    'CO':   (1.060e9, 1300.0, 'carbon monoxide'),
    'CH4':  (4.186e9, 1600.0, 'methane'),
    'N2O':  (2.543e8, 2300.0, 'nitrous oxide'),
    'C2H6': (2.063e10, 2000.0, 'ethane'),
    'C3H8': (1.520e10, 2500.0, 'propane'),
    'NO':   (6.580e9, 1480.0, 'nitric oxide'),
    'Ar':   (4.049e9, 1400.0, 'argon'),
    'He':   (2.705e10, 230.0, 'helium'),
    'Kr':   (2.398e9, 1700.0, 'krypton'),
    'Xe':   (4.510e8, 2100.0, 'xenon'),
    'Cl2':  (9.670e6, 2500.0, 'chlorine'),
}

_HENRY_ALIASES: Dict[str, str] = {
    'carbon dioxide': 'CO2', 'co2': 'CO2',
    'hydrogen sulfide': 'H2S', 'h2s': 'H2S',
    'sulfur dioxide': 'SO2', 'so2': 'SO2',
    'ammonia': 'NH3', 'nh3': 'NH3',
    'nitrogen': 'N2', 'n2': 'N2',
    'oxygen': 'O2', 'o2': 'O2',
    'hydrogen': 'H2', 'h2': 'H2',
    'carbon monoxide': 'CO', 'co': 'CO',
    'methane': 'CH4', 'ch4': 'CH4',
    'nitrous oxide': 'N2O', 'n2o': 'N2O',
    'ethane': 'C2H6', 'c2h6': 'C2H6',
    'propane': 'C3H8', 'c3h8': 'C3H8',
    'nitric oxide': 'NO', 'no': 'NO',
    'argon': 'Ar', 'ar': 'Ar',
    'helium': 'He', 'he': 'He',
    'krypton': 'Kr', 'kr': 'Kr',
    'xenon': 'Xe', 'xe': 'Xe',
    'chlorine': 'Cl2', 'cl2': 'Cl2',
}

_LIGHT_GAS_NAMES: frozenset = frozenset(_HENRY_WATER_DB.keys())
_WATER_NAMES: frozenset = frozenset({'water', 'Water', 'H2O', 'h2o'})


def _resolve_henry_name(name: str) -> Optional[str]:
    if name in _HENRY_WATER_DB:
        return name
    lower = name.lower().strip()
    mapped = _HENRY_ALIASES.get(lower)
    if mapped:
        return mapped
    chem = _get_chem(name)
    if chem:
        chem_name = getattr(chem, 'name', name)
        if chem_name in _HENRY_WATER_DB:
            return chem_name
        mapped2 = _HENRY_ALIASES.get(chem_name.lower().strip())
        if mapped2:
            return mapped2
    return None


def _is_gas_water_system(components: List[str]) -> bool:
    has_gas = False
    has_water = False
    for c in components:
        if _resolve_henry_name(c) is not None:
            has_gas = True
        if c in _WATER_NAMES or c.lower().strip() in ('water', 'h2o'):
            has_water = True
    return has_gas and has_water


def calc_henry_constant(name: str, T: float = 298.15,
                        solvent: str = 'water') -> Optional[float]:
    """
    气体在水中的亨利常数 H_x [Pa]（摩尔分数基准）。
    Sander (2015), Atmos. Chem. Phys., 15, 4399-4981.
    有效范围 273–373 K。
    """
    if solvent.lower().strip() not in ('water', 'h2o'):
        return None
    canon = _resolve_henry_name(name)
    if canon is None:
        return None
    H_ref, B, _ = _HENRY_WATER_DB[canon]
    lnH = math.log(H_ref) - B * (1.0 / T - 1.0 / _T_HENRY_REF)
    return math.exp(lnH)


def _has_supercritical_component(components: List[str], T: float) -> bool:
    for name in components:
        crit = get_critical_constants(name)
        if crit and T >= crit[0]:
            return True
    return False


def _auto_k_method(components: List[str], T: float, P: float) -> str:
    """自动选择 K 值计算方法。

    方法优先级（数据驱动分层策略）:
      ① henry — 气体-水体系（低压），Henry 定律处理溶解气体
      ② eos_pr — 高压 (>15 bar) 或超临界组分
      ③ unifac — 所有组分有 UNIFAC 基团数据（低压非理想体系默认方法）
      ④ wilson — 二元体系有实验回归参数（仅当显式指定时启用）
      ⑤ nrtl — 二元体系有实验回归参数（仅当显式指定时启用）
      ⑥ raoult — 兜底 Raoult 定律

    设计原则: UNIFAC 基于基团贡献法，无需实验拟合，通用精度最高；
    Wilson/NRTL 参数为拟合常数，仅在有实验数据回归时才优于 UNIFAC。
    """
    # ① 气体-水体系 + 低压 → Henry 定律
    if _is_gas_water_system(components) and P < 1.5e6:
        return "henry"
    # ② 超临界组分 → EOS
    if _has_supercritical_component(components, T):
        return "eos_pr"
    # ③ 高压 → EOS
    Pr_max = max((get_reduced_properties(n, T, P)[1] for n in components), default=0.0)
    if P > 1.5e6 or Pr_max > 0.15:
        return "eos_pr"
    # ④ UNIFAC 可用 → 基团贡献法（默认首选）
    if _can_use_unifac(components):
        return "unifac"
    # ⑤ Wilson 活度系数（二元体系，需实验回归参数）
    if _can_use_wilson(components):
        return "wilson"
    # ⑥ NRTL 活度系数（二元体系，需实验回归参数）
    if _can_use_nrtl(components):
        return "nrtl"
    # ⑦ 兜底 Raoult
    return "raoult"


def _calc_Ks(components: List[str], T: float, P: float,
             zs: List[float], method: str = "auto") -> List[float]:
    if method == "auto":
        method = _auto_k_method(components, T, P)

    if method == "henry":
        Ks = []
        for c in components:
            canon = _resolve_henry_name(c)
            if canon is not None:
                H = calc_henry_constant(c, T)
                Ks.append(H / P if H else 1.0)
            else:
                Ks.append((_safe_float(calc_vapor_pressure(c, T)) or P) / P)
        return Ks

    if method in ("eos_pr", "eos_srk"):
        EOS_cls = PRMIX if method == "eos_pr" else SRKMIX
        Tcs, Pcs, omegas, CASs, ok = [], [], [], [], True
        for c in components:
            crit = get_critical_constants(c)
            if not crit:
                ok = False; break
            Tcs.append(crit[0]); Pcs.append(crit[1]); omegas.append(crit[4])
            ch = _get_chem(c)
            CASs.append(ch.CAS if ch else None)
        if ok:
            try:
                kijs = (IPDB.get_ip_symmetric_matrix('ChemSep PR', CASs, 'kij')
                        if IPDB and all(CASs) else None)
                eos = EOS_cls(T=T, P=P, Tcs=Tcs, Pcs=Pcs, omegas=omegas,
                              zs=zs, kijs=kijs)
                # [FIX-8] 检查 lnphis_l / lnphis_g 是否有效，防止单相根崩溃
                lnphis_l = getattr(eos, 'lnphis_l', None)
                lnphis_g = getattr(eos, 'lnphis_g', None)
                if (lnphis_l is None or lnphis_g is None
                        or len(lnphis_l) != len(components)
                        or len(lnphis_g) != len(components)):
                    raise ValueError(
                        "EOS 只有单相根（lnphis_l 或 lnphis_g 为 None），"
                        "可能处于单相区，回退到 UNIFAC/Raoult。"
                    )
                return [math.exp(float(lnphis_l[i]) - float(lnphis_g[i]))
                        for i in range(len(components))]
            except Exception:
                pass  # [FIX-8] 回退到下方 UNIFAC/Raoult

    # ---- Wilson 活度系数法 ----
    if method == "wilson":
        gammas = _calc_gamma_wilson(components, T, zs)
        if gammas is not None:
            return [gammas[i] * (calc_vapor_pressure(components[i], T) or 0.0) / P
                    for i in range(len(components))]

    # ---- NRTL 活度系数法 ----
    if method == "nrtl":
        gammas = _calc_gamma_nrtl(components, T, zs)
        if gammas is not None:
            return [gammas[i] * (calc_vapor_pressure(components[i], T) or 0.0) / P
                    for i in range(len(components))]

    if method == "unifac" and UNIFAC_gammas_fn:
        try:
            cgs = [_get_chem(c).UNIFAC_groups for c in components]
            gammas = UNIFAC_gammas_fn(T=T, xs=zs, chemgroups=cgs)
            return [gammas[i] * (calc_vapor_pressure(components[i], T) or 0.0) / P
                    for i in range(len(components))]
        except Exception:
            pass

    return [(_safe_float(calc_vapor_pressure(c, T)) or P) / P for c in components]


def calc_activity_coefficients_unifac(components: List[str], T: float,
                                      xs: List[float]) -> Optional[List[float]]:
    """UNIFAC 活度系数 γi"""
    if not _can_use_unifac(components):
        return None
    try:
        cgs = [_get_chem(c).UNIFAC_groups for c in components]
        return list(UNIFAC_gammas_fn(T=T, xs=xs, chemgroups=cgs))
    except Exception:
        return None


def _wilson_k_initial(components: List[str], T: float, P: float) -> List[float]:
    Ks = []
    for name in components:
        crit = get_critical_constants(name)
        if crit:
            Tc, Pc, _, _, omega, _ = crit
            Ks.append(max((Pc / P) * math.exp(5.37 * (1 + omega) * (1 - Tc / T)), 1e-6))
        else:
            Ks.append(1.0)
    return Ks


def calc_bubble_point_T(components: List[str], x: List[float], P: float,
                        method: str = "auto", tol: float = 1e-4,
                        max_iter: int = 200,
                        precision: str = "normal") -> Dict:
    """泡点温度，返回 {T, T_C, y, K, method, converged, iterations, azeotrope_near}

    升级特性:
      - 初值猜测: 摩尔分率加权正常沸点
      - 收敛加速: Wegstein 割线法（超线性收敛）
      - 共沸点检测: 所有 K_i ≈ 1 时切换 scipy.optimize.brentq
      - K 值方法: 迭代中动态更新（温度变化导致对比压力跨阈值时自动切换）
      - 自适应精度: precision="normal"(1e-4) / "high"(1e-6) / "ultra"(1e-8)
    """
    n = len(components)
    # 自适应精度
    _tol_map = {"normal": 1e-4, "high": 1e-6, "ultra": 1e-8}
    tol = _tol_map.get(precision, tol)

    # ── 第一层: thermo Mixture.Tbubble 直接求解 ──
    try:
        m = Mixture(components, zs=x, P=P)
        Tb = float(m.Tbubble)
        return {"T": round(Tb, 4), "T_C": round(Tb-273.15, 4), "y": list(m.ys or []),
                "K": [m.ys[i]/x[i] if x[i] > 0 else 0 for i in range(n)],
                "method": "Mixture", "converged": True, "iterations": 0, "azeotrope_near": False}
    except Exception:
        pass

    # ── 第二层: 改进的 K 值迭代法 ──
    # ① 初值猜测: 摩尔分率加权正常沸点
    T_guess = 350.0
    Tb_sum = 0.0
    w_sum = 0.0
    for i in range(n):
        if x[i] > 1e-10:
            Tb_i = get_normal_boiling_point(components[i])
            if Tb_i and Tb_i > 0:
                Tb_sum += x[i] * Tb_i
                w_sum += x[i]
    if w_sum > 0:
        T_guess = Tb_sum / w_sum
    T = max(150.0, min(1100.0, T_guess))

    # ② 迭代求解 (Wegstein 割线法加速)
    method_used = _auto_k_method(components, T, P) if method == "auto" else method
    T_prev = None
    f_prev = None
    azeotrope_near = False
    iterations = 0

    for k in range(max_iter):
        iterations = k + 1
        # 动态更新 K 值方法
        if method == "auto":
            method_used = _auto_k_method(components, T, P)

        Ks = _calc_Ks(components, T, P, x, method_used)
        sum_Kx = sum(Ks[i]*x[i] for i in range(n))
        f_curr = sum_Kx - 1.0

        # 收敛判据
        if abs(f_curr) < tol:
            break

        # 共沸点检测: 所有 K_i ≈ 1 (y ≈ x)
        if all(abs(Ks[i] - 1.0) < 0.05 for i in range(n)):
            azeotrope_near = True
            # 切换 scipy.optimize.brentq 高精度求解
            try:
                from scipy.optimize import brentq
                def _obj(T_try):
                    Ks_try = _calc_Ks(components, T_try, P, x, method_used)
                    return sum(Ks_try[i]*x[i] for i in range(n)) - 1.0
                T_lo = max(150.0, T - 20.0)
                T_hi = min(1100.0, T + 20.0)
                if _obj(T_lo) * _obj(T_hi) < 0:
                    T = brentq(_obj, T_lo, T_hi, xtol=tol * 0.1, maxiter=100)
                    break
            except Exception:
                pass  # 回退到普通迭代

        # Wegstein 割线法加速
        if T_prev is not None and f_prev is not None and abs(f_curr - f_prev) > 1e-14:
            # 割线外推
            dT = -f_curr * (T - T_prev) / (f_curr - f_prev)
            # 限幅: 单步最大 ±50 K
            dT = max(-50.0, min(50.0, dT))
            T_new = T + dT
        else:
            # 首步: 阻尼牛顿法
            dT = (1.0 - sum_Kx) * T * 0.03
            dT = max(-30.0, min(30.0, dT))
            T_new = T + dT

        T_prev = T
        f_prev = f_curr
        T = max(100.0, min(1200.0, T_new))

    # 最终 K 值
    Ks = _calc_Ks(components, T, P, x, method_used)
    y = [Ks[i]*x[i] for i in range(n)]
    converged = abs(sum(Ks[i]*x[i] for i in range(n)) - 1.0) < tol
    return {"T": round(T, 4), "T_C": round(T-273.15, 4), "y": y, "K": Ks,
            "method": method_used, "converged": converged,
            "iterations": iterations, "azeotrope_near": azeotrope_near}


def calc_dew_point_T(components: List[str], y: List[float], P: float,
                     method: str = "auto", tol: float = 1e-4,
                     max_iter: int = 200,
                     precision: str = "normal") -> Dict:
    """露点温度，返回 {T, T_C, x, K, method, converged, iterations}

    升级特性: 同 calc_bubble_point_T（加权初值、割线法加速、共沸点检测、动态 K 值、自适应精度）
    """
    n = len(components)
    _tol_map = {"normal": 1e-4, "high": 1e-6, "ultra": 1e-8}
    tol = _tol_map.get(precision, tol)

    # ── 第一层: thermo Mixture.Tdew ──
    try:
        m = Mixture(components, zs=y, P=P)
        Td = float(m.Tdew)
        return {"T": round(Td,4), "T_C": round(Td-273.15,4), "x": list(m.xs or []),
                "K": [y[i]/m.xs[i] if m.xs[i] > 0 else 0 for i in range(n)],
                "method": "Mixture", "converged": True, "iterations": 0}
    except Exception:
        pass

    # ── 第二层: 改进的 K 值迭代法 ──
    # 初值: 摩尔分率加权正常沸点
    T_guess = 360.0
    Tb_sum, w_sum = 0.0, 0.0
    for i in range(n):
        if y[i] > 1e-10:
            Tb_i = get_normal_boiling_point(components[i])
            if Tb_i and Tb_i > 0:
                Tb_sum += y[i] * Tb_i
                w_sum += y[i]
    if w_sum > 0:
        T_guess = Tb_sum / w_sum
    T = max(150.0, min(1100.0, T_guess))

    method_used = _auto_k_method(components, T, P) if method == "auto" else method
    T_prev, f_prev = None, None
    iterations = 0

    for k in range(max_iter):
        iterations = k + 1
        if method == "auto":
            method_used = _auto_k_method(components, T, P)
        Ks = _calc_Ks(components, T, P, y, method_used)
        # 露点方程: Σ(y_i/K_i) = 1
        sum_yK = sum(y[i]/Ks[i] if Ks[i] > 1e-15 else 0 for i in range(n))
        f_curr = sum_yK - 1.0

        if abs(f_curr) < tol:
            break

        # Wegstein 割线法
        if T_prev is not None and f_prev is not None and abs(f_curr - f_prev) > 1e-14:
            dT = -f_curr * (T - T_prev) / (f_curr - f_prev)
            dT = max(-50.0, min(50.0, dT))
            T_new = T + dT
        else:
            dT = (sum_yK - 1.0) * T * 0.03
            dT = max(-30.0, min(30.0, dT))
            T_new = T + dT

        T_prev, f_prev = T, f_curr
        T = max(100.0, min(1200.0, T_new))

    Ks = _calc_Ks(components, T, P, y, method_used)
    x = [y[i]/Ks[i] if Ks[i] > 1e-15 else 0 for i in range(n)]
    converged = abs(sum(y[i]/Ks[i] if Ks[i] > 1e-15 else 0 for i in range(n)) - 1.0) < tol
    return {"T": round(T,4), "T_C": round(T-273.15,4), "x": x, "K": Ks,
            "method": method_used, "converged": converged, "iterations": iterations}


def calc_VLE_envelope(components: List[str], P: float, n_points: int = 20,
                      method: str = "auto", precision: str = "normal") -> Dict:
    """生成完整的等压 T-x-y 包络线（泡点线 + 露点线），含共沸点检测。

    Parameters
    ----------
    components : List[str]
        二元组分列表（仅支持二元体系）
    P : float
        系统压力 (Pa)
    n_points : int
        包络线采样点数（默认 20）
    method : str
        K 值计算方法（"auto" 自动选择）
    precision : str
        收敛精度 ("normal"/"high"/"ultra")

    Returns
    -------
    Dict with keys:
        - bubble_line: List[Dict] — 泡点线数据 [{x, T, y, method, iterations}, ...]
        - dew_line: List[Dict] — 露点线数据 [{y, T, x, method, iterations}, ...]
        - azeotrope: Dict or None — 共沸点信息 {T, T_C, z, method}
        - n_points: int
    """
    if len(components) != 2:
        return {"error": "calc_VLE_envelope 仅支持二元体系", "bubble_line": [], "dew_line": [], "azeotrope": None}

    # 生成 x 采样点（避开 0 和 1 的奇点）
    xs = [i / (n_points - 1) for i in range(n_points)]
    xs[0] = 1e-6   # 避免 x=0
    xs[-1] = 1 - 1e-6  # 避免 x=1

    bubble_line = []
    dew_line = []
    azeotrope = None

    # 泡点线
    for x1 in xs:
        x = [x1, 1.0 - x1]
        res = calc_bubble_point_T(components, x, P, method=method, precision=precision)
        bubble_line.append({
            "x": x, "T": res["T"], "T_C": res["T_C"],
            "y": res.get("y", []), "method": res.get("method", ""),
            "iterations": res.get("iterations", 0),
        })
        # 共沸点检测
        if res.get("azeotrope_near") and azeotrope is None:
            azeotrope = {"T": res["T"], "T_C": res["T_C"], "z": x, "method": res.get("method", "")}

    # 露点线
    for y1 in xs:
        y = [y1, 1.0 - y1]
        res = calc_dew_point_T(components, y, P, method=method, precision=precision)
        dew_line.append({
            "y": y, "T": res["T"], "T_C": res["T_C"],
            "x": res.get("x", []), "method": res.get("method", ""),
            "iterations": res.get("iterations", 0),
        })

    # 精细共沸点搜索: 在泡点线中找最低温度点
    if azeotrope is None and bubble_line:
        min_bubble = min(bubble_line, key=lambda b: b["T"])
        # 检查该点是否接近共沸 (y ≈ x)
        y_bubble = min_bubble["y"]
        x_bubble = min_bubble["x"]
        if len(y_bubble) == 2 and len(x_bubble) == 2:
            if abs(y_bubble[0] - x_bubble[0]) < 0.02:
                azeotrope = {
                    "T": min_bubble["T"], "T_C": min_bubble["T_C"],
                    "z": x_bubble, "method": min_bubble["method"],
                    "note": "detected by minimum boiling point"
                }

    return {
        "bubble_line": bubble_line,
        "dew_line": dew_line,
        "azeotrope": azeotrope,
        "n_points": n_points,
    }


def calc_bubble_point_P(components: List[str], x: List[float], T: float,
                        method: str = "auto") -> Dict:
    """泡点压力"""
    try:
        m = Mixture(components, zs=x, T=T)
        Pb = float(m.Pbubble)
        return {"P": Pb, "P_MPa": Pb/1e6,
                "K": [m.ys[i]/x[i] if x[i]>0 else 0 for i in range(len(x))],
                "method": "Mixture", "converged": True}
    except Exception:
        pass
    Pb = sum(x[i] * (calc_vapor_pressure(components[i], T) or 0) for i in range(len(components)))
    return {"P": Pb, "P_MPa": Pb/1e6, "K": [], "method": "Raoult", "converged": True}


def calc_dew_point_P(components: List[str], y: List[float], T: float,
                     method: str = "auto") -> Dict:
    """露点压力"""
    try:
        m = Mixture(components, zs=y, T=T)
        Pd = float(m.Pdew)
        return {"P": Pd, "P_MPa": Pd/1e6,
                "K": [y[i]/m.xs[i] if m.xs[i]>0 else 0 for i in range(len(y))],
                "method": "Mixture", "converged": True}
    except Exception:
        pass
    inv = sum(y[i]/(calc_vapor_pressure(components[i], T) or 1.0) for i in range(len(components)))
    Pd = 1.0/inv if inv > 0 else 0.0
    return {"P": Pd, "P_MPa": Pd/1e6, "K": [], "method": "Raoult", "converged": True}


def calc_flash_rachford_rice(components: List[str], z: List[float],
                             T: float, P: float,
                             method: str = "auto") -> Dict:
    """等温闪蒸 (Rachford-Rice), 返回 {beta, x, y, K, phase, method}"""
    n = len(components)
    method_used = _auto_k_method(components, T, P) if method == "auto" else method
    Ks = _calc_Ks(components, T, P, z, method_used)

    def rr_eq(b):
        return sum(z[i]*(Ks[i]-1)/(1+b*(Ks[i]-1)) for i in range(n))

    beta = 0.5
    K_min, K_max = min(Ks), max(Ks)
    if K_max <= 1.0:
        return {"beta": 0.0, "x": list(z), "y": list(z), "K": Ks,
                "phase": "liquid", "method": method_used}
    if K_min >= 1.0:
        return {"beta": 1.0, "x": list(z), "y": list(z), "K": Ks,
                "phase": "vapor", "method": method_used}
    try:
        from scipy.optimize import brentq
        lo = max(0.0, 1.0/(1.0-K_max)+1e-10)
        hi = min(1.0, 1.0/(1.0-K_min)-1e-10)
        beta = brentq(rr_eq, lo, hi, xtol=1e-10, maxiter=100)
    except Exception:
        for _ in range(200):
            f = rr_eq(beta)
            df = -sum(z[i]*(Ks[i]-1)**2/(1+beta*(Ks[i]-1))**2 for i in range(n))
            if abs(df) < 1e-30: break
            beta = max(0.0, min(1.0, beta - f/df))
            if abs(f) < 1e-10: break
    x = [z[i]/(1+beta*(Ks[i]-1)) for i in range(n)]
    y = [Ks[i]*x[i] for i in range(n)]
    sx, sy = sum(x), sum(y)
    if sx > 0: x = [v/sx for v in x]
    if sy > 0: y = [v/sy for v in y]
    phase = "liquid" if beta <= 0 else ("vapor" if beta >= 1 else "two_phase")
    return {"beta": round(beta, 6), "x": x, "y": y, "K": Ks,
            "phase": phase, "method": method_used}


# ============================================================
# 纯物质综合属性
# ============================================================
def get_chemical_properties(component_name: str, T: float = 298.15,
                            P: float = 101325.0) -> Dict:
    """纯物质全面物性字典 (thermo + CoolProp 双引擎)"""
    regime  = classify_regime(component_name, T, P)
    Tr, Pr  = get_reduced_properties(component_name, T, P)
    chem    = _get_chem(component_name, T=T, P=P)
    cp_name = _get_coolprop_name(component_name)

    MW   = get_molecular_weight(component_name)
    crit = get_critical_constants(component_name)
    Tc   = crit[0] if crit else None
    Pc   = crit[1] if crit else None
    Vc   = crit[2] if crit else None
    Zc   = crit[3] if crit else None
    omega= crit[4] if crit else None

    cp_props: Dict = {}
    if cp_name:
        keys_cp = ['D','Cpmolar','Cvmolar','speed_of_sound','Z','viscosity',
                   'conductivity','Prandtl','isothermal_compressibility',
                   'isobaric_expansion_coefficient','isentropic_expansion_coefficient',
                   'Hmolar','Smolar','Umolar','Gmolar','A','Dmolar']
        for k in keys_cp:
            cp_props[k] = _cp_prop(k, T, P, cp_name)
        f_cp, phi_cp = _cp_fugacity(T, P, cp_name)
        cp_props['fugacity'] = f_cp
        cp_props['phi']      = phi_cp

    if chem is None and not cp_name:
        blank = {k: None for k in [
            "MW","Tc","Pc","Vc","Zc","omega","Tb","Tm","Tt","Pt","Tr","Pr","regime",
            "rho_L","rho_G","rho_S","Vm_L","Vm_G","Z_L","Z_G",
            "mu_L","mu_G","nu_L","nu_G","k_L","k_G","alpha_L","alpha_G",
            "Pr_L","Pr_G","sigma","permittivity","Parachor",
            "Psat","Psub","Hvap","Hfus","Hsub",
            "Cp_L","Cp_G","Cp_S","Cv_G","Cv_L","gamma_G","gamma_L","c_G","c_L","kappa_s_L",
            "H_molar","S_molar","U_molar","G_molar","A_molar",
            "Hf_gas","Hf_liq","Gf_gas","Gf_liq","Hc_HHV","Hc_LHV","S0_gas","S0_liq",
            "Bvirial","JT_G","JT_L","beta_G","beta_L","kappa_G","kappa_L",
            "f_G","f_L","phi_G","phi_L","source",
        ]}
        blank.update({"MW": MW, "Tc": Tc, "Pc": Pc, "Vc": Vc, "Zc": Zc, "omega": omega,
                      "Tr": Tr, "Pr": Pr, "regime": regime, "source": "none"})
        return blank

    rho_L = calc_liquid_density(component_name, T, P)
    rho_G = calc_gas_density(component_name, T, P)
    rho_S = calc_solid_density(component_name, T)
    Vm_L  = calc_molar_volume(component_name, T, P, "liquid")
    Vm_G  = calc_molar_volume(component_name, T, P, "gas")

    Z_val = cp_props.get('Z')
    Z_L = Z_val if (Z_val and classify_regime(component_name, T, P) == PropertyRegime.LIQUID) else None
    Z_G = Z_val if (Z_val and classify_regime(component_name, T, P) != PropertyRegime.LIQUID) else None
    if Z_L is None: Z_L = calc_compressibility_factor(component_name, T, P, "liquid")
    if Z_G is None: Z_G = calc_compressibility_factor(component_name, T, P, "gas")

    mu_L = calc_viscosity_liquid(component_name, T, P)
    mu_G = calc_viscosity_gas(component_name, T, P)
    nu_L = calc_kinematic_viscosity(component_name, T, P, "liquid")
    nu_G = calc_kinematic_viscosity(component_name, T, P, "gas")
    k_L  = calc_thermal_conductivity_liquid(component_name, T, P)
    k_G  = calc_thermal_conductivity_gas(component_name, T, P)
    alpha_L = calc_thermal_diffusivity(component_name, T, P, "liquid")
    alpha_G = calc_thermal_diffusivity(component_name, T, P, "gas")
    Pr_L = calc_prandtl_number(component_name, T, P, "liquid")
    Pr_G = calc_prandtl_number(component_name, T, P, "gas")
    sigma = calc_surface_tension(component_name, T)

    Cp_L = calc_heat_capacity(component_name, T, P, "liquid")
    Cp_G = calc_heat_capacity(component_name, T, P, "gas")
    Cp_S = calc_heat_capacity(component_name, T, P, "solid")
    Cv_G = calc_Cv(component_name, T, P, "gas")
    Cv_L = calc_Cv(component_name, T, P, "liquid")
    gamma_G = calc_isentropic_exponent(component_name, T, P, "gas")
    gamma_L = calc_isentropic_exponent(component_name, T, P, "liquid")
    c_G = calc_speed_of_sound(component_name, T, P, "gas")
    c_L = calc_speed_of_sound(component_name, T, P, "liquid")
    kappa_s_L = calc_isentropic_compressibility(component_name, T, P, "liquid")

    H_mol = cp_props.get('Hmolar') or calc_enthalpy(component_name, T, P)
    S_mol = cp_props.get('Smolar') or calc_entropy(component_name, T, P)
    U_mol = cp_props.get('Umolar') or calc_internal_energy(component_name, T, P)
    G_mol = cp_props.get('Gmolar') or calc_gibbs_energy(component_name, T, P)
    A_mol = calc_helmholtz_energy(component_name, T, P)

    Hf_gas = get_formation_enthalpy(component_name, T, P, "gas")
    Hf_liq = get_formation_enthalpy(component_name, T, P, "liquid")
    Gf_gas = get_formation_gibbs(component_name, "gas")
    Gf_liq = get_formation_gibbs(component_name, "liquid")
    Hc_HHV = get_combustion_enthalpy(component_name, "gas", True)
    Hc_LHV = get_combustion_enthalpy(component_name, "gas", False)
    S0_gas = _safe_float(getattr(chem, 'S0gm', None)) if chem else None
    S0_liq = _safe_float(getattr(chem, 'S0m',  None)) if chem else None

    Bvirial  = calc_second_virial_coefficient(component_name, T)
    JT_G     = calc_joule_thomson_coefficient(component_name, T, P, "gas")
    JT_L     = calc_joule_thomson_coefficient(component_name, T, P, "liquid")
    beta_G   = calc_isobaric_expansion(component_name, T, P, "gas")
    beta_L   = calc_isobaric_expansion(component_name, T, P, "liquid")
    kappa_G  = calc_isothermal_compressibility(component_name, T, P, "gas")
    kappa_L  = calc_isothermal_compressibility(component_name, T, P, "liquid")

    f_cp   = cp_props.get('fugacity')
    phi_cp = cp_props.get('phi')

    Psat = calc_vapor_pressure(component_name, T)
    Psub = calc_sublimation_pressure(component_name, T)
    Hvap = calc_enthalpy_vaporization(component_name, T)
    Hfus = calc_enthalpy_fusion(component_name)
    Hsub = calc_enthalpy_sublimation(component_name)
    Tb   = get_normal_boiling_point(component_name)
    Tm   = get_normal_melting_point(component_name)
    tp   = get_triple_point(component_name)
    perm = calc_permittivity(component_name, T)
    parachor = calc_parachor(component_name, T)

    def _r(v, d=4): return round(v, d) if v is not None else None

    return {
        "MW": MW, "Tc": Tc, "Pc": Pc, "Vc": Vc, "Zc": Zc, "omega": omega,
        "Tb": Tb, "Tm": Tm,
        "Tt": tp[0] if tp else None, "Pt": tp[1] if tp else None,
        "Tr": _r(Tr, 5), "Pr": _r(Pr, 5), "regime": regime,
        "rho_L": _r(rho_L, 3), "rho_G": _r(rho_G, 5), "rho_S": _r(rho_S, 3),
        "Vm_L": Vm_L, "Vm_G": Vm_G, "Z_L": _r(Z_L, 6), "Z_G": _r(Z_G, 6),
        "mu_L": mu_L, "mu_G": mu_G, "nu_L": nu_L, "nu_G": nu_G,
        "k_L": k_L, "k_G": k_G, "alpha_L": alpha_L, "alpha_G": alpha_G,
        "Pr_L": _r(Pr_L, 4), "Pr_G": _r(Pr_G, 4),
        "sigma": sigma, "permittivity": perm, "Parachor": parachor,
        "Psat": Psat, "Psub": Psub, "Hvap": Hvap, "Hfus": Hfus, "Hsub": Hsub,
        "Cp_L": Cp_L, "Cp_G": Cp_G, "Cp_S": Cp_S,
        "Cv_G": Cv_G, "Cv_L": Cv_L,
        "gamma_G": _r(gamma_G, 4), "gamma_L": _r(gamma_L, 4),
        "c_G": _r(c_G, 2), "c_L": _r(c_L, 2), "kappa_s_L": kappa_s_L,
        "H_molar": H_mol, "S_molar": S_mol, "U_molar": U_mol,
        "G_molar": G_mol, "A_molar": A_mol,
        "Hf_gas": Hf_gas, "Hf_liq": Hf_liq,
        "Gf_gas": Gf_gas, "Gf_liq": Gf_liq,
        "Hc_HHV": Hc_HHV, "Hc_LHV": Hc_LHV,
        "S0_gas": S0_gas, "S0_liq": S0_liq,
        "Bvirial": Bvirial,
        "JT_G": JT_G, "JT_L": JT_L,
        "beta_G": beta_G, "beta_L": beta_L,
        "kappa_G": kappa_G, "kappa_L": kappa_L,
        "f_G": f_cp, "f_L": f_cp,
        "phi_G": phi_cp, "phi_L": phi_cp,
        "source": f"thermo+CoolProp[{cp_name}]" if cp_name else "thermo",
    }


# ============================================================
# 混合物综合属性
# ============================================================
def get_mixture_properties(components: List[str], z: List[float],
                           T: float = 298.15, P: float = 101325.0) -> Dict:
    """
    混合物物性字典

    混合规则:
      密度  : PRMIX EOS → 体积加权回退
      粘度  : 液相对数加权（质量分率）；气相 Graham 简化规则（摩尔分率×√MW 加权）
              [FIX-11] 气相粘度混合规则为 Wilke 简化版（Graham 规则），
              与完整 Wilke 公式相比误差约 ±5-15%；完整 Wilke 公式需要
              二元交互参数 Φ_ij，此处用 √(MW_i/MW_j) 简化，适合工程估算。
      导热  : 质量分率加权
      Cp    : 摩尔分率加权
      表面张力: Macleod-Sugden (Parachor) → 线性回退
               [FIX-9] 过滤 sigma=None 的组分（超临界组分无表面张力）
    """
    n = len(components)
    pure_props = [get_chemical_properties(c, T, P) for c in components]
    MWs = [p["MW"] for p in pure_props]
    if any(mw is None for mw in MWs):
        raise ValueError(f"无法获取完整分子量: {components}")

    MW_mix = sum(z[i]*MWs[i] for i in range(n))
    w = [z[i]*MWs[i]/MW_mix for i in range(n)] if MW_mix > 0 else list(z)

    k_method = _auto_k_method(components, T, P)
    Ks = _calc_Ks(components, T, P, z, k_method)

    # EOS 混合密度
    rho_L_mix = rho_V_mix = Z_L_mix = Z_V_mix = None
    try:
        Tcs = [p["Tc"] for p in pure_props]; Pcs = [p["Pc"] for p in pure_props]
        oms = [p["omega"] for p in pure_props]
        CASs = [_get_chem(c).CAS for c in components]
        if all(Tcs) and all(Pcs) and all(oms) and all(CASs):
            kijs = (IPDB.get_ip_symmetric_matrix('ChemSep PR', CASs, 'kij') if IPDB else None)
            em = PRMIX(T=T, P=P, Tcs=Tcs, Pcs=Pcs, omegas=oms, zs=z, kijs=kijs)
            if em.V_l:
                rho_L_mix = (MW_mix/1000.0)/float(em.V_l); Z_L_mix = float(em.Z_l)
            if em.V_g:
                rho_V_mix = (MW_mix/1000.0)/float(em.V_g); Z_V_mix = float(em.Z_g)
    except Exception:
        pass

    if rho_L_mix is None:
        try:
            inv = sum(w[i]/p["rho_L"] for i,p in enumerate(pure_props) if p["rho_L"])
            if inv > 0:
                rho_L_mix = 1.0/inv
                Z_L_mix = P*(MW_mix/1000.0)/(rho_L_mix*R_GAS*T)
        except Exception:
            pass

    if rho_V_mix is None:
        try:
            # [FIX-5] 补充 inv > 0 保护，防止全部 rho_G 为 None 时 ZeroDivisionError
            valid_items = [(w[i], p["rho_G"]) for i, p in enumerate(pure_props) if p["rho_G"]]
            if valid_items:
                inv = sum(wi / rho for wi, rho in valid_items)
                if inv > 0:
                    rho_V_mix = 1.0 / inv
                else:
                    rho_V_mix = P * MW_mix / 1000.0 / (R_GAS * T)
            else:
                # 全部 rho_G 为 None，回退理想气体
                rho_V_mix = P * MW_mix / 1000.0 / (R_GAS * T)
            Z_V_mix = P * (MW_mix / 1000.0) / (rho_V_mix * R_GAS * T)
        except Exception:
            pass

    # 粘度
    mu_L_mix = mu_V_mix = None
    try:
        if all(p["mu_L"] for p in pure_props):
            mu_L_mix = math.exp(sum(w[i]*math.log(max(p["mu_L"],1e-15))
                                     for i,p in enumerate(pure_props)))
    except Exception: pass
    try:
        # [FIX-11] Graham 简化规则（Wilke 简化版），误差 ±5-15%
        if all(p["mu_G"] for p in pure_props):
            denom = sum(z[i]*math.sqrt(max(MWs[i],1.0)) for i in range(n))
            if denom > 0:
                mu_V_mix = (sum(z[i]*pure_props[i]["mu_G"]*math.sqrt(max(MWs[i],1.0))
                               for i in range(n)) / denom)
    except Exception: pass

    k_L_mix = (sum(w[i]*p["k_L"] for i,p in enumerate(pure_props))
               if all(p["k_L"] for p in pure_props) else None)
    k_V_mix = (sum(w[i]*p["k_G"] for i,p in enumerate(pure_props))
               if all(p["k_G"] for p in pure_props) else None)
    Cp_L_mix = (sum(z[i]*p["Cp_L"] for i,p in enumerate(pure_props))
                if all(p["Cp_L"] for p in pure_props) else None)
    Cp_V_mix = (sum(z[i]*p["Cp_G"] for i,p in enumerate(pure_props))
                if all(p["Cp_G"] for p in pure_props) else None)

    # 表面张力 (Macleod-Sugden)
    # [FIX-9] 过滤 sigma=None（超临界组分），避免 None 参与计算
    sigma_mix = None
    try:
        parachors = [p["Parachor"] for p in pure_props]
        rhos_L    = [p["rho_L"]    for p in pure_props]
        rhos_G    = [p["rho_G"]    for p in pure_props]
        # 所有组分均需有效的 Parachor、rho_L、rho_G（超临界组分 rho_L 可能 None）
        if all(parachors) and all(rhos_L) and all(rhos_G):
            ys = [Ks[i]*z[i] for i in range(n)]
            sy = sum(ys); ys = [v/sy for v in ys] if sy > 0 else ys
            sq = sum(parachors[i]*(z[i]*rhos_L[i]/(MWs[i]/1000.0)
                                   - ys[i]*rhos_G[i]/(MWs[i]/1000.0)) for i in range(n))
            sigma_mix = sq**4 if sq > 0 else None
    except Exception:
        try:
            # [FIX-9] 线性回退时过滤 sigma=None 的组分
            valid_sigma = [(w[i], p["sigma"]) for i, p in enumerate(pure_props)
                           if p["sigma"] is not None]
            if valid_sigma:
                sigma_mix = sum(wi * s for wi, s in valid_sigma)
        except Exception:
            pass

    H_mix = (sum(z[i]*p["H_molar"] for i,p in enumerate(pure_props))
             if all(p["H_molar"] is not None for p in pure_props) else None)
    S_mix = (sum(z[i]*p["S_molar"] for i,p in enumerate(pure_props))
             if all(p["S_molar"] is not None for p in pure_props) else None)
    S_mix_w = None
    if S_mix is not None:
        try:
            ds = -R_GAS*sum(z[i]*math.log(z[i]) for i in range(n) if z[i] > 0)
            S_mix_w = S_mix + ds
        except Exception: pass

    Pr_L_mix = Pr_V_mix = None
    mk = MW_mix/1000.0
    if mu_L_mix and Cp_L_mix and k_L_mix and k_L_mix > 0:
        Pr_L_mix = mu_L_mix*(Cp_L_mix/mk)/k_L_mix
    if mu_V_mix and Cp_V_mix and k_V_mix and k_V_mix > 0:
        Pr_V_mix = mu_V_mix*(Cp_V_mix/mk)/k_V_mix

    gammas = calc_activity_coefficients_unifac(components, T, z) if _can_use_unifac(components) else None

    return {
        "MW_mix": round(MW_mix, 4),
        "rho_L": round(rho_L_mix, 3) if rho_L_mix else None,
        "rho_V": round(rho_V_mix, 5) if rho_V_mix else None,
        "Z_L": round(Z_L_mix, 5) if Z_L_mix else None,
        "Z_V": round(Z_V_mix, 5) if Z_V_mix else None,
        "mu_L": mu_L_mix, "mu_V": mu_V_mix,
        "k_L": k_L_mix,   "k_V": k_V_mix,
        "Cp_L": Cp_L_mix, "Cp_V": Cp_V_mix,
        "sigma_mix": sigma_mix,
        "Pr_L": round(Pr_L_mix, 4) if Pr_L_mix else None,
        "Pr_V": round(Pr_V_mix, 4) if Pr_V_mix else None,
        "H_mix": H_mix, "S_mix": S_mix, "S_mix_with_mixing": S_mix_w,
        "K": Ks, "K_method": k_method,
        "activity_coefficients": gammas,
        "per_component": pure_props,
    }


# ============================================================
# Agent 门面函数
# ============================================================
def get_fluid_density(fluid, T, P=101325.0, phase="liquid", z=None):
    """流体密度 (kg/m³)"""
    if isinstance(fluid, list):
        if z is None: raise ValueError("混合物必须提供组成 z")
        props = get_mixture_properties(fluid, z, T, P)
        val = props.get("rho_V" if phase in ("gas","vapor") else "rho_L")
        if val is None: raise ValueError(f"无法计算混合物 {fluid} 的密度")
        return val
    if phase in ("gas","vapor"):
        val = calc_gas_density(fluid, T, P)
    elif phase == "solid":
        val = calc_solid_density(fluid, T)
    else:
        val = calc_liquid_density(fluid, T, P)
    if val is None: raise ValueError(f"无法计算 '{fluid}' 的密度")
    return val

def get_fluid_viscosity(fluid, T, P=101325.0, phase="liquid", z=None):
    """流体粘度 (Pa·s)"""
    if isinstance(fluid, list):
        if z is None: raise ValueError("混合物必须提供组成 z")
        props = get_mixture_properties(fluid, z, T, P)
        val = props.get("mu_V" if phase in ("gas","vapor") else "mu_L")
        if val is None: raise ValueError(f"无法计算混合物 {fluid} 的粘度")
        return val
    val = calc_viscosity_gas(fluid, T, P) if phase in ("gas","vapor") else calc_viscosity_liquid(fluid, T, P)
    if val is None: raise ValueError(f"无法计算 '{fluid}' 的粘度")
    return val

def get_fluid_Cp(fluid, T, P=101325.0, phase="liquid", z=None):
    """流体定压热容 (J/mol/K)"""
    if isinstance(fluid, list):
        if z is None: raise ValueError("混合物必须提供组成 z")
        props = get_mixture_properties(fluid, z, T, P)
        val = props.get("Cp_V" if phase in ("gas","vapor") else "Cp_L")
        if val is None: raise ValueError(f"无法计算混合物 {fluid} 的热容")
        return val
    val = calc_heat_capacity(fluid, T, P, phase)
    if val is None: raise ValueError(f"无法计算 '{fluid}' 的热容")
    return val

def get_fluid_thermal_conductivity(fluid, T, P=101325.0, phase="liquid", z=None):
    """导热系数 (W/m/K)"""
    if isinstance(fluid, list):
        if z is None: raise ValueError("混合物必须提供组成 z")
        props = get_mixture_properties(fluid, z, T, P)
        val = props.get("k_V" if phase in ("gas","vapor") else "k_L")
        if val is None: raise ValueError(f"无法计算混合物 {fluid} 的导热系数")
        return val
    val = calc_thermal_conductivity_gas(fluid, T, P) if phase in ("gas","vapor") \
          else calc_thermal_conductivity_liquid(fluid, T, P)
    if val is None: raise ValueError(f"无法计算 '{fluid}' 的导热系数")
    return val

def get_fluid_surface_tension(fluid, T, z=None):
    """表面张力 (N/m)"""
    if isinstance(fluid, list):
        if z is None: raise ValueError("混合物必须提供组成 z")
        props = get_mixture_properties(fluid, z, T)
        val = props.get("sigma_mix")
        if val is None: raise ValueError(f"无法计算混合物 {fluid} 的表面张力")
        return val
    val = calc_surface_tension(fluid, T)
    if val is None: raise ValueError(f"无法计算 '{fluid}' 的表面张力")
    return val

def get_fluid_MW(fluid, z=None):
    """分子量 (g/mol)"""
    if isinstance(fluid, list):
        if z is None: raise ValueError("混合物必须提供组成 z")
        return get_mixture_properties(fluid, z)["MW_mix"]
    val = get_molecular_weight(fluid)
    if val is None: raise ValueError(f"无法获取 '{fluid}' 的分子量")
    return val

def get_fluid_vapor_pressure(fluid, T, z=None):
    """饱和蒸气压 / 泡点压力 (Pa)"""
    if isinstance(fluid, list):
        if z is None: raise ValueError("混合物必须提供组成 z")
        bp = calc_bubble_point_P(fluid, z, T)
        if bp and bp.get("P"): return bp["P"]
        raise ValueError(f"无法计算混合物 {fluid} 的泡点压力")
    val = calc_vapor_pressure(fluid, T)
    if val is None: raise ValueError(f"无法计算 '{fluid}' 在 {T}K 的饱和蒸气压")
    return val

def get_mixture_k_values(components, T, P, zs, method="auto"):
    """混合物 K 值 (Ki = yi/xi)"""
    return _calc_Ks(components, T, P, zs, method)


# ============================================================
# 测试
# ============================================================
if __name__ == "__main__":
    print(f"\n{'='*76}")
    print("thermo_helper — thermo + CoolProp 双引擎测试（修复版）")
    print(f"  CoolProp 可用: {HAS_COOLPROP}")
    print(f"{'='*76}")

    fmt = lambda v, f='', u='': (format(v, f) + (' ' + u if u else '')) if v is not None else 'N/A'

    print("\n[液态水 298.15K, 1atm]")
    for name in ['water', 'ethanol']:
        T, P = 298.15, 101325.0
        print(f"\n  {name}:")
        print(f"    ρ_L    = {fmt(calc_liquid_density(name,T,P), '.3f', 'kg/m³')}")
        print(f"    μ_L    = {fmt(calc_viscosity_liquid(name,T,P), '.4e', 'Pa·s')}")
        print(f"    k_L    = {fmt(calc_thermal_conductivity_liquid(name,T,P), '.4f', 'W/m/K')}")
        print(f"    Cp_L   = {fmt(calc_heat_capacity(name,T,P,'liquid'), '.2f', 'J/mol/K')}")
        print(f"    σ      = {fmt(calc_surface_tension(name,T), '.5f', 'N/m')}")
        print(f"    Psat   = {fmt(calc_vapor_pressure(name,T), '.1f', 'Pa')}")

    print("\n[CO2 超临界 320K, 10MPa — 相态检查测试]")
    T, P = 320.0, 1e7
    print(f"    σ_CO2 超临界 = {calc_surface_tension('CO2', T)} (应为 None)")
    print(f"    μ_L CO2 = {fmt(calc_viscosity_liquid('CO2',T,P), '.4e', 'Pa·s')}")
    print(f"    μ_G CO2 = {fmt(calc_viscosity_gas('CO2',T,P), '.4e', 'Pa·s')}")

    print("\n[气相扩散系数 — Chapman-Enskog + 高压 Takahashi]")
    print(f"    D(N2-O2) @ 298K, 1atm = {fmt(calc_diffusion_gas('nitrogen','oxygen',298.15,101325), '.3e', 'm²/s')} (ref ~1.8e-5)")
    print(f"    D(N2-O2) @ 298K, 50bar = {fmt(calc_diffusion_gas('nitrogen','oxygen',298.15,5e6), '.3e', 'm²/s')} (应低于低压值)")

    print("\n[液相扩散系数路由测试]")
    D_hl  = calc_diffusion_liquid_hayduk_laudie('oxygen', 'water', 298.15)
    D_hm  = calc_diffusion_liquid_hayduk_minhas('oxygen', 'water', 298.15)
    D_wc  = calc_diffusion_liquid_wilke_chang('oxygen', 'water', 298.15)
    print(f"    Hayduk-Laudie (水溶液优先) = {fmt(D_hl, '.3e', 'm²/s')} (ref ~2.4e-9)")
    print(f"    Hayduk-Minhas (非水优先)  = {fmt(D_hm, '.3e', 'm²/s')}")
    print(f"    Wilke-Chang (通用兜底)    = {fmt(D_wc, '.3e', 'm²/s')}")

    print("\n[Wilke-Chang phi 修复测试 — 非水/醇溶剂应返回 phi=1.0]")
    D_benzene = calc_diffusion_liquid_wilke_chang('naphthalene', 'benzene', 298.15)
    print(f"    D(萘/苯) @ 298K = {fmt(D_benzene, '.3e', 'm²/s')} (phi=1.0, 应非零)")

    print("\n[混合物气相密度 — FIX-5 保护测试]")
    try:
        mp = get_mixture_properties(['benzene','toluene'],[0.4,0.6],350,101325)
        print(f"    ρ_L={mp['rho_L']}  μ_L={fmt(mp['mu_L'],'.4e')}  Pr_L={mp['Pr_L']}")
    except Exception as e:
        print(f"    ❌ {e}")

    print("\n[异常测试]")
    try:
        get_fluid_density("unobtanium", 300)
    except ValueError as e:
        print(f"    ✓ 正确捕获: {e}")

    print(f"\n{'='*76}\n测试完成")