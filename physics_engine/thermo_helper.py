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
  │  蒸气压              │  CoolProp (精度最高)│  thermo (Wagner/DIPPR/Antoine) │
  │  升华压              │  thermo           │  —                   │
  │  汽化焓              │  CoolProp         │  thermo              │
  │  动力粘度 (液/气)    │  CoolProp         │  thermo              │
  │  导热系数 (液/气)    │  CoolProp         │  thermo              │
  │  表面张力            │  CoolProp (两相区)│  thermo              │
  │  Prandtl 数          │  CoolProp         │  thermo              │
  │  扩散系数 (气)       │  Chapman-Enskog   │  —                   │
  │  扩散系数 (液)       │  Wilke-Chang      │  Hayduk-Minhas       │
  │  安全/环境属性       │  thermo           │  —                   │
  │  分子描述符          │  thermo           │  —                   │
  └─────────────────────────────────────────────────────────────────┘

■ CoolProp 覆盖流体：~124 种纯流体 (NIST REFPROP 级别精度)
  不在 CoolProp 中的物质自动回退到 thermo

■ 扩散系数补充 (CoolProp 不提供):
  气相: Chapman-Enskog 理论 (精度约±5%)
  液相: Wilke-Chang 方程 (精度约±10%)  /  Hayduk-Minhas 方程

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

# thermo 库可用标志（thermo_helper 导入成功即表示 thermo 可用）
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
    """
    将任意化学品名称/CAS 转换为 CoolProp 可识别的流体名称。
    优先用 thermo 获取 CAS 号，再用 CAS 查 CoolProp（最稳健）。
    回退顺序: CAS → 原名 → 小写名
    """
    if not HAS_COOLPROP:
        return None
    # 方法1: 通过 thermo 获取 CAS，然后用 CAS 查 CoolProp
    chem = _get_chem_cached(name, 298.15, 101325.0)
    if chem is not None:
        cas = getattr(chem, 'CAS', None)
        if cas:
            try:
                _PropsSI('molar_mass', 'T', 300, 'P', 101325, cas)
                return cas
            except Exception:
                pass
    # 方法2: 直接用名称
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


# 相态分组 (CoolProp 7.x)
_CP_LIQUID_PHASES = frozenset({0, 3, 6})       # liquid, supercritical_liquid, twophase
_CP_GAS_PHASES    = frozenset({2, 5, 8})        # supercritical_gas, gas, not_imposed
# phase 1 (supercritical) 同时属于两组 — 气液不可区分


def _cp_fugacity(T: float, P: float, fluid_cp: str) -> Tuple[Optional[float], Optional[float]]:
    """
    用 CoolProp AbstractState 计算逸度 (Pa) 和逸度系数 φ。
    返回 (f, phi)，失败返回 (None, None)。
    """
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


def _cp_sat_liquid_prop(prop: str, T: float, fluid_cp: str) -> Optional[float]:
    """
    CoolProp 饱和液相 (Q=0) 物性, 按 T 求值。
    用于用户给定 (T,P) 处不存在液相 (P 低于饱和压, CoolProp 判定气相)
    但 T < Tc 的情形: 该温度下液相在饱和压力处存在,
    其物性仍有良定义 — 这是液相物性的正确求值状态。
    """
    if not HAS_COOLPROP or fluid_cp is None:
        return None
    try:
        v = _safe_float(_PropsSI(prop, 'T', T, 'Q', 0, fluid_cp))
        return v if v and v > 0 else None
    except Exception:
        return None


def _cp_sat_gas_prop(prop: str, T: float, fluid_cp: str) -> Optional[float]:
    """
    CoolProp 饱和气相 (Q=1) 物性, 按 T 求值。
    与 _cp_sat_liquid_prop 对称: 用于 (T,P) 处气相不存在 (CoolProp 判定液相)
    但 T < Tc 的情形, 取该温度饱和压力处的气相物性作为兑底。
    """
    if not HAS_COOLPROP or fluid_cp is None:
        return None
    try:
        v = _safe_float(_PropsSI(prop, 'T', T, 'Q', 1, fluid_cp))
        return v if v and v > 0 else None
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


def _is_supercritical_T(name: str, T: float) -> bool:
    """T 是否已达到/超过组分临界温度 (Tc 未知时返回 False, 不误降级)"""
    crit = get_critical_constants(name)
    return bool(crit and crit[0] and T >= crit[0])


def _select_method_in_range(prop_obj, T: float) -> Optional[str]:
    """
    在 thermo 物性对象上挑选有效范围覆盖 T 的方法。

    优先级:
      1. 当前激活方法 (若其 T_limits 覆盖 T, 保持不变避免频繁切换)
      2. 降级到范围最广 (Tmax 最高) 且覆盖 T 的方法
    无方法覆盖 T 时返回 None, 由调用方返回 None 而非静默外推出非物理值。
    """
    if prop_obj is None:
        return None
    T_limits = getattr(prop_obj, 'T_limits', None)
    method = getattr(prop_obj, 'method', None)
    if not T_limits:
        return method  # 无范围信息, 沿用当前方法
    if method and method in T_limits:
        lim = T_limits[method]
        if lim and lim[0] <= T <= lim[1]:
            return method
    covered = [(lim[1], m) for m, lim in T_limits.items()
               if lim and lim[0] <= T <= lim[1]]
    if covered:
        covered.sort(reverse=True)
        return covered[0][1]
    return None


def _eval_thermo_corr(prop_obj, T: float) -> Optional[float]:
    """
    带温度范围自动检查地求值 thermo 关联式。

    T 超出当前方法范围时自动切换到覆盖 T 的方法 (范围更广的方法兜底);
    无任何方法覆盖时返回 None, 不做静默外推。
    """
    if prop_obj is None or not callable(prop_obj):
        return None
    m = _select_method_in_range(prop_obj, T)
    if m is None:
        return None
    if getattr(prop_obj, 'method', None) != m:
        try:
            prop_obj.method = m
        except Exception:
            pass
    try:
        return _safe_float(prop_obj(T))
    except Exception:
        return None


def _supercritical_bulk_density(name: str, T: float, P: float) -> Optional[float]:
    """
    超临界流体体密度 (T ≥ Tc)。
    超临界后气液不可区分, 降级使用范围最广的方法:
      1. CoolProp HEOS (EOS 本身有效到 Tmax, 覆盖超临界区)
      2. thermo EOS 体密度
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_prop('D', T, P, cp_name)
        if val and val > 0:
            return val
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'rho', None))
        if val and val > 0:
            return val
    return None


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
# 蒸气压 & 升华压   (精度优先: CoolProp HEOS > thermo)
# ============================================================
def calc_vapor_pressure(name: str, T: float) -> Optional[float]:
    """
    蒸气压 Psat (Pa)

    路由 (精度优先):
      1. CoolProp HEOS  (NIST 级精度, 0.01~0.1%, 需 T < Tc)
      2. thermo VaporPressure callable  (Wagner/DIPPR/Antoine, 自动选覆盖 T 的方法)
      3. chem.Psat 属性                 (thermo 备用)
      4. 近临界区 (|T-Tc|<0.5K): Lee-Kesler 等宽范围方法兜底
      T ≥ Tc + 0.5K → 返回 None (远离临界区, 饱和蒸气压无物理定义)
    """
    # ── 0. 获取临界温度, 判断区域 ──
    crit = get_critical_constants(name)
    Tc = crit[0] if crit else None
    if Tc is None:
        pass  # Tc 未知, 继续尝试计算
    else:
        # T ≥ Tc + 0.5K: 远离临界区, 饱和蒸气压无物理定义
        if T >= Tc + 0.5:
            return None
        # 近临界区 (|T-Tc| < 0.5K): 用宽范围方法兜底 (如 Lee-Kesler)
        if abs(T - Tc) < 0.5:
            chem = _get_chem(name, T=T)
            if chem is not None:
                vp_obj = getattr(chem, 'VaporPressure', None)
                if vp_obj is not None:
                    # 优先用覆盖到 Tc 的宽范围方法
                    for method in ['LEE_KESLER_PSAT', 'AMBROSE_WALTON', 'DIPPR_PERRY_8E']:
                        if method in vp_obj.all_methods:
                            try:
                                vp_obj.method = method
                                val = _safe_float(vp_obj(T))
                                if val and val > 0:
                                    return val
                            except Exception:
                                pass
    # ── 1. CoolProp HEOS (精度最高) ──
    cp_name = _get_coolprop_name(name)
    if cp_name:
        try:
            # 获取临界温度，T >= Tc 时 CoolProp 无法计算饱和蒸气压
            Tc_cp = _safe_float(_PropsSI('Tcrit', '', 0, '', 0, cp_name))
            if Tc_cp is not None and T < Tc_cp:
                val = _safe_float(_PropsSI('P', 'T', T, 'Q', 0, cp_name))
                if val and val > 0:
                    return val
        except Exception:
            pass
    # ── 2. thermo VaporPressure (自动挑选覆盖 T 的方法: Wagner/DIPPR/Antoine) ──
    chem = _get_chem(name, T=T)
    if chem is not None:
        val = _eval_thermo_corr(getattr(chem, 'VaporPressure', None), T)
        if val and val > 0:
            return val
        val = _safe_float(getattr(chem, 'Psat', None))
        if val and val > 0:
            return val
    return None


def calc_sublimation_pressure(name: str, T: float) -> Optional[float]:
    """
    升华压 Psub (Pa), T < 三相点温度

    路由:
      1. thermo SublimationPressure callable
      2. Clausius-Clapeyron: ln(P/Pt) = -ΔHsub/R * (1/T - 1/Tt)
      T ≥ Tt (三相点) → 返回 None (升华压无物理定义, 不外推)
    """
    chem = _get_chem(name, T=T)
    if chem is None:
        return None
    # 温度范围检查: 高于三相点不存在升华压, 防止关联式外推
    Tt = _safe_float(getattr(chem, 'Tt', None))
    if Tt and T >= Tt:
        return None
    val = _eval_thermo_corr(getattr(chem, 'SublimationPressure', None), T)
    if val and val > 0:
        return val
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
# 策略: CoolProp 优先 (NIST 级精度), 但必须检查相态一致性
# 关键差异: 高压气体/超临界区 thermo PR EOS 误差可达 30%+
# ============================================================
def calc_liquid_density(name: str, T: float, P: float = 101325.0) -> Optional[float]:
    """
    液相质量密度 (kg/m³)

    路由:
      1. CoolProp HEOS (实验关联式，精度 <0.5%) — 仅当 CoolProp 判定为液相/两相时采用
      2. thermo VolumeLiquid + 压力修正 (精度约 1-3%)

    ⚠ 修复: CoolProp PropsSI('D','T',T,'P',P) 不区分气液，
       若 CoolProp 判定为气相则跳过，回退到 thermo 的 rhol
    ⚠ 超过临界温度: 液相不存在, 降级到 CoolProp HEOS 体密度
       (范围最广的方法, EOS 本身有效至 Tmax, 覆盖超临界区)
    """
    # T ≥ Tc → 降级到范围最广的方法: 超临界体密度
    if _is_supercritical_T(name, T):
        return _supercritical_bulk_density(name, T, P)
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        if phase is None or phase in _CP_LIQUID_PHASES or phase == 1:
            val = _cp_prop('D', T, P, cp_name)
            if val and val > 0:
                return val
        elif phase in _CP_GAS_PHASES:
            # P 低于饱和压被判定气相: 该温度下液相在饱和压处存在, 取饱和液相
            val = _cp_sat_liquid_prop('D', T, cp_name)
            if val:
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
      1. CoolProp HEOS (高压/超临界区精度显著优于 PR EOS) — 仅当 CoolProp 判定为气相/超临界气时采用
      2. thermo (含 EOS 路由)
      3. 理想气体 (仅 Pr < 0.01 时)

    ⚠ 修复: CoolProp PropsSI('D','T',T,'P',P) 不区分气液，
       若 CoolProp 判定为液相则跳过，回退到 thermo 的 rhog
    ⚠ 超过临界温度: 气液不可区分, 降级到 CoolProp HEOS 体密度
    """
    # T ≥ Tc → 降级到范围最广的方法: 超临界体密度
    if _is_supercritical_T(name, T):
        return _supercritical_bulk_density(name, T, P)
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        if phase is None or phase in _CP_GAS_PHASES or phase == 1 or phase == 3:
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
    """固相质量密度 (kg/m³)，仅 thermo"""
    chem = _get_chem(name, T=T)
    if chem is None:
        return None
    val = _safe_float(getattr(chem, 'rhos', None))
    return val if val and val > 0 else None


def calc_molar_volume(name: str, T: float, P: float,
                      phase: str = "liquid") -> Optional[float]:
    """
    摩尔体积 (m³/mol)

    路由: 从质量密度推导 (密度更可靠); 固相用 thermo
    """
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
        return (MW / 1000.0) / rho  # kg/mol / (kg/m³) = m³/mol
    return None


# ============================================================
# 压缩因子 Z
# ============================================================
def calc_compressibility_factor(name: str, T: float, P: float,
                                phase: str = "gas") -> Optional[float]:
    """
    压缩因子 Z = PV/(nRT)

    路由:
      1. CoolProp 直接输出 (基于 Helmholtz EOS，精度最高)
      2. thermo EOS 推导
      3. 维里方程 (低压气体)
      4. 理想气体 Z=1

    注: 液相 Z << 1 (水在 25°C Z ≈ 0.00074)
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
            z = _safe_float(getattr(eos, 'Z_g' if phase in ("gas","vapor") else 'Z_l', None))
            if z is not None:
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
    """
    逸度 f (Pa)

    路由:
      1. CoolProp AbstractState.fugacity(0) (最精确)
      2. thermo EOS
      3. 理想气体 f = P
    """
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
    """
    逸度系数 φ = f/P (无量纲)

    路由:
      1. CoolProp AbstractState.fugacity_coefficient(0)
      2. 从 calc_fugacity 推导
    """
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
      1. CoolProp Cpmolar (含相态检查，避免液相请求得到气相 Cp)
      2. thermo DIPPR 关联式

    ⚗ 近临界区 (Tr 0.9~1.1) Cp 急剧发散, CoolProp 处理更可靠
    """
    cp_name = _get_coolprop_name(name)
    sc = _is_supercritical_T(name, T)
    if cp_name:
        # 相态检查: 用户请求 liquid 但 CoolProp 判定为 gas 时回退 thermo
        # 超过临界温度后气液不可区分, 跳过相态检查 (范围更广的方法兜底)
        cp_phase = _cp_phase_index(T, P, cp_name)
        phase_ok = True
        if cp_phase is not None and not sc:
            if phase in ("liquid",) and cp_phase not in _CP_LIQUID_PHASES and cp_phase != 1:
                phase_ok = False  # 用户要液相，CoolProp 给气相
            elif phase in ("gas", "vapor") and cp_phase not in _CP_GAS_PHASES and cp_phase != 1:
                phase_ok = False  # 用户要气相，CoolProp 给液相
        if phase_ok:
            val = _cp_prop('Cpmolar', T, P, cp_name)
            if val and val > 0:
                return val
        elif phase == "liquid" and cp_phase in _CP_GAS_PHASES:
            # P 低于饱和压被判定气相: 该温度下液相在饱和压处存在, 取饱和液相
            val = _cp_sat_liquid_prop('Cpmolar', T, cp_name)
            if val:
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
    if (val is None or val <= 0) and sc and phase != "solid":
        # 超过临界温度: 液相关联式超出范围失效, 降级到气相 DIPPR 关联式 (范围更广)
        val = _safe_float(getattr(chem, 'Cpgm', None))
    return val if val and val > 0 else None


def calc_Cv(name: str, T: float, P: float = 101325.0,
            phase: str = "gas") -> Optional[float]:
    """
    定容热容 Cv (J/mol/K)

    路由:
      1. CoolProp Cvmolar (液相精确, 气相精确)
      2. 气相退化: Cv = Cp - R
      3. 液相热力学恒等式: Cv = Cp - T·V·β²/κ_T (需 EOS，精度有限)
    """
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
    # 超过临界温度: 液相关联式不再适用, 按气相处理 (降级到范围更广的方法)
    if phase in ("gas", "vapor") or _is_supercritical_T(name, T):
        cp = calc_heat_capacity(name, T, P, "gas")
        return (cp - R_GAS) if cp else None
    # 液相热力学恒等式
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
    """
    绝热指数 γ = Cp/Cv (无量纲)

    路由:
      气相:
        - CoolProp isentropic_expansion_coefficient (=Cp/Cv，仅气相有效)
        - 退化: Cp_g / Cv_g
      液相:
        - γ_l = Cp_l / Cv_l  (CoolProp Cpmolar / Cvmolar, 精确)
        - 注意: CoolProp 的 isentropic_expansion_coefficient 对液相
          是热力学广义定义 γ = -V/P·(∂P/∂V)_S，数值非常大，≠ Cp/Cv
          本函数液相统一用 Cp/Cv

    用途: 声速 (气相), 等熵压缩 (气相设备)
    """
    if phase in ("gas", "vapor"):
        cp_name = _get_coolprop_name(name)
        if cp_name:
            val = _cp_prop('isentropic_expansion_coefficient', T, P, cp_name)
            if val and 1.0 < val < 20.0:  # 合理气相范围
                return val
        cp = calc_heat_capacity(name, T, P, "gas")
        cv = calc_Cv(name, T, P, "gas")
        if cp and cv and cv > 0:
            return cp / cv
        return None
    else:
        # 液相: 始终用 Cp/Cv (CoolProp 给出的 isentropic_expansion_coefficient 液相不是 Cp/Cv)
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
    """
    相对焓 H (J/mol 或 J/kg)

    路由:
      1. CoolProp Hmolar / H  (基于 NIST 参考态，更一致)
      2. thermo Hm (参考态不同，仅用于差值计算)

    ⚠  CoolProp 与 thermo 参考态不同，绝对值不可交叉比较，
       但同一引擎内的差值 ΔH 是可靠的。
    """
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
    """
    内能 U (J/mol 或 J/kg)

    路由:
      1. CoolProp Umolar  (精确)
      2. 推导: U = H - PV = H - P/ρ (从 H 和密度推导)
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        key = 'Umolar' if molar else 'U'
        val = _cp_prop(key, T, P, cp_name)
        if val is not None:
            return val
    # 推导
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
    """
    Gibbs 自由能 G = H - TS (J/mol 或 J/kg)

    路由:
      1. CoolProp Gmolar
      2. G = H - T*S 推导
    """
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
    """
    Helmholtz 自由能 A = U - TS (J/mol 或 J/kg)

    路由:
      1. CoolProp A  (质量基)
      2. A = U - T*S 推导
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        # CoolProp 'A' 是质量基
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

    路由 (精度优先):
      1. CoolProp: ΔHvap = H_gas(T,Q=1) - H_liq(T,Q=0)  (NIST 级精度)
      2. thermo EnthalpyVaporization (Watson/DIPPR 106, 覆盖更广)
      T ≥ Tc → 返回 None
    """
    crit = get_critical_constants(name)
    if crit and T >= crit[0]:
        return None
    # ── 1. CoolProp (精度最高) ──
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
    # ── 2. thermo (Watson/DIPPR 106, 覆盖更广) ──
    chem = _get_chem(name, T=T)
    if chem is not None:
        val = _safe_float(getattr(chem, 'Hvapm', None))
        if val and val > 0:
            return val
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
    """标准摩尔生成焓 ΔHf° (J/mol)，由 thermo 提供

    ⚠ 不做跨相态回退：气相请求缺失 Hfgm 时不会静默用液相 Hfm 顶替
      （二者物理上相差一个汽化焓，量级可达数十 kJ/mol，不能互换）。
      缺失即返回 None，由调用方决定是否用 Hvap 做显式修正。
    """
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
                return None
            vals.append(hf)
        return sum(z[i] * vals[i] for i in range(len(fluid)))
    else:
        chem = _get_chem(fluid, T, P)
        if chem is None:
            return None
        hf = _safe_float(getattr(chem, 'Hfgm' if phase == "gas" else 'Hfm', None))
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
# 策略: CoolProp 优先 (实验关联式), 回退 thermo
# ============================================================
def calc_viscosity_liquid(name: str, T: float, P: float = 101325.0) -> Optional[float]:
    """
    液相动力粘度 μ (Pa·s)

    路由:
      1. CoolProp viscosity (NIST 关联式, 精度 <2%, 含相态一致性检查)
      2. thermo DIPPR 101 / Vogel (含高压 Lucas 修正)
      超过临界温度 → 降级到 CoolProp 体流体粘度 (范围最广)
    """
    sc = _is_supercritical_T(name, T)
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        # 超临界后气液不可区分直接接受; 否则要求 CoolProp 判定为液相
        if sc or phase is None or phase in _CP_LIQUID_PHASES or phase == 1:
            val = _cp_prop('viscosity', T, P, cp_name)
            if val and val > 0:
                return val
        elif phase in _CP_GAS_PHASES:
            # P 低于饱和压被判定气相: 该温度下液相在饱和压处存在, 取饱和液相
            val = _cp_sat_liquid_prop('viscosity', T, cp_name)
            if val:
                return val
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'mul', None))
        if val and val > 0:
            return val
        if sc:
            # 超过临界温度: 液相关联式超出范围失效, 降级到气相关联式 (范围更广)
            val = _safe_float(getattr(chem, 'mug', None))
            if val and val > 0:
                return val
    return None


def calc_viscosity_gas(name: str, T: float, P: float = 101325.0) -> Optional[float]:
    """
    气相动力粘度 μ (Pa·s)

    路由:
      1. CoolProp viscosity (含高压 Jossi-Stiel-Thodos 修正, 含相态一致性检查)
      2. thermo Chapman-Enskog / Reichenberg
      超过临界温度 → 降级到 CoolProp 体流体粘度 (范围最广)
    """
    sc = _is_supercritical_T(name, T)
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        if sc or phase is None or phase in _CP_GAS_PHASES or phase in (1, 3):
            val = _cp_prop('viscosity', T, P, cp_name)
            if val and val > 0:
                return val
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'mug', None))
        if val and val > 0:
            return val
    # thermo 也失效, 取该温度饱和气相粘度 (气相物性在该温度仍有良定义)
    if cp_name:
        val = _cp_sat_gas_prop('viscosity', T, cp_name)
        if val:
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
# ============================================================
def calc_thermal_conductivity_liquid(name: str, T: float,
                                     P: float = 101325.0) -> Optional[float]:
    """
    液相导热系数 λ (W/m/K)

    路由:
      1. CoolProp conductivity (NIST 关联式, 含近临界修正, 精度 <3%, 含相态一致性检查)
      2. thermo DIPPR 100 / Sato-Riedel
      超过临界温度 → 降级到 CoolProp 体流体导热 (范围最广)
    """
    sc = _is_supercritical_T(name, T)
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        if sc or phase is None or phase in _CP_LIQUID_PHASES or phase == 1:
            val = _cp_prop('conductivity', T, P, cp_name)
            if val and val > 0:
                return val
        elif phase in _CP_GAS_PHASES:
            # P 低于饱和压被判定气相: 该温度下液相在饱和压处存在, 取饱和液相
            val = _cp_sat_liquid_prop('conductivity', T, cp_name)
            if val:
                return val
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'kl', None))
        if val and val > 0:
            return val
        if sc:
            # 超过临界温度: 液相关联式超出范围失效, 降级到气相关联式 (范围更广)
            val = _safe_float(getattr(chem, 'kg', None))
            if val and val > 0:
                return val
    return None


def calc_thermal_conductivity_gas(name: str, T: float,
                                  P: float = 101325.0) -> Optional[float]:
    """
    气相导热系数 λ (W/m/K)

    路由:
      1. CoolProp conductivity (含高压 Stiel-Thodos 修正, 含相态一致性检查)
      2. thermo Eucken 法
      超过临界温度 → 降级到 CoolProp 体流体导热 (范围最广)
    """
    sc = _is_supercritical_T(name, T)
    cp_name = _get_coolprop_name(name)
    if cp_name:
        phase = _cp_phase_index(T, P, cp_name)
        if sc or phase is None or phase in _CP_GAS_PHASES or phase in (1, 3):
            val = _cp_prop('conductivity', T, P, cp_name)
            if val and val > 0:
                return val
    chem = _get_chem(name, T=T, P=P)
    if chem is not None:
        val = _safe_float(getattr(chem, 'kg', None))
        if val and val > 0:
            return val
    # thermo 也失效, 取该温度饱和气相导热 (气相物性在该温度仍有良定义)
    if cp_name:
        val = _cp_sat_gas_prop('conductivity', T, cp_name)
        if val:
            return val
    return None


# ============================================================
# 表面张力
# ============================================================
def calc_surface_tension(name: str, T: float) -> Optional[float]:
    """
    表面张力 σ (N/m)

    路由:
      1. CoolProp surface_tension (两相区界面, Q=0 调用)
      2. thermo DIPPR 106 / Macleod-Sugden
      T ≥ Tc → σ = 0
    """
    crit = get_critical_constants(name)
    if crit and T >= crit[0]:
        return 0.0
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
    """
    Prandtl 数 Pr = μ·Cp_mass/λ (无量纲)

    路由:
      1. CoolProp Prandtl (直接输出, 最一致)
      2. 从 μ, Cp, λ 手动计算
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_prop('Prandtl', T, P, cp_name)
        if val and val > 0:
            return val
    # 手动计算
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
    """
    等温压缩率 κ_T = -(1/V)(∂V/∂P)_T (1/Pa)

    路由:
      1. CoolProp isothermal_compressibility (精确, 液相误差 <1%)
      2. 理想气体: κ_T = 1/P
      3. thermo EOS 推导 (液相误差约 ±30%, 仅备用)
    """
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
    """
    等熵压缩率 κ_s = κ_T / γ (1/Pa)

    用途: 声速 c = 1/√(ρ·κ_s)
    注: 液相 κ_s 由 CoolProp 输出的 κ_T 和 Cp/Cv 计算，精度 <5%
    """
    kT = calc_isothermal_compressibility(name, T, P, phase)
    gamma = calc_isentropic_exponent(name, T, P, phase)
    if kT and gamma and gamma > 0:
        return kT / gamma
    return None


def calc_isobaric_expansion(name: str, T: float, P: float,
                            phase: str = "liquid") -> Optional[float]:
    """
    等压热膨胀系数 β = (1/V)(∂V/∂T)_P (1/K)

    路由:
      1. CoolProp isobaric_expansion_coefficient
      2. 理想气体: β = 1/T
      3. thermo EOS
    """
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
    """
    Joule-Thomson 系数 μ_JT = (∂T/∂P)_H (K/Pa)

    路由:
      1. CoolProp joule_thomson_coefficient (精确, 全范围有效含超临界)
      2. 理想气体: μ_JT = 0
      3. thermo EOS JTg / JTl
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_prop('joule_thomson_coefficient', T, P, cp_name)
        if val is not None:
            return val
    if classify_regime(name, T, P) == PropertyRegime.IDEAL_GAS:
        return 0.0
    chem = _get_chem(name, T=T, P=P)
    if chem is None:
        return None
    return _safe_float(getattr(chem, 'JTg' if phase in ("gas","vapor") else 'JTl', None))


def calc_speed_of_sound(name: str, T: float, P: float = 101325.0,
                        phase: str = "gas") -> Optional[float]:
    """
    声速 c (m/s)

    路由:
      1. CoolProp speed_of_sound (精确, 液相误差 <0.1%, 气相 <0.5%)
      2. 气相退化: c = √(γ·Z·R·T/M)
      3. 液相退化: c = √(1/(ρ·κ_s))  (κ_s 由 κ_T/γ 推导)

    ⚠ 原 thermo PR EOS 液相路线误差约 60%, 已完全由 CoolProp 替代；
      本函数的液相兜底仅在 CoolProp 缺该流体时启用，精度弱于 CoolProp。
    """
    cp_name = _get_coolprop_name(name)
    if cp_name:
        val = _cp_prop('speed_of_sound', T, P, cp_name)
        if val and val > 0:
            return val
    if phase in ("gas", "vapor"):
        # 气相退化: c = sqrt(gamma * Z * R * T / M)
        MW = get_molecular_weight(name)
        if MW:
            gamma = calc_isentropic_exponent(name, T, P, "gas")
            Z = calc_compressibility_factor(name, T, P, "gas") or 1.0
            if gamma:
                return math.sqrt(gamma * Z * R_GAS * T / (MW / 1000.0))
        return None
    else:
        # 液相退化: c = sqrt(1/(rho * kappa_s))
        rho = calc_liquid_density(name, T, P)
        kappa_s = calc_isentropic_compressibility(name, T, P, "liquid")
        if rho and rho > 0 and kappa_s and kappa_s > 0:
            return math.sqrt(1.0 / (rho * kappa_s))
        return None


def calc_second_virial_coefficient(name: str, T: float) -> Optional[float]:
    """第二维里系数 B (m³/mol)"""
    chem = _get_chem(name, T=T)
    return _safe_float(getattr(chem, 'Bvirial', None)) if chem else None


def calc_parachor(name: str, T: float) -> Optional[float]:
    """Parachor [N^0.25·m^2.75/mol]，用于表面张力混合规则"""
    chem = _get_chem(name, T=T)
    return _safe_float(getattr(chem, 'Parachor', None)) if chem else None


def calc_permittivity(name: str, T: float) -> Optional[float]:
    """相对介电常数 (无量纲)，液相"""
    chem = _get_chem(name, T=T)
    return _safe_float(getattr(chem, 'permittivity', None)) if chem else None


# ============================================================
# 扩散系数  (CoolProp 不提供, 用经典方程实现)
# ============================================================

def calc_diffusion_gas(
    name_A: str, name_B: str, T: float, P: float = 101325.0
) -> Optional[float]:
    """
    气相二元扩散系数 D_AB (m²/s) — Chapman-Enskog 理论

    公式 (Reid et al., Properties of Gases and Liquids, 5th Ed. Eq.11-3.2):
        D_AB [cm²/s] = 1.8583e-3 * T^1.5 * √(1/MA + 1/MB)
                       / (P[atm] * σ_AB[Å]² * Ω_D)

    碰撞积分 Ω_D (Neufeld-Janzen-Aziz 关联式):
        T* = k_B*T / ε_AB
        Ω_D = A/(T*^B) + C/exp(D*T*) + E/exp(F*T*) + G/exp(H*T*)

    参数来源: Stockmayer 势 (ε/k_B, σ) 从 thermo 读取
    精度: ±5% (非极性气体), ±10% (极性气体)

    ⚠ 修复: 原实现里先按错误单位算了一遍 D_AB (σ 用米、P 用 Pa 直接
      套 CGS 系数) 但从未使用该结果，属于死代码，已删除；同时把
      P_bar 改为标准公式要求的 P_atm，避免 ~1.3% 的系统偏差。
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
    sigma_AB = (sigma_A + sigma_B) / 2.0          # Å
    eps_AB   = math.sqrt(eps_A * eps_B)            # K
    T_star   = T / eps_AB
    # Neufeld-Janzen-Aziz Ω_D
    A, B, C, D, E, F = 1.06036, 0.15610, 0.19300, 0.47635, 1.03587, 1.52996
    G, H = 1.76474, 3.89411
    try:
        Omega_D = (A / T_star**B + C / math.exp(D * T_star)
                   + E / math.exp(F * T_star) + G / math.exp(H * T_star))
    except Exception:
        return None
    # Chapman-Enskog 标准公式: P 用 atm, σ 用 Å, M 用 g/mol → cm²/s
    P_atm = P / 101325.0
    D_cm2_s = (1.8583e-3 * T**1.5 * math.sqrt(1.0/MA + 1.0/MB)
               / (P_atm * sigma_AB**2 * Omega_D))
    return D_cm2_s * 1e-4  # cm²/s → m²/s


def _get_Vb_cm3_mol(solute: str) -> Optional[float]:
    """
    获取溶质在正常沸点下的液相摩尔体积 V_b (cm³/mol)，用于扩散系数计算。

    路由:
      1. CoolProp: PropsSI('D', 'T', Tb, 'Q', 0, fluid) → 两相边界液相密度
      2. thermo: Chemical(T=Tb-1K).rhol (在沸点略低处取值)
    注: thermo 在精确 Tb 点判定为气相导致 rhol=None，需偏移 1K
    """
    chem_A = _get_chem(solute)
    if chem_A is None:
        return None
    MW_A = _safe_float(getattr(chem_A, 'MW', None))
    Tb_A = _safe_float(getattr(chem_A, 'Tb', None))
    if not MW_A or not Tb_A:
        return None
    # 方法1: CoolProp 两相边界液相密度 (最可靠)
    cp_name = _get_coolprop_name(solute)
    if cp_name and HAS_COOLPROP:
        try:
            rho_tb = _safe_float(_PropsSI('D', 'T', Tb_A, 'Q', 0, cp_name))
            if rho_tb and rho_tb > 0:
                return MW_A / rho_tb * 1000.0  # cm³/mol
        except Exception:
            pass
    # 方法2: thermo 在 Tb-1K 处取液相密度 (避免边界判定为气相)
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
    液相无限稀释扩散系数 D°_AB (m²/s) — Wilke-Chang 方程

    公式 (Wilke-Chang, AIChE J. 1955):
        D°_AB = 7.4e-8 * (φ·M_B)^0.5 * T / (μ_B · V_A^0.6)

    参数:
        φ   = 溶剂缔合因子: 水=2.6, 甲醇=1.9, 乙醇=1.5, 苯/非极性=1.0
        M_B = 溶剂摩尔质量 (g/mol)
        μ_B = 溶剂粘度 (cP = mPa·s)
        V_A = 溶质在正常沸点下的液相摩尔体积 (cm³/mol)
            → 通过 _get_Vb_cm3_mol 获取 (CoolProp 两相边界优先)

    精度: ±10% (水为溶剂), ±25% (有机溶剂)
    结果单位: m²/s
    超过溶剂临界温度 → 返回 None (液相不存在, 关联式不适用)
    """
    if _is_supercritical_T(solvent, T):
        return None
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

    # 溶剂缔合因子 φ
    cas_B = getattr(chem_B, 'CAS', None)
    solvent_assoc = {'7732-18-5': 2.6, '67-56-1': 1.9, '64-17-5': 1.5}
    phi = solvent_assoc.get(cas_B or '', 0) or           {'water': 2.6, 'methanol': 1.9, 'ethanol': 1.5}.get(solvent.lower(), 1.0)

    D_cm2_s = 7.4e-8 * math.sqrt(phi * M_B) * T / (mu_B_cP * V_A**0.6)
    return D_cm2_s * 1e-4  # cm²/s → m²/s


def calc_diffusion_liquid_hayduk_minhas(
    solute: str, solvent: str, T: float
) -> Optional[float]:
    """
    液相无限稀释扩散系数 D°_AB (m²/s) — Hayduk-Minhas 方程

    Hayduk & Minhas (Can. J. Chem. Eng. 1982) 给出两套系数：
      · 水溶液体系:
            D°_AB = 1.25e-8 * (V_A^-0.19 - 0.292) * T^1.52 * μ_B^(9.58/V_A - 1.12)
      · 非水（正常烷烃等）有机溶剂体系:
            D°_AB = 13.3e-8 * T^1.47 * μ_B^(10.2/V_A - 0.791) / V_A^0.71

    ⚠ 修复: 原实现只有非水一套系数，却在测试用例里把它用在
      水溶剂体系上，二者不是同一关联式，精度会打折扣。现在按
      solvent 是否为水自动切换公式。

    V_A 获取路由同 Wilke-Chang (见 _get_Vb_cm3_mol)。
    精度: ±10% (对口场景：水溶液或非极性有机溶剂)
    超过溶剂临界温度 → 返回 None (液相不存在, 关联式不适用)
    """
    if _is_supercritical_T(solvent, T):
        return None
    mu_B = calc_viscosity_liquid(solvent, T)
    if not mu_B or mu_B <= 0:
        return None
    mu_B_cP = mu_B * 1000.0

    V_A = _get_Vb_cm3_mol(solute)
    if V_A is None or V_A <= 0:
        return None

    chem_B = _get_chem(solvent, T=T)
    cas_B = getattr(chem_B, 'CAS', None) if chem_B else None
    is_water = (cas_B == '7732-18-5') or (solvent.lower().strip() in ('water', 'h2o'))

    try:
        if is_water:
            exp_mu = 9.58 / V_A - 1.12
            D_cm2_s = 1.25e-8 * (V_A**-0.19 - 0.292) * T**1.52 * mu_B_cP**exp_mu
        else:
            exp_mu = 10.2 / V_A - 0.791
            D_cm2_s = 13.3e-8 * T**1.47 * mu_B_cP**exp_mu / V_A**0.71
    except Exception:
        return None
    return D_cm2_s * 1e-4  # cm²/s → m²/s


def calc_diffusion_self(name: str, T: float, P: float = 101325.0) -> Optional[float]:
    """
    自扩散系数 D_self (m²/s) — 气相近似

    气相: 用 A-A 自扩散，Chapman-Enskog 简化形式
    精度: 定性参考
    """
    return calc_diffusion_gas(name, name, T, P)


# ============================================================
# 独立属性函数 (相变温度 / 分子常数)
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
    """折射率 nD (589 nm, 液相 20°C)"""
    chem = _get_chem(name)
    return _safe_float(getattr(chem, 'RI', None)) if chem else None


def get_standard_entropy(name: str, phase: str = "gas") -> Optional[float]:
    """标准摩尔熵 S° (J/mol/K)，298.15 K, 101325 Pa"""
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
# Wilson / NRTL 二元交互参数数据库
# ============================================================
_WILSON_NRTL_DB: Dict[Tuple[str, str], Dict] = {
    ("ethanol", "water"): {
        "wilson": {"Lambda_12": 526.06, "Lambda_21": 1447.86},
        "nrtl":   {"tau_12": 631.05, "tau_21": 1197.41, "alpha": 0.3},
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
        "wilson": {"Lambda_12": 0.0, "Lambda_21": 0.0},
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
# Henry 常数数据库 (Sander 2015, Atmos. Chem. Phys.)
# 适用于稀薄气体溶解于水，摩尔分数基准
# 关联式: ln(H_x / Pa) = A - B * (1/T - 1/T_ref)
#   T_ref = 298.15 K, A = ln(H_ref)
#   B 直接取自 Sander(2015) 表中 d ln(H^cp)/d(1/T) 系数 (均为正值)
#
# ⚠ 修复说明:
#   Sander(2015) 表里的 B 是对 H^cp [mol/(m³·Pa)]（溶解度型亨利常数，
#   数值越大代表越易溶）定义的: d ln(H^cp)/d(1/T) = B (>0)，即温度越低
#   （1/T 越大）该气体越易溶、H^cp 越大。
#   本模块用的是 H_x [Pa]（分压型: p = H_x·x, K = H_x/P），它与 H^cp
#   互为倒数关系，因此温度依赖的符号要反过来：
#       ln(H_x) = ln(H_x,ref) - B * (1/T - 1/T_ref)
#   原代码写成 "+ B*(...)"，导致算出的 H_x（进而 K=H_x/P）随温度升高
#   反而下降，方向与实际相反（CO2/O2/N2/CH4 等常见气体在 0~150°C 范围
#   内都是温度越高溶解度越低、H_x 应越大）。现已改为减号。
# ============================================================
_T_HENRY_REF = 298.15  # K

# {canonical_name: (H_ref [Pa] @298.15K, B [K], description)}
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

# 别名映射 → canonical name
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

# 轻气体集合（用于检测气体-水体系）
_LIGHT_GAS_NAMES: frozenset = frozenset(_HENRY_WATER_DB.keys())
_WATER_NAMES: frozenset = frozenset({'water', 'Water', 'H2O', 'h2o'})


def _resolve_henry_name(name: str) -> Optional[str]:
    """将物质名称映射到 Henry 数据库的 canonical name。"""
    if name in _HENRY_WATER_DB:
        return name
    lower = name.lower().strip()
    mapped = _HENRY_ALIASES.get(lower)
    if mapped:
        return mapped
    # 尝试 thermo 名称归一化
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
    """检测体系是否为气体-水混合物（存在轻气体 + 水）。"""
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
    计算气体在水中的亨利常数 H_x [Pa]（摩尔分数基准, p = H_x·x）。

    基于 Sander (2015) 文献关联式（已修正符号，见模块级注释）：
        ln(H_x / Pa) = A - B * (1/T - 1/T_ref)
        T_ref = 298.15 K

    Parameters
    ----------
    name : str
        气体名称（如 'CO2', 'H2S', 'N2', 'O2', 'CH4' 等）
    T : float
        温度 [K]，有效范围 273-373 K
    solvent : str
        溶剂名称，当前仅支持 'water'

    Returns
    -------
    float or None
        亨利常数 H_x [Pa]，无数据时返回 None

    Notes
    -----
    K_i = H_x / P（当溶剂为水、稀溶液条件下）
    对于非理想性强的气体-水体系，Henry 定律远优于 Raoult 定律。
    H_x 随温度升高而增大（多数常见气体在 0~150°C 范围内溶解度随温度
    升高而降低），例如 CO2: 25°C≈1670atm → 50°C≈2960atm → 100°C≈5710atm。

    References
    ----------
    Sander, R. (2015). Atmos. Chem. Phys., 15, 4399-4981.
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

    # ---- Henry 定律: 气体-水体系 ----
    if method == "henry":
        Ks = []
        for c in components:
            canon = _resolve_henry_name(c)
            if canon is not None:
                H = calc_henry_constant(c, T)
                Ks.append(H / P if H else 1.0)
            else:
                # 溶剂 (水) 或其他非气体组分 → Raoult
                Ks.append((_safe_float(calc_vapor_pressure(c, T)) or P) / P)
        return Ks

    # ---- EOS (Peng-Robinson / SRK) ----
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
                return [math.exp(float(eos.lnphis_l[i]) - float(eos.lnphis_g[i]))
                        for i in range(len(components))]
            except Exception:
                pass

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

    # ---- UNIFAC 活度系数法 ----
    if method == "unifac" and UNIFAC_gammas_fn:
        try:
            cgs = [_get_chem(c).UNIFAC_groups for c in components]
            gammas = UNIFAC_gammas_fn(T=T, xs=zs, chemgroups=cgs)
            return [gammas[i] * (calc_vapor_pressure(components[i], T) or 0.0) / P
                    for i in range(len(components))]
        except Exception:
            pass

    # ---- 兜底 Raoult ----
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
    """泡点温度，返回 {T, T_C, y, K, method, converged, iterations}"""
    n = len(components)
    _tol_map = {"normal": 1e-4, "high": 1e-6, "ultra": 1e-8}
    _tol = _tol_map.get(precision, tol)
    method_used = _auto_k_method(components, 400.0, P) if method == "auto" else method
    try:
        m = Mixture(components, zs=x, P=P)
        Tb = float(m.Tbubble)
        return {"T": round(Tb, 4), "T_C": round(Tb-273.15, 4), "y": list(m.ys or []),
                "K": [m.ys[i]/x[i] if x[i] > 0 else 0 for i in range(n)],
                "method": "Mixture", "converged": True, "iterations": 0}
    except Exception:
        pass
    T_init = sum(x[i] * (get_normal_boiling_point(c) or 350.0) for i, c in enumerate(components))
    T_init = max(100.0, min(1200.0, T_init))
    T_prev, T_curr = T_init, T_init
    f_prev = None
    iter_count = 0
    for iter_i in range(max_iter):
        iter_count = iter_i + 1
        Ks = _calc_Ks(components, T_curr, P, x, method_used)
        sum_Kx = sum(Ks[i]*x[i] for i in range(n))
        f_curr = sum_Kx - 1.0
        if abs(f_curr) < _tol:
            break
        if f_prev is not None and abs(f_curr - f_prev) > 1e-12:
            T_next = T_curr - f_curr * (T_curr - T_prev) / (f_curr - f_prev)
        else:
            T_next = T_curr - f_curr * T_curr * 0.03
        T_next = max(100.0, min(1200.0, T_next))
        T_prev, T_curr = T_curr, T_next
        f_prev = f_curr
        if iter_i > 30 and method_used in ("wilson", "nrtl", "unifac", "raoult"):
            all_near_1 = all(abs(Ks[i] - 1.0) < 0.3 for i in range(n))
            if all_near_1:
                try:
                    from scipy.optimize import brentq
                    def _bubble_residual(T_val):
                        Ks_v = _calc_Ks(components, T_val, P, x, method_used)
                        return sum(Ks_v[i]*x[i] for i in range(n)) - 1.0
                    T_lo, T_hi = min(T_curr, T_curr - 50), max(T_curr, T_curr + 50)
                    T_lo = max(100.0, T_lo)
                    T_hi = min(1200.0, T_hi)
                    if _bubble_residual(T_lo) * _bubble_residual(T_hi) < 0:
                        T_root = brentq(_bubble_residual, T_lo, T_hi, xtol=_tol)
                        T_curr = T_root
                        iter_count += 1
                        break
                except Exception:
                    pass
    Ks = _calc_Ks(components, T_curr, P, x, method_used)
    y = [Ks[i]*x[i] for i in range(n)]
    sum_y = sum(y)
    if sum_y > 0:
        y = [yi / sum_y for yi in y]
    aze_near = all(abs(Ks[i] - 1.0) < 0.1 for i in range(n)) if n > 1 else False
    return {"T": round(T_curr, 6), "T_C": round(T_curr-273.15, 6), "y": y, "K": Ks,
            "method": method_used,
            "converged": abs(sum(Ks[i]*x[i] for i in range(n)) - 1.0) < _tol,
            "iterations": iter_count, "azeotrope_near": aze_near}


def calc_dew_point_T(components: List[str], y: List[float], P: float,
                     method: str = "auto", tol: float = 1e-4,
                     max_iter: int = 200,
                     precision: str = "normal") -> Dict:
    """露点温度，返回 {T, T_C, x, K, method, converged, iterations}"""
    n = len(components)
    _tol_map = {"normal": 1e-4, "high": 1e-6, "ultra": 1e-8}
    _tol = _tol_map.get(precision, tol)
    method_used = _auto_k_method(components, 400.0, P) if method == "auto" else method
    try:
        m = Mixture(components, zs=y, P=P)
        Td = float(m.Tdew)
        return {"T": round(Td,4), "T_C": round(Td-273.15,4), "x": list(m.xs or []),
                "K": [y[i]/m.xs[i] if m.xs[i] > 0 else 0 for i in range(n)],
                "method": "Mixture", "converged": True, "iterations": 0}
    except Exception:
        pass
    T_init = sum(y[i] * (get_normal_boiling_point(c) or 360.0) for i, c in enumerate(components))
    T_init = max(100.0, min(1200.0, T_init))
    T_prev, T_curr = T_init, T_init
    f_prev = None
    iter_count = 0
    for iter_i in range(max_iter):
        iter_count = iter_i + 1
        Ks = _calc_Ks(components, T_curr, P, y, method_used)
        sum_yK = sum(y[i]/Ks[i] for i in range(n))
        f_curr = sum_yK - 1.0
        if abs(f_curr) < _tol:
            break
        if f_prev is not None and abs(f_curr - f_prev) > 1e-12:
            T_next = T_curr - f_curr * (T_curr - T_prev) / (f_curr - f_prev)
        else:
            T_next = T_curr - f_curr * T_curr * 0.03
        T_next = max(100.0, min(1200.0, T_next))
        T_prev, T_curr = T_curr, T_next
        f_prev = f_curr
    Ks = _calc_Ks(components, T_curr, P, y, method_used)
    x = [y[i]/Ks[i] for i in range(n)]
    sum_x = sum(x)
    if sum_x > 0:
        x = [xi / sum_x for xi in x]
    return {"T": round(T_curr, 6), "T_C": round(T_curr-273.15, 6), "x": x, "K": Ks,
            "method": method_used,
            "converged": abs(sum(y[i]/Ks[i] for i in range(n)) - 1.0) < _tol,
            "iterations": iter_count}


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


def calc_VLE_envelope(components: List[str], P: float,
                      n_points: int = 20, method: str = "auto",
                      precision: str = "normal") -> Dict:
    """生成等压 VLE 包络线（泡点线 + 露点线）。"""
    if len(components) != 2:
        return {"error": "仅支持二元体系", "bubble_line": [], "dew_line": []}
    c1, c2 = components
    bubble_line, dew_line = [], []
    x_vals = [i / (n_points - 1) for i in range(n_points)]
    for xi in x_vals:
        x = [xi, 1.0 - xi]
        res_b = calc_bubble_point_T(components, x, P, method=method, precision=precision)
        if res_b["converged"]:
            bubble_line.append({"x": x, "T": res_b["T"], "T_C": res_b["T_C"],
                                "y": res_b["y"], "K": res_b["K"]})
    for yi in x_vals:
        y = [yi, 1.0 - yi]
        res_d = calc_dew_point_T(components, y, P, method=method, precision=precision)
        if res_d["converged"]:
            dew_line.append({"y": y, "T": res_d["T"], "T_C": res_d["T_C"],
                             "x": res_d["x"], "K": res_d["K"]})
    azeotrope = None
    if bubble_line:
        min_T_point = min(bubble_line, key=lambda p: p["T"])
        min_T_C = min_T_point["T_C"]
        T_pure = []
        for xi in [0.0, 1.0]:
            x = [xi, 1.0 - xi]
            res = calc_bubble_point_T(components, x, P, method=method, precision=precision)
            T_pure.append(res["T_C"])
        if min_T_C < min(T_pure) - 0.5:
            z_aze = min_T_point["x"]
            azeotrope = {"T_C": min_T_C, "T": min_T_point["T"],
                         "z": z_aze, "type": "minimum"}
        else:
            max_T_point = max(bubble_line, key=lambda p: p["T"])
            max_T_C = max_T_point["T_C"]
            if max_T_C > max(T_pure) + 0.5:
                z_aze = max_T_point["x"]
                azeotrope = {"T_C": max_T_C, "T": max_T_point["T"],
                             "z": z_aze, "type": "maximum"}
            else:
                azeotrope = {"note": "未检测到明显共沸点"}
    return {"bubble_line": bubble_line, "dew_line": dew_line,
            "azeotrope": azeotrope, "P": P}


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
        return {"beta": 0.0, "x": list(z), "y": list(z), "K": Ks, "phase": "liquid", "method": method_used}
    if K_min >= 1.0:
        return {"beta": 1.0, "x": list(z), "y": list(z), "K": Ks, "phase": "vapor", "method": method_used}
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
    return {"beta": round(beta, 6), "x": x, "y": y, "K": Ks, "phase": phase, "method": method_used}


# ============================================================
# 纯物质综合属性
# ============================================================
def get_chemical_properties(component_name: str, T: float = 298.15,
                            P: float = 101325.0) -> Dict:
    """
    纯物质全面物性字典 (thermo + CoolProp 双引擎)

    各字段来源说明见模块 docstring 中的引擎分工表。
    标注 [CP] 表示 CoolProp 主导，[TH] 表示 thermo 主导。

    主要字段单位:
      MW g/mol | Tc/Tb/Tm/Tt K | Pc/Pt Pa | Vc m³/mol | omega/Zc -
      rho_L/G/S kg/m³ | Vm_L/G m³/mol | Z_L/G -
      mu_L/G Pa·s | k_L/G W/m/K | nu_L/G m²/s | alpha_L/G m²/s
      Pr_L/G - | sigma N/m
      Psat/Psub Pa | Hvap/Hfus/Hsub J/mol
      Cp_L/G/S J/mol/K | Cv_G/L J/mol/K | gamma_G/L -
      c_G/L m/s | kappa_s_L 1/Pa
      H/S/U/G/A_molar J/mol | Hf_gas/liq J/mol | Gf_gas/liq J/mol
      Hc_HHV/LHV J/mol | S0_gas/liq J/mol/K
      Bvirial m³/mol | JT_G/L K/Pa | beta_G/L 1/K | kappa_G/L 1/Pa
      f_G/L Pa | phi_G/L -
    """
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

    # CoolProp 批量读取 (单次构建 AbstractState 更高效)
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

    # --- 密度 & 体积 ---
    rho_L = calc_liquid_density(component_name, T, P)
    rho_G = calc_gas_density(component_name, T, P)
    rho_S = calc_solid_density(component_name, T)
    Vm_L  = calc_molar_volume(component_name, T, P, "liquid")
    Vm_G  = calc_molar_volume(component_name, T, P, "gas")

    # --- 压缩因子 ---
    Z_val = cp_props.get('Z')
    Z_L = Z_val if (Z_val and classify_regime(component_name, T, P) == PropertyRegime.LIQUID) else None
    Z_G = Z_val if (Z_val and classify_regime(component_name, T, P) != PropertyRegime.LIQUID) else None
    if Z_L is None: Z_L = calc_compressibility_factor(component_name, T, P, "liquid")
    if Z_G is None: Z_G = calc_compressibility_factor(component_name, T, P, "gas")

    # --- 传输 ---
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

    # --- 热容 & 声速 ---
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

    # --- 热力学势 (用 is None 判断, 防止合法值 0 被 or 吞掉) ---
    H_mol = cp_props.get('Hmolar')
    if H_mol is None: H_mol = calc_enthalpy(component_name, T, P)
    S_mol = cp_props.get('Smolar')
    if S_mol is None: S_mol = calc_entropy(component_name, T, P)
    U_mol = cp_props.get('Umolar')
    if U_mol is None: U_mol = calc_internal_energy(component_name, T, P)
    G_mol = cp_props.get('Gmolar')
    if G_mol is None: G_mol = calc_gibbs_energy(component_name, T, P)
    A_mol = calc_helmholtz_energy(component_name, T, P)

    # --- 生成/燃烧 (thermo) ---
    Hf_gas = get_formation_enthalpy(component_name, T, P, "gas")
    Hf_liq = get_formation_enthalpy(component_name, T, P, "liquid")
    Gf_gas = get_formation_gibbs(component_name, "gas")
    Gf_liq = get_formation_gibbs(component_name, "liquid")
    Hc_HHV = get_combustion_enthalpy(component_name, "gas", True)
    Hc_LHV = get_combustion_enthalpy(component_name, "gas", False)
    S0_gas = _safe_float(getattr(chem, 'S0gm', None)) if chem else None
    S0_liq = _safe_float(getattr(chem, 'S0m',  None)) if chem else None

    # --- 推导量 ---
    Bvirial  = calc_second_virial_coefficient(component_name, T)
    JT_G     = calc_joule_thomson_coefficient(component_name, T, P, "gas")
    JT_L     = calc_joule_thomson_coefficient(component_name, T, P, "liquid")
    beta_G   = calc_isobaric_expansion(component_name, T, P, "gas")
    beta_L   = calc_isobaric_expansion(component_name, T, P, "liquid")
    kappa_G  = calc_isothermal_compressibility(component_name, T, P, "gas")
    kappa_L  = calc_isothermal_compressibility(component_name, T, P, "liquid")

    # --- 逸度 ---
    f_cp   = cp_props.get('fugacity')
    phi_cp = cp_props.get('phi')

    # --- 相变 ---
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
        # 密度 & 体积
        "rho_L": _r(rho_L, 3), "rho_G": _r(rho_G, 5), "rho_S": _r(rho_S, 3),
        "Vm_L": Vm_L, "Vm_G": Vm_G, "Z_L": _r(Z_L, 6), "Z_G": _r(Z_G, 6),
        # 传输
        "mu_L": mu_L, "mu_G": mu_G, "nu_L": nu_L, "nu_G": nu_G,
        "k_L": k_L, "k_G": k_G, "alpha_L": alpha_L, "alpha_G": alpha_G,
        "Pr_L": _r(Pr_L, 4), "Pr_G": _r(Pr_G, 4),
        "sigma": sigma, "permittivity": perm, "Parachor": parachor,
        # 相变
        "Psat": Psat, "Psub": Psub, "Hvap": Hvap, "Hfus": Hfus, "Hsub": Hsub,
        # 热容 & 速度
        "Cp_L": Cp_L, "Cp_G": Cp_G, "Cp_S": Cp_S,
        "Cv_G": Cv_G, "Cv_L": Cv_L,
        "gamma_G": _r(gamma_G, 4), "gamma_L": _r(gamma_L, 4),
        "c_G": _r(c_G, 2), "c_L": _r(c_L, 2), "kappa_s_L": kappa_s_L,
        # 热力学势
        "H_molar": H_mol, "S_molar": S_mol, "U_molar": U_mol,
        "G_molar": G_mol, "A_molar": A_mol,
        # 生成/燃烧
        "Hf_gas": Hf_gas, "Hf_liq": Hf_liq,
        "Gf_gas": Gf_gas, "Gf_liq": Gf_liq,
        "Hc_HHV": Hc_HHV, "Hc_LHV": Hc_LHV,
        "S0_gas": S0_gas, "S0_liq": S0_liq,
        # 推导量
        "Bvirial": Bvirial,
        "JT_G": JT_G, "JT_L": JT_L,
        "beta_G": beta_G, "beta_L": beta_L,
        "kappa_G": kappa_G, "kappa_L": kappa_L,
        # 逸度
        "f_G": f_cp, "f_L": f_cp,
        "phi_G": phi_cp, "phi_L": phi_cp,
        "source": f"thermo+CoolProp[{cp_name}]" if cp_name else "thermo",
    }


# ============================================================
# 气相混合物粘度 — Wilke 方法 (含 phi_ij 相互作用因子)
# ============================================================
def _wilke_gas_viscosity(mu_i: List[float], M_i: List[float],
                         z: List[float]) -> Optional[float]:
    """
    Wilke (1950) 气体混合物粘度混合规则:
        mu_mix = Σ_i [ z_i * mu_i / Σ_j (z_j * phi_ij) ]
        phi_ij = [1 + (mu_i/mu_j)^0.5 * (M_j/M_i)^0.25]^2
                 / [8*(1 + M_i/M_j)]^0.5
    比简单的 sqrt(M) 加权平均更严格，是气体粘度混合的标准工程方法。
    """
    n = len(mu_i)
    if n == 0 or any(v is None or v <= 0 for v in mu_i) or any(v is None or v <= 0 for v in M_i):
        return None
    try:
        phi = [[0.0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                num = (1.0 + math.sqrt(mu_i[i]/mu_i[j]) * (M_i[j]/M_i[i])**0.25) ** 2
                den = math.sqrt(8.0 * (1.0 + M_i[i]/M_i[j]))
                phi[i][j] = num / den
        mu_mix = 0.0
        for i in range(n):
            denom = sum(z[j]*phi[i][j] for j in range(n))
            if denom <= 0:
                return None
            mu_mix += z[i]*mu_i[i]/denom
        return mu_mix
    except Exception:
        return None


# ============================================================
# 混合物综合属性
# ============================================================
def get_mixture_properties(components: List[str], z: List[float],
                           T: float = 298.15, P: float = 101325.0) -> Dict:
    """
    混合物物性字典

    混合规则:
      密度  : PRMIX EOS → 体积加权回退
      粘度  : 液相对数加权 (Arrhenius型); 气相 Wilke 方法 (含 phi_ij)
      导热  : 质量加权
      Cp    : 摩尔分数加权
      表面张力: Macleod-Sugden (Parachor, 用混合物摩尔密度) → 线性回退
      K值   : _auto_k_method 路由
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
            inv = sum(w[i]/p["rho_G"] for i,p in enumerate(pure_props) if p["rho_G"])
            rho_V_mix = 1.0/inv if inv > 0 else P*MW_mix/1000.0/(R_GAS*T)
            Z_V_mix = P*(MW_mix/1000.0)/(rho_V_mix*R_GAS*T)
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
        if all(p["mu_G"] for p in pure_props):
            mu_V_mix = _wilke_gas_viscosity(
                [p["mu_G"] for p in pure_props], MWs, z)
    except Exception: pass

    k_L_mix = (sum(w[i]*p["k_L"] for i,p in enumerate(pure_props))
               if all(p["k_L"] for p in pure_props) else None)
    k_V_mix = (sum(w[i]*p["k_G"] for i,p in enumerate(pure_props))
               if all(p["k_G"] for p in pure_props) else None)
    Cp_L_mix = (sum(z[i]*p["Cp_L"] for i,p in enumerate(pure_props))
                if all(p["Cp_L"] for p in pure_props) else None)
    Cp_V_mix = (sum(z[i]*p["Cp_G"] for i,p in enumerate(pure_props))
                if all(p["Cp_G"] for p in pure_props) else None)

    # 表面张力 (Macleod-Sugden, 用混合物摩尔密度而非纯组分密度)
    #   sigma_mix^(1/4) = Σ_i Pi * (x_i*rho_L_mix/M_mix - y_i*rho_V_mix/M_mix)
    sigma_mix = None
    try:
        parachors = [p["Parachor"] for p in pure_props]
        if all(parachors) and rho_L_mix and rho_V_mix and MW_mix > 0:
            ys = [Ks[i]*z[i] for i in range(n)]
            sy = sum(ys); ys = [v/sy for v in ys] if sy > 0 else ys
            molar_rho_L = rho_L_mix / (MW_mix/1000.0)   # mol/m³
            molar_rho_V = rho_V_mix / (MW_mix/1000.0)   # mol/m³
            sq = sum(parachors[i]*(z[i]*molar_rho_L - ys[i]*molar_rho_V) for i in range(n))
            sigma_mix = sq**4 if sq > 0 else None
    except Exception:
        pass
    if sigma_mix is None:
        try:
            sigma_mix = sum(w[i]*p["sigma"] for i,p in enumerate(pure_props)
                           if p["sigma"] is not None)
        except Exception: pass

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
    print("thermo_helper — thermo + CoolProp 双引擎测试 (修复版)")
    print(f"  CoolProp 可用: {HAS_COOLPROP}")
    print(f"{'='*76}")

    fmt = lambda v, f='', u='': (format(v, f) + (' ' + u if u else '')) if v is not None else 'N/A'

    # 0. 亨利常数符号修复验证
    print("\n[Henry 常数修复验证] CO2 / O2 随温度应升高")
    for gname in ['CO2', 'O2']:
        print(f"  {gname}:")
        for T in [273.15, 298.15, 323.15, 373.15]:
            H = calc_henry_constant(gname, T)
            print(f"    T={T:.2f}K  H_x={fmt(H,'.3e','Pa')}  ({fmt(H/101325 if H else None,'.0f')} atm)")

    # 1. 精度对比: liquid water
    print("\n[液态水 298.15K, 1atm] — 与 NIST 对照")
    for name in ['water', 'ethanol']:
        T, P = 298.15, 101325.0
        print(f"\n  {name}:")
        print(f"    ρ_L    = {fmt(calc_liquid_density(name,T,P), '.3f', 'kg/m³')}")
        print(f"    μ_L    = {fmt(calc_viscosity_liquid(name,T,P), '.4e', 'Pa·s')}")
        print(f"    k_L    = {fmt(calc_thermal_conductivity_liquid(name,T,P), '.4f', 'W/m/K')}")
        print(f"    Cp_L   = {fmt(calc_heat_capacity(name,T,P,'liquid'), '.2f', 'J/mol/K')}")
        print(f"    Cv_L   = {fmt(calc_Cv(name,T,P,'liquid'), '.2f', 'J/mol/K')}")
        print(f"    γ_L    = {fmt(calc_isentropic_exponent(name,T,P,'liquid'), '.4f')}")
        print(f"    c_L    = {fmt(calc_speed_of_sound(name,T,P,'liquid'), '.1f', 'm/s')} (water ref=1497)")
        print(f"    σ      = {fmt(calc_surface_tension(name,T), '.5f', 'N/m')}")
        print(f"    Pr_L   = {fmt(calc_prandtl_number(name,T,P,'liquid'), '.3f')}")
        print(f"    κ_T    = {fmt(calc_isothermal_compressibility(name,T,P,'liquid'), '.3e', '1/Pa')}")
        print(f"    φ      = {fmt(calc_fugacity_coefficient(name,T,P,'liquid'), '.5f')}")
        print(f"    Psat   = {fmt(calc_vapor_pressure(name,T), '.1f', 'Pa')}")

    # 2. 高压气体: CO2 超临界
    print("\n[CO2 超临界 280, 1bar] — 高压区精度对比")
    T, P = 280.0, 1e6
    name = 'CO2'
    print(f"    ρ_G    = {fmt(calc_gas_density(name,T,P), '.4f', 'kg/m³')} ")
    print(f"    μ_G    = {fmt(calc_viscosity_gas(name,T,P), '.4e', 'Pa·s')}")
    print(f"    Cp_G   = {fmt(calc_heat_capacity(name,T,P,'gas'), '.2f', 'J/mol/K')}")
    print(f"    c_G    = {fmt(calc_speed_of_sound(name,T,P,'gas'), '.1f', 'm/s')}")
    print(f"    Z      = {fmt(calc_compressibility_factor(name,T,P,'gas'), '.5f')}")
    print(f"    H      = {fmt(calc_enthalpy(name,T,P), '.1f', 'J/mol')}")
    print(f"    S      = {fmt(calc_entropy(name,T,P), '.3f', 'J/mol/K')}")
    print(f"    G      = {fmt(calc_gibbs_energy(name,T,P), '.1f', 'J/mol')}")
    print(f"    A      = {fmt(calc_helmholtz_energy(name,T,P), '.1f', 'J/mol')}")
    print(f"    U      = {fmt(calc_internal_energy(name,T,P), '.1f', 'J/mol')}")

    # 3. 扩散系数
    print("\n[扩散系数]")
    print("  气相 N2-O2 @ 298K, 1atm (ref ~1.8e-5 m²/s):")
    D_gas = calc_diffusion_gas('nitrogen', 'oxygen', 298.15, 101325)
    print(f"    D_AB (Chapman-Enskog) = {fmt(D_gas, '.3e', 'm²/s')}")
    print("  液相 O2 溶于水 @ 298K (ref ~2.4e-9 m²/s):")
    D_wc = calc_diffusion_liquid_wilke_chang('oxygen', 'water', 298.15)
    D_hm = calc_diffusion_liquid_hayduk_minhas('oxygen', 'water', 298.15)
    print(f"    D (Wilke-Chang)      = {fmt(D_wc, '.3e', 'm²/s')}")
    print(f"    D (Hayduk-Minhas, 水系数) = {fmt(D_hm, '.3e', 'm²/s')}")

    # 4. 内能 / Helmholtz 能 (新增)
    print("\n[新增热力学势: methane 200K 50bar]")
    T, P, name = 200.0, 5e6, 'methane'
    print(f"    U = {fmt(calc_internal_energy(name,T,P), '.1f', 'J/mol')}")
    print(f"    A = {fmt(calc_helmholtz_energy(name,T,P), '.1f', 'J/mol')}")

    # 5. 混合物
    print("\n[混合物: benzene(40%)+toluene(60%) 350K 1atm]")
    try:
        mp = get_mixture_properties(['benzene','toluene'],[0.4,0.6],350,101325)
        print(f"    ρ_L={mp['rho_L']}  μ_L={fmt(mp['mu_L'],'.4e')}  Pr_L={mp['Pr_L']}")
        print(f"    K={[round(k,3) for k in mp['K']]}  γ_UNIFAC={mp['activity_coefficients']}")
    except Exception as e:
        print(f"    ❌ {e}")

    # 6. 升华压
    print("\n[升华压: CO2 @ 200K (< Tt=216.6K, ref~155kPa)]")
    print(f"    Psub = {fmt(calc_sublimation_pressure('CO2',200), '.0f', 'Pa')}")

    # 7. 异常
    print("\n[异常测试: unobtanium]")
    try:
        get_fluid_density("unobtanium", 300)
    except ValueError as e:
        print(f"    ✓ 正确捕获: {e}")

    print(f"\n{'='*76}\n测试完成")