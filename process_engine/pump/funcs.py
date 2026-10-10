"""泵/压缩机：绝热、无反应、忽略动能和势能变化，指定出口压力与等熵效率。

设计原则：每个 public calc_* 接收原始参数，不互相调用，只返回自身结果。
共享私有 _solve 完成 TP -> PS -> PH；独立重复调用时可能重复闪蒸。

安装: pip install thermo==0.6.0
本文件独立使用 thermo（与 thermo_helper 的基础库相同），无需修改 helper；
不用 helper 的纯组分焓/熵加权、固定 K 闪蒸或混合焓基准。

模型必须在相连设备中保持一致：
- PR (默认)：适合许多烃类气体；水/醇等强非理想液体精度有限。
- UNIFAC：原始 UNIFAC 液相 + 理想气相，包含超额焓和 Poynting 项；
  限定 <=1.5 MPa 且各组分亚临界；液体体积仍为混合近似。
- IAPWS95：仅纯水/水蒸气，用于水泵或蒸汽压缩。
未知模型或数据缺失不静默改用另一套模型；PR 缺 kij 的零值假设会警告。

流量 mol/s，温度 K，压力绝对 Pa，摩尔焓 J/mol，功率 W。
组成必须非负、总和为 1；别名按 CAS 合并，零组分删除。
机器入口和出口限单相；泵限液体，压缩机限气相。
纯物质超临界状态在此按单相致密流体处理，设备适用性由用户核对。
不支持湿压缩、液液/固相、电解质、反应、等温或多变压缩、冷却器及 NPSH 校核。
入口须由 TP 唯一确定；饱和纯物流须增加气化率规格，此版本明确拒绝歧义。

效率定义：eta_s=(h2s-h1)/(h2-h1)，不是多变效率或电机总效率。
流体功率=n*(h2-h1)，轴输入功率=流体功率/机械效率，
电输入功率=轴输入功率/电机效率。
这里额外机械/电机损失均假设不回流到流体；若给定效率已包含这些损失，
不可重复除以它们。提供的 PR/UNIFAC/IAPWS95 必须与流程其他模块统一。
"""
from __future__ import annotations
import math
import os
import warnings
from functools import lru_cache
from thermo import (Chemical, ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX,
                    IdealGas, GibbsExcessLiquid, FlashVL, FlashPureVLS,
                    IAPWS95Gas, IAPWS95Liquid, iapws_constants, iapws_correlations)
from thermo.interaction_parameters import IPDB
from thermo.unifac import UNIFAC, UFIP, UFSG

_VERBOSE = os.getenv('FUNCS_VERBOSE', '0') == '1'
def _vprint(*args, **kwargs):
    if _VERBOSE:
        print(*args, **kwargs)




def _number(value, name, positive=False, nonnegative=False):
    if isinstance(value, (bool, str, bytes)):
        raise ValueError(f'{name} 必须是数值')
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f'{name} 必须是有限数值') from exc
    if not math.isfinite(value) or (positive and value <= 0) or (nonnegative and value < 0):
        raise ValueError(f'{name} 超出允许范围: {value}')
    return value


@lru_cache(maxsize=256)
def _cas(name):
    if not isinstance(name, str) or not name.strip():
        raise ValueError('组分名必须是非空字符串')
    try:
        return Chemical(name).CAS
    except Exception as exc:
        raise ValueError(f'无法识别组分 {name!r}') from exc


def _package(cases, model):
    model = str(model).upper()
    if model == 'IAPWS95':
        if cases != ['7732-18-5']:
            raise ValueError('IAPWS95 仅支持纯水/水蒸气')
        f = FlashPureVLS(iapws_constants, iapws_correlations,
                         gas=IAPWS95Gas(T=300., P=101325., zs=[1.]),
                         liquids=[IAPWS95Liquid(T=300., P=101325., zs=[1.])], solids=[])
        return f, iapws_constants, model, []
    if model not in ('PR', 'UNIFAC'):
        raise ValueError("thermodynamic_model 仅支持 'PR'、'UNIFAC' 或 'IAPWS95'")
    constants, props = ChemicalConstantsPackage.from_IDs(cases)
    n = len(cases)
    zs = [1.0/n]*n
    notes = ['仅模拟一个气相和至多一个液相；不检查液液分相或固相稳定性。']
    common = dict(HeatCapacityGases=props.HeatCapacityGases, T=298.15, P=101325., zs=zs)
    if model == 'PR':
        if any(v is None or not math.isfinite(v) for seq in
               (constants.Tcs, constants.Pcs, constants.omegas) for v in seq):
            raise ValueError('PR 缺少完整临界常数或偏心因子')
        kijs = IPDB.get_ip_symmetric_matrix('ChemSep PR', cases, 'kij')
        missing = [(cases[i], cases[j]) for i in range(n) for j in range(i+1, n)
                   if not IPDB.has_ip_specific('ChemSep PR', [cases[i], cases[j]], 'kij')]
        if missing:
            notes.append(f'以下组分对缺少 ChemSep PR kij，采用 0: {missing}')
        if any(c in cases for c in ('7732-18-5', '64-17-5', '67-56-1')):
            notes.append('PR 对水/醇等缔合或强非理想液相的精度有限；低压体系可显式选择 UNIFAC。')
        kwargs = dict(Tcs=constants.Tcs, Pcs=constants.Pcs,
                      omegas=constants.omegas, kijs=kijs)
        gas = CEOSGas(PRMIX, eos_kwargs=kwargs, **common)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=kwargs, **common)
    else:
        groups = constants.UNIFAC_groups
        if not all(groups):
            raise ValueError('UNIFAC 基团不完整；请明确选择适合体系的模型，未自动回退')
        # 拒绝缺失的组间参数，避免库将未知相互作用隐式当作零。
        mains = {UFSG[g].main_group_id for group in groups for g in group}
        missing = [(a, b) for a in mains for b in mains
                   if a != b and b not in UFIP.get(a, {})]
        if missing:
            raise ValueError(f'UNIFAC 缺少主基团相互作用参数: {missing}')
        ge = UNIFAC.from_subgroups(T=298.15, xs=zs, chemgroups=groups,
                                  version=0, subgroups=UFSG, interaction_data=UFIP)
        gas = IdealGas(**common)
        liquid = GibbsExcessLiquid(VaporPressures=props.VaporPressures,
                                  VolumeLiquids=props.VolumeLiquids,
                                  GibbsExcessModel=ge,
                                  equilibrium_basis='Poynting', caloric_basis='Poynting',
                                  **common)
        notes.append('原始 UNIFAC 预测模型 + 理想气相；液体体积采用纯液体体积混合规则，非实验回归保证。')
    if n == 1:
        flasher = FlashPureVLS(constants, props, gas=gas, liquids=[liquid], solids=[])
    else:
        flasher = FlashVL(constants, props, gas=gas, liquid=liquid)
    return flasher, constants, model, notes


def _domain(constants, model, T, P):
    _number(T, '温度', positive=True)
    _number(P, '绝对压力', positive=True)
    if model == 'UNIFAC':
        if P > 1.5e6 or any(tc is None or T >= tc for tc in constants.Tcs):
            raise ValueError('本 UNIFAC 路线限定 P<=1.5 MPa 且所有组分亚临界；请选其他物性模型')


def _check_state(state, zs):
    if not state.phases or len(state.phases) != len(state.betas):
        raise ValueError('闪蒸未返回有效相状态')
    betas = state.betas
    if any(not math.isfinite(b) or b < 0 or b > 1 for b in betas):
        raise ValueError('闪蒸相分率无效')
    if abs(sum(betas)-1) > 1e-8:
        raise ValueError('相分率不闭合')
    for phase in state.phases:
        if any(not math.isfinite(z) or z < 0 for z in phase.zs) or abs(sum(phase.zs)-1) > 1e-7:
            raise ValueError('相组成无效')
        _number(phase.V(), '相摩尔体积', positive=True)
    err = max(abs(sum(b*p.zs[j] for b, p in zip(betas, state.phases))-z)
              for j, z in enumerate(zs))
    if err > 1e-7:
        raise ValueError(f'闪蒸组分分配不闭合: {err}')
    if len(state.phases) == 2:
        a, b = [p.fugacities() for p in state.phases]
        ferr = max(abs(x-y)/max(abs(x), abs(y), 1e-8)
                   for x, y, z in zip(a, b, zs) if z > 1e-12)
        if ferr > 1e-5:
            raise ValueError(f'闪蒸相平衡未收敛，最大相对逸度残差 {ferr}')
    _number(state.H(), '摩尔焓')
    return err


def _flash(flasher, constants, model, zs, P, T=None, H=None, VF=None, S=None):
    if T is not None:
        _domain(constants, model, T, P)
    specs = dict(P=P, zs=zs)
    if VF is not None:
        vf = _number(VF, '摩尔气化率', nonnegative=True)
        if vf > 1:
            raise ValueError('气化率必须在 [0,1]')
        specs['VF'] = vf
    elif S is not None:
        specs['S'] = _number(S, '摩尔熵')
    elif H is not None:
        specs['H'] = _number(H, '入口摩尔焓')
    else:
        # 纯组分在饱和点 TP 不能唯一指定焓。
        if len(zs) == 1 and T is not None and T < constants.Tcs[0] and P < constants.Pcs[0]:
            sat = flasher.flash(P=P, VF=0., zs=zs)
            if abs(sat.T-T) <= 1e-4:
                raise ValueError('纯物质入口/出口位于饱和点；请指定气化率或同模型焓')
        specs['T'] = T
    state = flasher.flash(**specs)
    # thermo 0.6.0 的纯物质 PH 路径可能将同一个液体 EOS 根包装为气相。
    # 对单相且根类型与标签冲突的结果，用已求出的 TP 重新稳定判相；
    # 随后的焓校核保证重新判相没有破坏能量衡算。
    if model == 'PR' and len(state.phases) == 1:
        root_kind = getattr(state.phases[0], 'phase', None)
        reported_kind = 'g' if state.gas is not None else 'l'
        if root_kind in ('g', 'l') and root_kind != reported_kind:
            state = flasher.flash(T=state.T, P=P, zs=zs)
    _domain(constants, model, state.T, state.P)
    if T is not None and (VF is not None or H is not None) and abs(state.T-T) > max(1e-3, 1e-6*T):
        raise ValueError(f'给定 T={T} 与 P/VF 或 P/H 状态不一致，闪蒸 T={state.T}')
    _check_state(state, zs)
    if H is not None and abs(state.H()-H) > max(1e-3, 1e-7*max(abs(H), 1.)):
        raise ValueError('PH 闪蒸焓残差超限')
    if S is not None and abs(state.S()-S) > max(1e-7, 1e-8*max(abs(S), 1.)):
        raise ValueError('PS 闪蒸熵残差超限')
    return state


def _efficiency(value, name):
    eta = _number(value, name, positive=True)
    if eta > 1:
        raise ValueError(f'{name} 必须在 (0,1]，请用 0.75 而不是 75')
    return eta


def _pressures(inlet_pressure_Pa, target_pressure_Pa):
    p1 = _number(inlet_pressure_Pa, '入口绝对压力', positive=True)
    p2 = _number(target_pressure_Pa, '出口绝对压力', positive=True)
    if p2 <= p1:
        raise ValueError('此升压模型要求目标压力严格大于入口压力；降压应使用阀或膨胀机模型')
    return p1, p2


def _composition(components, mole_fractions):
    if not isinstance(components, (list, tuple)) or not components:
        raise ValueError('components 必须是非空列表')
    if not isinstance(mole_fractions, (list, tuple)) or len(components) != len(mole_fractions):
        raise ValueError('组分和摩尔分率列表必须等长')
    zs = [_number(z, '摩尔分率', nonnegative=True) for z in mole_fractions]
    total = sum(zs)
    if not math.isfinite(total) or abs(total-1) > 1e-8:
        raise ValueError('摩尔分率之和必须为 1')
    merged = {}
    for name, z in zip(components, zs):
        cas = _cas(name)
        if z > 0:
            merged[cas] = merged.get(cas, 0.) + z/total
    return list(merged), list(merged.values())


def _machine_phase(state, constants, device_type, where):
    active = [p for b, p in zip(state.betas, state.phases) if b > 1e-8]
    if len(active) != 1:
        raise ValueError(f'{where} 为两相：不支持两相泵送或湿压缩')
    # 超临界纯物质没有汽液界面；不用 VF 强行区分设备类别。
    if len(state.zs) == 1 and state.T >= constants.Tcs[0] and state.P >= constants.Pcs[0]:
        warnings.warn(f'{where} 为纯物质超临界流体，请另行确认设备适用性', RuntimeWarning, stacklevel=3)
        return
    kind = getattr(active[0], 'phase', None)
    if kind not in ('g', 'l'):
        kind = 'g' if active[0] is state.gas else 'l'
    expected = 'l' if device_type == 'pump' else 'g'
    if kind != expected:
        raise ValueError(f'{where} 相态 {kind} 不适用于 {device_type}')


def _solve(components, mole_fractions, inlet_temperature_K, inlet_pressure_Pa,
           target_pressure_Pa, isentropic_efficiency, *, device_type='pump',
           thermodynamic_model='PR'):
    if device_type not in ('pump', 'compressor'):
        raise ValueError("device_type 必须为 'pump' 或 'compressor'")
    p1, p2 = _pressures(inlet_pressure_Pa, target_pressure_Pa)
    t1 = _number(inlet_temperature_K, '入口温度', positive=True)
    eta = _efficiency(isentropic_efficiency, '等熵效率')
    cases, zs = _composition(components, mole_fractions)
    f, constants, model, notes = _package(cases, thermodynamic_model)
    for note in notes:
        if 'kij' in note or '精度有限' in note:
            warnings.warn(note, RuntimeWarning, stacklevel=3)
    inlet = _flash(f, constants, model, zs, p1, T=t1)
    _machine_phase(inlet, constants, device_type, '入口')
    h1, s1 = float(inlet.H()), float(inlet.S())
    ideal = _flash(f, constants, model, zs, p2, S=s1)
    dhs = ideal.H()-h1
    if not math.isfinite(dhs) or dhs <= 0:
        raise ValueError('升压的等熵焓升非正，检查物性模型或闪蒸收敛')
    dh = dhs/eta
    target_h = h1+dh
    outlet = _flash(f, constants, model, zs, p2, H=target_h)
    _machine_phase(outlet, constants, device_type, '出口')
    # 用焓升而非任意焓参考值确定残差容限，避免小压升被大基准值掩盖。
    herr = outlet.H()-target_h
    if abs(herr) > max(1e-6, 1e-6*abs(dh)):
        raise ValueError(f'实际出口焓残差超限: {herr} J/mol')
    serr = outlet.S()-s1
    if not math.isfinite(serr) or serr < -1e-6:
        raise ValueError('绝热压缩熵降低，检查求解或物性模型')
    _vprint(f'{device_type}: T1={t1:.6f} K, T2={outlet.T:.6f} K, model={model}')
    return dict(temperature=float(outlet.T), enthalpy_rise=float(outlet.H()-h1),
                isentropic_enthalpy_rise=float(dhs), enthalpy_residual=float(herr),
                entropy_residual=float(ideal.S()-s1), entropy_generation=float(serr),
                inlet=inlet, outlet=outlet)


def calc_outlet_pressure(inlet_pressure_Pa: float, target_pressure_Pa: float) -> dict:
    """仅校验并返回指定出口压力及压升；此函数无需温度或效率。"""
    p1, p2 = _pressures(inlet_pressure_Pa, target_pressure_Pa)
    return {'outlet_pressure_Pa': p2, 'pressure_rise_Pa': p2-p1}



def calc_outlet_temperature(
    components: list, mole_fractions: list,
    inlet_temperature_K: float, inlet_pressure_Pa: float, target_pressure_Pa: float,
    isentropic_efficiency: float, *, device_type: str = 'pump', thermodynamic_model: str = 'PR',
) -> dict:
    """出口温度 K；输入等熵效率，不能传多变效率。"""
    result = _solve(components, mole_fractions, inlet_temperature_K, inlet_pressure_Pa,
                    target_pressure_Pa, isentropic_efficiency, device_type=device_type,
                    thermodynamic_model=thermodynamic_model)
    return {'outlet_temperature_K': result['temperature']}


def calc_enthalpy_rise(
    components: list, mole_fractions: list,
    inlet_temperature_K: float, inlet_pressure_Pa: float, target_pressure_Pa: float,
    isentropic_efficiency: float, *, device_type: str = 'pump', thermodynamic_model: str = 'PR',
) -> dict:
    """实际摩尔焓升 J/mol；可作为流程能量接口。"""
    result = _solve(components, mole_fractions, inlet_temperature_K, inlet_pressure_Pa,
                    target_pressure_Pa, isentropic_efficiency, device_type=device_type,
                    thermodynamic_model=thermodynamic_model)
    return {'enthalpy_rise_J_per_mol': result['enthalpy_rise']}


def calc_fluid_power(
    components: list, mole_fractions: list,
    inlet_temperature_K: float, inlet_pressure_Pa: float, target_pressure_Pa: float,
    isentropic_efficiency: float, inlet_total_molar_flow_mol_per_s: float,
    *, device_type: str = 'pump', thermodynamic_model: str = 'PR',
) -> dict:
    """传给流体的功率 W，不含额外机械或电机损失。"""
    flow = _number(inlet_total_molar_flow_mol_per_s, '总摩尔流量', positive=True)
    result = _solve(components, mole_fractions, inlet_temperature_K, inlet_pressure_Pa,
                    target_pressure_Pa, isentropic_efficiency, device_type=device_type,
                    thermodynamic_model=thermodynamic_model)
    power = flow*result['enthalpy_rise']
    if not math.isfinite(power):
        raise ValueError('功率溢出')
    return {'fluid_power_W': power}


def calc_shaft_power(
    components: list, mole_fractions: list,
    inlet_temperature_K: float, inlet_pressure_Pa: float, target_pressure_Pa: float,
    isentropic_efficiency: float, inlet_total_molar_flow_mol_per_s: float, mechanical_efficiency: float,
    *, device_type: str = 'pump', thermodynamic_model: str = 'PR',
) -> dict:
    """轴输入功率 W；机械损失按不回到流体处理。"""
    flow = _number(inlet_total_molar_flow_mol_per_s, '总摩尔流量', positive=True)
    em = _efficiency(mechanical_efficiency, '机械效率')
    result = _solve(components, mole_fractions, inlet_temperature_K, inlet_pressure_Pa,
                    target_pressure_Pa, isentropic_efficiency, device_type=device_type,
                    thermodynamic_model=thermodynamic_model)
    power = flow*result['enthalpy_rise']/em
    if not math.isfinite(power):
        raise ValueError('功率溢出')
    return {'shaft_power_W': power}


def calc_electrical_power(
    components: list, mole_fractions: list,
    inlet_temperature_K: float, inlet_pressure_Pa: float, target_pressure_Pa: float,
    isentropic_efficiency: float, inlet_total_molar_flow_mol_per_s: float, mechanical_efficiency: float, motor_efficiency: float,
    *, device_type: str = 'pump', thermodynamic_model: str = 'PR',
) -> dict:
    """电输入功率 W；机械/电机损失均按不回到流体处理。"""
    flow = _number(inlet_total_molar_flow_mol_per_s, '总摩尔流量', positive=True)
    em = _efficiency(mechanical_efficiency, '机械效率')
    ee = _efficiency(motor_efficiency, '电机效率')
    result = _solve(components, mole_fractions, inlet_temperature_K, inlet_pressure_Pa,
                    target_pressure_Pa, isentropic_efficiency, device_type=device_type,
                    thermodynamic_model=thermodynamic_model)
    power = flow*result['enthalpy_rise']/em/ee
    if not math.isfinite(power):
        raise ValueError('功率溢出')
    return {'electrical_power_W': power}


if __name__ == '__main__':
    water = dict(components=['water'], mole_fractions=[1.], inlet_temperature_K=300.,
                 inlet_pressure_Pa=1e5, target_pressure_Pa=5e6,
                 isentropic_efficiency=.75, device_type='pump', thermodynamic_model='IAPWS95')
    print('水泵出口温度:', calc_outlet_temperature(**water))
    print('水泵流体功率:', calc_fluid_power(**water, inlet_total_molar_flow_mol_per_s=100.))
    gas = dict(components=['methane', 'ethane'], mole_fractions=[.9, .1],
               inlet_temperature_K=300., inlet_pressure_Pa=1e5, target_pressure_Pa=1e6,
               isentropic_efficiency=.75, device_type='compressor', thermodynamic_model='PR')
    print('压缩机出口温度:', calc_outlet_temperature(**gas))
    print('压缩机轴功率:', calc_shaft_power(**gas, inlet_total_molar_flow_mol_per_s=10.,
                                          mechanical_efficiency=.95))

