"""非反应、绝热、无轴功混合器；流量 mol/s、温度 K、绝对压力 Pa。

依赖: thermo>=0.6.0, scipy，以及项目 property/thermo_helper.py。
也支持本文件同目录的 thermo_helper.py。已在 thermo 0.6.0 验证。
保留原六个 public 接口；它们不互相调用，只共享私有实现。
新增 calc_mixer_state 一次获得相互一致的所有结果，推荐流程计算使用。

重要物性约定:
* 使用提供的 thermo_helper.Chemical 解析物质，直接使用其依赖 thermo
  的 FlashVL/FlashPureVLS；不调用 helper 的纯组分焓加权 H_mix 或固定 K 闪蒸。
* thermodynamic_model='PR' (默认): PR 气液状态方程及剩余焓；数据库无 kij
  时明确报告零值假设。PR 并非水/醇等强非理想液体的推荐精度模型。
* thermodynamic_model='UNIFAC': 原始 UNIFAC 液相 + 理想气相，液相焓含
  超额焓及 Poynting 校正。仅低压、亚临界、最多一液相体系；不静默回退。
* 整次计算固定同一模型和组分顺序，不混用 CoolProp/thermo 的焓基准。
* 仅允许 V、L、VL；不计算液液分相、固相、电解质反应或化学反应。
  指定模型仍需按体系验证；模型内闭合不等于 Aspen 或实验精度一致。
* 入口默认 T/P 平衡态。纯物质在饱和点的 T/P 不确定气化率，须指定
  inlet_vapor_fractions (摩尔气化率)；也可指定同模型基准的 inlet_enthalpies_J_per_mol。
  两个可选列表均按入口索引对齐，不指定的元素为 None。
* 出口压力默认最低有效入口压力；可给 outlet_pressure_Pa 或 pressure_drop_Pa。
* 零流量入口不参与热力计算及最低压力选择，负数/非有限数拒绝。
* 组分别名合并，以首次出现的名称为输出键。

示例:
    feeds = [{'Water': 5.0, 'Ethanol': 2.0}, {'Water': 3.0}]
    result = calc_mixer_state(feeds, [350., 300.], [2e5, 1.5e5],
                              thermodynamic_model='UNIFAC')
    # 体积接口沿用原签名；独立 TP 调用不能推断纯物质两相气化率，需传
    # outlet_vapor_fraction=result['vapor_fraction']，并保持同一物性模型。
"""
from __future__ import annotations

import math
import sys
from functools import lru_cache
from pathlib import Path

# 保持原项目目录结构；同目录模块仅在项目模块不存在时使用。
_parent = str(Path(__file__).resolve().parent.parent)
if _parent not in sys.path:
    sys.path.insert(0, _parent)
try:
    from property import thermo_helper as _th
except ModuleNotFoundError as exc:
    if exc.name not in ('property', 'property.thermo_helper'):
        raise
    import thermo_helper as _th

from thermo import (ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX,
                    IdealGas, GibbsExcessLiquid, FlashVL, FlashPureVLS)
from thermo.interaction_parameters import IPDB
from thermo.unifac import UNIFAC, UFIP, UFSG


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
        return _th.Chemical(name).CAS
    except Exception as exc:
        raise ValueError(f'无法识别组分 {name!r}') from exc


def _feeds(streams):
    if not isinstance(streams, (list, tuple)) or not streams:
        raise ValueError('至少需要一股进料')
    labels, clean = {}, []
    for i, stream in enumerate(streams):
        if not isinstance(stream, dict):
            raise ValueError(f'入口 {i} 必须是组分流量字典')
        row = {}
        for name, value in stream.items():
            flow = _number(value, f'入口 {i} 的 {name} 流量', nonnegative=True)
            cas = _cas(name)
            if flow == 0:
                continue
            labels.setdefault(cas, name)
            row[cas] = row.get(cas, 0.0) + flow
        clean.append(row)
    total = {cas: sum(row.get(cas, 0.0) for row in clean) for cas in labels}
    if not math.isfinite(sum(total.values())):
        raise ValueError('总流量溢出')
    return clean, labels, total


def _list(values, n, name, optional=False):
    if values is None and optional:
        return [None]*n
    if not isinstance(values, (list, tuple)) or len(values) != n:
        raise ValueError(f'{name} 必须是长度 {n} 的列表')
    return values


def _pressure(pressures, outlet, drop):
    if not pressures:
        raise ValueError('没有有效入口压力')
    minimum = min(_number(p, '入口绝对压力', positive=True) for p in pressures)
    drop = _number(drop, '压降', nonnegative=True)
    if outlet is not None and drop != 0:
        raise ValueError('出口压力和非零压降不能同时指定')
    p = _number(outlet, '出口绝对压力', positive=True) if outlet is not None else minimum-drop
    if p <= 0 or p > minimum:
        raise ValueError('无轴功混合器出口压力必须 >0 且不高于最低有效入口压力')
    return p, minimum


def _package(cases, model):
    model = str(model).upper()
    if model not in ('PR', 'UNIFAC'):
        raise ValueError("thermodynamic_model 仅支持 'PR' 或 'UNIFAC'")
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


def _flash(flasher, constants, model, zs, P, T=None, H=None, VF=None):
    if T is not None:
        _domain(constants, model, T, P)
    specs = dict(P=P, zs=zs)
    if VF is not None:
        vf = _number(VF, '摩尔气化率', nonnegative=True)
        if vf > 1:
            raise ValueError('气化率必须在 [0,1]')
        specs['VF'] = vf
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
    return state


def _state_dict(state, total, labels, model, notes):
    phase_data = []
    for beta, phase in zip(state.betas, state.phases):
        kind = 'g' if phase is state.gas else 'l'
        vm = float(phase.V())
        phase_data.append(dict(phase=kind, mole_fraction=float(beta),
                               molar_flow_mol_per_s=total*beta,
                               composition=dict(zip(labels, map(float, phase.zs))),
                               molar_volume_m3_per_mol=vm,
                               volumetric_flow_m3_per_s=total*beta*vm))
    phase_name = 'two_phase' if len(phase_data) > 1 else phase_data[0]['phase']
    return dict(outlet_temperature_K=float(state.T), outlet_pressure_Pa=float(state.P),
                outlet_molar_enthalpy_J_per_mol=float(state.H()),
                vapor_fraction=float(state.VF), phase=phase_name,
                outlet_volumetric_flow_m3_per_s=sum(p['volumetric_flow_m3_per_s'] for p in phase_data),
                phases=phase_data, thermodynamic_model=model, model_notes=notes)


def _solve(streams, temperatures, pressures, *, thermodynamic_model='PR',
           outlet_pressure_Pa=None, pressure_drop_Pa=0.,
           inlet_vapor_fractions=None, inlet_enthalpies_J_per_mol=None):
    rows, labels, totals = _feeds(streams)
    n = len(rows)
    ts = _list(temperatures, n, '入口温度')
    ps = _list(pressures, n, '入口压力')
    vfs = _list(inlet_vapor_fractions, n, '入口气化率', True)
    hs = _list(inlet_enthalpies_J_per_mol, n, '入口焓', True)
    active = [i for i, row in enumerate(rows) if sum(row.values()) > 0]
    if not active:
        raise ValueError('所有入口流量为零，无法定义出口热力状态')
    p_out, p_min = _pressure([ps[i] for i in active], outlet_pressure_Pa, pressure_drop_Pa)
    cases = list(totals)
    flasher, constants, model, notes = _package(cases, thermodynamic_model)
    hin, inlet_data = 0., []
    for i in active:
        T = _number(ts[i], f'入口 {i} 温度', positive=True)
        P = _number(ps[i], f'入口 {i} 压力', positive=True)
        if vfs[i] is not None and hs[i] is not None:
            raise ValueError('同一入口不能同时指定气化率和焓')
        flow = sum(rows[i].values())
        zs = [rows[i].get(c, 0.)/flow for c in cases]
        # union 组分中只有一个非零组分时，同样检查纯物流饱和歧义。
        present = [j for j, z in enumerate(zs) if z > 0]
        if len(cases) > 1 and len(present) == 1 and vfs[i] is None and hs[i] is None:
            j = present[0]
            if T < constants.Tcs[j] and P < constants.Pcs[j]:
                sat = flasher.flash(P=P, VF=0., zs=zs)
                if abs(sat.T-T) <= 1e-4:
                    raise ValueError(f'入口 {i} 是饱和纯物流，必须指定气化率或焓')
        st = _flash(flasher, constants, model, zs, P, T=T, H=hs[i], VF=vfs[i])
        hin += flow*st.H()
        inlet_data.append(dict(index=i, molar_enthalpy_J_per_mol=float(st.H()),
                               vapor_fraction=float(st.VF)))
    flow_out = sum(totals.values())
    zs = [totals[c]/flow_out for c in cases]
    out = _flash(flasher, constants, model, zs, p_out, H=hin/flow_out)
    result = _state_dict(out, flow_out, list(labels.values()), model, notes)
    hout = flow_out*out.H()
    scale = max(abs(hin), abs(hout), flow_out, 1.)
    if abs(hout-hin) > max(1e-3*flow_out, 1e-7*scale):
        raise ValueError(f'总能量衡算不闭合: {hout-hin} W')
    result.update(outlet_molar_flows_mol_per_s={labels[c]: v for c, v in totals.items()},
                  outlet_total_molar_flow_mol_per_s=flow_out,
                  outlet_mole_fractions=dict(zip(labels.values(), zs)),
                  component_CAS=dict(zip(labels.values(), cases)),
                  min_inlet_pressure_Pa=p_min, inlet_states=inlet_data,
                  inlet_enthalpy_flow_W=hin, outlet_enthalpy_flow_W=hout,
                  energy_residual_W=hout-hin,
                  energy_residual_relative=(hout-hin)/scale,
                  phase_component_balance_residual=_check_state(out, zs))
    return result


def calc_outlet_molar_flows(inlet_molar_flows_list: list) -> dict:
    """出口逐组分摩尔流量；别名合并，零组分省略。"""
    _, labels, totals = _feeds(inlet_molar_flows_list)
    return {'outlet_molar_flows_mol_per_s': {labels[c]: v for c, v in totals.items()}}


def calc_outlet_total_molar_flow(inlet_molar_flows_list: list) -> dict:
    """出口总摩尔流量；所有入口为零时返回零。"""
    _, _, totals = _feeds(inlet_molar_flows_list)
    return {'outlet_total_molar_flow_mol_per_s': sum(totals.values())}


def calc_outlet_composition(inlet_molar_flows_list: list) -> dict:
    """出口总组成，不是两相中任一相的组成。"""
    _, labels, totals = _feeds(inlet_molar_flows_list)
    total = sum(totals.values())
    if total <= 0:
        raise ValueError('出口总流量为零，组成未定义')
    return {'outlet_mole_fractions': {labels[c]: v/total for c, v in totals.items()}}


def calc_outlet_pressure(inlet_pressures_Pa: list, *, outlet_pressure_Pa=None,
                         pressure_drop_Pa=0., inlet_molar_flows_list=None) -> dict:
    """默认最低入口压力；可传流量列表以排除零流量入口。"""
    if not isinstance(inlet_pressures_Pa, (list, tuple)):
        raise ValueError('入口压力必须为列表')
    pressures = inlet_pressures_Pa
    if inlet_molar_flows_list is not None:
        rows, _, _ = _feeds(inlet_molar_flows_list)
        _list(pressures, len(rows), '入口压力')
        pressures = [p for p, row in zip(pressures, rows) if sum(row.values()) > 0]
    p, pmin = _pressure(pressures, outlet_pressure_Pa, pressure_drop_Pa)
    return {'outlet_pressure_Pa': p, 'min_inlet_pressure_Pa': pmin}


def calc_outlet_temperature(inlet_molar_flows_list: list, inlet_temperatures_K: list,
                            inlet_pressures_Pa: list, **options) -> dict:
    """PH 闪蒸出口温度；options 同 calc_mixer_state，额外返回气化率和残差。"""
    state = _solve(inlet_molar_flows_list, inlet_temperatures_K, inlet_pressures_Pa, **options)
    keys = ('outlet_temperature_K', 'outlet_pressure_Pa', 'vapor_fraction', 'phase',
            'thermodynamic_model', 'model_notes', 'energy_residual_W', 'energy_residual_relative')
    return {key: state[key] for key in keys}


def calc_outlet_volumetric_flow(inlet_molar_flows_list: list,
                                outlet_temperature_K: float, outlet_pressure_Pa: float,
                                *, thermodynamic_model='PR', outlet_vapor_fraction=None) -> dict:
    """给定出口 TP 的体积计算；两相按各相流量及摩尔体积求和。

    本接口不校验入口能量衡算。必须使用温度计算时同一物性模型；纯物质
    两相需另给 outlet_vapor_fraction。流程计算优先使用 calc_mixer_state。
    """
    _, labels, totals = _feeds(inlet_molar_flows_list)
    flow = sum(totals.values())
    if flow <= 0:
        raise ValueError('总流量为零，出口状态未定义')
    zs = [v/flow for v in totals.values()]
    flasher, constants, model, notes = _package(list(totals), thermodynamic_model)
    st = _flash(flasher, constants, model, zs,
                _number(outlet_pressure_Pa, '出口压力', positive=True),
                T=_number(outlet_temperature_K, '出口温度', positive=True), VF=outlet_vapor_fraction)
    return _state_dict(st, flow, list(labels.values()), model, notes)


def calc_mixer_state(inlet_molar_flows_list: list, inlet_temperatures_K: list,
                     inlet_pressures_Pa: list, *, thermodynamic_model='PR',
                     outlet_pressure_Pa=None, pressure_drop_Pa=0.,
                     inlet_vapor_fractions=None, inlet_enthalpies_J_per_mol=None) -> dict:
    """一次返回物料、温压、相分率、各相组成/体积、焓及闭合残差。

    inlet_enthalpies_J_per_mol 只能使用本模型同基准的焓，不能直接传入
    Aspen/CoolProp 的绝对焓；提供的入口温度仍会与 PH/VF 状态核对。
    所有 optional 列表必须与入口等长。公开函数之间无互相调用。
    """
    return _solve(inlet_molar_flows_list, inlet_temperatures_K, inlet_pressures_Pa,
                  thermodynamic_model=thermodynamic_model, outlet_pressure_Pa=outlet_pressure_Pa,
                  pressure_drop_Pa=pressure_drop_Pa, inlet_vapor_fractions=inlet_vapor_fractions,
                  inlet_enthalpies_J_per_mol=inlet_enthalpies_J_per_mol)


if __name__ == '__main__':
    from pprint import pprint
    pprint(calc_mixer_state([{'Water': 5., 'Ethanol': 2.}, {'Water': 3.}],
                            [350., 300.], [2e5, 1.5e5], thermodynamic_model='UNIFAC'))
