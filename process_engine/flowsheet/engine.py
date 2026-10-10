"""回流迭代引擎 — 通用带循环化工流程求解器

三层架构:
  用户层  — flowsheet 配置字典（描述拓扑 + 标记循环流）
  核心层  — 本模块：Tear Stream 迭代 + Wegstein 加速 + 统一收敛判据
  设备层  — 适配 process_engine 现有 run_* 工具

用法:
    from process_engine.flowsheet.engine import FlowsheetEngine
    eng = FlowsheetEngine(config)
    result = eng.solve()
"""
from __future__ import annotations
import math
import copy
from typing import Dict, List, Any, Optional, Tuple


# ════════════════════════════════════════════════════════════════
#  设备适配器 — 将 flowsheet 流股模型翻译为 run_* 函数调用
# ════════════════════════════════════════════════════════════════

def _adapt_mixer(inlets: List[Dict], config: Dict) -> Dict:
    """混合器: N 股进料 → 1 股出料"""
    from multi_equipment_series.funcs import run_mixer
    flows_list = [_component_flows(s) for s in inlets]
    temps = [s.get("_T", 298.15) for s in inlets]
    pressures = [s.get("_P", 101325.0) for s in inlets]
    result = run_mixer(
        inlet_molar_flows_list=flows_list,
        inlet_temperatures_K=temps,
        inlet_pressures_Pa=pressures,
    )
    out = dict(result["outlet_molar_flows_mol_per_s"])
    out["_T"] = result["outlet_temperature_K"]
    out["_P"] = result["outlet_pressure_Pa"]
    return out


def _adapt_pump(inlets: List[Dict], config: Dict) -> Dict:
    """泵: 只改变压力，流量/温度不变"""
    from multi_equipment_series.funcs import run_pump
    inlet = inlets[0]
    p_in = inlet.get("_P", 101325.0)
    p_out = config["params"]["target_pressure_Pa"]
    result = run_pump(inlet_pressure_Pa=p_in, target_pressure_Pa=p_out)
    out = _copy_flows(inlet)
    out["_P"] = result["outlet_pressure_Pa"]
    out["_T"] = inlet.get("_T", 298.15)
    return out


def _adapt_valve(inlets: List[Dict], config: Dict) -> Dict:
    """减压阀: 只改变压力，流量/温度不变"""
    from multi_equipment_series.funcs import run_pressure_reducing_valve
    inlet = inlets[0]
    p_in = inlet.get("_P", 101325.0)
    p_out = config["params"]["target_pressure_Pa"]
    result = run_pressure_reducing_valve(
        inlet_pressure_Pa=p_in, target_pressure_Pa=p_out)
    out = _copy_flows(inlet)
    out["_P"] = result["outlet_pressure_Pa"]
    out["_T"] = inlet.get("_T", 298.15)
    return out


def _adapt_heater(inlets: List[Dict], config: Dict) -> Dict:
    """换热器（简化）: 设定出口温度，压力不变"""
    inlet = inlets[0]
    params = config["params"]
    out = _copy_flows(inlet)
    out["_T"] = params["outlet_temp_K"]
    out["_P"] = params.get("outlet_pressure_Pa", inlet.get("_P", 101325.0))
    return out


def _adapt_reactor(inlets: List[Dict], config: Dict) -> Dict:
    """化学计量反应器: 按化学计量比改变流量"""
    from multi_equipment_series.funcs import run_stoichiometric_reactor
    inlet = inlets[0]
    params = config["params"]
    feed = _copy_flows(inlet)
    # 确保产物组分存在（初始为 0）
    for comp in params.get("all_components", []):
        if comp not in feed:
            feed[comp] = 0.0
    result = run_stoichiometric_reactor(
        inlet_molar_flows_mol_per_s=feed,
        stoichiometric_coefficients=params["stoichiometric_coefficients"],
        key_component=params["key_component"],
        key_component_conversion=params["key_component_conversion"],
        reactor_temp_K=inlet.get("_T", 298.15),
        reactor_pressure_Pa=inlet.get("_P", 101325.0),
    )
    out = dict(result["outlet_molar_flows_mol_per_s"])
    out["_T"] = inlet.get("_T", 298.15)
    out["_P"] = inlet.get("_P", 101325.0)
    return out


def _adapt_distillation(inlets: List[Dict], config: Dict) -> Tuple[Dict, Dict]:
    """精馏塔（简化质量衡算）: 1 股进料 → 2 股出料（塔顶 + 塔底）"""
    from tray_distillation.funcs import calc_mass_balance
    inlet = inlets[0]
    params = config["params"]
    feed = _copy_flows(inlet)
    result = calc_mass_balance(
        feed_molar_flows=feed,
        distillate_purity=params["distillate_purity"],
        bottoms_purity=params["bottoms_purity"],
        light_key_component=params["light_key_component"],
        heavy_key_component=params["heavy_key_component"],
        light_components=params["light_components"],
        heavy_components=params["heavy_components"],
    )
    distillate = dict(result["distillate_flows_mol_per_s"])
    bottoms = dict(result["bottoms_flows_mol_per_s"])
    # 温度/压力继承进料（简化模型，后处理时再精确计算）
    for d in (distillate, bottoms):
        d["_T"] = inlet.get("_T", 298.15)
        d["_P"] = inlet.get("_P", 101325.0)
    return distillate, bottoms


def _adapt_column_recovery(inlets: List[Dict], config: Dict) -> Tuple[Dict, Dict]:
    """回收率 + wt% 规格分离塔: 1 股进料 → 2 股出料（塔顶 + 塔底）。

    params: recover_to_top, top_mass_fraction, bottom_mass_fraction,
            swing_component, MW(可选)。温/ pres 继承进料（简化质衡算）。
    """
    from multi_equipment_series.funcs import run_column_recovery
    inlet = inlets[0]
    params = config["params"]
    feed = _copy_flows(inlet)
    result = run_column_recovery(
        feed_molar_flows_mol_per_s=feed,
        recover_to_top=params["recover_to_top"],
        top_mass_fraction=params.get("top_mass_fraction"),
        bottom_mass_fraction=params.get("bottom_mass_fraction"),
        swing_component=params.get("swing_component"),
        MW=params.get("MW"),
    )
    top = dict(result["top_flows_mol_per_s"])
    bottom = dict(result["bottom_flows_mol_per_s"])
    for d in (top, bottom):
        d["_T"] = inlet.get("_T", 298.15)
        d["_P"] = inlet.get("_P", 101325.0)
    return top, bottom


def _adapt_flash(inlets: List[Dict], config: Dict) -> Tuple[Dict, Dict]:
    """闪蒸罐: 给定 T, P → 气液两相"""
    from multi_equipment_series.funcs import run_flash_drum
    inlet = inlets[0]
    params = config["params"]
    feed = _copy_flows(inlet)
    result = run_flash_drum(
        flash_pressure_Pa=params.get("flash_pressure_Pa", inlet.get("_P", 101325.0)),
        inlet_molar_flows_mol_per_s=feed,
        key_component=params["key_component"],
        key_component_recovery=params["key_component_recovery"],
        recovery_phase=params.get("recovery_phase", "vapor"),
    )
    vapor = dict(result["vapor_molar_flows_mol_per_s"])
    liquid = dict(result["liquid_molar_flows_mol_per_s"])
    T_flash = result["flash_temperature_K"]
    P_flash = params.get("flash_pressure_Pa", inlet.get("_P", 101325.0))
    for d in (vapor, liquid):
        d["_T"] = T_flash
        d["_P"] = P_flash
    return vapor, liquid


# ── 适配器注册表 ──
_EQUIPMENT_ADAPTERS = {
    "mixer":         {"fn": _adapt_mixer,         "outlets": 1},
    "pump":          {"fn": _adapt_pump,          "outlets": 1},
    "compressor":    {"fn": _adapt_pump,          "outlets": 1},
    "valve":         {"fn": _adapt_valve,         "outlets": 1},
    "heater":        {"fn": _adapt_heater,        "outlets": 1},
    "cooler":        {"fn": _adapt_heater,        "outlets": 1},
    "reactor":       {"fn": _adapt_reactor,       "outlets": 1},
    "stoichiometric_reactor": {"fn": _adapt_reactor, "outlets": 1},
    "distillation":  {"fn": _adapt_distillation,  "outlets": 2},
    "column_recovery": {"fn": _adapt_column_recovery, "outlets": 2},
    "flash":         {"fn": _adapt_flash,         "outlets": 2},
}


# ── 辅助函数 ──
def _component_flows(stream: Dict) -> Dict[str, float]:
    """从流股中提取纯组分流量（去掉 _T/_P 元数据）"""
    return {k: v for k, v in stream.items() if not k.startswith("_")}


def _copy_flows(stream: Dict) -> Dict[str, float]:
    """复制组分流量"""
    return {k: v for k, v in stream.items() if not k.startswith("_")}


# ════════════════════════════════════════════════════════════════
#  FlowsheetEngine — 核心引擎
# ════════════════════════════════════════════════════════════════

class FlowsheetEngine:
    """带循环收敛的化工流程求解引擎

    Parameters
    ----------
    config : dict
        流程配置字典，包含:
        - equipment : List[dict]  — 设备列表（按执行顺序）
        - recycles  : List[dict]  — 循环流定义 [{from, to, tear}]
        - feeds     : dict        — 新鲜进料 {stream_name: {component: flow}}
        - defaults  : dict(可选)  — 流股默认温度/压力

    Example
    -------
    >>> config = {
    ...     "feeds": {"FRESH": {"A": 100}},
    ...     "equipment": [
    ...         {"id": "MX", "type": "mixer",
    ...          "inlet": ["FRESH", "RECYCLE"], "outlet": "S1"},
    ...         {"id": "R",  "type": "reactor",
    ...          "inlet": "S1", "outlet": "S2",
    ...          "params": {...}},
    ...         {"id": "T",  "type": "distillation",
    ...          "inlet": "S2", "outlet": ["D", "W"],
    ...          "params": {...}},
    ...     ],
    ...     "recycles": [{"from": "D", "to": "RECYCLE"}],
    ... }
    >>> eng = FlowsheetEngine(config)
    >>> result = eng.solve()
    """

    def __init__(self, config: Dict):
        self.config = config
        self.feeds: Dict[str, Dict] = {}
        self.equipment: List[Dict] = []
        self.recycles: List[Dict] = []
        self.streams: Dict[str, Dict] = {}
        self._equipment_log: List[Dict] = []
        self._parse_config()

    # ────────────────────────────────────────────
    #  配置解析
    # ────────────────────────────────────────────
    def _parse_config(self):
        """解析配置字典"""
        cfg = self.config

        # 默认 T/P
        defaults = cfg.get("defaults", {})
        self._default_T = defaults.get("temperature_K", 298.15)
        self._default_P = defaults.get("pressure_Pa", 101325.0)

        # 新鲜进料
        for name, flows in cfg.get("feeds", {}).items():
            s = dict(flows)
            s.setdefault("_T", self._default_T)
            s.setdefault("_P", self._default_P)
            self.feeds[name] = s

        # 设备列表
        self.equipment = list(cfg.get("equipment", []))

        # 循环定义
        self.recycles = list(cfg.get("recycles", []))

    # ────────────────────────────────────────────
    #  Tear Stream 管理
    # ────────────────────────────────────────────
    def _get_tear_info(self) -> List[Dict]:
        """提取所有循环流的信息"""
        info = []
        for r in self.recycles:
            info.append({
                "from_stream": r["from"],
                "to_stream": r["to"],
                "initial": r.get("initial", None),
            })
        return info

    def _init_tear_streams(self, tear_info: List[Dict]) -> Dict[str, Dict]:
        """初始化循环流（全零或用户猜测值）"""
        tears = {}
        for t in tear_info:
            if t["initial"] is not None:
                s = dict(t["initial"])
                s.setdefault("_T", self._default_T)
                s.setdefault("_P", self._default_P)
                tears[t["to_stream"]] = s
            else:
                tears[t["to_stream"]] = {"_T": self._default_T, "_P": self._default_P}
        return tears

    def _extract_tear_values(self, streams: Dict, tear_info: List[Dict]) -> Dict[str, Dict]:
        """从流股字典中提取循环流的纯组分流量"""
        values = {}
        for t in tear_info:
            src = t["from_stream"]
            if src in streams:
                values[src] = _component_flows(streams[src])
            else:
                values[src] = {}
        return values

    def _tear_residual(self, new: Dict, old: Dict, tear_info: List[Dict]) -> float:
        """计算循环流总残差 sum(|new_i - old_i|)"""
        total = 0.0
        for t in tear_info:
            src = t["from_stream"]
            n = new.get(src, {})
            o = old.get(src, {})
            all_keys = set(list(n.keys()) + list(o.keys()))
            for k in all_keys:
                total += abs(n.get(k, 0.0) - o.get(k, 0.0))
        return total

    # ────────────────────────────────────────────
    #  流程模拟（单次遍历）
    # ────────────────────────────────────────────
    def _simulate_once(self, tear_streams: Dict) -> Tuple[Dict, Dict, List]:
        """按设备顺序执行一次完整遍历

        Returns: (outlet_streams, tear_new_values, equipment_log)
        """
        streams: Dict[str, Dict] = {}

        # ① 注入新鲜进料
        for name, flow in self.feeds.items():
            streams[name] = dict(flow)

        # ② 注入循环流
        for name, flow in tear_streams.items():
            streams[name] = dict(flow)

        log = []

        # ③ 按序执行各设备
        for eq in self.equipment:
            eq_id = eq["id"]
            eq_type = eq["type"].lower().replace(" ", "_")
            adapter_info = _EQUIPMENT_ADAPTERS.get(eq_type)
            if adapter_info is None:
                raise ValueError(f"未知设备类型: {eq['type']}")

            adapter_fn = adapter_info["fn"]
            n_outlets = adapter_info["outlets"]

            # 收集入口流股
            inlet_names = eq["inlet"]
            if isinstance(inlet_names, str):
                inlet_names = [inlet_names]
            inlets = []
            for name in inlet_names:
                if name not in streams:
                    raise ValueError(
                        f"设备 {eq_id}: 入口流股 '{name}' 未定义。"
                        f"可用流股: {list(streams.keys())}")
                inlets.append(streams[name])

            # 调用适配器
            if n_outlets == 1:
                outlet = adapter_fn(inlets, eq)
                outlet_name = eq["outlet"]
                if isinstance(outlet_name, list):
                    outlet_name = outlet_name[0]
                streams[outlet_name] = outlet
                log.append({
                    "id": eq_id, "type": eq["type"],
                    "inlet_streams": inlet_names,
                    "outlet_streams": [outlet_name],
                    "outlet": outlet,
                })
            else:
                out1, out2 = adapter_fn(inlets, eq)
                out_names = eq["outlet"]
                if not isinstance(out_names, list) or len(out_names) != 2:
                    raise ValueError(
                        f"设备 {eq_id}: 双出口设备 outlet 必须是长度为 2 的列表")
                streams[out_names[0]] = out1
                streams[out_names[1]] = out2
                log.append({
                    "id": eq_id, "type": eq["type"],
                    "inlet_streams": inlet_names,
                    "outlet_streams": out_names,
                    "outlet": {out_names[0]: out1, out_names[1]: out2},
                })

        # ④ 提取新的循环流值
        tear_info = self._get_tear_info()
        tear_new = self._extract_tear_values(streams, tear_info)

        return streams, tear_new, log

    # ────────────────────────────────────────────
    #  Wegstein 加速
    # ────────────────────────────────────────────
    @staticmethod
    def _wegstein_update(
        x_old: Dict[str, float],
        x_new: Dict[str, float],
        x_prev: Optional[Dict[str, float]],
        f_prev: Optional[Dict[str, float]],
        omega: float = 0.5,
    ) -> Tuple[Dict[str, float], Dict[str, float]]:
        """Wegstein 割线法加速循环收敛

        Returns: (x_next, f_curr)  其中 f_curr = x_new - x_old
        """
        all_keys = set(list(x_old.keys()) + list(x_new.keys()))
        f_curr = {}
        x_next = {}

        for k in all_keys:
            old_v = x_old.get(k, 0.0)
            new_v = x_new.get(k, 0.0)
            f_curr[k] = new_v - old_v

            if (x_prev is not None and f_prev is not None
                    and k in x_prev and k in f_prev):
                # Wegstein 割线外推
                delta_f = f_curr[k] - f_prev[k]
                if abs(delta_f) > 1e-14:
                    w = f_curr[k] / delta_f
                    # 限制加速因子，防止发散
                    w = max(-3.0, min(3.0, w))
                    next_v = new_v - w * f_curr[k]
                else:
                    next_v = old_v + omega * f_curr[k]
            else:
                # 首次迭代：阻尼逐次代入
                next_v = old_v + omega * f_curr[k]

            # 物理约束：流量不能为负
            x_next[k] = max(0.0, next_v)

        return x_next, f_curr

    # ────────────────────────────────────────────
    #  主求解方法
    # ────────────────────────────────────────────
    def solve(
        self,
        max_iter: int = 500,
        tol: float = 1e-4,
        acceleration: str = "wegstein",
        omega: float = 0.5,
        verbose: bool = False,
    ) -> Dict[str, Any]:
        """求解带循环的化工流程

        Parameters
        ----------
        max_iter : int
            最大迭代次数
        tol : float
            收敛容差（循环流总残差）
        acceleration : str
            加速方法: "wegstein"（默认）/ "successive"（逐次代入）/ "none"
        omega : float
            阻尼因子（逐次代入法用，默认 0.5）
        verbose : bool
            是否打印迭代过程

        Returns
        -------
        dict with keys:
            - converged : bool
            - iterations : int
            - residual : float
            - streams : Dict[str, Dict]  — 所有流股
            - equipment_log : List[Dict] — 设备执行记录
            - feeds : Dict — 新鲜进料
        """
        tear_info = self._get_tear_info()
        has_recycle = len(tear_info) > 0

        if not has_recycle:
            # 无循环 → 单次模拟
            streams, _, log = self._simulate_once({})
            return {
                "converged": True,
                "iterations": 0,
                "residual": 0.0,
                "streams": streams,
                "equipment_log": log,
                "feeds": self.feeds,
            }

        # ── 初始化循环流 ──
        tear_streams = self._init_tear_streams(tear_info)
        tear_old = self._extract_tear_values(tear_streams, tear_info)

        x_prev: Optional[Dict] = None
        f_prev: Optional[Dict] = None
        converged = False
        iterations = 0
        final_residual = float("inf")

        if verbose:
            print("=" * 60)
            print("回流迭代引擎 — 开始收敛")
            print(f"  循环流: {[t['to_stream'] for t in tear_info]}")
            print(f"  加速方法: {acceleration}, ω={omega}")
            print(f"  容差: {tol}, 最大迭代: {max_iter}")
            print("=" * 60)

        for it in range(max_iter):
            iterations = it + 1

            # ① 用当前循环流模拟一次
            streams, tear_new, log = self._simulate_once(tear_streams)

            # ② 计算残差
            residual = self._tear_residual(tear_new, tear_old, tear_info)
            final_residual = residual

            if verbose:
                detail = ""
                for t in tear_info:
                    src = t["from_stream"]
                    vals = tear_new.get(src, {})
                    total = sum(v for v in vals.values())
                    detail += f"  {src}: Σ={total:.2f}"
                print(f"  iter {iterations:3d} | residual={residual:.2e} | {detail}")

            # ③ 收敛判断
            if residual < tol:
                converged = True
                break

            # ④ 更新循环流
            if acceleration == "wegstein":
                # 合并所有循环流为单字典进行 Wegstein
                flat_old = {}
                flat_new = {}
                for t in tear_info:
                    src = t["from_stream"]
                    tgt = t["to_stream"]
                    for k, v in tear_old.get(src, {}).items():
                        flat_old[f"{tgt}::{k}"] = v
                    for k, v in tear_new.get(src, {}).items():
                        flat_new[f"{tgt}::{k}"] = v

                flat_next, _ = self._wegstein_update(
                    flat_old, flat_new, x_prev, f_prev, omega)

                # 拆分回各循环流
                tear_streams = {}
                for t in tear_info:
                    tgt = t["to_stream"]
                    s = {}
                    for key, val in flat_next.items():
                        if key.startswith(f"{tgt}::"):
                            comp = key.split("::", 1)[1]
                            s[comp] = val
                    s.setdefault("_T", self._default_T)
                    s.setdefault("_P", self._default_P)
                    tear_streams[tgt] = s

                x_prev = dict(flat_old)
                f_prev = {k: flat_new.get(k, 0) - flat_old.get(k, 0)
                          for k in set(list(flat_old.keys()) + list(flat_new.keys()))}

            elif acceleration == "successive":
                # 阻尼逐次代入
                for t in tear_info:
                    src = t["from_stream"]
                    tgt = t["to_stream"]
                    new_s = {}
                    for k, v in tear_new.get(src, {}).items():
                        old_v = tear_old.get(src, {}).get(k, 0.0)
                        new_s[k] = max(0.0, old_v + omega * (v - old_v))
                    new_s.setdefault("_T", self._default_T)
                    new_s.setdefault("_P", self._default_P)
                    tear_streams[tgt] = new_s

            else:
                # 无加速（纯逐次代入）
                for t in tear_info:
                    src = t["from_stream"]
                    tgt = t["to_stream"]
                    new_s = dict(tear_new.get(src, {}))
                    new_s.setdefault("_T", self._default_T)
                    new_s.setdefault("_P", self._default_P)
                    tear_streams[tgt] = new_s

            tear_old = tear_new

        # ── 最终模拟（确保所有流股一致） ──
        streams, _, log = self._simulate_once(tear_streams)

        return {
            "converged": converged,
            "iterations": iterations,
            "residual": final_residual,
            "streams": streams,
            "equipment_log": log,
            "feeds": self.feeds,
        }

    # ────────────────────────────────────────────
    #  后处理 — 用收敛后的流股调用完整设备计算
    # ────────────────────────────────────────────
    def post_process(self, streams: Dict, equipment_log: List[Dict]) -> Dict[str, Any]:
        """用收敛后的流股调用完整设备函数，获取设计参数

        对每种设备执行完整计算（换热器面积、精馏塔板数等），
        而不仅仅是质量衡算。

        Returns
        -------
        dict : {equipment_id: full_result_dict}
        """
        results = {}

        for eq in self.equipment:
            eq_id = eq["id"]
            eq_type = eq["type"].lower().replace(" ", "_")
            params = eq.get("params", {})

            # 获取入口流股
            inlet_names = eq["inlet"]
            if isinstance(inlet_names, str):
                inlet_names = [inlet_names]
            inlets = {n: streams[n] for n in inlet_names if n in streams}

            try:
                if eq_type in ("heater", "cooler"):
                    results[eq_id] = self._post_heater(inlets, eq)
                elif eq_type in ("reactor", "stoichiometric_reactor"):
                    results[eq_id] = self._post_reactor(inlets, eq)
                elif eq_type == "distillation":
                    results[eq_id] = self._post_distillation(inlets, eq)
                elif eq_type in ("pump", "compressor"):
                    results[eq_id] = self._post_pump(inlets, eq)
                elif eq_type == "valve":
                    results[eq_id] = self._post_valve(inlets, eq)
                elif eq_type == "mixer":
                    results[eq_id] = self._post_mixer(inlets, eq)
                else:
                    results[eq_id] = {"note": f"后处理不支持设备类型 {eq_type}"}
            except Exception as e:
                results[eq_id] = {"error": str(e)}

        return results

    def _post_heater(self, inlets: Dict, eq: Dict) -> Dict:
        from multi_equipment_series.funcs import run_heat_exchanger
        inlet_name = eq["inlet"] if isinstance(eq["inlet"], str) else eq["inlet"][0]
        s = inlets[inlet_name]
        flows = _component_flows(s)
        params = eq["params"]
        utility = params.get("utility", {})
        return run_heat_exchanger(
            process_fluid_molar_flows_mol_per_s=flows,
            process_fluid_temp_in_K=params.get("inlet_temp_K", s.get("_T", 298.15)),
            process_fluid_temp_out_K=params["outlet_temp_K"],
            process_fluid_pressure_Pa=params.get("outlet_pressure_Pa", s.get("_P", 101325.0)),
            utility_fluid_temp_in_K=utility.get("temp_in_K", 298.15),
            utility_fluid_temp_out_K=utility.get("temp_out_K", 323.15),
            utility_type=utility.get("type", "cooling_water"),
        )

    def _post_reactor(self, inlets: Dict, eq: Dict) -> Dict:
        from multi_equipment_series.funcs import run_stoichiometric_reactor
        inlet_name = eq["inlet"] if isinstance(eq["inlet"], str) else eq["inlet"][0]
        s = inlets[inlet_name]
        feed = _component_flows(s)
        params = eq["params"]
        for comp in params.get("all_components", []):
            feed.setdefault(comp, 0.0)
        return run_stoichiometric_reactor(
            inlet_molar_flows_mol_per_s=feed,
            stoichiometric_coefficients=params["stoichiometric_coefficients"],
            key_component=params["key_component"],
            key_component_conversion=params["key_component_conversion"],
            reactor_temp_K=s.get("_T", 298.15),
            reactor_pressure_Pa=s.get("_P", 101325.0),
        )

    def _post_distillation(self, inlets: Dict, eq: Dict) -> Dict:
        from multi_equipment_series.funcs import run_distillation_column
        inlet_name = eq["inlet"] if isinstance(eq["inlet"], str) else eq["inlet"][0]
        s = inlets[inlet_name]
        feed = _component_flows(s)
        params = eq["params"]
        return run_distillation_column(
            feed_molar_flows_mol_per_s=feed,
            distillate_purity=params["distillate_purity"],
            bottoms_purity=params["bottoms_purity"],
            light_key_component=params["light_key_component"],
            heavy_key_component=params["heavy_key_component"],
            light_components=params["light_components"],
            heavy_components=params["heavy_components"],
            feed_temp_K=s.get("_T", 298.15),
            feed_pressure_Pa=s.get("_P", 101325.0),
            reflux_factor=params.get("reflux_factor", 1.2),
        )

    def _post_pump(self, inlets: Dict, eq: Dict) -> Dict:
        from multi_equipment_series.funcs import run_pump
        inlet_name = eq["inlet"] if isinstance(eq["inlet"], str) else eq["inlet"][0]
        s = inlets[inlet_name]
        return run_pump(
            inlet_pressure_Pa=s.get("_P", 101325.0),
            target_pressure_Pa=eq["params"]["target_pressure_Pa"],
        )

    def _post_valve(self, inlets: Dict, eq: Dict) -> Dict:
        from multi_equipment_series.funcs import run_pressure_reducing_valve
        inlet_name = eq["inlet"] if isinstance(eq["inlet"], str) else eq["inlet"][0]
        s = inlets[inlet_name]
        return run_pressure_reducing_valve(
            inlet_pressure_Pa=s.get("_P", 101325.0),
            target_pressure_Pa=eq["params"]["target_pressure_Pa"],
        )

    def _post_mixer(self, inlets: Dict, eq: Dict) -> Dict:
        from multi_equipment_series.funcs import run_mixer
        inlet_names = eq["inlet"] if isinstance(eq["inlet"], list) else [eq["inlet"]]
        flows_list = [_component_flows(inlets[n]) for n in inlet_names if n in inlets]
        temps = [inlets[n].get("_T", 298.15) for n in inlet_names if n in inlets]
        pressures = [inlets[n].get("_P", 101325.0) for n in inlet_names if n in inlets]
        return run_mixer(
            inlet_molar_flows_list=flows_list,
            inlet_temperatures_K=temps,
            inlet_pressures_Pa=pressures,
        )

    # ────────────────────────────────────────────
    #  便捷方法
    # ────────────────────────────────────────────
    def solve_full(
        self,
        max_iter: int = 500,
        tol: float = 1e-4,
        acceleration: str = "wegstein",
        omega: float = 0.5,
        verbose: bool = False,
        do_post_process: bool = True,
    ) -> Dict[str, Any]:
        """完整求解：收敛循环 + 后处理

        Returns
        -------
        dict with keys:
            - converged, iterations, residual
            - streams : 所有流股
            - feeds   : 新鲜进料
            - equipment_summary : 各设备简要结果
            - post_process : 完整设备设计参数（do_post_process=True 时）
        """
        # ① 收敛循环
        solve_result = self.solve(
            max_iter=max_iter, tol=tol,
            acceleration=acceleration, omega=omega,
            verbose=verbose,
        )

        result = {
            "converged": solve_result["converged"],
            "iterations": solve_result["iterations"],
            "residual": solve_result["residual"],
            "streams": solve_result["streams"],
            "feeds": solve_result["feeds"],
            "equipment_summary": self._build_equipment_summary(
                solve_result["equipment_log"]),
        }

        # ② 后处理
        if do_post_process and solve_result["converged"]:
            result["post_process"] = self.post_process(
                solve_result["streams"],
                solve_result["equipment_log"],
            )

        return result

    def _build_equipment_summary(self, log: List[Dict]) -> List[Dict]:
        """构建设备结果摘要"""
        summary = []
        for entry in log:
            s = {"id": entry["id"], "type": entry["type"],
                 "inlet_streams": entry["inlet_streams"],
                 "outlet_streams": entry["outlet_streams"]}
            # 提取关键数值
            outlet = entry.get("outlet", {})
            if isinstance(outlet, dict) and "_T" in outlet:
                s["outlet_T_K"] = outlet.get("_T")
                s["outlet_P_Pa"] = outlet.get("_P")
                total = sum(v for k, v in outlet.items() if not k.startswith("_"))
                s["outlet_total_mol_per_s"] = total
            summary.append(s)
        return summary
