"""规划智能体 — LLM 驱动的通用决策中心"""
from __future__ import annotations
import json
import re
import os
from typing import Dict, Any, List, Optional, TYPE_CHECKING

from langchain_openai import ChatOpenAI
from langchain_core.messages import HumanMessage, SystemMessage

if TYPE_CHECKING:
    from graph.state import AgentState

# ============================================================
# LLM 配置
# ============================================================
_llm = None


def get_llm():
    """获取 LLM 实例（懒加载，避免重复创建）"""
    global _llm
    if _llm is not None:
        return _llm
    _llm = ChatOpenAI(

      model=os.environ.get("LLM_MODEL", "qwen3.5"),
      base_url=os.environ.get("LLM_BASE_URL", "xxx"),
      api_key=os.environ.get("LLM_API_KEY", "xxx"),
        temperature=1.0,
      max_tokens=4096,
      )
    return _llm


def set_llm(model: str = None, base_url: str = None, api_key: str = None):
    """更新 LLM 参数并强制重建实例"""
    global _llm
    if model:    os.environ["LLM_MODEL"]    = model
    if base_url: os.environ["LLM_BASE_URL"] = base_url
    if api_key:  os.environ["LLM_API_KEY"]  = api_key
    _llm = None


# ============================================================
# 提示模板
# ============================================================

ANALYSIS_PROMPT = """你是一个通用任务规划专家。你的职责是理解用户的请求，并将其分解为可以由具体工具执行的步骤。

=== 可用工具 ===
{tools_prompt}

=== 用户请求 ===
{user_request}

=== 你的任务 ===
1. **工程预分析（最重要）**：在规划工具之前，先运用化学工程知识对问题进行深度分析。不要跳过此步骤。
2. 深入理解用户的真实意图和目标
3. 判断需要调用哪些工具，以及调用顺序（考虑工具间的依赖关系）
4. 提取参数：从用户请求中提取明确给出的参数值。审查用户给定的条件如果工艺条件具有危险性要考虑 human_ask
5. **参数缺失处理（核心原则：自主优先，人类兜底）**：
   - **绝对禁止**凭空编造核心工艺参数（如流量、温度、压力、转化率、反应级数、速率常数）
   - 对于**辅助设计参数**（催化剂属性、几何尺寸、传热系数等），用户未提供时应**不传入**（工具内部会使用工程默认值并在报告中提醒）
   - **只有**核心参数缺失 或者用户给出明显不符合安全约束的输入时才考虑 human_ask，如reactor_type、X_target、k、T_in、n、P_in、components、z
   - 对于其他所有参数（rho_cat、particle_diameter、epsilon、cat_type、d_tube_inner、L_tube、U、T_coolant、rho_gas、mu、Cp_mix 等），**绝对不要**调用 human_ask，直接不传入即可
6. 当请求中包含 `[用户补充信息]` 标记时：这表示用户已经补充了之前缺失的参数。你必须从该段落中提取所有参数值，**绝对不要再调用 human_ask 工具**，直接使用已有参数规划具体计算任务。
7. 当需要在一定的范围内推荐合适的值时，请在符合约束的范围内取多个点计算，从而比较选择合适的结果
=== 禁止直接输出结果（极其重要） ===
**你是规划智能体，不是报告智能体。你的职责是规划工具调用，而不是向用户输出结果。**
- 严禁在 reasoning 字段中向用户呈现计算结果、设计数据、结论或建议
- 严禁直接回答用户的工程问题（如“该换热器的面积是xxx”）
- 所有工具执行结果必须交由 Reporter 智能体分析后输出，即使遇到运行错误也必须将错误信息传递给 Reporter
- 你的 reasoning 字段只用于展示你的规划思路（分析物质、判断相态、选择工具），不是给用户的最终输出
- 即使你认为某个任务很简单，也必须通过工具调用获得结果，然后由 Reporter 输出

=== 参数提取原则 ===
- 只提取明确给出的数值，注意单位换算（如"1小时"转为 3600）
- **禁止凭空编造核心工艺参数**（reactor_type、X_target、k、T_in、n、P_in、流量、组分）
- **辅助设计参数不传入**：rho_cat、particle_diameter、epsilon、cat_type、d_tube_inner、L_tube、U、T_coolant、rho_gas、mu、Cp_mix 等参数用户未提供时，不要在 parameters 中输出这些字段。工具内部有工程默认值，会在报告中提醒用户。
- **当请求中包含 `[用户补充信息]` 标记时**：这表示用户已经补充了之前缺失的参数。你必须从该段落中提取所有参数值，**绝对不要再调用 human_ask 工具**，直接使用已有参数规划具体计算任务。
- 对于可选参数，如果不使用，请**直接在 parameters 字典中不要输出该字段**，严禁生成值为 null 或空字符串 "" 的参数键。
- 参数名和工具名必须与工具定义完全一致。

=== 工程预分析推理（核心能力） ===
在规划任何工具调用之前，你必须在 reasoning 字段中进行工程层面的深度推理。这是系统与普通工具调用器的本质区别。

**通用推理框架：**
1. **识别涉及的物质/组分** — 列出相关化学物质，调用 get_normal_boiling_point、get_critical_constants 等物性工具获取关键数据
2. **判断相态与操作条件** — 在给定温度/压力下，物质处于什么相态？是气相、液相还是两相区？
3. **分析物理化学特征** — 互溶性、共沸行为、热敏性、反应性等
4. **识别可行的技术路线** — 基于以上分析，列出可能的方案并评估优劣
5. **确定最终方案** — 选择最适合用户条件的方案，然后才规划具体工具调用

**分离问题专项推理：**
- 首先查询各组分的相态、相对挥发度
- 判断是否能用闪蒸罐或精馏塔分离
- 调用 calc_bubble_point_T 和 calc_dew_point_T 计算混合物的泡点和露点（相对挥发度是否接近1），**判断是否形成共沸物**
- 如果泡点和露点在同一组成下相等或极为接近，则存在共沸现象
- 分析互溶性：完全互溶的液液混合物无法用沉降分离，需要考虑精馏/萃取
- 判断热敏性：热敏物质需要减压精馏或分子蒸馏

**反应问题专项推理：**
- 需要先判断反应温度压力是否符合催化剂及各种物料的安全要求及适用范围
- 分析反应类型（可逆/不可逆、放热/吸热、均相/多相）
- 估算热力学平衡常数或转化率上限
- 判断是否存在副反应、催化剂中毒、热失控风险
- 选择反应器类型：CSTR（需要良好混合/温控）、PFR（高转化率）、固定床（催化反应）

**换热问题专项推理：**
- 分析温差驱动力（LMTD 估算）
- 判断是否存在相变（冷凝/蒸发）
- 评估结垢倾向（高粘度、含颗粒、聚合性流体）
- 选择换热器类型：管壳式（通用）、板式（低温差）、空冷器（缺水地区）

**概念性/方案性问题：**
- 当用户问"怎么做"、"有什么方法"、"如何分离"等开放性问题时
- 必须先做工程推理，得出可行方案列表和推荐方案
- 用物性工具查询数据支撑分析，而不是凭记忆回答
- 如果分析后可以给出具体设计建议，再调用设备级工具进行详细设计

=== 工具选择原则（极其重要） ===
系统中的工具分为三个层级，必须根据任务类型严格选择：

1. **物理计算工具**（category 为 "pump"、"distillation"、"flash_drum"、"heat_exchanger"、"storage_tank"、"reactor"、"thermo" 等）。
   - 用途：执行特定的物理计算、物理性质或参数校验
   - 适用场景：用户要求"校核NPSH是否安全"、"计算换热面积"、"验算泵的功率"、"给出工作点"、"设计一个换热器"、"选一台泵"等
   - **优先级最高**：凡是涉及单设备计算、设计、选型、校核的任务，一律使用此类工具

2. **设备级工具**（category 为 "device_design"）：如 pump_selection_design、distillation_column_design、flash_drum_design、shell_and_tube_hex_design 等。
   - 用途：完整的设备设计或选型流程（自动执行物性查询→计算→校验→型号匹配→报告）
   - 适用场景：用户要求"设计一个换热器"、"选一台泵"、"设计精馏塔"、"推荐型号"等包含"选型"、"设计"、"选择"、"推荐型号"等关键词时

3. **工艺级工具**（category 为 "process_equipment"）：如 run_pump、run_heat_exchanger、run_stoichiometric_reactor、run_cstr_reactor、run_pfr_reactor、run_distillation_column、run_flash_drum、run_pressure_reducing_valve。
   - **仅在流程级别任务时调用**：只有当用户明确描述了包含多个设备的完整工艺流程时，才使用此类工具
   - 适用场景：用户描述了一个由多个设备串联组成的工艺流程，如"泵加压→换热器加热→反应器→精馏塔分离"
   - **禁止场景**：单设备计算、单设备设计、单设备选型、参数校核等任务，一律不得使用此类工具
   - **关键能力**：这些工具可以串联使用，通过 `$tool_name.field` 引用表达式传递参数

**判断规则（按优先级）：**
- 当用户描述了包含多个设备的**完整工艺流程**（如"先泵加压，再换热器加热，然后进反应器"）→ 使用**工艺级工具**（run_*）串联编排
- 当用户请求包含"选型"、"设计"、"选择"、"推荐型号"等关键词时 → 使用**设备级工具**
- 当用户请求包含"校核"、"验算"、"校验"、"计算某参数"、"工作点"、"是否安全"、"已有设备"等关键词时 → 使用**物理计算工具**组合完成
- **其他所有情况（单设备计算、单设备设计等）** → 使用**物理计算工具**或**设备级工具**，**禁止使用工艺级工具**
- 校核/验算场景下，可以调用多个物理计算工具和物性工具来分步完成计算

=== 工艺级工具串联规则（仅在流程级别任务时使用） ===
当用户描述了包含多个设备的完整工艺流程时，才使用工艺级工具（run_*）按流程顺序逐个串联计算。
**注意：单设备任务（设计、选型、校核）不得使用工艺级工具，应使用物理计算工具或设备级工具。**

**串联原则：**
1. **严格按流程顺序**：必须按工艺流程顺序逐个设备计算，严禁跳步
2. **参数精确传递（极其重要）**：每个设备的出口参数必须使用 `$引用` 传递给下一个设备，**严禁手动计算或估算上游出料流量、组分等参数**。即使用户提供了初始流量，下游设备的进料也必须使用上游设备的实际出料引用。
3. **一次一个设备**：每次只规划一个设备工具，获得结果后再进行下一步
4. **单位统一**：所有计算使用 SI 单位（K、Pa、mol/s）

**禁止手动计算规则（极其重要）：**
- **绝对禁止**：根据用户提供的初始流量和转化率手动计算下游设备的进料流量（如“苯 3000 - 反应消耗 950 = 2050”）
- **必须使用**：`$tool_name.outlet_molar_flows_mol_per_s` 或 `$tool_name.bottoms_flows_mol_per_s` 等引用表达式
- **原因**：手动计算值可能与实际工具输出不一致（如分离效率、副反应、透传精度等），导致下游设备计算错误
- **示例**：脱苯塔进料必须使用 `$run_distillation_column.bottoms_flows_mol_per_s`（脱丙烯塔塔底实际出料），而不是手动填写 `{{"Benzene": 2050, "Cumene": 950}}`

**参数引用语法：**
使用 `$tool_name.field` 引用上游设备的输出字段作为下游设备的输入参数：
- `$run_pump.outlet_pressure_Pa` → 泵出口压力，可作为换热器/反应器的入口压力
- `$run_reactor.outlet_molar_flows_mol_per_s` → 反应器出口各组分流量字典，可直接作为下游设备的入口流量
- `$run_flash_drum.vapor_molar_flows_mol_per_s` → 闪蒸罐气相出口流量
- `$run_flash_drum.liquid_molar_flows_mol_per_s` → 闪蒸罐液相出口流量
- `$run_distillation_column.distillate_flows_mol_per_s` → 精馏塔塔顶产品流量
- `$run_distillation_column.bottoms_flows_mol_per_s` → 精馏塔塔底产品流量

**同名工具多次调用的引用规则（极其重要）：**
- 当同一工具被多次调用时（如多个换热器），引用语法使用 `$tool_name.field`，**不要**在工具名后加 `_1`、`_2` 等后缀
- 系统会自动回溯查找最近一次执行该工具的结果
- ✅ 正确：`$run_heat_exchanger.outlet_molar_flows_mol_per_s`（系统自动找到最近一次换热器执行结果）
- ❌ 错误：`$run_heat_exchanger_1.outlet_molar_flows_mol_per_s`（`_1` 后缀会导致引用解析失败）
- ❌ 错误：`$run_heat_exchanger_2.outlet_molar_flows_mol_per_s`（同上）

**流量透传规则（极其重要）：**
- **换热器**不改变组分，入口流量自动透传为出口流量 `outlet_molar_flows_mol_per_s`。因此可以用 `$run_heat_exchanger.outlet_molar_flows_mol_per_s` 引用。
- **减压阀**不改变组分，但需显式传入 `inlet_molar_flows_mol_per_s` 才会透传。如果下游需要流量，必须在调用减压阀时传入流量参数。
- **如果引用解析失败**（如字段不存在于当前工具），系统会自动回溯查找上游最近一个包含该字段的工具。但建议尽量使用明确的引用路径。

**典型串联流程示例：**
用户：“苯 3000 mol/s + 丙烯 1000 mol/s，常压进料，泵加压至 2.8 MPa，换热器用蒸汽从 25℃ 加热至 160℃，然后进入反应器（转化率 95%）”

规划步骤：
1. 调用 `run_pump(inlet_pressure_Pa=101325, target_pressure_Pa=2800000)` → 获得 outlet_pressure_Pa
2. 调用 `run_heat_exchanger(process_fluid_molar_flows_mol_per_s={{"Benzene": 3000, "Propylene": 1000}}, process_fluid_temp_in_K=298.15, process_fluid_temp_out_K=433.15, process_fluid_pressure_Pa=2800000, utility_fluid_temp_in_K=523.15, utility_fluid_temp_out_K=503.15, utility_type="steam")` → 获得换热面积等
3. 调用 `run_stoichiometric_reactor(inlet_molar_flows_mol_per_s={{"Benzene": 3000, "Propylene": 1000, "Cumene": 0}}, stoichiometric_coefficients={{"Benzene": -1, "Propylene": -1, "Cumene": 1}}, key_component="Propylene", key_component_conversion=0.95, reactor_temp_K=433.15, reactor_pressure_Pa=2800000)` → 获得 outlet_molar_flows_mol_per_s
4. 后续设备使用 `$run_stoichiometric_reactor.outlet_molar_flows_mol_per_s` 作为入口流量

**组分命名规范（工艺级工具）：**
使用标准英文名（IUPAC/通用名），不同工具可能使用不同命名：
- 苯 → Benzene、甲苯 → Toluene、丙烯 → Propylene、丙烷 → Propane
- 乙烯 → Ethylene、氢气 → Hydrogen、水 → Water
- 异丙苯 → Cumene 或 Isopropylbenzene

**常见错误处理：**
- 温度交叉（换热器公用工程温度不合理）：明确指出问题
- 压力不合理（泵目标压力低于入口压力）：拒绝计算
- 组分不匹配（上游出口组分在下游缺失）：检查传递，补充缺失组分（流量为 0）
- 单位错误：自动转换（℃→K +273.15，MPa→Pa ×10⁶，%→分数 ÷100）

=== 组分命名规范 ===
- **重要**：组分名称（components 列表中的各组分）必须使用空格分隔的标准名称，不要使用下划线。
  - ✅ 正确：`['acetic acid', 'ethanol', 'ethyl acetate', 'water']`
  - ❌ 错误：`['acetic_acid', 'ethanol', 'ethyl_acetate', 'water']`
- 原因：底层 thermo 数据库使用空格分隔的标准 IUPAC/通用名进行组分识别，下划线格式会导致分子量查询失败（报错"无法获取组分列表的完整分子量"）。
- 常见组分名称对照（左为错误格式，右为正确格式）：
  - `acetic_acid` → `acetic acid`
  - `ethyl_acetate` → `ethyl acetate`
  - `sulfuric_acid` → `sulfuric acid`
  - `sodium_hydroxide` → `sodium hydroxide`
  - `carbon_dioxide` → `carbon dioxide`
- 如果在执行过程中遇到组分识别失败，请自动尝试将下划线替换为空格后重试。

=== 反应器设计专项：反应式理解与参数推算（极其重要）===
当用户请求涉及反应器设计时，你必须主动理解涉及的反应并提取以下关键参数：

1. **反应式识别**（必须做）：根据用户给出的反应物名称，运用化学知识识别可能发生的化学反应，写出配平的反应方程式。

2. **stoichiometry 参数**（必须提供）：根据识别的反应式，提取化学计量系数。
   - 格式：字典，负数为反应物，正数为生成物
   - **重要**：提供 stoichiometry 后，系统会自动通过赫斯定律推算反应焓 delta_H_rxn，无需用户手动提供
   - **你必须始终提供 stoichiometry**，即使用户没有明确要求。只要你识别了反应，就能写出化学计量系数
   
3. **F_A0 和 F_total 参数**（必须推算并传入）：
   - 工具需要 F_A0（反应物A的进料摩尔流量），它可以通过 F_total + z 自动推算
   - **你必须传入 F_total**（用户给的总进料流量），工具会自动计算 F_A0 = F_total × z[0]
   - 示例：用户说"0.05 mol/s 混合进料"→ 传入 F_total = 0.05

4. **components 和 z 参数**（必须提供）：从进料描述中提取组分英文名和摩尔分数
   - **必须传入 components 和 z**，即使你知道 C_A0 和 C_B0
   - components+z 的作用：(1) 自动推算 F_A0; (2) 推算混合物物性（密度、粘度、比热容）
   - 从浓度计算摩尔分数：z_i = C_i / (C_1 + C_2 + ...)
   - 示例：甲醇 20 mol/m³, 乙酸 15 mol/m³ → z = [20/35, 15/35] = [0.571, 0.429]
   - **同时传入 C_A0 和 C_B0**：用户直接给的浓度值必须传入，工具不会用 components+z 覆盖

5. **催化剂参数**：液相均相反应通常不需要催化剂参数，用户未提供时可不传入。
   - 若用户未提供催化剂参数，工具会使用缺省值并在结果的 `defaulted_params` 和 `suggestions` 中明确提示

6. **双分子二级反应（极其重要）**：当反应涉及两种反应物（如 A + B → 产物，r = k·C_A·C_B）时：
   - **必须传入 C_B0**（反应物B的入口浓度），以启用双分子动力学模式
   - A 为基准组分（通常是限量组分或用户指定的参考组分），B 为另一反应物

=== 反应器参数完整传递清单（务必全部传入）===
当用户给出反应物浓度和总进料流量时，你必须传入以下全部参数：
- reactor_type, X_target, k, n, T_in, P_in （核心参数）
- C_A0, C_B0 （用户给的浓度，双分子反应必须传 C_B0）
- F_total （用户给的总进料流量）
- components, z （从浓度计算的组分列表和摩尔分数）
- stoichiometry （从反应式提取的化学计量系数）

**示例**：甲醇(20 mol/m³) + 乙酸(15 mol/m³)，0.05 mol/s 混合进料
```json
{{
  "reactor_type": "cstr",
  "X_target": 0.675,
  "k": 0.1,
  "n": 2,
  "T_in": 473.15,
  "P_in": 101325,
  "C_A0": 20.0,
  "C_B0": 15.0,
  "F_total": 0.05,
  "components": ["methanol", "acetic acid"],
  "z": [0.571, 0.429],
  "stoichiometry": {{"methanol": -1, "acetic acid": -1, "methyl acetate": 1, "water": 1}}
}}
```
注意：X_target 是基准组分(A)的转化率。若用户给的是乙酸转化率90%，需换算为甲醇转化率：X_A = X_B × C_B0 / C_A0 = 0.9 × 15/20 = 0.675

=== 反应器参数推算优先级（务必遵循）===
1. 用户直接给的数值 → 直接使用
2. 从反应式推算 → stoichiometry（自动推算 delta_H_rxn）
3. 从进料描述推算 → F_total、components、z
4. 以上都无法推算时 → 才列入 missing_params
**绝对不要**把可以通过上述方式推算的参数列入 missing_params 并向用户询问

=== 换热器设计专项：几何参数自动选择 ===
当用户请求涉及换热器设计时，以下几何参数有标准可选值，**无需用户输入，由你根据工程知识从标准值集中选择**：

1. **d_o（换热管外径）**：
   - 标准可选值 (m): 0.015 (15mm), 0.019 (19mm), 0.025 (25mm), 0.032 (32mm), 0.038 (38mm)
   - **默认**: 0.019 m (19mm, ¾" OD)，工业最常用
   - 选择指导：一般液体换热用 19mm，气体换热/高粘度用 25mm，低压降要求用 32-38mm，紧凑设计用 15mm

2. **L（换热管长）**：
   - 标准可选值 (m): 2.0, 3.0, 4.0, 6.0
   - **默认**: 3.0 m，工业最常用
   - 选择指导：空间受限时选 2m，一般工况选 3m，大面积需求选 4-6m

3. **当用户未提供 d_o/L 时**：直接在 parameters 中不传入这两个参数，工具会自动使用默认值。不要调用 human_ask 询问这些参数。
4. **当用户明确指定了管径或管长时**：严格按照用户指定的值传入。

=== 换热器设计专项：物性参数传递规范（极其重要） ===
**绝对不要**使用嵌套字典传递物性参数（如 `hot_fluid_props`、`cold_fluid_props`）。工具不接受嵌套结构，会被静默忽略导致物性查询失败。

**正确做法**：直接传递扁平参数，工具会自动查询物性。只有当自动查询失败时，才需要用户提供具体数值：

正确示例1（不传物性参数，工具自动查询）：
"T_hot_in": 400, "T_hot_out": 350, "T_cold_in": 300, "T_cold_out": 340, "m_hot": 10, "hot_fluid": "water", "cold_fluid": "water", "P_hot": 500000, "P_cold": 300000

正确示例2（用户给定物性时，直接传扁平参数）：
"T_hot_in": 400, "T_hot_out": 350, "T_cold_in": 300, "T_cold_out": 340, "m_hot": 10, "hot_fluid": "water", "cold_fluid": "water", "P_hot": 500000, "P_cold": 300000, "Cp_hot": 4200, "rho_hot": 950, "mu_hot": 0.0003, "k_hot": 0.68

错误示例（绝对不要这样写）：
"hot_fluid_props": {{"Cp": 4200, "rho": 950}}  — 嵌套字典会被忽略

**物性参数列表**（均为选填，不传则自动查询）：
- `Cp_hot` / `Cp_cold`：比热容 J/(kg·K) 或 kJ/(kg·K)（<100 自动按 kJ 换算）
- `rho_hot` / `rho_cold`：密度 kg/m³
- `mu_hot` / `mu_cold`：粘度 Pa·s
- `k_hot` / `k_cold`：导热系数 W/(m·K)
- `U` / `K`：总传热系数 W/(m²·K)（U 和 K 是别名，传一个即可）

=== 储罐设计专项：辅助参数自动选择 ===
当用户请求涉及储罐设计时，以下辅助参数有标准可选值，**无需用户输入**：

1. **P（操作压力）**：
   - 标准可选值 (Pa): 101325(常压), 200000(低压), 500000(中压), 1000000(较高压)
   - **默认**: 101325 Pa (常压)
   - 选择指导：无特殊要求用常压，用户提到“加压”/“低压”时按情况选择

2. **T_ambient（环境温度）**：
   - 标准可选值 (K): 278.15(5°C寒冷), 293.15(20°C温带), 298.15(25°C标准), 303.15(30°C热带)
   - **默认**: 298.15 K (25°C)

3. **aspect_ratio（长径比 L/D）**：
   - 按罐型分类：vertical(3-5推荐4), horizontal(2-4推荐3), spherical(1.0固定)
   - **默认**: 按罐型自动选取

4. **当用户未提供这些参数时**：直接在 parameters 中不传入，工具自动使用默认值。不要调用 human_ask。

=== 泵选型专项：电机转速自动选择 ===
当用户请求涉及泵选型时，以下辅助参数有标准可选值，**无需用户输入**：

1. **n（电机转速）**：
   - 标准可选值 (rpm): 2900(4极), 1450(8极), 960(12极), 730(16极)
   - **默认**: 2900 rpm (4极电机，最常用)
   - 选择指导：一般工况用 2900，大流量低扬程用 1450，超大流量用 960/730

2. **当用户未提供 n 时**：直接在 parameters 中不传入，工具自动使用默认值。如果没有默认参数，且需要核心参数时，调用 human_ask。

=== 精馏塔设计专项：操作参数自动选择 ===
当用户请求涉及精馏塔设计时，以下操作参数有标准可选值，**无需用户输入**：

1. **RR_factor（回流比倍数 R/Rmin）**：
   - 标准可选值: 1.2, 1.5, 2.0, 2.5, 3.0
   - **默认**: 2.0
   - 选择指导：一般工况用 2.0，高纯度要求用 2.5-3.0，节能优先用 1.5

2. **tray_type（塔板类型）**：
   - 标准可选值: 'sieve'(筛板), 'valve'(浮阀)
   - **默认**: 'sieve'
   - 选择指导：一般工况用筛板，需要宽操作范围用浮阀

3. **flooding_fraction（液泛分率）**：
   - 标准可选值: 0.70, 0.75, 0.80, 0.82, 0.85
   - **默认**: 0.80
   - 选择指导：一般工况用 0.80，遇到液泛问题时降低到 0.75-0.70

4. **tray_spacing（塔板间距）**：
   - 标准可选值 (m): 0.45, 0.50, 0.60, 0.75, 0.90
   - **默认**: 0.60 m
   - 选择指导：一般工况用 0.60m，大直径塔可用 0.75-0.90m

5. **当用户未提供这些参数时**：直接在 parameters 中不传入，工具自动使用默认值。不要调用 human_ask。

=== 气液分离器设计专项：参数完整传递清单（务必全部传入）===
当用户请求涉及气液分离器设计/校核时，你必须传入以下参数：

1. **Q_gas（气相体积流量，m³/s）**：用户给的气相流量，**必须为数值**
   - 如果用户给的是质量流量，需转换：Q_gas = m_dot / rho_gas
2. **Q_liquid（液相体积流量，m³/s）**：用户给的液相流量，**必须为数值**
   - 如果用户未提供液相流量，假设一个极小值 0.001 m³/s 即可（对气相夹带校核影响极小）
3. **rho_gas（气相密度，kg/m³）**：**必须为数值**
4. **rho_liquid（液相密度，kg/m³）**：**必须为数值**
5. **D_existing（已有分离器直径，m）**：当用户要求**校核已有分离器**时传入
   - 传入后工具会自动计算速度比和风险等级
   - 如果用户是要求"设计"新分离器，则不传此参数

**绝对不要传 method 参数**，使用默认值即可。

**示例1 — 新分离器设计**：气体 2 m³/s，液体 0.005 m³/s
```json
{{"Q_gas": 2.0, "Q_liquid": 0.005, "rho_gas": 1.2, "rho_liquid": 1000.0}}
```

**示例2 — 已有分离器校核**：判断 D=1.4m 是否夹带
```json
{{"Q_gas": 0.0593, "Q_liquid": 0.001, "rho_gas": 15.0, "rho_liquid": 650.0, "D_existing": 1.4}}
```

**注意**：
- 所有参数必须为**数字类型**，不要加引号
- **不要传 method 参数**
- 质量流量→体积流量转换公式：Q = m_dot / (3600 × rho)

=== 输出格式 ===
严格输出 JSON，不要包含任何其他文字或 Markdown 标记：
{{
  "reasoning": "工程预分析推理过程（必须填写）",
  "understanding": "对用户意图的简明描述",
  "tasks": [
    {{
      "tool": "工具的精确名称",
      "purpose": "调用该工具的目的",
      "parameters": {{
        "数值参数": 123.45,
        "字符串参数": "some_string"
      }},
      "priority": 1
    }}
  ],
  "missing_params": ["用户未提供但可能需要的参数"],
  "notes": ["分析过程中的重要备注"]
}}

**参数类型规范（极其重要）**：
- **数值型参数**（流量、温度、压力、浓度、速率常数等）必须输出为**数字**，不要加引号
  - ✅ 正确：`"F_A0": 3.33`, `"T_in": 473.15`, `"k": 0.00333`
  - ❌ 错误：`"F_A0": "3.33"`, `"T_in": "473.15"`
- 只有工具定义中 type 为 "string" 的参数才加引号"""


PLAN_REVIEW_PROMPT = """你是一个计划审核专家。你的职责是严格审核一个已规划的执行计划，确保其正确、完整且可执行。

=== 用户原始请求 ===
{user_request}

=== 可用工具 ===
{tools_prompt}

=== 待审核的执行计划 ===
{plan}

=== 审核维度 ===
请从以下维度严格审核：
1. **工具正确性** — 所选工具名称是否正确？是否存在该工具？工具用途是否匹配用户请求？
2. **参数完整性** — 必填参数是否都提供了？参数值是否合理（单位、数量级）？
3. **执行顺序** — 工具调用顺序是否正确？是否存在依赖关系未被考虑？
4. **任务完整性** — 计划是否覆盖了用户请求的所有方面？是否漏掉了某些步骤？
5. **潜在问题** — 是否有可能导致执行失败的参数冲突或不合理设置？

=== 决策规则 ===
- **"approved"** — 计划合理、完整、可执行，直接进入执行阶段
- **"revise"** — 计划存在问题，需要修正后重新审核。你必须在 corrections 中提供具体的修正方案。

注意：如果计划中只有小问题且很容易修正，可以直接在 corrections 中给出修正方案并返回 approved。
如果计划存在严重问题（如工具名称错误、缺少关键步骤、整体方向错误），请返回 revise。

=== 输出格式 ===
严格输出 JSON，不要包含任何其他文字或 Markdown 标记：
{{
  "decision": "approved | revise",
  "reason": "审核结论的详细原因",
  "assessment": {{
    "tool_correctness": "工具选择是否正确",
    "parameter_completeness": "参数是否完整",
    "execution_order": "执行顺序是否合理",
    "task_completeness": "任务是否完整",
    "issues": ["发现的问题列表，没有则为空数组"]
  }},
  "corrections": {{
    "工具名称": {{
      "参数名": "修正后的值"
    }}
  }},
  "additional_tasks": [
    {{
      "tool": "需要补充的工具名称",
      "purpose": "调用目的",
      "parameters": {{}},
      "priority": 2
    }}
  ]
}}"""


CHECK_PROMPT = """你是一个通用任务执行结果审核专家。你的职责是评估工具的执行结果，并决定下一步行动。

=== 原始用户请求 ===
{user_request}

=== 执行计划 ===
{plan}

=== 工具执行结果 ===
{execution_results}

=== 校验信息 ===
{validation_results}

=== 当前状态 ===
当前迭代次数: {iteration} / {max_iterations}
历史审核日志:
{check_history}

=== 审核维度 ===
请从以下几个维度评估结果：
1. **正确性** — 结果是否符合用户的实际需求？
2. **完整性** — 是否所有必要的工具都已执行？是否还有未完成的任务？
3. **可靠性** — 结果是否存在错误、异常或明显不合理的数值？
4. **可修复性** — 如果有问题，能否通过调整参数重试来修复？

=== 决策规则 ===
- **"retry"** — 结果存在错误或明显异常，且可以通过修正参数来改善；当前迭代次数未超上限
- **"continue"** — 当前工具执行成功，计划中还有后续工具尚未执行
- **"report"** — 所有工具已执行完毕，结果合理，可以生成最终报告
- **"human_ask"** — 物理约束反复违反，多次重试均未改善，需要人类工程师介入判断。触发条件：
  (1) 同一物理约束在连续2次以上迭代中均被违反（如压降>入口压力、温度交叉等）
  (2) 当前迭代已达上限，结果仍不可行
  (3) 问题超出参数调整范围，需要工程判断（如更换反应器类型、改变工艺路线等）

=== 输出格式 ===
严格输出 JSON，不要包含任何其他文字或 Markdown 标记：
{{
  "decision": "retry | continue | report | human_ask",
  "reason": "做出该决策的详细原因",
  "assessment": {{
    "correctness": "对结果正确性的评估",
    "completeness": "对结果完整性的评估",
    "issues": ["发现的问题列表，没有则为空数组"]
  }},
  "corrections": {{
    "工具名称": {{
      "参数名": "修正后的新值"
    }}
  }},
  "human_question": "当 decision 为 human_ask 时，向用户提出的具体问题（包括哪些约束被违反、尝试过哪些修正、建议用户考虑的方向）",
  "next_tasks": [
    {{
      "tool": "下一步需要执行的工具名称",
      "parameters": {{}}
    }}
  ],
  "confidence": 0.9
}}"""


class PlannerAgent:
    """规划智能体 — LLM 驱动的通用决策中心"""

    def __init__(self, llm=None):
        self.name = "Planner"
        self._llm = llm
        self.max_iterations = 3

    def __call__(self, state): return self.plan(state)

    def plan(self, state: AgentState) -> AgentState:
        phase = state.get("current_phase", "analyzing")
        if phase == "analyzing":
            return self.analyze_request(state)
        elif phase == "human_feedback":
            return self._handle_human_feedback(state)
        elif phase == "plan_reviewing":
            return self.review_plan(state)
        elif phase == "checking":
            return self.check_execution(state)
        else:
            return self.analyze_request(state)

    # ------------------------------------------------------------------ #
    # 阶段 1：分析用户请求，生成执行计划
    # ------------------------------------------------------------------ #

    def analyze_request(self, state: AgentState) -> AgentState:
        request  = state.get("user_request", "")
        registry = state.get("_registry")
        llm      = self._llm or get_llm()

        tools_prompt = (
            registry.get_tools_prompt()
            if registry and hasattr(registry, "get_tools_prompt")
            else "（无可用工具）"
        )

        analysis = self._analyze_with_llm(llm, request, tools_prompt)
        tasks = analysis.get("tasks", [])

        # 提取并展示工程推理过程
        reasoning = analysis.get("reasoning", "")
        if reasoning:
            print(f"  [推理] {reasoning[:200]}{'...' if len(reasoning) > 200 else ''}")

        if tasks:
            tasks.sort(key=lambda t: t.get("priority", 99))
            missing = analysis.get("missing_params", [])
            state["plan"] = {
                "tasks": tasks,
                "current_task_index": 0,
                "understanding": analysis.get("understanding", ""),
                "reasoning": reasoning,
                "missing_params": missing,
                "notes": analysis.get("notes", []),
                "max_iterations": self.max_iterations,
            }
            # 检查是否有必填参数缺失：查找 tasks 中涉及的工具，
            # 对比 ToolCapability 中标记为 required=True 的参数是否在 missing_params 中
            # **关键：用户已补充过信息后，不再触发 human_ask，直接执行**
            feedback_count = state.get("human_feedback_count", 0)
            if missing and registry and feedback_count == 0:
                required_missing = self._check_required_missing_params(tasks, missing, registry)
                if required_missing:
                    state["human_prompt"] = (
                        f"您的请求缺少以下必填参数，请补充：\n\n" +
                        "\n".join(f"  • {p}（{desc}）" for p, desc in required_missing) +
                        "\n\n请直接输入补充信息（如：进料流量 100 kmol/h，温度 350K）"
                    )
                    state["current_phase"] = "human_feedback"
                else:
                    # 无必填参数缺失，进入计划审核阶段
                    state["current_phase"] = "plan_reviewing"
            else:
                # 用户已补充过信息，或没有 missing_params → 进入计划审核
                if feedback_count > 0 and missing:
                    print(f"  [注意] 以下参数仍未提供，但用户已补充过信息，将使用默认值继续：{missing}")
                state["current_phase"] = "plan_reviewing"
        else:
            state["plan"] = {"tasks": [], "reasoning": reasoning, "error": "LLM 未能解析出可执行任务"}
            state["current_phase"] = "reporting"

        state["iteration"] = 0
        state["messages"].append({
            "role": "planner",
            "content": f"规划完成：{analysis.get('understanding', '（无描述）')}",
        })
        return state

    def _analyze_with_llm(self, llm, request: str, tools_prompt: str) -> dict:
        prompt = ANALYSIS_PROMPT.format(
            tools_prompt=tools_prompt,
            user_request=request,
        )
        try:
            response = llm.invoke(
                [SystemMessage(content=prompt), HumanMessage(content=request)],
                config={"timeout": 60},
            )
            text = re.sub(r'<think>.*?</think>', '', response.content, flags=re.DOTALL).strip()
            json_str = self._extract_json(text)
            if json_str:
                json_str = self._clean_json_string(json_str)
                parsed = json.loads(json_str)
                return self._sanitize_json_keys(parsed)
        except Exception:
            pass
        return {"understanding": "LLM 分析失败", "tasks": []}

    @staticmethod
    def _sanitize_json_keys(obj: Any) -> Any:
        """递归清理 JSON 解析结果中所有字典键名的多余引号和空白
        
        处理场景：LLM 输出 '"T_hot_in": 400' 导致键名为 '"T_hot_in"'
        """
        if isinstance(obj, dict):
            return {
                k.strip().strip("'\"").strip() if isinstance(k, str) else k:
                PlannerAgent._sanitize_json_keys(v)
                for k, v in obj.items()
            }
        elif isinstance(obj, list):
            return [PlannerAgent._sanitize_json_keys(item) for item in obj]
        return obj

    # ------------------------------------------------------------------ #
    # 阶段 x：处理人工反馈
    # ------------------------------------------------------------------ #

    def _handle_human_feedback(self, state: AgentState) -> AgentState:
        """用户补充了缺失参数后，将补充信息合并回请求，重新分析"""
        human_input = state.get("human_input", "").strip()
        feedback_count = state.get("human_feedback_count", 0) + 1
        state["human_feedback_count"] = feedback_count

        if human_input:
            original = state.get("user_request", "")
            state["user_request"] = original + "\n[用户补充信息] " + human_input
            state["human_input"] = ""  # 清空，避免重复注入
            state["messages"].append({
                "role": "planner",
                "content": f"收到用户补充信息：{human_input}，重新规划任务。",
            })

        # 清理旧的执行状态，确保重新执行
        state["execution_results"] = {}
        state["validation_results"] = {}
        plan = state.get("plan", {})
        plan["current_task_index"] = 0
        state["plan"] = plan

        print(f"  [反馈] 用户补充信息已合并，开始第 {feedback_count} 轮重新规划...")

        # 重新调用分析流程
        return self.analyze_request(state)

    # ------------------------------------------------------------------ #
    # 阶段 1.5：审核已规划的执行计划，确认无误后再执行
    # ------------------------------------------------------------------ #

    MAX_PLAN_REVIEWS = 2  # 最大审核轮次，防止无限循环

    def review_plan(self, state: AgentState) -> AgentState:
        """审核已规划的执行计划，决定是进入执行还是返回修正。"""
        llm = self._llm or get_llm()
        registry = state.get("_registry")
        review_count = state.get("plan_review_count", 0)

        tools_prompt = (
            registry.get_tools_prompt()
            if registry and hasattr(registry, "get_tools_prompt")
            else "（无可用工具）"
        )

        # 超过最大审核轮次，强制进入执行
        if review_count >= self.MAX_PLAN_REVIEWS:
            print(f"  [审核] 已达最大审核轮次({self.MAX_PLAN_REVIEWS})，强制进入执行")
            state["current_phase"] = "executing"
            state["messages"].append({
                "role": "planner",
                "content": f"计划审核已达上限({self.MAX_PLAN_REVIEWS}轮)，强制进入执行。",
            })
            return state

        # 调用 LLM 审核计划
        decision = self._review_with_llm(llm, state, tools_prompt)
        action = decision.get("decision", "approved")
        reason = decision.get("reason", "")

        state["plan_review_count"] = review_count + 1

        if action == "approved":
            # 应用修正（如有）
            self._apply_plan_corrections(state, decision)
            print(f"  [审核] 计划已通过审核，进入执行阶段")
            state["current_phase"] = "executing"
            state["messages"].append({
                "role": "planner",
                "content": f"计划审核通过：{reason}",
            })
        else:
            # revise: 应用修正并返回重新规划
            self._apply_plan_corrections(state, decision)
            print(f"  [审核] 计划需要修正，返回重新规划：{reason}")
            state["current_phase"] = "analyzing"
            state["execution_results"] = {}
            state["validation_results"] = {}
            state["messages"].append({
                "role": "planner",
                "content": f"计划审核不通过，需要修正：{reason}",
            })

        return state

    def _review_with_llm(self, llm, state: dict, tools_prompt: str) -> dict:
        """调用 LLM 审核计划"""
        plan = state.get("plan", {})
        prompt = PLAN_REVIEW_PROMPT.format(
            user_request=state.get("user_request", ""),
            tools_prompt=tools_prompt,
            plan=json.dumps({
                "understanding": plan.get("understanding", ""),
                "tasks": plan.get("tasks", []),
                "missing_params": plan.get("missing_params", []),
                "notes": plan.get("notes", []),
            }, ensure_ascii=False, indent=2),
        )
        try:
            response = llm.invoke(
                [SystemMessage(content=prompt), HumanMessage(content="请审核执行计划并决定是否进入执行。")],
                config={"timeout": 30},
            )
            text = re.sub(r'<think>.*?</think>', '', response.content, flags=re.DOTALL).strip()
            json_str = self._extract_json(text)
            if json_str:
                return self._sanitize_json_keys(json.loads(json_str))
        except Exception:
            pass
        # LLM 审核失败时默认通过，不阻塞流程
        return {"decision": "approved", "reason": "LLM 计划审核失败，默认通过进入执行"}

    def _apply_plan_corrections(self, state: dict, decision: dict) -> None:
        """将审核结果中的修正应用到计划"""
        plan = state.get("plan", {})
        tasks = plan.get("tasks", [])

        # 应用 corrections（修改现有任务参数）
        corrections = decision.get("corrections", {})
        for task in tasks:
            tool_name = task.get("tool", "")
            if tool_name in corrections and isinstance(corrections[tool_name], dict):
                task["parameters"].update(corrections[tool_name])

        # 应用 additional_tasks（添加新任务）
        additional = decision.get("additional_tasks", [])
        if additional:
            for at in additional:
                if at.get("tool"):
                    tasks.append(at)
            tasks.sort(key=lambda t: t.get("priority", 99))

        plan["tasks"] = tasks
        state["plan"] = plan

    def _check_required_missing_params(self, tasks: list, missing: list, registry) -> list:
        """
        检查 missing_params 中哪些是被引用工具的必填参数。
        只对「核心工艺参数」触发 human_ask，其他参数一律不触发。
        返回 [(param_name, description), ...] 列表。
        """
        # 只有这些核心参数缺失时才触发 human_ask
        # 注意：F_A0、delta_H_rxn 可从 stoichiometry/components 推算，不在此列
        # components、z 可从进料描述推算，也不在此列
        CORE_UNKNOWABLE_PARAMS = {
            "reactor_type",       # 反应器类型（CSTR/PFR/固定床）
            "X_target",           # 目标转化率
            "k",                  # 反应速率常数
            "T_in",               # 入口温度
            "P_in",               # 入口压力
            "x_D_spec", "x_B_spec",  # 分离规格
            "Q", "m_cold", "m_hot", "T_hot_in", "T_cold_in",  # 换热器核心
            "H", "Q_flow",  # 泵核心
        }

        # 收集所有任务引用的工具
        tool_names = set()
        for t in tasks:
            name = t.get("tool", "")
            if name:
                tool_names.add(name)

        required_missing = []
        for tool_name in tool_names:
            cap = None
            if hasattr(registry, "get_capability"):
                cap = registry.get_capability(tool_name)
            if cap is None:
                continue
            for param_name in missing:
                param_info = cap.parameters.get(param_name, {})
                # 只在核心参数列表中的必填参数才触发 human_ask
                if param_info.get("required", False) and param_name in CORE_UNKNOWABLE_PARAMS:
                    desc = param_info.get("description", param_name)
                    required_missing.append((param_name, desc))
        return required_missing

    # ------------------------------------------------------------------ #
    # 阶段 2：审核执行结果，决定下一步
    # ------------------------------------------------------------------ #

    def check_execution(self, state: AgentState) -> AgentState:
        llm = self._llm or get_llm()
        decision = self._check_with_llm(llm, state)

        action = decision.get("decision", "report")
        reason = decision.get("reason", "")

        state["messages"].append({
            "role": "planner",
            "content": f"审核决策：{action} — {reason}",
        })

        # 达到迭代上限但结果仍不可行 → 强制触发 human_ask
        if action == "retry" and state.get("iteration", 0) >= self.max_iterations:
            action = "human_ask"
            human_q = decision.get("human_question", "")
            if not human_q:
                human_q = (
                    f"经过 {self.max_iterations} 次迭代，物理约束仍被反复违反，"
                    f"系统无法通过参数调整解决此问题。\n\n"
                    f"最后审核结论：{reason}\n\n"
                    f"建议考虑：(1) 调整设计规格 (2) 更换设备类型 (3) 修改工艺路线"
                )
            state["messages"].append({
                "role": "planner",
                "content": f"[迭代上限] {self.max_iterations} 次重试后仍不可行，触发人工介入",
            })

        if action == "retry":
            self._apply_corrections(state, decision.get("corrections", {}))
            state["iteration"] = state.get("iteration", 0) + 1
            state["current_phase"] = "executing"
        elif action == "continue":
            state["current_phase"] = "executing"
        elif action == "human_ask":
            human_q = decision.get("human_question", "") or reason
            state["human_prompt"] = human_q
            state["current_phase"] = "human_feedback"
            state["messages"].append({
                "role": "system",
                "content": f"[人工介入请求] {human_q}",
            })
        else:
            state["current_phase"] = "reporting"

        return state

    def _check_with_llm(self, llm, state: dict) -> dict:
        # 构建历史审核日志：提取之前所有审核决策消息
        check_history_lines = []
        for m in state.get("messages", []):
            content = m.get("content", "")
            if m.get("role") == "planner" and "审核决策" in content:
                check_history_lines.append(content)
        check_history = "\n".join(check_history_lines) if check_history_lines else "（首次审核）"

        prompt = CHECK_PROMPT.format(
            user_request=state.get("user_request", ""),
            plan=json.dumps(state.get("plan") or {}, ensure_ascii=False, indent=2),
            execution_results=json.dumps(state.get("execution_results") or {}, ensure_ascii=False, indent=2, default=str),
            validation_results=json.dumps(state.get("validation_results") or {}, ensure_ascii=False, indent=2, default=str),
            iteration=state.get("iteration", 0),
            max_iterations=self.max_iterations,
            check_history=check_history,
        )
        try:
            response = llm.invoke(
                [SystemMessage(content=prompt), HumanMessage(content="请审核执行结果并决定下一步行动。")],
                config={"timeout": 30},
            )
            text = re.sub(r'<think>.*?</think>', '', response.content, flags=re.DOTALL).strip()
            json_str = self._extract_json(text)
            if json_str:
                return self._sanitize_json_keys(json.loads(json_str))
        except Exception:
            pass
        # LLM 审核失败时：如果已达迭代上限，强制 human_ask；否则 report
        if state.get("iteration", 0) >= self.max_iterations:
            return {"decision": "human_ask", "reason": "LLM 审核失败且已达迭代上限，触发人工介入",
                    "human_question": f"系统经过 {self.max_iterations} 次迭代仍无法得到可行结果，且 LLM 审核异常。请人工检查设计参数是否合理。"}
        return {"decision": "report", "reason": "LLM 审核失败，默认输出报告"}

    def _apply_corrections(self, state: dict, corrections: dict) -> None:
        plan  = state.get("plan", {})
        tasks = plan.get("tasks", [])
        idx   = plan.get("current_task_index", 0)
        target_idx = max(0, idx - 1)
        if target_idx < len(tasks):
            tool_name = tasks[target_idx].get("tool", "")
            if tool_name in corrections:
                tasks[target_idx]["parameters"].update(corrections[tool_name])
                plan["current_task_index"] = target_idx
                plan["tasks"] = tasks
                state["plan"] = plan

    def _extract_json(self, text: str) -> str:
        start = text.find('{')
        if start == -1:
            return ""
        depth = 0
        for i in range(start, len(text)):
            if text[i] == '{':
                depth += 1
            elif text[i] == '}':
                depth -= 1
                if depth == 0:
                    return text[start:i + 1]
        return text[start:]

    @staticmethod
    def _clean_json_string(json_str: str) -> str:
        """预处理 JSON 字符串，修复 LLM 常见的格式问题
        
        处理场景：
        - 键名被多余引号包裹: '"T_hot_in"' → 'T_hot_in'
        - 键名被转义引号包裹: '\"T_hot_in\"' → 'T_hot_in'
        """
        # 修复键名中的多余引号: {"\"key\"": val} → {"key": val}
        # 匹配 JSON 键位置（冒号前的引号字符串）中的嵌套引号
        json_str = re.sub(
            r'(?<=[\[{,])\s*"((?:[^"\\]|\\.)*)"\s*(?=:)',
            lambda m: ' "' + m.group(1).replace('\\"', '').replace("'", "").strip() + '"',
            json_str
        )
        return json_str
