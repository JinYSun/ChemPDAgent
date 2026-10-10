# ChemPDAgent — Intelligent Chemical Process Equipment Design MCP Agent

> LLM-Driven Chemical Process Equipment Design System
> A chemical equipment design and process simulation system built on a **LangGraph multi-agent architecture + dual CoolProp/thermo property engines**

Describe a design task in plain natural language (e.g., "Design a benzene–toluene distillation column, feed 100 kmol/h, atmospheric pressure, saturated liquid feed"), and the system automatically completes the entire workflow: **requirement understanding → engineering planning → property calculation → equipment design → result verification → report generation**. It can also expose its tools through the **MCP protocol (stdio / HTTP)**.

---

## 1. Key Features

- **Fully LLM-driven workflow**: three collaborating agents — Planner (planning) → Executor (execution) → Reporter (reporting) — with support for retry, continue, plan review, and human-in-the-loop parameter supplementation
- **Dual property engines**: [thermo_helper.py](physics_engine/thermo_helper.py) calls **CoolProp HEOS** first (~124 pure fluids, NIST REFPROP-grade accuracy) and automatically falls back to the **thermo** library (PR/SRK EOS, UNIFAC, Henry, DIPPR)
- **Three-layer computation hierarchy**:
  1. `physics_engine` — low-level properties and physics calculations
  2. `device_tools` — six complete equipment-level design pipelines
  3. `process_engine` — ten chainable process units plus a recycle-converging flowsheet solver
- **Unified tool registry**: adding a new tool only requires "write the function + define a ToolCapability + register it" — no changes to agent code
- **SI units throughout**: K, Pa, mol/s, m³; every result is fully traceable to its basis (flows / density / area)

---

## 2. Project Structure

```
chempdagent/
├── main.py                     # Main entry: interactive / single task / tests / MCP server
├── server.py                   # MCP Server: stdio JSON-RPC and HTTP REST API
│
├── agents/                     # Agent definitions
│   ├── planner.py              #   PlannerAgent — requirement analysis, tool selection, result review
│   ├── executor.py             #   ExecutorAgent — dynamic tool invocation and chained execution
│   └── reporter.py             #   ReporterAgent — natural-language design reports
│
├── graph/                      # LangGraph workflow
│   ├── workflow.py             #   State-graph orchestration (plan→execute→review→report)
│   └── state.py                #   AgentState definition
│
├── tools/                      # Unified tool registry
│   ├── registry.py             #   UniversalToolRegistry / ToolCapability
│   ├── base.py                 #   Base registry and MCP JSON-RPC handling
│   ├── physics_tools.py        #   22 physics calculation tools
│   └── extended_physics_tools.py  # Extended property / phase-equilibrium / equipment-parameter tools
│
├── device_tools/               # Equipment-level design tools (6 complete pipelines)
│   ├── distillation_device.py      # Distillation column (FUG + hydraulics)
│   ├── flash_drum_device.py        # Flash drum
│   ├── heatexchanger_device.py     # Shell-and-tube heat exchanger
│   ├── reactor_device.py           # Reactor (zero/1st/2nd-order kinetics, multi-stage)
│   ├── pump_device.py              # Pump selection
│   └── storage_tank_device.py      # Storage tank
│
├── physics_engine/             # Physics calculation engine (low level)
│   ├── thermo_helper.py        #   Unified property engine (CoolProp + thermo, ~2800 lines)
│   ├── distillation.py         #   FUG method, column diameter/tray hydraulics, Fair flooding
│   ├── reactor.py              #   Fixed bed / fluidized bed / CSTR
│   ├── heatexchanger.py        #   LMTD, ε-NTU, area, pressure drop
│   ├── flash_drum.py           #   Rachford-Rice flash, separator sizing
│   ├── pump.py                 #   Pump power, NPSH, system curve
│   ├── storage_tank.py         #   Vertical/horizontal tanks, spheres, breathing losses
│   ├── pipe.py                 #   Darcy-Weisbach, Friedel two-phase flow
│   ├── valve.py                #   Control valve Cv (IEC 60534), safety valve (API 520)
│   ├── packed_bed.py           #   Ergun equation + 13 pressure-drop correlations
│   ├── instrument.py           #   Orifice flow meter (ISO 5167)
│   └── common.py               #   Material/energy balances, Antoine, Raoult
│
├── process_engine/             # Process unit library (chainable, 10 units)
│   ├── cstr/                   #   CSTR (12 independent calc_* tools)
│   ├── pfr/                    #   PFR / fixed bed
│   ├── stor/                   #   Stoichiometric reactor
│   ├── tray_distillation/      #   Tray distillation column
│   ├── heat_exchangers/        #   Heat exchanger
│   ├── flash_drum/             #   Flash drum
│   ├── pump/                   #   Pump
│   ├── pressure_reducing_valve/ #  Pressure reducing valve
│   ├── mixer/                  #   Mixer (multiple feeds, adiabatic mixing)
│   ├── multi_equipment_series/ #   run_* equipment-level wrappers
│   └── flowsheet/engine.py     #   FlowsheetEngine: Tear Stream + Wegstein convergence
│
├── process_tools/              # Registers process_engine as generic tools
└──  test/                       # 35+ design examples and device/engine comparison scripts

```

---

## 3. Architecture

### 3.1 Multi-Agent Workflow

```
Natural-language user input
            ↓
[Planner] Engineering pre-analysis + task decomposition + parameter extraction
            ↓
[Executor] Dynamically invokes registered tools per the plan (supports chained multi-tool execution)
            ↓
[Planner] Reviews execution results
      ├─ retry    → re-execute
      ├─ continue → execute the next task
      └─ report   → [Reporter] generate the design report
                          ↓
                     Final report
```

- The workflow is a **LangGraph** state graph; all decisions are made by the LLM, with no hardcoded equipment rules
- Supports the `human_ask` tool: when core parameters are missing, the flow pauses to ask the user; auxiliary parameters automatically take engineering defaults
- Infinite-loop protection: both human supplementation and plan review have round-count limits

### 3.2 Dual Property-Engine Routing

| Property category | Primary engine | Fallback engine |
|---|---|---|
| Density / molar volume | CoolProp (high-pressure preferred) | thermo (EOS / ideal gas) |
| Enthalpy / entropy / internal energy / Gibbs / Helmholtz | CoolProp | thermo |
| Cp / Cv / γ, speed of sound, Z, fugacity | CoolProp | thermo EOS |
| Vapor pressure / enthalpy of vaporization | CoolProp | thermo Wagner / DIPPR / Antoine |
| Viscosity / thermal conductivity / surface tension | CoolProp | thermo |
| Gas-phase diffusion | Chapman-Enskog (±5%) | — |
| Liquid-phase diffusion | Wilke-Chang (±10%) | Hayduk-Minhas |
| Safety / environmental properties, molecular descriptors | thermo | — |

K-values for mixture VLE are routed automatically based on operating conditions: **Henry (gas + water, <1.5 MPa) → PR/SRK EOS (supercritical components / high pressure) → UNIFAC → Wilson → NRTL → Raoult**.

### 3.3 Process Simulation Engine

The [FlowsheetEngine](process_engine/flowsheet/engine.py) describes the process topology through a configuration dictionary (feeds / equipment / recycles), automatically identifies the **Tear Stream**, and accelerates recycle-stream iteration to convergence using the **Wegstein secant method**. It can also perform post-processing to output complete equipment design parameters (heat-transfer area, actual tray count, etc.).

---

## 4. Capability Overview

### Equipment-Level Design (device_tools, 6 types)

| Equipment | Core entry point | Key outputs |
|---|---|---|
| Distillation column | `distillation_device._full_design()` | N_min, R_min, R, theoretical/actual stages, feed location, column diameter, height, tray pressure drop, flooding check |
| Flash drum | `flash_drum_device._full_design()` | Vapor fraction, vapor/liquid flows and compositions, drum dimensions, residence time |
| Shell-and-tube heat exchanger | `heatexchanger_device._full_he_design()` | Heat duty, LMTD, U, heat-transfer area, tube-/shell-side pressure drop |
| Reactor | `reactor_device._full_reactor_design()` | Reaction extent, outlet composition, reactor volume, heat-transfer area (supports side reactions and multi-stage reactions) |
| Pump | `pump_device._full_design()` | Head, hydraulic/shaft/motor power, efficiency, NPSH |
| Storage tank | `storage_tank_device._full_design()` | Tank dimensions, aspect ratio, safe fill height, evaporation loss |

### Process Units (process_engine, chainable)

`run_pump`, `run_pressure_reducing_valve`, `run_heat_exchanger`, `run_stoichiometric_reactor`, `run_cstr_reactor`, `run_pfr_reactor`, `run_distillation_column`, `run_flash_drum`, `run_mixer`, `run_flowsheet`.

The outlet dictionaries of upstream units (`outlet_molar_flows_mol_per_s`, `outlet_pressure_Pa`, etc.) can be fed directly as inputs to downstream units, enabling chained multi-unit process orchestration.

### Distillation Design Method (FUG)

- **Fenske**: minimum number of theoretical stages (including distribution of multicomponent non-key components)
- **Underwood**: minimum reflux ratio (including calculation of the feed thermal-condition parameter q)
- **Gilliland (Molokanov continuous correlation)**: theoretical stages at the operating reflux ratio
- **Kirkbride**: feed-stage location
- **Column diameter**: Souders-Brown / Fair (1961) capacity-factor correlation (using the flow parameter F_LV and tray spacing as independent variables)
- **Tray hydraulics**: sieve-tray hole velocity, dry/liquid pressure drop, flooding factor, weeping check

---

## 5. Installation

### Requirements

- Python 3.10+ (developed on Python 3.13 / Windows)
- Required libraries:

```bash
pip install thermo CoolProp numpy scipy fluids chemicals
pip install langgraph langchain-openai langchain-core
```

> `procee_engine/requirements.txt` additionally lists verified minimum versions (thermo≥0.2.27, scipy≥1.10, numpy≥1.23, fluids≥1.0.25, chemicals≥1.1.5).

### LLM Configuration (Environment Variables)

| Variable | Description |
|---|---|
| `LLM_API_KEY` | LLM API key (**must be set before use**) |
| `LLM_MODEL` | Model name |
| `LLM_BASE_URL` | Custom API endpoint (OpenAI-compatible interface) |

Windows PowerShell example:

```powershell
$env:LLM_API_KEY="your-api-key"
$env:LLM_MODEL="gpt-4o"
$env:LLM_BASE_URL="https://your-endpoint/v1"
```

---

## 6. Quick Start

### 1. Command-Line Interactive Mode

```bash
python main.py
```

```
>>> Design a benzene–toluene distillation column, feed 100 kmol/h, molar composition 40% benzene / 60% toluene, atmospheric pressure, saturated liquid feed
```

### 2. Single Task

```bash
python main.py --task "Select a centrifugal pump for 25°C clean water, flow rate 36 m3/h, head 30 m"
```

### 3. Run Tests

```bash
python main.py --test          # Full suite (9 scenarios)
python main.py --test pump     # Single test: imports/registry/pump/dist/hex/reactor/flash/tank/e2e
```

### 4. MCP / HTTP Server

```bash
python main.py --stdio         # MCP stdio (JSON-RPC, line-based)
python main.py --http 8080     # HTTP REST API
```

HTTP endpoints:

| Path | Function |
|---|---|
| `POST /task` | Natural-language task: `{"request": "..."}` → `{"report": "..."}` |
| `POST /mcp` | Standard MCP JSON-RPC request |
| `POST /tools` | List all registered tools |

### 5. Direct Invocation in Python

**Property queries:**

```python
from physics_engine import calc_liquid_density, calc_vapor_pressure

rho = calc_liquid_density(name="water", T=373.15, P=101325.0)   # kg/m³
psat = calc_vapor_pressure(name="ethanol", T=351.52)            # Pa
```

**Equipment design pipeline:**

```python
from device_tools.distillation_device import _full_design

result = _full_design(
    components=["benzene", "toluene"],
    z=[0.40, 0.60],
    F=100.0,            # mol/s
    T_feed=298.15,      # K
    P=101325.0,         # Pa
    x_D=0.95, x_B=0.05,
)
```

**Chained process units:**

```python
from process_tools import run_pump, run_heat_exchanger

p1 = run_pump(inlet_pressure_Pa=101325, target_pressure_Pa=2_800_000)
hx = run_heat_exchanger(
    process_fluid_molar_flows_mol_per_s={"benzene": 10.0},
    process_fluid_temp_in_K=298.15,
    process_fluid_temp_out_K=350.0,
    process_fluid_pressure_Pa=p1["outlet_pressure_Pa"],
    utility_fluid_temp_in_K=523.15,
    utility_fluid_temp_out_K=503.15,
)
```

Debug mode: set `FUNCS_VERBOSE=1` to print intermediate values from each `calc_*` function.

---

## 7. Engineering Conventions

The following unified conventions are built into the system and should be maintained when using or extending it:

- **SI units** throughout; `calc_*` tool functions have a single responsibility, must not call each other, and return dictionaries with unit suffixes
- The **net area for tray hydraulics is 90% of the total column cross-sectional area** (A_net = 0.90·A_col)
- Actual tray count: **N_a = ⌈(N_theoretical − 1) / E₀⌉** (the reboiler stage is deducted before applying tray efficiency)
- Total column pressure drop: N_a × per-tray pressure drop; bottom pressure = top pressure + total column pressure drop
- Bubble-point calculations for strongly non-ideal systems such as ethanol–water use the **UNIFAC** activity-coefficient model, falling back to Raoult's law on failure
- Column diameter is determined by the governing section (rectifying or stripping) with the larger calculated diameter, then rounded to the standard diameter series
- Single-pass tray design is assumed for D ≤ 2 m

---

## 8. Tests and Examples

The [test/](test/) directory provides 35+ directly runnable scripts covering design examples for all equipment types and "device-level vs. physics-engine" result comparisons, for example:

- `example_distillation_benzene_toluene.py`, `example_distillation_ethylene_ethane.py`
- `example_flash_methanol_water.py`, `example_heatexchanger_design.py`
- `example_pump_npsh.py`, `example_storage_tank_150m3_ethanol.py`
- `compare_flash_device_vs_engine.py`, `compare_pump_device_vs_engine.py`

The repository root also contains the VCM (vinyl chloride monomer) full-process consistency test and VLE benchmarking scripts, such as `test_vcm_crosslayer_consistency.py` and `benchmark_vle.py`.

---


