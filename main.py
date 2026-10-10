"""多智能体系统 — 主入口

基于 LangGraph + LangChain 的通用多智能体系统。
所有测试均通过自然语言请求驱动完整工作流，由 LLM 完成参数理解与提取。

使用方式:
    # 1. 配置 LLM（通过环境变量）
    set LLM_API_KEY=your-api-key
    set LLM_MODEL=gpt-4o               # 可选，默认 gpt-4o
    set LLM_BASE_URL=https://...       # 可选，自定义 API 地址

    # 2. 交互模式
    python main.py

    # 3. 单次任务
    python main.py --task "选型一台流量 36m³/h、扬程 30m 的水泵"

    # 4. 运行所有测试
    python main.py --test

    # 5. 运行单个测试
    python main.py --test pump

    # 6. MCP Server
    python main.py --stdio
    python main.py --http 8080
"""
from __future__ import annotations
import os
import sys
import time

# 确保直接运行时能找到项目根目录
_PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PARENT_DIR not in sys.path:
    sys.path.insert(0, _PARENT_DIR)


# ============================================================
# LLM 配置
# ============================================================

def setup_llm():
    """从环境变量配置 LLM"""
    from agents.planner import set_llm
    set_llm(
   model=os.environ.get("LLM_MODEL", "qwen3.5"),
      base_url=os.environ.get("LLM_BASE_URL", "xxx"),
      api_key=os.environ.get("LLM_API_KEY", "xxx"),
       )


# ============================================================
# 注册中心构建
# ============================================================

def build_registry():
    """创建并注册所有工具（物理计算 + 设备级 + process_engine 流程设备）"""
    from tools.registry import UniversalToolRegistry
    from tools import register_all_tools_to_universal

    registry = UniversalToolRegistry()
    # 1. 注册原有工具（device_tools + physics_tools）
    register_all_tools_to_universal(registry)
    # 2. 注册 process_engine 流程设备工具（8 种可串联设备）
    from process_tools import register_process_tools
    register_process_tools(registry)
    return registry


# ============================================================
# 核心执行函数
# ============================================================

def run_task(request: str, registry) -> str:
    """通过多智能体工作流执行自然语言任务，返回最终报告"""
    from graph.workflow import run_task as _run_task
    return _run_task(request, registry, interactive=True)


def run_task_state(request: str, registry) -> dict:
    """执行任务并返回完整 state（用于测试断言，不阻塞人工交互）"""
    from graph.workflow import create_workflow
    workflow = create_workflow(registry, interactive=False)
    return workflow(request)


# ============================================================
# 交互模式
# ============================================================

def run_interactive(registry):
    """命令行交互模式"""
    setup_llm()

    print("=" * 70)
    print("  通用多智能体系统（LangGraph + LLM）")
    print(f"  已注册工具: {registry.tool_count} 个  |  类别: {registry.categories}")
    print("=" * 70)
    print("\n输入自然语言描述你的任务，智能体将自动理解并执行。")
    print("命令: ls=列出工具  test=运行测试  quit=退出\n")

    while True:
        try:
            user_input = input(">>> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n再见!")
            break

        if not user_input:
            continue
        if user_input.lower() in ("quit", "exit", "q"):
            print("再见!")
            break
        if user_input.lower() in ("ls", "list", "tools"):
            for t in registry.list_tools():
                print(f"  [{t.category or '?'}] {t.name}: {t.description[:60]}...")
            continue
        if user_input.lower() == "test":
            run_all_tests(registry)
            continue

        t0 = time.time()
        report = run_task(user_input, registry)
        print(report)
        print(f"\n（耗时 {time.time() - t0:.1f}s）\n")


# ============================================================
# 各场景测试函数（全部通过自然语言驱动工作流）
# ============================================================

def _assert_success(state: dict, name: str) -> bool:
    """检查是否至少有一个工具执行成功"""
    results = state.get("execution_results", {})
    ok = any(isinstance(v, dict) and v.get("success") for v in results.values())
    if ok:
        print(f"  ✓ [{name}] 通过")
    else:
        print(f"  ✗ [{name}] 无工具执行成功。errors={state.get('errors', [])}")
    return ok


def test_imports() -> bool:
    """验证核心依赖是否可以正常导入"""
    print("\n--- 测试: 核心依赖导入 ---")
    results = {}

    for module in ("langgraph.graph", "langchain_openai", "langchain_core.messages"):
        try:
            __import__(module)
            results[module] = True
            print(f"  ✓ {module}")
        except ImportError as e:
            results[module] = False
            print(f"  ✗ {module}: {e}")

    ok = all(results.values())
    if ok:
        print("  ✓ [依赖导入] 全部通过")
    return ok


def test_registry(registry) -> bool:
    """验证工具注册中心"""
    print("\n--- 测试: 工具注册中心 ---")
    count = registry.tool_count
    print(f"  注册工具数: {count}，类别: {registry.categories}")
    ok = count > 0
    if ok:
        print(f"  ✓ [工具注册中心] 通过")
    else:
        print(f"  ✗ [工具注册中心] 无工具注册")
    return ok


def test_pump(registry) -> bool:
    """泵选型 — LLM 从自然语言中提取流量、扬程、介质等"""
    print("\n--- 测试: 泵选型 ---")
    request = (
        "请帮我选型一台离心泵，输送 25°C 清水，"
        "设计流量 36 m³/h，所需扬程 30 m。"
    )
    print(f"  请求: {request}")
    state = run_task_state(request, registry)
    print(f"  理解: {state.get('plan', {}).get('understanding', '—')}")
    return _assert_success(state, "泵选型")


def test_distillation(registry) -> bool:
    """精馏塔设计 — LLM 解析组成、进料量、操作条件"""
    print("\n--- 测试: 精馏塔设计 ---")
    request = (
        "设计苯-甲苯精馏塔，进料 100 kmol/h，"
        "摩尔组成苯 40%、甲苯 60%，常压操作，饱和液体进料。"
    )
    print(f"  请求: {request}")
    state = run_task_state(request, registry)
    print(f"  理解: {state.get('plan', {}).get('understanding', '—')}")
    return _assert_success(state, "精馏塔设计")


def test_heat_exchanger(registry) -> bool:
    """换热器设计 — LLM 解析冷热侧参数"""
    print("\n--- 测试: 换热器设计 ---")
    request = (
        "设计一台管壳式换热器，热侧热水从 400 K 冷却至 350 K，流量 1 kg/s；"
        "冷侧冷水从 300 K 加热至 340 K，流量 1 kg/s。"
    )
    print(f"  请求: {request}")
    state = run_task_state(request, registry)
    print(f"  理解: {state.get('plan', {}).get('understanding', '—')}")
    return _assert_success(state, "换热器设计")


def test_reactor(registry) -> bool:
    """反应器设计 — LLM 解析反应动力学和操作参数"""
    print("\n--- 测试: 反应器设计 ---")
    request = (
        "设计全混流反应器，甲醇和乙酸0.05mol/s混合进料，甲醇占摩尔比0.4，二级不可逆反应，反应常数为0.1m³/(mol·s)，目标转化率 90%，温度 200 ℃ ，压力 1 atm"
    )
    print(f"  请求: {request}")
    state = run_task_state(request, registry)
    print(f"  理解: {state.get('plan', {}).get('understanding', '—')}")
    return _assert_success(state, "反应器设计")


def test_flash_drum(registry) -> bool:
    """闪蒸罐设计 — LLM 解析多组分体系和操作条件"""
    print("\n--- 测试: 闪蒸罐设计 ---")
    request = (
        "设计甲醇-水混合物闪蒸罐，"
        "进料量 100 kmol/h，摩尔组成甲醇 40%、水 60%，"
        "操作温度 350 K，常压。"
    )
    print(f"  请求: {request}")
    state = run_task_state(request, registry)
    print(f"  理解: {state.get('plan', {}).get('understanding', '—')}")
    return _assert_success(state, "闪蒸罐设计")


def test_storage_tank(registry) -> bool:
    """储罐设计 — LLM 解析容积和罐型"""
    print("\n--- 测试: 储罐设计 ---")
    request = "设计一座 500 m³ 的立式储罐，用于常温储存甲醇。"
    print(f"  请求: {request}")
    state = run_task_state(request, registry)
    print(f"  理解: {state.get('plan', {}).get('understanding', '—')}")
    return _assert_success(state, "储罐设计")


def test_workflow_integrity(registry) -> bool:
    """验证 Planner→Executor→Reporter 三阶段均正确参与"""
    print("\n--- 测试: 工作流完整性 ---")
    request = "选型一台输送清水的泵，流量 20 m³/h，扬程 20 m。"
    state = run_task_state(request, registry)

    messages = state.get("messages", [])
    roles = {m.get("role") for m in messages}

    checks = {
        "planner 参与":   "planner"  in roles,
        "executor 参与":  "executor" in roles,
        "reporter 参与":  "reporter" in roles,
        "有最终报告":     len(state.get("final_report", "")) > 50,
        "有执行结果":     bool(state.get("execution_results")),
    }

    all_pass = True
    for name, passed in checks.items():
        print(f"  {'✓' if passed else '✗'} {name}")
        if not passed:
            all_pass = False

    if all_pass:
        print("  ✓ [工作流完整性] 通过")
    return all_pass


# ============================================================
# 测试套件
# ============================================================

# 注册所有可选测试（key 用于命令行指定单测）
_TEST_SUITE = [
    ("imports",     "核心依赖导入",       lambda r: test_imports()),
    ("registry",    "工具注册中心",       test_registry),
    ("pump",        "泵选型",             test_pump),
    ("dist",        "精馏塔设计",         test_distillation),
    ("hex",         "换热器设计",         test_heat_exchanger),
    ("reactor",     "反应器设计",         test_reactor),
    ("flash",       "闪蒸罐设计",         test_flash_drum),
    ("tank",        "储罐设计",           test_storage_tank),
    ("e2e",         "工作流完整性",       test_workflow_integrity),
]


def run_all_tests(registry=None, filter_key: str = None):
    """运行全量（或指定）测试"""
    if registry is None:
        registry = build_registry()

    print("\n" + "=" * 70)
    print("  多智能体系统全量测试")
    print("  所有测试均通过自然语言请求驱动，LLM 负责理解和参数提取")
    print("=" * 70)

    suite = _TEST_SUITE
    if filter_key:
        suite = [(k, n, f) for k, n, f in _TEST_SUITE if k == filter_key]
        if not suite:
            print(f"  未找到测试 '{filter_key}'，可选: {[k for k,_,_ in _TEST_SUITE]}")
            return False

    passed, failed = 0, []
    t_total = time.time()

    for key, name, fn in suite:
        t0 = time.time()
        try:
            ok = fn(registry)
            elapsed = time.time() - t0
            if ok is not False:
                passed += 1
            else:
                failed.append(name)
            print(f"  耗时: {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - t0
            print(f"\n  ✗ [{name}] 异常: {e}")
            import traceback
            traceback.print_exc()
            failed.append(name)
            print(f"  耗时: {elapsed:.1f}s（异常）")

    total_elapsed = time.time() - t_total
    print("\n" + "=" * 70)
    print(f"  结果: {passed}/{len(suite)} 通过  |  总耗时: {total_elapsed:.1f}s")
    if failed:
        print(f"  失败: {failed}")
    else:
        print("  全部通过 ✓")
    print("=" * 70)

    return len(failed) == 0


# ============================================================
# 主入口
# ============================================================

def main():
    # 修复 Windows 控制台编码
    if sys.stdout.encoding and sys.stdout.encoding.lower().replace("-", "") != "utf8":
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except Exception:
            pass

    args = sys.argv[1:]
    setup_llm()

    if "--test" in args:
        registry = build_registry()
        # 支持指定单测：--test pump
        idx = args.index("--test")
        filter_key = args[idx + 1] if idx + 1 < len(args) and not args[idx + 1].startswith("--") else None
        success = run_all_tests(registry, filter_key=filter_key)
        sys.exit(0 if success else 1)

    elif "--task" in args:
        idx = args.index("--task")
        request = args[idx + 1] if idx + 1 < len(args) else ""
        if not request:
            print("用法: python main.py --task \"你的任务描述\"")
            sys.exit(1)
        registry = build_registry()
        print(f"\n执行任务: {request}\n")
        report = run_task(request, registry)
        print(report)

    elif "--stdio" in args:
        from server import run_stdio_server
        run_stdio_server()

    elif "--http" in args:
        from server import run_http_server
        idx = args.index("--http")
        port = int(args[idx + 1]) if idx + 1 < len(args) else 8080
        run_http_server(port)

    else:
        registry = build_registry()
        run_interactive(registry)


if __name__ == "__main__":
    main()
