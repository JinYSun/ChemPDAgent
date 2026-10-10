"""MCP Server 入口

启动方式:
    python server.py              # 命令行交互模式
    python server.py --stdio      # MCP stdio 协议
    python server.py --http 8080  # HTTP REST API
"""
from __future__ import annotations
import sys
import json

from tools.base import ToolRegistry, handle_mcp_request
from tools import register_all_physics_tools
from device_tools import register_all_device_tools
from graph.workflow import run_task


def build_registry() -> ToolRegistry:
    """构建并返回工具注册中心"""
    registry = ToolRegistry("multi-agent-system")
    register_all_physics_tools(registry)
    register_all_device_tools(registry)
    return registry


registry = build_registry()


# ============================================================
# Server 运行模式
# ============================================================

def run_stdio_server():
    """MCP stdio 协议模式"""
    print(f"[MCP Server] 已就绪，工具数: {registry.tool_count}", file=sys.stderr)
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            resp = handle_mcp_request(json.loads(line), registry)
        except json.JSONDecodeError as e:
            resp = {"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": str(e)}}
        print(json.dumps(resp, ensure_ascii=False))
        sys.stdout.flush()


def run_http_server(port: int = 8080):
    """HTTP REST API 模式"""
    from http.server import HTTPServer, BaseHTTPRequestHandler

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):
            length = int(self.headers.get("Content-Length", 0))
            body = self.rfile.read(length).decode()
            try:
                req = json.loads(body)

                if self.path == "/mcp":
                    # 标准 MCP JSON-RPC 请求
                    resp = handle_mcp_request(req, registry)

                elif self.path == "/task":
                    # 自然语言任务接口：{"request": "..."}
                    request_text = req.get("request", "")
                    report = run_task(request_text, registry)
                    resp = {
                        "jsonrpc": "2.0",
                        "id": req.get("id"),
                        "result": {"report": report},
                    }

                elif self.path == "/tools":
                    # 列出所有已注册工具
                    resp = {
                        "jsonrpc": "2.0",
                        "id": req.get("id"),
                        "result": {"tools": registry.list_tools()},
                    }

                else:
                    resp = {
                        "jsonrpc": "2.0",
                        "id": req.get("id"),
                        "error": {"code": -32601, "message": f"未知路径: {self.path}"},
                    }

                self.send_response(200)
                self.send_header("Content-Type", "application/json; charset=utf-8")
                self.send_header("Access-Control-Allow-Origin", "*")
                self.end_headers()
                self.wfile.write(json.dumps(resp, ensure_ascii=False).encode())

            except Exception as e:
                self.send_response(500)
                self.end_headers()
                self.wfile.write(json.dumps({"error": str(e)}).encode())

        def do_OPTIONS(self):
            self.send_response(200)
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
            self.send_header("Access-Control-Allow-Headers", "Content-Type")
            self.end_headers()

        def log_message(self, fmt, *args):
            # 只打印错误，减少日志噪音
            if args and str(args[1]) not in ("200", "204"):
                super().log_message(fmt, *args)

    server = HTTPServer(("0.0.0.0", port), Handler)
    print(f"[MCP Server] HTTP 服务已启动: http://localhost:{port}", file=sys.stderr)
    print(f"  接口: POST /task  POST /mcp  POST /tools", file=sys.stderr)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n[MCP Server] 已停止", file=sys.stderr)


def run_interactive():
    """命令行交互模式"""
    print("=" * 70)
    print("  多智能体系统 — 交互模式")
    print(f"  注册工具: {registry.tool_count} 个  |  类别: {registry.categories}")
    print("=" * 70)
    print("\n输入自然语言描述你的任务，智能体将自动理解并执行。")
    print("命令: ls=列出工具  quit=退出\n")

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
                print(f"  [{t.get('category', '?')}] {t.name}: {t.description[:60]}...")
            continue

        print("\n正在执行多智能体工作流...\n")
        report = run_task(user_input, registry)
        print(report)
        print()


# ============================================================
# 主入口
# ============================================================

def main():
    args = sys.argv[1:]
    if "--stdio" in args:
        run_stdio_server()
    elif "--http" in args:
        idx = args.index("--http")
        port = int(args[idx + 1]) if idx + 1 < len(args) else 8080
        run_http_server(port)
    else:
        run_interactive()


if __name__ == "__main__":
    main()