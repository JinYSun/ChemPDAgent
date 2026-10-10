# -*- coding: utf-8 -*-
"""检查物料衡算"""
import sys, io, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.distillation_device import _full_design

result = _full_design(
    components=["benzene", "toluene"],
    z=[0.4, 0.6],
    F=100.0,
    x_D_spec=0.99,
    x_B_spec=0.01,
    P=101325.0,
    q=1.0,
    RR_factor=2.0,
)

print("=== 物料衡算检查 ===\n")

mb = result.get("material_balance", {})
print(f"进料 F = {mb.get('F_kmolh')} kmol/h")
print(f"塔顶 D = {mb.get('D_kmolh')} kmol/h")
print(f"塔底 B = {mb.get('B_kmolh')} kmol/h")
print(f"x_D_LK = {mb.get('x_D_LK')}")
print(f"x_B_LK = {mb.get('x_B_LK')}")

print("\n组分分配:")
for c in mb.get("component_distribution", []):
    print(f"  {c['component']}: 进料={c['feed_kmolh']}, 塔顶={c['distillate_kmolh']}, 塔底={c['bottoms_kmolh']}")

# 验证
F = mb.get("F_kmolh", 0)
D = mb.get("D_kmolh", 0)
B = mb.get("B_kmolh", 0)
print(f"\n--- 验证 ---")
print(f"总物料平衡: F = D + B => {F} = {D} + {B} = {D + B:.2f} {'OK' if abs(F - D - B) < 0.1 else 'ERROR'}")

comp = mb.get("component_distribution", [])
if len(comp) >= 2:
    benzene_feed = comp[0]["feed_kmolh"]
    benzene_dist = comp[0]["distillate_kmolh"]
    benzene_bot = comp[0]["bottoms_kmolh"]
    print(f"苯平衡: {benzene_feed} = {benzene_dist} + {benzene_bot} = {benzene_dist + benzene_bot:.2f} {'OK' if abs(benzene_feed - benzene_dist - benzene_bot) < 0.1 else 'ERROR'}")
    
    # 计算实际组成
    x_D_actual = benzene_dist / D if D > 0 else 0
    x_B_actual = benzene_bot / B if B > 0 else 0
    print(f"\n实际组成:")
    print(f"  x_D (苯) = {x_D_actual:.4f} (要求 >= 0.99)")
    print(f"  x_B (苯) = {x_B_actual:.4f} (要求 <= 0.01)")
