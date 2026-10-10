# -*- coding: utf-8 -*-
"""甲醇-水闪蒸罐设计示例"""
import sys, io, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, r'c:\Users\12002\.openclaw\workspace\chem-design-agents\equipment\mcp_agents')

from device_tools.flash_drum_device import _full_design

result = _full_design(
    components=["methanol", "water"],   # 甲醇-水体系
    z=[0.4, 0.6],                       # 进料摩尔组成：甲醇40%、水60%
    F=100.0,                            # 进料流量 100 kmol/h
    T=350.0,                            # 闪蒸温度 350 K (77°C)
    P=101325.0,                         # 常压操作 101.325 kPa
)

print("=" * 70)
print("  甲醇-水闪蒸罐设计结果")
print("=" * 70)

# 打印完整结果
print("\n--- 完整结果 ---")
print(json.dumps(result, indent=2, default=str, ensure_ascii=False))
