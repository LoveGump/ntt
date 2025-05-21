#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# 读取CSV文件
df = pd.read_csv('/home/gump/桌面/并行/ntt/ntt/performance_stats_pivot.csv')

# 创建一个新的DataFrame来存储计算的指标
metrics_df = pd.DataFrame()
metrics_df['Program'] = df['Program']
metrics_df['Run'] = df['Run']

# 计算CPI (Cycles Per Instruction)
# CPI = cycles / instructions
metrics_df['core_CPI'] = df['cpu_core_cycles'] / df['cpu_core_instructions']
metrics_df['atom_CPI'] = df['cpu_atom_cycles'] / df['cpu_atom_instructions']

# 计算L1缓存命中率
# L1 hit rate = 1 - (L1 misses / L1 loads)
# atom核心不支持L1-dcache-load-misses，所以只计算core核心
metrics_df['core_L1_hit_rate'] = 1 - (df['cpu_core_L1_dcache_load_misses'] / df['cpu_core_L1_dcache_loads'])

# 计算L2缓存命中率 
# L2 hit rate = 1 - (L2 misses / L2 references)
metrics_df['core_L2_hit_rate'] = 1 - (df['cpu_core_l2_rqsts_miss'] / df['cpu_core_l2_rqsts_references'])

# 计算L3 (LLC)缓存命中率
# LLC hit rate = 1 - (LLC misses / LLC loads)
metrics_df['core_L3_hit_rate'] = 1 - (df['cpu_core_LLC_load_misses'] / df['cpu_core_LLC_loads'])
metrics_df['atom_L3_hit_rate'] = 1 - (df['cpu_atom_LLC_load_misses'] / df['cpu_atom_LLC_loads'])

# 计算每个程序的平均值
avg_metrics = metrics_df.groupby('Program').mean().reset_index()

# 保存计算结果
metrics_df.to_csv('perf_results/calculated_metrics.csv', index=False)
avg_metrics.to_csv('perf_results/avg_calculated_metrics.csv', index=False)

print("计算完成并保存到 calculated_metrics.csv 和 avg_calculated_metrics.csv")

# 输出计算结果摘要
print("\n各程序平均指标：")
print(avg_metrics)