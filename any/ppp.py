import pandas as pd
import numpy as np

# 读取CSV文件
df = pd.read_csv('performance_stats_pivot.csv')

# 处理<not值，替换为NaN
df = df.replace('<not', np.nan)

# 将数值列转换为float类型
for col in df.columns[2:]:  # 跳过Program和Run列
    if col not in ['sys_time', 'time_elapsed', 'user_time']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

# 指定要排除的列
excluded_cols = ['Run', 'sys_time', 'time_elapsed', 'user_time', 'cpu_atom_L1_dcache_load_misses']

# 按Program分组并计算平均值，排除指定列
result = df.drop(columns=excluded_cols).groupby('Program').mean()

# 保存结果到CSV文件
result.to_csv('avg_performance_metrics.csv')

print("已按程序类型计算平均值，并排除时间相关列和cpu_atom_L1_dcache_load_misses列")
print(f"结果已保存到 avg_performance_metrics.csv")
print("\n数据预览:")
print(result.head())