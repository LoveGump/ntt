#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 设置中文字体支持
try:
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体为黑体
    plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像时负号'-'显示为方块的问题
except:
    print("警告: 无法设置中文字体，图表中的中文可能无法正确显示")

# 读取数据
metrics_df = pd.read_csv('perf_results/calculated_metrics.csv')

# 读取原始性能数据以获取运行时间
try:
    perf_df = pd.read_csv('perf_results/performance_stats_pivot.csv')
    has_time_data = True
except:
    print("警告: 找不到原始性能数据文件，将不包括运行时间分析")
    has_time_data = False

# 计算每个程序的平均指标
avg_metrics = metrics_df.groupby('Program').mean().reset_index()

# 创建结果目录
import os
os.makedirs('perf_results/visualizations', exist_ok=True)

#---------- 1. CPI分析 ----------#
plt.figure(figsize=(12, 6))
ax = plt.subplot(111)

# 设置x位置
x = np.arange(len(avg_metrics['Program']))
width = 0.35

# 使用双柱状图比较core和atom的CPI
core_bars = ax.bar(x - width/2, avg_metrics['core_CPI'], width, label='Core CPI', color='royalblue')
atom_bars = ax.bar(x + width/2, avg_metrics['atom_CPI'], width, label='Atom CPI', color='tomato')

# 添加文本标签
ax.set_xlabel('程序实现')
ax.set_ylabel('CPI值 (每指令周期数)')
ax.set_title('不同NTT实现的Core和Atom CPI比较')
ax.set_xticks(x)
ax.set_xticklabels(avg_metrics['Program'])
ax.legend()

# 在柱状图上添加数值标签
def add_labels(bars):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3点垂直偏移
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=9)

add_labels(core_bars)
add_labels(atom_bars)

plt.tight_layout()
plt.savefig('perf_results/visualizations/cpi_comparison.png', dpi=300)
plt.close()

#---------- 2. 缓存命中率分析 ----------#
plt.figure(figsize=(15, 10))

# L1缓存命中率 (仅core)
plt.subplot(2, 2, 1)
sns.barplot(x='Program', y='core_L1_hit_rate', data=avg_metrics, palette='Blues_d')
plt.title('Core L1缓存命中率')
plt.xlabel('程序实现')
plt.ylabel('命中率')
plt.ylim(0.98, 1.0)  # 调整Y轴以更好地显示差异

# 添加数值标签
for i, v in enumerate(avg_metrics['core_L1_hit_rate']):
    plt.text(i, v, f'{v:.4f}', ha='center', va='bottom', fontsize=9)

# L2缓存命中率 (仅core)
plt.subplot(2, 2, 2)
sns.barplot(x='Program', y='core_L2_hit_rate', data=avg_metrics, palette='Greens_d')
plt.title('Core L2缓存命中率')
plt.xlabel('程序实现')
plt.ylabel('命中率')
plt.ylim(0.65, 0.71)  # 调整Y轴以更好地显示差异

# 添加数值标签
for i, v in enumerate(avg_metrics['core_L2_hit_rate']):
    plt.text(i, v, f'{v:.4f}', ha='center', va='bottom', fontsize=9)

# L3缓存命中率比较
plt.subplot(2, 2, 3)
x = np.arange(len(avg_metrics['Program']))
width = 0.35

ax = plt.gca()
core_bars = ax.bar(x - width/2, avg_metrics['core_L3_hit_rate'], width, label='Core', color='purple')
atom_bars = ax.bar(x + width/2, avg_metrics['atom_L3_hit_rate'], width, label='Atom', color='orchid')

plt.title('L3 (LLC)缓存命中率比较')
plt.xlabel('程序实现')
plt.ylabel('命中率')
plt.xticks(x, avg_metrics['Program'])
plt.legend()
plt.ylim(0.7, 1.0)  # 调整Y轴以更好地显示差异

# 添加数值标签
def add_labels(bars):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.4f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=8)

add_labels(core_bars)
add_labels(atom_bars)

# 所有缓存层级的命中率汇总柱状图
plt.subplot(2, 2, 4)
df_melt = pd.melt(avg_metrics, 
                  id_vars=['Program'], 
                  value_vars=['core_L1_hit_rate', 'core_L2_hit_rate', 'core_L3_hit_rate'],
                  var_name='Cache Level', 
                  value_name='Hit Rate')

# 替换名称使其更易读
df_melt['Cache Level'] = df_melt['Cache Level'].replace({
    'core_L1_hit_rate': 'L1 缓存',
    'core_L2_hit_rate': 'L2 缓存',
    'core_L3_hit_rate': 'L3 缓存'
})

sns.barplot(x='Program', y='Hit Rate', hue='Cache Level', data=df_melt, palette='YlOrRd')
plt.title('Core各级缓存命中率对比')
plt.xlabel('程序实现')
plt.ylabel('命中率')
plt.ylim(0.6, 1.01)  # 调整Y轴以适应所有缓存层级的值

plt.tight_layout()
plt.savefig('perf_results/visualizations/cache_hit_rates.png', dpi=300)
plt.close()

#---------- 3. 箱线图分析变异性 ----------#
plt.figure(figsize=(16, 10))

# CPI箱线图
plt.subplot(2, 2, 1)
df_melt = pd.melt(metrics_df, 
                  id_vars=['Program'], 
                  value_vars=['core_CPI', 'atom_CPI'],
                  var_name='CPU Type', 
                  value_name='CPI')
df_melt['CPU Type'] = df_melt['CPU Type'].replace({'core_CPI': 'Core', 'atom_CPI': 'Atom'})
sns.boxplot(x='Program', y='CPI', hue='CPU Type', data=df_melt, palette='Set2')
plt.title('CPI分布箱线图')
plt.xlabel('程序实现')
plt.ylabel('CPI值')

# Core缓存命中率箱线图
plt.subplot(2, 2, 2)
df_melt = pd.melt(metrics_df, 
                  id_vars=['Program'], 
                  value_vars=['core_L1_hit_rate', 'core_L2_hit_rate', 'core_L3_hit_rate'],
                  var_name='Cache Level', 
                  value_name='Hit Rate')
df_melt['Cache Level'] = df_melt['Cache Level'].replace({
    'core_L1_hit_rate': 'L1 Cache', 
    'core_L2_hit_rate': 'L2 Cache', 
    'core_L3_hit_rate': 'L3 Cache'
})
sns.boxplot(x='Program', y='Hit Rate', hue='Cache Level', data=df_melt, palette='Set3')
plt.title('Core各级缓存命中率分布')
plt.xlabel('程序实现')
plt.ylabel('命中率')

# 如果有运行时间数据，添加运行时间分析
if has_time_data:
    plt.subplot(2, 2, 3)
    # 提取时间相关列
    time_cols = [col for col in perf_df.columns if 'time' in col]
    time_df = perf_df[['Program', 'Run'] + time_cols]
    
    time_melt = pd.melt(time_df, 
                         id_vars=['Program', 'Run'], 
                         value_vars=time_cols,
                         var_name='Time Type', 
                         value_name='Seconds')
    
    # 绘制运行时间箱线图
    sns.boxplot(x='Program', y='Seconds', hue='Time Type', data=time_melt, palette='Paired')
    plt.title('运行时间分布')
    plt.xlabel('程序实现')
    plt.ylabel('时间 (秒)')
    plt.legend(title='时间类型')

# 处理器利用率对比 - IPC (Instruction Per Cycle) - CPI的倒数
plt.subplot(2, 2, 4)
avg_metrics['core_IPC'] = 1 / avg_metrics['core_CPI']
avg_metrics['atom_IPC'] = 1 / avg_metrics['atom_CPI']

x = np.arange(len(avg_metrics['Program']))
width = 0.35

ax = plt.gca()
core_bars = ax.bar(x - width/2, avg_metrics['core_IPC'], width, label='Core', color='teal')
atom_bars = ax.bar(x + width/2, avg_metrics['atom_IPC'], width, label='Atom', color='coral')

plt.title('处理器利用率 (IPC)')
plt.xlabel('程序实现')
plt.ylabel('每周期指令数 (IPC)')
plt.xticks(x, avg_metrics['Program'])
plt.legend()

# 添加数值标签
for i, v in enumerate(avg_metrics['core_IPC']):
    plt.text(i - width/2, v, f'{v:.2f}', ha='center', va='bottom', fontsize=9)
for i, v in enumerate(avg_metrics['atom_IPC']):
    plt.text(i + width/2, v, f'{v:.2f}', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig('perf_results/visualizations/performance_distribution.png', dpi=300)
plt.close()

#---------- 4. 热图分析不同程序间的指标相关性 ----------#
# 根据程序准备数据
programs = metrics_df['Program'].unique()
correlation_dfs = {}

for program in programs:
    prog_data = metrics_df[metrics_df['Program'] == program].drop(['Program', 'Run'], axis=1)
    correlation_dfs[program] = prog_data.corr()

# 创建热图
plt.figure(figsize=(18, 6))

for i, program in enumerate(programs):
    plt.subplot(1, 3, i+1)
    sns.heatmap(correlation_dfs[program], annot=True, cmap='coolwarm', fmt='.2f', linewidths=.5, vmin=-1, vmax=1)
    plt.title(f'{program}实现的指标相关性')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig('perf_results/visualizations/correlation_heatmap.png', dpi=300)
plt.close()

#---------- 5. 雷达图比较综合性能 ----------#
# 准备雷达图数据 - 需要对数据进行标准化以便在同一个尺度上比较
from sklearn.preprocessing import MinMaxScaler

# 选择要在雷达图上显示的指标
radar_metrics = ['core_CPI', 'core_L1_hit_rate', 'core_L2_hit_rate', 'core_L3_hit_rate']

# 对于CPI，我们需要取倒数，因为较低的CPI表示更好的性能
avg_metrics_radar = avg_metrics.copy()
avg_metrics_radar['core_CPI'] = 1 / avg_metrics_radar['core_CPI']  # 转换为IPC

# 对选择的指标进行标准化，使值在0到1之间
scaler = MinMaxScaler()
avg_metrics_radar[radar_metrics] = scaler.fit_transform(avg_metrics_radar[radar_metrics])

# 设置雷达图
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, polar=True)

# 角度设置
angles = np.linspace(0, 2*np.pi, len(radar_metrics), endpoint=False).tolist()
angles += angles[:1]  # 闭合图形

# 添加指标标签
labels = ['IPC (更高更好)', 'L1命中率', 'L2命中率', 'L3命中率']
plt.xticks(angles[:-1], labels)

# 绘制每个程序的雷达图
for i, program in enumerate(avg_metrics_radar['Program']):
    values = avg_metrics_radar.loc[i, radar_metrics].values.tolist()
    values += values[:1]  # 闭合数据
    ax.plot(angles, values, linewidth=2, linestyle='solid', label=program)
    ax.fill(angles, values, alpha=0.1)

plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
plt.title('NTT实现的综合性能对比', size=15)
plt.tight_layout()
plt.savefig('perf_results/visualizations/radar_chart.png', dpi=300)
plt.close()

print("可视化分析完成！所有图表已保存在 perf_results/visualizations/ 目录下")