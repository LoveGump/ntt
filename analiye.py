#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.font_manager as fm
import os

# 首先检查系统中可用的字体
fonts = [f.name for f in fm.fontManager.ttflist]
print("系统中的所有字体:")
for font in fonts:
    print(font)

# 使用跨平台兼容性更好的方式设置字体
# 方法一: 使用 fontproperties 直接指定字体文件
try:
    # 寻找系统中的字体文件
    font_paths = [
        # 在多个可能的位置查找中文字体
        '/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf',  # Ubuntu常见
        '/usr/share/fonts/wqy-microhei/wqy-microhei.ttc',            # WQY微米黑
        '/usr/share/fonts/wqy-zenhei/wqy-zenhei.ttc',                # WQY正黑
        '/usr/share/fonts/truetype/arphic/uming.ttc',                # AR PL UMing
        '/System/Library/Fonts/PingFang.ttc'                         # macOS
    ]
    
    font_file = None
    for path in font_paths:
        if os.path.exists(path):
            font_file = path
            print(f"找到可用字体文件: {path}")
            break
    
    if font_file is None:
        print("未找到预定义的字体文件，使用不显示字体或英文替代")
        # 使用英文标签作为备选
        use_english = True
    else:
        chinese_font = fm.FontProperties(fname=font_file)
        use_english = False
        
except Exception as e:
    print(f"设置字体时出错: {e}")
    use_english = True  # 出错时使用英文

# 读取数据
df = pd.read_csv('context-Switches.csv')

# 计算CPU利用率
df['CPU_Utilization'] = (df['User-Time(sec)'] + df['Sys-Time(sec)']) / df['Elapsed-Time(sec)']

# 定义标签文本（中英文双版本）
labels = {
    'title1': '线程调度开销对比 (对数刻度)' if not use_english else 'Thread Scheduling Overhead (Log Scale)',
    'xlabel1': '程序实现' if not use_english else 'Implementation',
    'ylabel1': '次数 (对数刻度)' if not use_english else 'Count (Log Scale)',
    'label1_1': '上下文切换' if not use_english else 'Context Switches',
    'label1_2': 'CPU迁移' if not use_english else 'CPU Migrations',
    'title2': '执行时间对比' if not use_english else 'Execution Time Comparison',
    'xlabel2': '程序实现' if not use_english else 'Implementation',
    'ylabel2': '时间 (秒)' if not use_english else 'Time (sec)',
    'label2_1': '实际时间' if not use_english else 'Elapsed Time',
    'label2_2': '用户时间' if not use_english else 'User Time',
    'label2_3': '系统时间' if not use_english else 'System Time',
    'title3': 'CPU利用率对比' if not use_english else 'CPU Utilization Comparison',
    'xlabel3': '程序实现' if not use_english else 'Implementation',
    'ylabel3': 'CPU利用率' if not use_english else 'CPU Utilization',
    'label3_1': '8核心理想利用率' if not use_english else 'Ideal 8-core Utilization',
    'title4': '并行效率关键指标热力图' if not use_english else 'Key Parallel Efficiency Metrics Heatmap',
    'label4_1': '归一化值 (高=更好)' if not use_english else 'Normalized Value (High = Better)'
}

# 创建可视化图表
plt.figure(figsize=(15, 15))

# 1. 线程调度开销（上下文切换和CPU迁移）- 使用对数刻度
ax1 = plt.subplot(3, 1, 1)
width = 0.35
x = np.arange(len(df['Program']))

# 使用对数刻度以便清晰显示数量级的差异
ax1.set_yscale('log')
bars1 = ax1.bar(x - width/2, df['Context-Switches'], width, label=labels['label1_1'], color='skyblue')
bars2 = ax1.bar(x + width/2, df['CPU-Migrations'], width, label=labels['label1_2'], color='salmon')

# 设置标签（使用fontproperties如果找到了字体文件）
if not use_english and 'chinese_font' in locals():
    ax1.set_xlabel(labels['xlabel1'], fontproperties=chinese_font)
    ax1.set_ylabel(labels['ylabel1'], fontproperties=chinese_font)
    ax1.set_title(labels['title1'], fontproperties=chinese_font)
    legend = ax1.legend(prop=chinese_font)
else:
    ax1.set_xlabel(labels['xlabel1'])
    ax1.set_ylabel(labels['ylabel1'])
    ax1.set_title(labels['title1'])
    ax1.legend()

ax1.set_xticks(x)
ax1.set_xticklabels(df['Program'])

# 在柱状图上添加数值标签
def add_labels(bars):
    for bar in bars:
        height = bar.get_height()
        if height > 0:  # 只为非零值添加标签
            ax1.annotate(f'{int(height)}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),  # 3点垂直偏移
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)

add_labels(bars1)
add_labels(bars2)

# 2. 执行时间对比
ax2 = plt.subplot(3, 1, 2)
width = 0.25
x = np.arange(len(df['Program']))

bars1 = ax2.bar(x - width, df['Elapsed-Time(sec)'], width, label=labels['label2_1'], color='cornflowerblue')
bars2 = ax2.bar(x, df['User-Time(sec)'], width, label=labels['label2_2'], color='lightcoral')
bars3 = ax2.bar(x + width, df['Sys-Time(sec)'], width, label=labels['label2_3'], color='mediumseagreen')

# 设置标签
if not use_english and 'chinese_font' in locals():
    ax2.set_xlabel(labels['xlabel2'], fontproperties=chinese_font)
    ax2.set_ylabel(labels['ylabel2'], fontproperties=chinese_font)
    ax2.set_title(labels['title2'], fontproperties=chinese_font)
    legend = ax2.legend(prop=chinese_font)
else:
    ax2.set_xlabel(labels['xlabel2'])
    ax2.set_ylabel(labels['ylabel2'])
    ax2.set_title(labels['title2'])
    ax2.legend()

ax2.set_xticks(x)
ax2.set_xticklabels(df['Program'])

# 添加数值标签
for i, v in enumerate(df['Elapsed-Time(sec)']):
    ax2.text(i - width, v + 0.02, f'{v:.3f}', ha='center', va='bottom', fontsize=9)
for i, v in enumerate(df['User-Time(sec)']):
    ax2.text(i, v + 0.02, f'{v:.3f}', ha='center', va='bottom', fontsize=9)
for i, v in enumerate(df['Sys-Time(sec)']):
    ax2.text(i + width, v + 0.02, f'{v:.3f}', ha='center', va='bottom', fontsize=9)

# 3. CPU利用率对比
ax3 = plt.subplot(3, 1, 3)
bars = ax3.bar(df['Program'], df['CPU_Utilization'], color=sns.color_palette("viridis", len(df)))

# 设置标签
if not use_english and 'chinese_font' in locals():
    ax3.set_xlabel(labels['xlabel3'], fontproperties=chinese_font)
    ax3.set_ylabel(labels['ylabel3'], fontproperties=chinese_font)
    ax3.set_title(labels['title3'], fontproperties=chinese_font)
else:
    ax3.set_xlabel(labels['xlabel3'])
    ax3.set_ylabel(labels['ylabel3'])
    ax3.set_title(labels['title3'])

# 添加数值标签
for i, v in enumerate(df['CPU_Utilization']):
    ax3.text(i, v + 0.1, f'{v:.2f}x', ha='center', va='bottom')
    
# 添加理想线 - 8核心系统的理想利用率为8
ax3.axhline(y=8, color='r', linestyle='--', alpha=0.7, label=labels['label3_1'])

if not use_english and 'chinese_font' in locals():
    legend = ax3.legend(prop=chinese_font)
else:
    ax3.legend()

plt.tight_layout()
plt.savefig('thread_scheduling_analysis.png', dpi=300)
plt.show()

# 创建一个组合的热力图，显示调度开销与性能的关系
plt.figure(figsize=(10, 8))

# 准备热力图数据
heatmap_data = df[['Context-Switches', 'CPU-Migrations', 'Elapsed-Time(sec)', 'CPU_Utilization']]
heatmap_data = heatmap_data.set_index(df['Program'])

# 对数据进行归一化，使其在0-1之间，便于比较
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
heatmap_normalized = pd.DataFrame(scaler.fit_transform(heatmap_data), 
                                 columns=heatmap_data.columns,
                                 index=heatmap_data.index)

# 为了让更低的时间对应更高的性能，反转elapsed time
heatmap_normalized['Elapsed-Time(sec)'] = 1 - heatmap_normalized['Elapsed-Time(sec)']
heatmap_normalized.rename(columns={'Elapsed-Time(sec)': 'Speed(normalized)'}, inplace=True)

# 绘制热力图
sns.heatmap(heatmap_normalized, annot=heatmap_data.values, fmt='.2f', cmap='YlGnBu', 
            linewidths=.5, cbar_kws={'label': labels['label4_1']})

if not use_english and 'chinese_font' in locals():
    plt.title(labels['title4'], fontproperties=chinese_font)
else:
    plt.title(labels['title4'])

plt.savefig('parallel_efficiency_heatmap.png', dpi=300)
plt.show()

print("可视化分析完成！图表已保存。")