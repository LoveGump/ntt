#!/bin/bash

# 分析不同实现的并行效率
for prog in naive openMP pthread; do
    echo "分析 $prog..."
    
    # 收集基本性能数据
    perf stat -e cycles,instructions,context-switches,cpu-migrations \
        -o perf_results/$prog-stat.txt ./$prog
    
    # 收集核心级别数据
    perf stat --per-core -e cycles,instructions \
        -o perf_results/$prog-cores.txt ./$prog
done

# 提取并计算并行效率指标
echo "计算并行效率..."
python3 - <<EOF
import re
import pandas as pd

# 读取时间数据
progs = ['naive', 'openMP', 'pthread']
times = {}
for prog in progs:
    with open(f'perf_results/{prog}-stat.txt', 'r') as f:
        content = f.read()
        # 提取elapsed时间
        elapsed = re.search(r'([0-9.]+) seconds time elapsed', content)
        if elapsed:
            times[prog] = float(elapsed.group(1))

# 计算加速比
if 'naive' in times:
    base_time = times['naive']
    print(f"基准时间(naive): {base_time:.4f}秒")
    
    for prog in progs:
        if prog != 'naive' and prog in times:
            speedup = base_time / times[prog]
            print(f"{prog}加速比: {speedup:.2f}x")
EOF