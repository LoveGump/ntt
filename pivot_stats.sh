#!/bin/bash

INPUT_CSV="perf_results/performance_stats.csv"
OUTPUT_CSV="perf_results/performance_stats_pivot.csv"

# 获取所有唯一的指标名称
metrics=$(tail -n +2 $INPUT_CSV | cut -d',' -f3 | sort | uniq)

# 创建表头行
echo -n "Program,Run" > $OUTPUT_CSV
for metric in $metrics; do
    echo -n ",$metric" >> $OUTPUT_CSV
done
echo "" >> $OUTPUT_CSV

# 为每个程序和运行次数创建一行
for program in naive openMP pthread; do
    for run in {1..5}; do
        echo -n "$program,$run" >> $OUTPUT_CSV
        
        # 对于每个指标，查找并添加对应的值
        for metric in $metrics; do
            value=$(grep "^$program,$run,$metric," $INPUT_CSV | cut -d',' -f4)
            # 如果找不到值，使用NA
            if [ -z "$value" ]; then
                value="NA"
            fi
            echo -n ",$value" >> $OUTPUT_CSV
        done
        echo "" >> $OUTPUT_CSV
    done
done

echo "转换后的CSV文件已创建: $OUTPUT_CSV"