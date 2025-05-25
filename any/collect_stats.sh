#!/bin/bash

# 定义CSV输出文件
CSV_OUTPUT="perf_results/performance_stats.csv"

# 定义要处理的程序
PROGRAM_NAMES=("naive" "openMP" "pthread")

# 定义要提取的指标
METRICS=(
    "cpu_atom/cycles/"
    "cpu_core/cycles/"
    "cpu_atom/instructions/"
    "cpu_core/instructions/"
    "cpu_atom/L1-dcache-loads/"
    "cpu_core/L1-dcache-loads/"
    "cpu_atom/L1-dcache-load-misses/"
    "cpu_core/L1-dcache-load-misses/"
    "cpu_core/l2_rqsts.references/"
    "cpu_core/l2_rqsts.miss/"
    "cpu_atom/LLC-loads/"
    "cpu_core/LLC-loads/"
    "cpu_atom/LLC-load-misses/"
    "cpu_core/LLC-load-misses/"
    "cpu_atom/branch-misses/"
    "cpu_core/branch-misses/"
    "seconds time elapsed"
    "seconds user"
    "seconds sys"
)

# 创建CSV表头
echo "Program,Run,Metric,Value" > $CSV_OUTPUT

# 处理每个程序的每个运行
for program in "${PROGRAM_NAMES[@]}"; do
    for run in {1..5}; do
        stat_file="perf_results/${program}/stat_run${run}.txt"
        
        # 检查文件是否存在
        if [ ! -f "$stat_file" ]; then
            echo "警告: 找不到文件 $stat_file"
            continue
        fi
        
        # 处理每个性能指标
        for metric in "${METRICS[@]}"; do
            # 根据指标类型选择不同的提取方法
            if [[ "$metric" == "seconds time elapsed" ]]; then
                # 提取时间经过的秒数
                value=$(grep "seconds time elapsed" "$stat_file" | awk '{print $1}')
                clean_metric="time_elapsed"
            elif [[ "$metric" == "seconds user" ]]; then
                # 提取用户时间
                value=$(grep "seconds user" "$stat_file" | awk '{print $1}')
                clean_metric="user_time"
            elif [[ "$metric" == "seconds sys" ]]; then
                # 提取系统时间
                value=$(grep "seconds sys" "$stat_file" | awk '{print $1}')
                clean_metric="sys_time"
            else
                # 提取性能计数器值
                value=$(grep "$metric" "$stat_file" | awk '{gsub(/,/, "", $1); print $1}')
                # 清理指标名称以适合CSV
                clean_metric=$(echo "$metric" | sed 's/\//_/g' | sed 's/\./_/g' | sed 's/-/_/g' | sed 's/ /_/g' | sed 's/_$//')
            fi
            
            # 如果找到值，则添加到CSV
            if [ ! -z "$value" ]; then
                echo "${program},${run},${clean_metric},${value}" >> $CSV_OUTPUT
            else
                # 如果找不到值，检查是否为不支持的指标
                not_supported=$(grep "<not supported>" "$stat_file" | grep -F "$metric")
                if [ ! -z "$not_supported" ]; then
                    echo "${program},${run},${clean_metric},not_supported" >> $CSV_OUTPUT
                else
                    # 否则记录为NA
                    echo "${program},${run},${clean_metric},NA" >> $CSV_OUTPUT
                fi
            fi
        done
    done
done

echo "CSV统计文件已创建: $CSV_OUTPUT"