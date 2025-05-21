#!/bin/bash

# 创建结果目录
mkdir -p perf_results

# 定义要分析的程序文件
PROGRAMS=("main_naive.cc" "main_openMP.cc" "main_pthread.cc")
PROGRAM_NAMES=("naive" "openMP" "pthread")

# 编译选项
COMPILER="g++"
OPTIMIZATION="-O2"
ARCH="-march=native"
THREAD_FLAGS="-pthread"
OMP_FLAGS="-fopenmp"
DEBUG="-g"  # 确保perf可以获取符号信息

# perf事件列表
PERF_EVENTS="cycles,instructions,L1-dcache-loads,L1-dcache-load-misses,l2_rqsts.references,l2_rqsts.miss,LLC-loads,LLC-load-misses,branch-misses"

# 重复次数
REPEAT_COUNT=5

# 创建汇总结果文件
OUTPUT_SUMMARY="perf_results/program_outputs.txt"
PERF_SUMMARY="perf_results/perf_results.txt"

# 初始化汇总文件
echo "# 程序输出汇总" > $OUTPUT_SUMMARY
echo "===================" >> $OUTPUT_SUMMARY
echo "" >> $OUTPUT_SUMMARY

echo "# 性能计数器汇总" > $PERF_SUMMARY
echo "===================" >> $PERF_SUMMARY
echo "" >> $PERF_SUMMARY

# 编译和运行所有程序
for i in "${!PROGRAMS[@]}"; do
    PROGRAM="${PROGRAMS[$i]}"
    NAME="${PROGRAM_NAMES[$i]}"
    
    echo "===== 处理程序: $NAME ====="
    
    # 编译程序（根据文件名添加合适的编译标志）
    if [[ "$PROGRAM" == *"openMP"* ]]; then
        # OpenMP程序需要额外的标志
        echo "编译 OpenMP 程序..."
        $COMPILER $PROGRAM -o "${NAME}" $OPTIMIZATION $ARCH $THREAD_FLAGS $OMP_FLAGS $DEBUG
    else
        # 其他程序
        echo "编译程序..."
        $COMPILER $PROGRAM -o "${NAME}" $OPTIMIZATION $ARCH $THREAD_FLAGS $DEBUG
    fi
    
    if [ $? -ne 0 ]; then
        echo "编译失败，跳过此程序"
        continue
    fi
    
    # 创建程序专用结果目录
    mkdir -p "perf_results/${NAME}"
    
    # 运行程序并记录输出
    echo "运行程序并记录基本输出..."
    ./"${NAME}" > "perf_results/${NAME}/output.txt" 2>&1
    
    # 添加程序输出到汇总文件
    {
        echo "## $NAME 程序输出"
        echo ""
        cat "perf_results/${NAME}/output.txt"
        echo ""
        echo "-------------------------"
        echo ""
    } >> $OUTPUT_SUMMARY
    
    # 使用perf进行重复测试并记录到汇总文件
    echo "## $NAME 性能分析 ($REPEAT_COUNT 次重复)" >> $PERF_SUMMARY
    
    for ((run=1; run<=$REPEAT_COUNT; run++)); do
        echo "  运行 $run/$REPEAT_COUNT"
        
        # 使用perf stat收集性能计数器数据
        echo "### 运行 $run" >> $PERF_SUMMARY
        perf stat -e $PERF_EVENTS -o "perf_results/${NAME}/stat_run${run}.txt" ./"${NAME}" > /dev/null 2>&1
        cat "perf_results/${NAME}/stat_run${run}.txt" >> $PERF_SUMMARY
        echo "" >> $PERF_SUMMARY
        
        # 获取core级别统计
        echo "#### Core级别统计 - 运行 $run" >> $PERF_SUMMARY
        perf stat --per-core -e $PERF_EVENTS -o "perf_results/${NAME}/core_stat_run${run}.txt" ./"${NAME}" > /dev/null 2>&1
        cat "perf_results/${NAME}/core_stat_run${run}.txt" >> $PERF_SUMMARY
        echo "" >> $PERF_SUMMARY
        
        # 获取更详细的性能数据记录
        perf record -F 99 -e $PERF_EVENTS -g -o "perf_results/${NAME}/perf_run${run}.data" ./"${NAME}" > /dev/null 2>&1
        
        # 生成报告
        perf report --stdio -i "perf_results/${NAME}/perf_run${run}.data" > "perf_results/${NAME}/report_run${run}.txt"
    done
    
    # 添加性能统计摘要
    echo "### $NAME 性能统计摘要" >> $PERF_SUMMARY
    echo "" >> $PERF_SUMMARY
    
    # 计算每个指标的平均值
    for metric in cycles instructions L1-dcache-loads L1-dcache-load-misses l2_rqsts.references l2_rqsts.miss LLC-loads LLC-load-misses branch-misses; do
        echo "#### $metric" >> $PERF_SUMMARY
        grep "$metric" "perf_results/${NAME}/stat_run"*.txt | awk -v metric="$metric" '
        BEGIN { sum = 0; count = 0; }
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[0-9,]+$/) {
                    gsub(/,/, "", $i);
                    sum += $i;
                    count++;
                    break;
                }
            }
        }
        END {
            if (count > 0) printf "平均值: %.2f\n", sum/count;
            else print "无数据";
        }' >> $PERF_SUMMARY
        echo "" >> $PERF_SUMMARY
    done
    
    echo "完成 ${NAME} 的分析" >> $PERF_SUMMARY
    echo "-------------------------" >> $PERF_SUMMARY
    echo "" >> $PERF_SUMMARY
    
    echo "完成 ${NAME} 的分析"
    echo "-------------------------"
done

# 添加比较总结部分到性能汇总文件
echo "# 性能比较总结" >> $PERF_SUMMARY
echo "===================" >> $PERF_SUMMARY
echo "" >> $PERF_SUMMARY

echo "## 运行时间比较" >> $PERF_SUMMARY
echo "" >> $PERF_SUMMARY

for NAME in "${PROGRAM_NAMES[@]}"; do
    echo "### $NAME" >> $PERF_SUMMARY
    grep "average latency" "perf_results/${NAME}/output.txt" >> $PERF_SUMMARY
    echo "" >> $PERF_SUMMARY
done

echo "分析完成！结果保存在 perf_results/ 目录"
echo "程序输出汇总: $OUTPUT_SUMMARY"
echo "性能数据汇总: $PERF_SUMMARY"