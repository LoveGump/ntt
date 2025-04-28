#!/bin/bash
# filepath: /root/ntt/perf_analysis.sh

# 创建结果目录
mkdir -p perf_results

# 定义要分析的算法数组
algorithms=(
    "FFT_multiply(a, b, ab, n_, p_);"
    "NTT_multiply(a, b, ab, n_, p_);"
    "NTT_multiply_Montgomery(a, b, ab, n_, p_);"
    "NTT_multiply_Montgomerybase2(a, b, ab, n_, p_);"
    "NTT_multiply_base4_Naive(a, b, ab, n_, p_);"
    "NTT_multiply_base4_m(a, b, ab, n_, p_);"
    "NTT_multiply_base4_Montgomery_domain(a, b, ab, n_, p_);"
    "NTT_multiply_base4_Montgomery32neon(a, b, ab, n_, p_);"
)

# 定义算法名称（用于结果文件）
algo_names=(
    "fft"
    "ntt_basic"
    "ntt_montgomery"
    "ntt_montgomery_domain"
    "ntt_radix4"
    "ntt_radix4_m"
    "ntt_radix4_montgomery"
    "ntt_radix4_montgomery_neon"
)

# 编译和性能分析选项
COMPILER="g++"
SOURCE="main.cc"
OUTPUT="ntt_perf"
OPTIMIZATION="-O2"
ARCH="-march=native"
DEBUG="-g"  # 使perf可以获取符号信息

# 备份原始main.cc
cp main.cc main.cc.bak

# 设置测试范围（建议使用较小的测试集进行性能分析）
sed -i 's/int test_begin = [0-9]\+;/int test_begin = 1;/g' main.cc
sed -i 's/int test_end = [0-9]\+;/int test_end = 1;/g' main.cc

# 对每个算法进行分析
for i in "${!algorithms[@]}"; do
    algo="${algorithms[$i]}"
    name="${algo_names[$i]}"
    
    echo "===== 分析算法: $name ====="
    
    # 修改main.cc启用当前算法，禁用其他算法
    cp main.cc.bak main.cc
    sed -i 's|^[ \t]*\([A-Za-z_]\+\)(a, b, ab, n_, p_);|        // \1(a, b, ab, n_, p_);|g' main.cc
    sed -i "s|^[ \t]*// *${algo%;}|        ${algo%;}|g" main.cc
    
    # 编译
    echo "编译程序..."
    $COMPILER -o $OUTPUT $SOURCE $OPTIMIZATION $ARCH $DEBUG
    
    if [ $? -ne 0 ]; then
        echo "编译失败，跳过此算法"
        continue
    fi
    
    # 使用perf进行性能分析
    echo "运行perf分析..."
    
    # CPU周期分析
    perf stat -e cycles,instructions,cache-references,cache-misses -o "perf_results/${name}_stat.txt" ./$OUTPUT
    
    # 详细性能记录
    perf record -F 99 -g -o "perf_results/${name}_perf.data" ./$OUTPUT
    
    # 生成报告
    perf report --stdio -i "perf_results/${name}_perf.data" > "perf_results/${name}_report.txt"
    
    echo "性能分析完成，结果保存在 perf_results/${name}_*.txt"
    echo ""
done

# 恢复原始main.cc
mv main.cc.bak main.cc

echo "所有算法性能分析完成！"
echo "结果保存在 perf_results/ 目录"

# 生成汇总报告
echo "生成性能汇总报告..."
echo "算法性能比较" > perf_results/summary.txt
echo "================" >> perf_results/summary.txt

for name in "${algo_names[@]}"; do
    echo "" >> perf_results/summary.txt
    echo "## $name" >> perf_results/summary.txt
    grep -A 5 "Performance counter stats" "perf_results/${name}_stat.txt" >> perf_results/summary.txt
done

echo "汇总报告已生成: perf_results/summary.txt"