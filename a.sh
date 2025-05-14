#!/bin/bash
# filepath: /root/ntt/a.sh

# 定义变量
OUTPUT_NAME="radix4_neon"
COMPILER="g++"
SOURCE="main.cc"
OPTIMIZATION="-O2"
ARCH="-march=native"

# 使用变量
$COMPILER -o $OUTPUT_NAME $SOURCE $OPTIMIZATION $ARCH
./$OUTPUT_NAME