# NTT多项式乘法算法实现

本项目实现了使用数论变换(Number Theoretic Transform, NTT)的多项式乘法算法，提供了三种不同的实现版本：串行版本和两种并行版本。同时支持普通模数和大模数(通过中国剩余定理CRT)情况下的多项式乘法。

## 项目文件说明

- **main_naive.cc**: 串行实现版本，提供基础的NTT和CRT-NTT实现，不使用任何并行技术。
- **main_pthread.cc**: 使用POSIX线程(pthread)实现的并行版本，采用线程池机制。
- **main_openMP.cc**: 使用OpenMP实现的并行版本，通过编译指令自动并行化。

## 算法介绍

### 普通NTT多项式乘法

将多项式转换为点值表示形式，进行点值乘法，再转换回系数表示：
1. 使用NTT将两个多项式系数表示转换为点值表示
2. 点值表示下进行系数乘法
3. 使用逆NTT将结果转换回系数表示

### CRT-NTT大模数多项式乘法

当模数过大时(超过uint64_t范围或不满足NTT友好条件)，采用中国剩余定理(CRT)：
1. 选取多个较小的NTT友好模数
2. 在每个模数下计算NTT多项式乘法
3. 使用CRT合并多个模数下的结果为最终结果

## 实现版本对比

### 1. 串行版本 (main_naive.cc)
- 基础的单线程实现
- 包含直接NTT和基于CRT的NTT实现
- 适合作为基准进行性能比较

### 2. Pthread版本 (main_pthread.cc)
- 使用POSIX线程库实现并行计算
- 采用线程池机制避免频繁创建/销毁线程
- 通过任务队列和barrier同步机制控制并行任务
- 适合细粒度控制线程行为和资源分配

### 3. OpenMP版本 (main_openMP.cc)
- 使用OpenMP编译指令简化并行实现
- 通过简单的`#pragma omp parallel for`指令实现循环并行
- 使用`schedule(static)`实现静态负载均衡
- 代码简洁，易于维护，同时保持高性能

## 性能优化

所有版本都实现了以下优化：
- 位逆序置换优化
- 模运算优化
- 内存预分配和重用

并行版本特有的优化：
- 多层次并行：NTT蝶形操作并行、点值乘法并行、CRT并行等
- 局部性优化：提高缓存命中率
- 动态任务分配：根据问题规模调整并行度

## 使用方法

编译程序：
```bash
g++ main_openMP.cc -o main -O2 -fopenmp -std=c++11
g++ main_pthread.cc -o main -O2 -lpthread -std=c++11
g++ main_naive.cc -o main -O2  -std=c++11
```

运行程序：
```bash
bash test.sh 2 1 1 # ntt_naive
bash test.sh 2 2 8 # ntt_pthread
bash test.sh 3 2 8 # ntt_openMP
```

程序将自动读取测试数据，执行多项式乘法，并验证结果正确性。

