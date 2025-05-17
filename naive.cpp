#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <pthread.h>
#include <string>

// 可以自行添加需要的头文件

void fRead(uint64_t *a, uint64_t *b, int *n, uint64_t *p, int input_id) {
    // 数据输入函数
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin >> *n >> *p;
    for (uint64_t i = 0; i < *n; i++) {
        fin >> a[i];
    }
    for (uint64_t i = 0; i < *n; i++) {
        fin >> b[i];
    }
}

void fCheck(uint64_t *ab, int n, int input_id) {
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++) {
        uint64_t x;
        fin >> x;
        if (x != ab[i]) {
            std::cout << "多项式乘法结果错误" << std::endl;
            return;
        }
    }
    std::cout << "多项式乘法结果正确" << std::endl;
    return;
}

void fWrite(uint64_t *ab, int n, int input_id) {
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++) {
        fout << ab[i] << '\n';
    }
}

#define NUM_THREADS 8

// 全局变量，用于传递给线程的参数
struct NTTThreadArgs {
    uint64_t *a;   // 输入数组
    int n;         // 数组长度
    int len;       // 当前蝶形长度
    int thread_id; // 线程ID
    uint64_t p;    // 模数
    uint64_t g_n;  // 当前单位根
};

struct MulThreadArgs {
    uint64_t *fa;  // 输入数组A
    uint64_t *fb;  // 输入数组B
    int len;       // 数组长度
    int thread_id; // 线程ID
    uint64_t p;
};

// 幂取模
uint64_t pow(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) {
            result = (uint64_t)((uint64_t)result * base % mod);
        }
        base = (uint64_t)((uint64_t)base * base % mod);
        exp >>= 1;
    }
    return result;
}

void reverse(uint64_t *a, int n, int bit) {
    // a 是输入数组，n 是数组长度，bit 是二进制位数
    int *rev = new int[n];
    // i 的反转 = (i/2) 的反转去掉最高位后，再补上 i 的最低位作为最高位。
    rev[0] = 0;
    for (int i = 0; i < n; i++) {
        // 进制反转, rev[i] = 0;
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1));
    }
    for (int i = 0; i < n; i++) {
        if (i < rev[i]) {
            // 使用异或交换，避免临时变量
            a[i] ^= a[rev[i]];
            a[rev[i]] ^= a[i];
            a[i] ^= a[rev[i]];
        }
    }
    delete[] rev;
}

// 处理单个块的线程函数
void *ntt_thread_function(void *arg) {
    NTTThreadArgs *args = (NTTThreadArgs *)arg;
    uint64_t *a = args->a;
    int n = args->n;
    int len = args->len;
    int thread_id = args->thread_id;

    uint64_t p = args->p;
    uint64_t g_n = args->g_n;

    // 每个线程处理一部分块
    for (int i = thread_id * len; i < n; i += len * NUM_THREADS) {
        uint64_t g = 1; // 初始单位根为1

        // 处理每个蝶形单元
        for (int j = 0; j < len / 2; j++) {
            uint64_t u = a[i + j];
            uint64_t v = (a[i + j + len / 2] * g) % p;
            // F(g^k*n) = G(g^k_{n/2}) + g^k_n*H(g^k_{n/2})
            a[i + j] = (u + v) % p;
            // F(g^{k+n/2}*n) = G(g^k*{n/2}) - g^k_n*H(g^k_{n/2})
            a[i + j + len / 2] =
                (u - v + p) % p; // 注意要加上p再取模，保证结果非负

            g = (g * g_n) % p; // 更新单位根
        }
    }

    return NULL;
}

// 点值乘法的线程函数
void *multiply_thread_function(void *arg) {
    MulThreadArgs *args = (MulThreadArgs *)arg;
    uint64_t *fa = args->fa;
    uint64_t *fb = args->fb;
    int len = args->len;
    int thread_id = args->thread_id;

    uint64_t p = args->p;

    // 每个线程处理一部分点值乘法
    for (int i = thread_id; i < len; i += NUM_THREADS) {
        fa[i] = (fa[i] * fb[i]) % p;
    }

    return NULL;
}

void NTT_parallel(uint64_t *a, uint64_t n, bool invert, uint64_t p, int g = 3) {
    // 算二进制位数
    int bit = 0;
    while ((1 << bit) < n)
        bit++;

    // 逆序置换（这部分不并行，因为有依赖关系）
    reverse(a, n, bit);

    // 准备线程
    pthread_t threads[NUM_THREADS];
    NTTThreadArgs thread_args[NUM_THREADS];

    // 形操作
    for (int len = 2; len <= n; len <<= 1) {
        // len 是当前的蝶形长度，len = 2^k
        // 算单位根
        uint64_t g_n = invert ? pow(g, (p - 1) - (p - 1) / len, p)
                              : pow(g, (p - 1) / len, p);

        // 初始化线程参数
        for (int i = 0; i < NUM_THREADS; i++) {
            thread_args[i].a = a;
            thread_args[i].n = n;
            thread_args[i].len = len;
            thread_args[i].thread_id = i;

            thread_args[i].p = p;
            thread_args[i].g_n = g_n;

            // 创建线程
            pthread_create(&threads[i], NULL, ntt_thread_function,
                           &thread_args[i]);
        }

        // 等待所有线程完成
        for (int i = 0; i < NUM_THREADS; i++) {
            pthread_join(threads[i], NULL);
        }
    }

    // 如果是逆变换，需要除以n（即乘以n的模p逆元）
    if (invert) {
        uint64_t inv_n = pow(n, p - 2, p); // 使用费马小定理计算逆元
        for (int i = 0; i < n; i++) {
            a[i] = a[i] * inv_n % p;
        }
    }
}

// 使用并行NTT的多项式乘法
void NTT_multiply_parallel(uint64_t *a, uint64_t *b, uint64_t *result, int n,
                           uint64_t p) {
    // 计算可以容纳结果的2的幂次长度
    int len = 1;
    while (len < 2 * n)
        len <<= 1;

    // 创建临时数组
    uint64_t *fa = new uint64_t[len];
    uint64_t *fb = new uint64_t[len];

    // 复制输入数组并填充0
    for (int i = 0; i < n; i++) {
        fa[i] = a[i];
        fb[i] = b[i];
    }
    for (int i = n; i < len; i++) {
        fa[i] = fb[i] = 0;
    }

    // 确定原根
    int g = 3;

    // 正向NTT
    NTT_parallel(fa, len, false, p, g);
    NTT_parallel(fb, len, false, p, g);

    // 点值表示下的乘法（并行）
    pthread_t threads[NUM_THREADS];
    MulThreadArgs mul_args[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; i++) {
        mul_args[i].fa = fa;
        mul_args[i].fb = fb;
        mul_args[i].len = len;
        mul_args[i].thread_id = i;
        mul_args[i].p = p;
        pthread_create(&threads[i], NULL, multiply_thread_function,
                       &mul_args[i]);
    }

    // 等待所有线程完成
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    // 逆向NTT得到结果
    NTT_parallel(fa, len, true, p, g);

    // 复制结果
    for (int i = 0; i < 2 * n - 1; i++) {
        result[i] = fa[i];
    }

    delete[] fa;
    delete[] fb;
}

// 示例主函数

uint64_t a[300000], b[300000], ab[300000];
int main(int argc, char *argv[]) {

    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1
    // 的形式 输入模数分别为 7340033 104857601 469762049 1337006139375617
    // 167772161 998244353 1004535809 469762049
    // 第四个模数超过了整型表示范围,
    // 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT,
    // 请在完成前三个模数的基础代码及优化后实现大模数 NTT 输入文件共五个,
    // 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久,
    // 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 4;
    for (int i = test_begin; i <= test_end; ++i) {
        long double ans = 0;
        int n_;
        uint64_t p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        auto Start = std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
        // poly_multiply(a, b, ab, n_, p_);
        // NTT_multiply(a, b, ab, n_, p_);
        // NTT_multiply_base4_Naive(a, b, ab, n_, p_);
        NTT_multiply_parallel(a, b, ab, n_, p_);

        auto End = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::ratio<1, 1000>> elapsed =
            End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);

        std::cout << "average latency for n = " << n_ << " p = " << p_ << " : "
                  << ans << " (us) " << std::endl;
        // 结果输出
        fWrite(ab, n_, i);
    }
    return 0;
}
