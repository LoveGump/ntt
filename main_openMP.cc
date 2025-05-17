#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/time.h>

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
        int x;
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

void poly_multiply(uint64_t *a, uint64_t *b, uint64_t *ab, int n, uint64_t p) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            ab[i + j] = (1LL * a[i] * b[j] % p + ab[i + j]) % p;
        }
    }
}

// 快速幂
uint64_t pow(int base, int exp, uint64_t mod) {
    uint64_t res = 1;
    while (exp > 0) {
        if (exp & 1) {
            res = 1LL * res * base % mod;
        }
        base = 1LL * base * base % mod;
        exp >>= 1;
    }
    return res;
}

void reverse(uint64_t *a, int n, int bit) {
    // a 是输入数组，n 是数组长度，bit 是二进制位数
    int *rev = new int[n];
    // i 的反转 = (i/2) 的反转去掉最高位后，再补上 i 的最低位作为最高位。
    rev[0] = 0;
    for (int i = 0; i < n; i++) {
        // 二进制反转, rev[i] = 0;
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1));
    }
    for (int i = 0; i < n; i++) {
        if (i < rev[i]) {
            std::swap(a[i], a[rev[i]]); // 位逆序置换
        }
    }
}

void NTT(uint64_t *a, uint64_t n, bool invert, uint64_t p, int g = 3) {
    // 计算二进制位数
    int bit = 0;
    while ((1 << bit) < n)
        bit++;

    // 位逆序置换
    reverse(a, n, bit);

    // 蝶形操作
    for (int len = 2; len <= n; len <<= 1) {
        // len 是当前的蝶形长度，len = 2^k
        // wlen 是当前的单位根，wlen = g^(p-1)/len
        // 计算单位根，原根为g，G^((p-1)/len)
        // 如果invert为true，则单位根为g^((p-1)/len) 的逆元
        uint64_t g_n = invert ? pow(g, (p - 1) - (p - 1) / len, p)
                              : pow(g, (p - 1) / len, p);

        // 处理每个块
        for (int i = 0; i < n; i += len) {
            uint64_t g = 1; // 初始单位根为1 k = 0 ，g = 1

            // 处理每个蝶形单元
            for (int j = 0; j < len / 2; j++) {
                uint64_t u = a[i + j];
                uint64_t v = (1LL * a[i + j + len / 2] * g) % p;
                // F(omega^k*n) = G(omega^k\_{n/2}) +
                // omega^k_n*H(omega^k\_{n/2})
                a[i + j] = (u + v) % p;
                // F(omega^{k+n/2}*n) = G(omega^k*{n/2}) -
                // omega^k_n*H(omega^k\_{n/2})
                a[i + j + len / 2] =
                    (u - v + p) % p; // 注意要加上p再取模，保证结果非负

                g = (1LL * g * g_n) % p; // 更新单位根
            }
        }
    }

    // 如果是逆变换，需要除以n（即乘以n的模p逆元）
    if (invert) {
        uint64_t inv_n = pow(n, p - 2, p); // 使用费马小定理计算逆元
        for (int i = 0; i < n; i++) {
            a[i] = (1LL * a[i] * inv_n) % p;
        }
    }
}

// 使用NTT的多项式乘法
void NTT_multiply(uint64_t *a, uint64_t *b, uint64_t *result, int n,
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
    NTT(fa, len, false, p, g);
    NTT(fb, len, false, p, g);

    // 点值表示下的乘法
    for (int i = 0; i < len; i++) {
        fa[i] = (1LL * fa[i] * fb[i]) % p;
    }

    // 逆向NTT得到结果
    NTT(fa, len, true, p, g);

    // 复制结果
    for (int i = 0; i < 2 * n - 1; i++) {
        result[i] = fa[i];
    }

    delete[] fa;
    delete[] fb;
}

// 基4NTT
// 这里的基4NTT是基于基2NTT的优化版本
// 通过将输入数据分成4个部分来减少蝶形操作的次数
// 基4的位逆序置换函数
void reverse_base4(uint64_t *a, int n) {

    int log4n = 0; // log4n是n的基4对数
    int temp = n;  // 临时变量
    while (temp > 1) {
        temp >>= 2;
        log4n++;
    }
    // 计算n的基4对数
    int *rev = new int[n];
    for (int i = 0; i < n; i++) {
        int reversed = 0;
        int num = i;
        for (int j = 0; j < log4n; j++) {
            reversed = (reversed << 2) | (num & 3); // 每次取最低两位并左移
            num >>= 2;
        }
        rev[i] = reversed;
    }

    // 位逆序置换
    for (int i = 0; i < n; i++) {
        if (i < rev[i]) {
            std::swap(a[i], a[rev[i]]);
        }
    }
    delete[] rev;
}

// ------------------------------------------------------------------------------------

// 基4 NTT/INTT核心函数
void NTT_base4_naive(uint64_t *a, int n, bool invert, uint64_t p, int g = 3) {
    reverse_base4(a, n);

    // 蝶形操作
    for (int len = 4; len <= n; len <<= 2) {
        // len 是当前的蝶形长度，len = 4^k
        // wlen 是当前的单位根，wlen = g^(p-1)/len
        // 计算单位根，原根为g，G^((p-1)/len)
        // 如果invert为true，则单位根为g^((p-1)/len) 的逆元
        // 使用Montgomery计算单位根
        // g_n^1
        uint64_t g_n = invert ? pow(g, (p - 1) - (p - 1) / len, p)
                              : pow(g, (p - 1) / len, p);

        // 预计算单位根幂
        uint64_t g_n2 = (1LL * g_n * g_n) % p;
        uint64_t g_n3 = (1LL * g_n2 * g_n) % p;

        int step = len >> 2;
        uint64_t g_pow_step = pow(g_n, step, p);

        for (int i = 0; i < n; i += len) {
            // 处理每个蝴蝶形单元
            // 这里的len是4的幂次

            uint64_t w[4] = {1, 1, 1, 1}; // 当前单位根幂次

            // 当前单位根幂次
            uint64_t u[4];
            for (int j = 0; j < len / 4; j++) {
                if (j == 0) {
                    for (int k = 0; k < 4; k++) {
                        u[k] = a[i + j + k * step];
                    }
                } else {
                    for (int k = 0; k < 4; k++) {
                        u[k] = 1LL * a[i + j + k * step] * w[k] % p;
                    }
                }
                uint64_t j_1 = 1LL * u[1] * g_pow_step % p;
                uint64_t j_3 = 1LL * u[3] * g_pow_step % p;

                a[i + j] = (u[0] + u[1] + u[2] + u[3]) % p;
                a[i + j + step] = (u[0] + j_1 + p - u[2] + p - j_3) % p;
                a[i + j + 2 * step] = (u[0] + p - u[1] + u[2] + p - u[3]) % p;
                a[i + j + 3 * step] = (u[0] + p - j_1 + p - u[2] + j_3) % p;
                w[1] = 1LL * w[1] * g_n % p;
                w[2] = 1LL * w[2] * g_n2 % p;
                w[3] = 1LL * w[3] * g_n3 % p;
            }
        }
    }
    // 处理逆变换的系数
    if (invert) {
        int inv_n = pow(n, p - 2, p);
        for (int i = 0; i < n; i++) {
            a[i] = 1LL * a[i] * inv_n % p;
        }
    }
}

// 基4多项式乘法
void NTT_multiply_base4_Naive(uint64_t *a, uint64_t *b, uint64_t *result, int n,
                              uint64_t p) {
    // 准备工作不变
    // len 为4的幂次
    int len = 1;
    while (len < 2 * n)
        len <<= 2;

    uint64_t *fa = new uint64_t[len];
    uint64_t *fb = new uint64_t[len];

    for (int i = 0; i < n; i++) {
        fa[i] = a[i];
        fb[i] = b[i];
    }
    for (int i = n; i < len; i++) {
        fa[i] = fb[i] = 0;
    }

    // 使用基4 NTT
    NTT_base4_naive(fa, len, false, p);
    NTT_base4_naive(fb, len, false, p);

    for (int i = 0; i < len; i++) {
        fa[i] = (1LL * fa[i] * fb[i]) % p;
    }

    // 逆NTT
    NTT_base4_naive(fa, len, true, p);

    // 复制结果
    for (int i = 0; i < 2 * n - 1; i++) {
        result[i] = fa[i];
    }

    delete[] fa;
    delete[] fb;
}

uint64_t a[300000], b[300000], ab[300000];
//
int main(int argc, char *argv[]) {

    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1
    // 的形式 输入模数分别为 7340033 104857601 469762049 1337006139375617
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
        NTT_multiply_base4_Naive(a, b, ab, n_, p_);

        auto End = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::ratio<1, 1000>> elapsed =
            End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);

        std::cout << "average latency for n = " << n_ << " p = " << p_ << " : "
                  << ans << " (us) " << std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab, n_, i);
    }
    return 0;
}
