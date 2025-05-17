#pragma once
#include <stdint.h>

#include <stdexcept>

// 针对中等大小模数优化的蒙哥马利乘法类(最大支持469762049)
class Montgomery32 {
  public:
    uint32_t N;         // 模数
    uint64_t R;         // 2^32
    uint32_t logR;      // 32
    uint32_t N_inv_neg; // -N^(-1) mod R
    uint32_t R2;        // R² mod N
    uint32_t R_minus_1; // R - 1
    uint32_t R_mod_N;   // R

  public:
    // 构造函数 - 针对 N <= 469762049 进行优化
    Montgomery32(uint32_t N) : N(N) {

        if (N == 0 || (N & 1) == 0) {
            throw std::runtime_error("N 必须是正奇数");
        }
        // 计算N的有效位数

        this->logR = 32; // logR = bits
        this->R = (1ull << logR);

        //  this->R = (1u << logR);
        this->R_minus_1 = 0xffffffff;

        // 计算 N⁻¹ mod R
        uint32_t N_inv = modinv(N, R);
        this->N_inv_neg = R - N_inv; // -N⁻¹ mod R

        // 预计算 R² mod N (避免溢出)
        uint64_t R_mod_N = R % N;
        this->R2 = (uint32_t)((R_mod_N * R_mod_N) % N);
        this->R_mod_N = this->REDC(this->R2);
    }

    // 改进的REDC，使用32位整数，减少溢出检查
    uint32_t REDC(uint64_t T) const {

        // uint32_t T_low = T & (R - 1);
        uint32_t T_low = T & R_minus_1;

        // m = (T_low * N_inv_neg) & R_minus_1;32位不会溢出
        uint32_t m = (T_low * N_inv_neg);

        uint32_t t = (uint32_t)((T + (uint64_t)m * N) >> logR);
        return (t >= N) ? t - N : t;
    }

  private:
    // 扩展欧几里得算法求模逆元
    static int64_t extendedGCD(int64_t a, int64_t b, int64_t &x, int64_t &y) {
        if (a == 0) {
            x = 0;
            y = 1;
            return b;
        }
        int64_t x1, y1;
        int64_t gcd = extendedGCD(b % a, a, x1, y1);
        x = y1 - (b / a) * x1;
        y = x1;
        return gcd;
    }

    // 计算a⁻¹ mod m
    static uint32_t modinv(uint64_t a, uint64_t m) {
        int64_t x, y;
        int64_t gcd = extendedGCD(a, m, x, y);

        if (gcd != 1) {
            throw std::runtime_error("不存在模逆元");
        }
        return (uint32_t)(x % uint64_t(m) + m) % m; // 确保结果为正
    }
};