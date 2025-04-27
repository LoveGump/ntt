#pragma once
#include <stdint.h>
#include <arm_neon.h>
#include <stdexcept>


// 针对中等大小模数优化的蒙哥马利乘法类(最大支持469762049)
class Montgomery32 {
public:
    uint32_t N;        // 模数
    uint32_t R;        // 
    uint32_t logR;     // R的位数
    uint32_t N_inv_neg;// -N^(-1) mod R
    uint32_t R2;       // R² mod N

public:
    // 构造函数 - 针对 N <= 469762049 进行优化
    Montgomery32(uint32_t N) : N(N) {
        if (N == 0 || (N & 1) == 0) {
            throw std::runtime_error("N 必须是正奇数");
        }
        // 计算N的有效位数
        uint32_t bits = 0;
        uint32_t temp = N;
        while (temp > 0) {
            temp >>= 1;
            bits++;
        }
        this->logR = bits ; // logR = bits
        this->R = (1u << logR);       
        if (R <= N) {
            throw std::runtime_error("N 必须小于 2^30");
        }
        
        // 计算 N⁻¹ mod R
        uint32_t N_inv = modinv(N, R);
        this->N_inv_neg = R - N_inv; // -N⁻¹ mod R
        
        // 预计算 R² mod N (避免溢出)
        uint64_t R_mod_N = R % N;
        this->R2 = (uint32_t)((R_mod_N * R_mod_N) % N);

 
    }
    
    // 改进的REDC，使用32位整数，减少溢出检查
    uint32_t REDC(uint64_t T) const {
        uint32_t T_low = T & (R - 1);
        uint32_t m = (T_low * N_inv_neg) & (R - 1);
        uint32_t t = (uint32_t)((T + (uint64_t)m * N) >> logR);
        return (t >= N) ? t - N : t;
    }
private:
    // 扩展欧几里得算法
    static int32_t extendedGCD(int32_t a, int32_t b, int32_t& x, int32_t& y) {
        if (a == 0) {
            x = 0;
            y = 1;
            return b;
        }
        int32_t x1, y1;
        int32_t gcd = extendedGCD(b % a, a, x1, y1);
        x = y1 - (b / a) * x1;
        y = x1;
        return gcd;
    }
    
    // 计算模逆元
    static uint32_t modinv(uint32_t a, uint32_t m) {
        int32_t x, y;
        int32_t gcd = extendedGCD(a, m, x, y);
        
        if (gcd != 1) {
            throw std::runtime_error("不存在模逆元");
        }
        return (x % (int32_t)m + m) % m;
    }
};