#pragma once
#include <stdint.h>

#include <arm_neon.h>
#include <stdexcept>


// 针对中等大小模数优化的蒙哥马利乘法类(最大支持469762049)
class Montgomery32 {
public:
    uint32_t N;        // 模数
    uint64_t R;        // 2^32
    uint32_t logR;     // 32
    uint32_t N_inv_neg;// -N^(-1) mod R
    uint32_t R2;       // R² mod N
    uint32_t R_minus_1; // R - 1
    uint32_t R_mod_N;  // R

    uint32x4_t N_vec;    // 模数
    uint32x4_t N_inv_neg_vec;  // -N^(-1) mod R
    uint32x4_t R2_vec;  // R² mod N
    uint32x4_t R_minus_1_vec; // R - 1

public:
    // 构造函数 - 针对 N <= 469762049 进行优化
    Montgomery32(uint32_t N) : N(N) {
        
        if (N == 0 || (N & 1) == 0) {
            throw std::runtime_error("N 必须是正奇数");
        }
        // 计算N的有效位数

        this->logR = 32 ; // logR = bits
        this->R =  (1ll << logR);
             
        //  this->R = (1u << logR);       
        this->R_minus_1 = 0xffffffff;
       
        // 计算 N⁻¹ mod R
        uint32_t N_inv = modinv(N, R);
        this->N_inv_neg = R - N_inv; // -N⁻¹ mod R
        
        // 预计算 R² mod N (避免溢出)
        uint64_t R_mod_N = R % N;
        this->R2 = (uint32_t)((R_mod_N * R_mod_N) % N);
        this->R_mod_N = this->REDC(this->R2);

        // 预计算向量
        this->N_vec = vdupq_n_u32(N);
        this->N_inv_neg_vec = vdupq_n_u32(N_inv_neg);
        this->R2_vec = vdupq_n_u32(R2);
        this->R_minus_1_vec = vdupq_n_u32(R - 1);
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
    // 接下来就是利用NEON指令集进行优化
    // NEON优化的REDC函数 - 一次处理4个数 数组处理
    uint32x4_t REDC_neon(const uint64_t T[4]) const {

        uint32_t result[4];
        
        // 步骤1: 提取每个64位T % R
        uint32_t T_low[4];
        for (int i = 0; i < 4; i++) {
            T_low[i] = T[i] & (R - 1);
        }
        
        // 加载低32位到NEON寄存器
        uint32x4_t T_low_vec = vld1q_u32(T_low);
        
        // 步骤2: 计算 m = (T_low * N_inv_neg) & (R - 1)
        uint32x4_t m_vec = vmulq_u32(T_low_vec, N_inv_neg_vec);
        m_vec = vandq_u32(m_vec, R_minus_1_vec);
        
        // 提取m值
        uint32_t m[4];
        vst1q_u32(m, m_vec);
        
        // 步骤3: 计算 t = (T + m*N) >> logR
        // 这一步涉及64位运算，无法直接用NEON指令实现
        for (int i = 0; i < 4; i++) {
            uint64_t t_full = T[i] + (uint64_t)m[i] * N;
            uint32_t t = (uint32_t)(t_full >> logR);
            result[i] = (t >= N) ? t - N : t;
        }
        
        // 加载结果到NEON寄存器并返回
        return vld1q_u32(result);
    }
    // NEON优化的REDC函数 - 一次处理4个数  寄存器处理（T 的低32位 和 T的高32位 ）
    uint32x4_t REDC_neon(uint32x4_t T_high  , uint32x4_t T_low )const {
        // 计算 m = (T_low * N_inv_neg) & R_minus_1
        uint32x4_t m = vmulq_u32(T_low, N_inv_neg_vec);
        
        // 由于需要64位完整结果，我们分解为低位和高位分别计算
        uint32x2_t m_low = vget_low_u32(m);   // m0 m1
        uint32x2_t m_high = vget_high_u32(m); // m2 m3
        uint32x2_t N_low = vget_low_u32(N_vec);  // N0 N1
        uint32x2_t N_high = vget_high_u32(N_vec);// N2 N3
        
        // 完整的乘法结果
        uint64x2_t mul_low = vmull_u32(m_low, N_low);
        uint64x2_t mul_high = vmull_u32(m_high, N_high);
        
        // 组合高 32 位
        uint32x4_t mul_result_high = vcombine_u32(
            vshrn_n_u64(mul_low, 32), 
            vshrn_n_u64(mul_high, 32)
        );
        
        // 右移32位
        uint32x4_t t = vaddq_u32(T_high, mul_result_high);
        
        // 最后处理减法
        return vbslq_u32(vcgeq_u32(t, N_vec), vsubq_u32(t, N_vec), t);
    }

    private:
            // 扩展欧几里得算法求模逆元
    static int64_t extendedGCD(int64_t a, int64_t b, int64_t& x, int64_t& y) {
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