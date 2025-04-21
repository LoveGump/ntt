#pragma once

// Montgomery乘法类
class Montgomery {
    public:
        uint64_t N;    // 模数
        uint64_t R;    // 通常选择2^k且满足 R > N, gcd(R,N)=1
        uint64_t logR; // R的二进制位数（如R=2^64 → logR=64）
        uint64_t N_inv_neg; // -N⁻¹ mod R
        uint64_t R2;   // R² mod N
    

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
        static uint64_t modinv(uint64_t a, uint64_t m) {
            int64_t x, y;
            int64_t gcd = extendedGCD(a, m, x, y);

            if (gcd != 1) {
                throw std::runtime_error("不存在模逆元");
            }
            return (x % uint64_t(m) + m) % m; // 确保结果为正
        }
    
    public:
        // 构造函数
        Montgomery(uint64_t N ) : N(N){
            if (N == 0 || (N & 1) == 0) {
                throw std::runtime_error("N 必须是正奇数");
            }
    
            // 计算R为大于N的最小2^k
            this->logR = 32; // 使用32位整数
            this->R = (1ULL << logR); // R = 2^32
            if (R <= N) {
                throw std::runtime_error("R 必须大于 N");
            }
    
            // 计算N⁻¹ mod R
            uint64_t N_inv = modinv(N, R);
            this->N_inv_neg = R - N_inv; // -N⁻¹ mod R
    
            // 预计算R² mod N
            this->R2 = (R % N) * (R % N) % N;
        }
    
        // Montgomery约简算法
        // T 是待约简的数
        // 返回值是约简后的结果
        // 结果范围在[0, N-1]之间
        uint64_t REDC(uint64_t T) const {
            uint64_t m = ((T & (R - 1)) * N_inv_neg) & (R - 1);
             // m = (T mod R) * (-N⁻¹) mod R
            uint64_t t = (T + m * N) >> logR;
             // t = (T + mN)/R
    
            // 结果规约到[0, N-1]
            return (t >= N) ? t - N : t;
        }
    
        // Montgomery乘法
        uint64_t multiply(uint64_t a, uint64_t b) const {
            if (a >= N || b >= N) {
                throw std::runtime_error("输入值必须小于模数 N");
            }
    
            // 转换为Montgomery形式
            uint64_t aR = REDC(a * R2);
            uint64_t bR = REDC(b * R2);
    
            // 标准乘法
            uint64_t T = aR * bR;
    
            // 约简并转换回普通形式
            return REDC(REDC(T));
        }
    
        // 辅助函数：快速模幂（利用Montgomery乘法）
        uint64_t pow(uint64_t x, uint64_t power) const {
            // 转换为Montgomery形式
            uint64_t xR = REDC(x * R2);
            uint64_t resultR = REDC(1 * R2); // 1 in Montgomery form
    
            while (power > 0) {
                if (power & 1) {
                    resultR = REDC(resultR * xR);
                }
                xR = REDC(xR * xR);
                power >>= 1;
            }
    
            // 转换回普通形式
            return REDC(resultR);
        }
};
  
