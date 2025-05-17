void reverse(int *a, int n, int bit) {
    // a 是输入数组，n 是数组长度，bit 是二进制位数
    int *rev = new int[n];
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
    delete[] rev;
}