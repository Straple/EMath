using u64 = unsigned long long;

// x, y < mod
// mod < LLONG_MAX(~1e18)
// returns (x * y) % mod
// non overflow
u64 safe_multip(u64 x, u64 y, u64 mod) {
    u64 res = 0;
    while (y) {
        if (y & 1) {
            res += x;
            if (res >= mod) {
                res -= mod;
            }
            y--;
        }
        else {
            x <<= 1;
            if (x >= mod) {
                x -= mod;
            }

            y >>= 1;
        }
    }
    return res;
}
