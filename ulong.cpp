// ulong
// Developed by Mob

#include<iostream>

// длинная арфиметика на беззнаковых числах длиной 64бит
template<size_t len>
class ulong {

    typedef uint64_t u64;
    typedef uint32_t u32;

    u64 words[len];

    // O(len)
    void clear() {
        for (int i = 0; i < len; i++) {
            words[i] = 0;
        }
    }

public:

    ulong() {
        clear();
    }
    ulong(u64 n) {
        clear();

        *words = n;
    }

    u64 operator [](int index) const {
        return words[index];
    }
    u64& operator [](int index) {
        return words[index];
    }

private:

    // O(len)
    int find_not_equall_ind(const ulong& Rhs) const {
        int i;
        for (i = len - 1; i >= 0 && words[i] == Rhs[i]; i--) {}
        return i;
    }

public:

    bool operator == (const ulong& Rhs) const {
        return find_not_equall_ind(Rhs) == -1;
    }
    bool operator != (const ulong& Rhs) const {
        return !(*this == Rhs);
    }

    bool operator < (const ulong& Rhs) const {
        int i = find_not_equall_ind(Rhs);
        return i != -1 && words[i] < Rhs[i];
    }
    bool operator > (const ulong& Rhs) const {
        int i = find_not_equall_ind(Rhs);
        return i != -1 && words[i] > Rhs[i];
    }

    bool operator <= (const ulong& Rhs) const {
        return !(*this > Rhs);
    }
    bool operator >= (const ulong& Rhs) const {
        return !(*this < Rhs);
    }

    // O(len)
    ulong operator + (const ulong& add) const {
        ulong res;

        for (int i = 0; i < len; i++) {
            res[i] += words[i] + add[i];
            if (res[i] < std::max(words[i], add[i]) && i + 1 < len) {
                res[i + 1]++;
            }
        }

        return res;
    }
    ulong& operator += (const ulong& add) {
        return *this = *this + add;
    }

    // O(len)
    ulong operator - (const ulong& sub) const {
        ulong res;

        bool c = 0;
        for (int i = 0; i < len; i++) {
            res[i] = words[i] - sub[i] - c;

            c = words[i] < sub[i] + c && i + 1 < len;
        }

        return res;
    }
    ulong& operator -= (const ulong& sub) {
        return *this = *this - sub;
    }

    // O((2 * len)^2)
    ulong operator * (const ulong& mult) const {
        auto mult1 = reinterpret_cast<const u32*>(words);
        auto mult2 = reinterpret_cast<const u32*>(mult.words);
        
        ulong result;
        u32* res = reinterpret_cast<u32*>(result.words);

        for (int i = 0; i < len * 2; i++) {
            u64 c = 0;
            for (int j = 0; i + j < len * 2; j++) {
                u64 k = static_cast<u64>(mult1[i]) * mult2[j] + c + res[i + j];

                c = k >> 32;

                res[i + j] = k - (c << 32);
            }
        }

        return result;
    }
    ulong& operator *= (const ulong& mult) {
        return *this = *this * mult;
    }

    // O(len)
    ulong& operator >>= (size_t pos) {
        pos %= len * 64;

        const int Wordshift = pos / 64;
        if (Wordshift != 0) {
            for (int Wpos = 0; Wpos < len; Wpos++) {
                words[Wpos] = Wordshift < len - Wpos ? words[Wpos + Wordshift] : 0;
            }
        }

        if ((pos %= 64) != 0) {
            for (int Wpos = 0; Wpos < len - 1; Wpos++) {
                words[Wpos] = (words[Wpos] >> pos) | (words[Wpos + 1] << (64 - pos));
            }

            words[len - 1] >>= pos;
        }
        return *this;
    }
    ulong operator >> (size_t pos) const {
        return ulong(*this) >>= pos;
    }

    // O(len)
    ulong& operator <<= (size_t pos) noexcept {
        pos %= len * 64;

        const int Wordshift = pos / 64;
        if (Wordshift != 0) {
            for (int Wpos = len - 1; 0 <= Wpos; Wpos--) {
                words[Wpos] = Wordshift <= Wpos ? words[Wpos - Wordshift] : 0;
            }
        }

        if ((pos %= 64) != 0) {
            for (int Wpos = len - 1; 0 < Wpos; Wpos--) {
                words[Wpos] = (words[Wpos] << pos) | (words[Wpos - 1] >> (64 - pos));
            }

            words[0] <<= pos;
        }
        return *this;
    }
    ulong operator << (size_t pos) const {
        return ulong(*this) <<= pos;
    }

    // O(len * 64)
    ulong operator / (const ulong& div) const {
        ulong res;

        int bit = 0; // сдвиг не нарушающий целостность div

        for (int i = 64 * len - 1; i >= 0; i--) {
            if (div[i / 64] & (static_cast<u64>(1) << i % 64)) {
                bit = 64 * len - i - 1;
                break;
            }
        }

        auto num = *this;
        for (; bit >= 0; bit--) {
            if ((div << bit) <= num) {
                num -= (div << bit);
                res[bit / 64] |= (static_cast<u64>(1) << (bit % 64));
            }
        }

        return res;
    }
    ulong& operator /= (const ulong& div) {
        return *this = *this / div;
    }

    // O(len * 64)
    ulong operator % (const ulong& div) const {
        int bit = 0; // сдвиг не нарушающий целостность div

        for (int i = 64 * len - 1; i >= 0; i--) {
            if (div[i / 64] & (static_cast<u64>(1) << i % 64)) {
                bit = 64 * len - i - 1;
                break;
            }
        }

        auto num = *this;
        for (; bit >= 0; bit--) {
            if ((div << bit) <= num) {
                num -= (div << bit);
            }
        }

        return num;
    }
    ulong& operator %= (const ulong& div) {
        return *this = *this % div;
    }


    ulong& operator &= (const ulong& val) {
        for (int i = 0; i < len; i++) {
            words[i] &= val[i];
        }
        return *this;
    }
    ulong operator & (const ulong& val) const {
        return ulong(*this) &= val;
    }

    ulong& operator |= (const ulong& val) {
        for (int i = 0; i < len; i++) {
            words[i] |= val[i];
        }
        return *this;
    }
    ulong operator | (const ulong& val) const {
        return ulong(*this) |= val;
    }

    ulong& operator ^= (const ulong& val) {
        for (int i = 0; i < len; i++) {
            words[i] ^= val[i];
        }
        return *this;
    }
    ulong operator ^ (const ulong& val) const {
        return ulong(*this) ^= val;
    }

    ulong operator ~ () const {
        ulong res;
        for (int i = 0; i < len; i++) {
            res[i] = ~words[i];
        }
        return res;
    }

    ulong& operator ++ () {
        *this += 1;
        return *this;
    }
    ulong operator ++ (int) {
        ulong temp = *this;
        ++(*this);
        return temp;
    }

    ulong& operator -- () {
        *this -= 1;
        return *this;
    }
    ulong operator -- (int) {
        ulong temp = *this;
        --(*this);
        return temp;
    }
};
