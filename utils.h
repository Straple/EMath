#pragma once

// utils
// Developed by Mob

#include<cstdint>

#include<string>

#include<algorithm>

typedef int8_t    s8;
typedef uint8_t   u8;
typedef int16_t   s16;
typedef uint16_t  u16;
typedef int32_t   s32;
typedef uint32_t  u32;
typedef int64_t   s64;
typedef uint64_t  u64;

using point_t = double;

// utils

template<typename T>
T clamp(const T& min, const T& val, const T& max) {
    _STL_VERIFY(min <= max, "WRONG BORDER");

    if (val < min) {
        return min;
    }
    else if (val > max) {
        return max;
    }
    else {
        return val;
    }
}

template<typename T>
bool is_between(const T& min, const T& val, const T& max) {
    _STL_VERIFY(min <= max, "WRONG BORDER");

    return min <= val && val <= max;
}

// u64 -> string
std::string cast(u64 num) {
    std::string str;

    if (num == 0) {
        str = "0";
    }
    else {
        while (num) {
            str += num % 10 + '0';
            num /= 10;
        }
        std::reverse(str.begin(), str.end());
    }
    return str;
}

// s64 -> string
std::string cast(s64 num) {

    std::string str = cast(static_cast<u64>(abs(num)));

    if (num < 0) {
        str.insert(str.begin(), '-');
    }
    return str;
}

// point_t -> string
std::string cast(point_t num, u32 precision) {
    std::string str;

    str = cast(static_cast<s64>(num));

    {
        num = abs(num);

        long double tmp;
        num = std::modfl(num, &tmp);
    }
    

    str += ".";

    while (num > 0 && precision) {
        precision--;
        num *= 10;

        str += trunc(num) + '0';

        num = num - trunc(num);
    }

    while (str.back() == '0') {
        str.pop_back();
    }

    if (str.back() == '.') {
        str.pop_back();
    }

    return str;
}


// algorithms


// округление вверъ до степени двойки
u64 roundTwo(u64 v) {
    v--;

    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;

    v++;
    return v;
}

template<typename T1, typename T2>
T1 epow(const T1& a, const T2& n) {
    if (n == 0) {
        return 1;
    }
    else {
        T1 z = epow(a, n / 2);
        z *= z;

        if (n % 2 == 1) {
            z *= a;
        }
        return z;
    }
}

template<typename T1, typename T2>
T1 epow(const T1& a, const T2& n, const T1& mod) {
    if (n == 0) {
        return 1;
    }
    else {
        T1 z = epow(a, n / 2);
        z = (z * z) % mod;

        if (n % 2 == 1) {
            z = (z * a) % mod;
        }
        return z;
    }
}

template<typename T>
T gcd(T a, T b) {
    if (a > b) {
        swap(a, b);
    }
    while (a != 0) {
        b %= a;
        swap(a, b);
    }
    return b;
}

template<typename T>
T lcm(const T& a, const T& b) {
    return (a / gcd(a, b)) * b;
}

// Расширенный Алгоритм Евклида
template<typename T>
T gcd(T a, T b, T& x, T& y) {
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }
    else {
        T x1, y1;
        T res = gcd(b % a, a, x1, y1);
        x = y1 - (b / a) * x1;
        y = x1;
        return res;
    }
}


// Линейные диофантовы уравнения с двумя переменными
// a * x + b * y = c
// a * xg + b * yg = g
template<typename T>
bool solveLDE(T a, T b, T c, T& xg, T& yg, T& g) {
    g = gcd(abs(a), abs(b), xg, yg);
    if (c % g != 0) {
        return false;
    }
    xg *= c / g;
    yg *= c / g;
    xg *= a < 0 ? -1 : 1;
    yg *= b < 0 ? -1 : 1;
    return true;
}

// Обратный элемент в кольце по модулю
// (a ^ -1) mod m = x
template<typename T>
bool reverseElement(T a, T m, T& x) {
    T y;
    T g = gcd(a, m, x, y);
    if (g != 1) {
        return false;
    }
    else {
        x = (x % m + m) % m;
        return true;
    }
}

// Метод Ньютона для поиска целочисленных корней
template<typename T>
T sqrtn(const T& n) {
    T x = 1;
    bool decreased = false;
    while (true) {
        T nx = (x + n / x) / 2;
        if (x == nx || nx > x && decreased) {
            break;
        }
        decreased = nx < x;
        x = nx;
    }
    return x;
}
