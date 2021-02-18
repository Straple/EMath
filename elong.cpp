// Long Arithmetic
// Developed by Mob


#include <iostream>
#include <iomanip>


#include "utils.cpp"
#include "edeque.cpp"

const s64 long_length = 18,
    long_base = epow<s64>(10, long_length),
    long_base_expansion = epow<s64>(10, long_length / 2);



// signed long integer
class elong {
    // unsigned basic long
    struct basicLong {
        edeque<s64> digits; // цифры

        basicLong() {}
        basicLong(const edeque<s64>& digits) {
            this->digits = digits;
        }
        basicLong(const std::string& str) {
            int i;
            for (i = str.size(); i >= long_length; i -= long_length) {
                digits.push_back(atoll(str.substr(i - long_length, long_length).c_str()));
            }
            if (i > 0) {
                digits.push_back(atoll(str.substr(0, i).c_str()));
            }
        }

        // удаляет ведущие нули
        void remove_leading_zeros() {
            while (!digits.empty() && digits.back() == 0) {
                digits.pop_back();
            }
        }

    private:
        enum class OP {
            less,
            equally,
            more,
        };

#define compareValues(a, b) (a < b ? OP::less : OP::more)

        // сравнивает 2 числа
        OP compare(const basicLong& Rhs) const {
            // если у них разные длины
            if (digits.size() != Rhs.digits.size()) {
                return compareValues(digits.size(), Rhs.digits.size());
            }
            else {
                int i = digits.size() - 1;
                while (i >= 0 && digits[i] == Rhs.digits[i]) {
                    i--;
                }
                return i >= 0 ? compareValues(digits[i], Rhs.digits[i]) : OP::equally;
            }
        }

    public:

        bool operator < (const basicLong& Rhs) const {
            return compare(Rhs) == OP::less;
        }
        bool operator == (const basicLong& Rhs) const {
            return compare(Rhs) == OP::equally;
        }
        bool operator > (const basicLong& Rhs) const {
            return compare(Rhs) == OP::more;
        }

        // операция сложения
        basicLong operator + (const basicLong& added) const {
            basicLong result = *this;
            bool k = 0;
            s64 i = 0;
            s64 len = std::max(result.digits.size(), added.digits.size());
            for (i = 0; i < len || k != 0; i++) {
                if (i == result.digits.size()) {
                    result.digits.push_back(0);
                }
                result.digits[i] += k + (i < added.digits.size() ? added.digits[i] : 0);
                k = result.digits[i] >= long_base;
                result.digits[i] -= k * long_base;
            }
            return result;
        }

        // операция вычитания
        basicLong operator - (const basicLong& subtrahend) const {
            basicLong result = *this;
            bool k = 0;
            for (int i = 0; i < subtrahend.digits.size() || k != 0; i++) {
                result.digits[i] -= k + (i < subtrahend.digits.size() ? subtrahend.digits[i] : 0);
                k = result.digits[i] < 0;
                result.digits[i] += k * long_base;
            }

            result.remove_leading_zeros();
            return result;
        }

        // расширяет длину числа, уменьшая его модуль
        basicLong expansion() const {
            basicLong result;
            result.digits.resize(digits.size() << 1);
            for (int i = 0; i < digits.size(); i++) {
                result.digits[i << 1] = digits[i] % long_base_expansion;
                result.digits[(i << 1) + 1] = digits[i] / long_base_expansion;
            }
            return result;
        }

        // сужает длину числа, увеличивая его модуль
        basicLong reduction() const {
            basicLong result;
            result.digits.resize(digits.size() >> 1);
            for (int i = 0; i < result.digits.size(); i++) {
                result.digits[i] = digits[i << 1] + digits[(i << 1) + 1] * long_base_expansion;
            }
            return result;
        }

        // Простое умножение двух чисел за O(n^2)
        basicLong naiveMul(const basicLong& a, const basicLong& b) const {
            basicLong result, mult1 = a.expansion(), mult2 = b.expansion();
            result.digits.resize(mult1.digits.size() + mult2.digits.size() + 1);

            u64 c, k;
            int i, j;
            for (i = 0; i < mult1.digits.size(); i++) {
                c = 0;
                for (j = 0; j < mult2.digits.size() || c != 0; j++) {
                    // mul1[i] < 1e9 && mult2[j] < 1e9
                    // mult1[i] * mult2[j] < 1e18
                    
                    k = result.digits[i + j];
                    k += c;
                    if (j < mult2.digits.size()) {
                        k += mult1.digits[i] * 1ULL * mult2.digits[j];
                    }

                    c = k / long_base_expansion;
                    result.digits[i + j] = k - c * long_base_expansion;
                }
            }

            result = result.reduction();
            result.remove_leading_zeros();
            return result;
        }

        // O(n ^ 1.54)
        basicLong KaratsubaMul(const basicLong& x, const basicLong& y, int n) const {
            if (n <= 128) {
                return naiveMul(x, y);
            }
            else {
                int half = n >> 1;
                // firstHalf - первая половина числа
                basicLong Xl = x.firstHalf(), Xr = x.secondHalf(), Yl = y.firstHalf(), Yr = y.secondHalf();
                basicLong sumX = Xl + Xr, sumY = Yl + Yr;

                basicLong P1 = KaratsubaMul(Xl, Yl, half);
                basicLong P2 = KaratsubaMul(Xr, Yr, half);
                basicLong P3;
                if (sumX.digits.size() != half || sumY.digits.size() != half) {
                    sumX.digits.resize(n);
                    sumY.digits.resize(n);
                    P3 = KaratsubaMul(sumX, sumY, n);
                }
                else {
                    P3 = KaratsubaMul(sumX, sumY, half);
                }
                // P1 * (base ^ n) + (P3 - P1 - P2) * (base ^ half) + P2

                auto result = (P1 + ((P3 - P1 - P2) << half) + (P2 << n));
                result.remove_leading_zeros();
                return result;
            }
        }

        basicLong operator * (const basicLong& mult) const {
            int k = std::max(this->digits.size(), mult.digits.size());
            if (k <= 128) {
                return naiveMul(*this, mult);
            }
            else {
                basicLong mult1 = *this, mult2 = mult;
                int len = roundTwo(k);
                mult1.digits.resize(len);
                mult2.digits.resize(len);

                return KaratsubaMul(mult1, mult2, len);
            }
        }

        // ( /, % )
        std::pair<basicLong, basicLong> div(const basicLong& divider) const {
            basicLong result, value, temp, t;

            int i = digits.size() - 1;
            while (i >= 0 && value < divider) {
                value.digits.push_front(digits[i]);
                i--;
            }
            value.digits.pop_front();
            i++;
            for (; i >= 0; i--) {
                value.remove_leading_zeros();
                value.digits.push_front(digits[i]);

                long long l = 0, r = long_base;
                while (l < r - 1) {
                    long long m = (l + r) >> 1;

                    temp.digits.push_back(m);
                    t = divider * temp;
                    temp.digits.pop_back();

                    if (value < t) {
                        r = m;
                    }
                    else {
                        l = m;
                    }
                }

                result.digits.push_front(l);
                temp.digits.push_back(l);
                value = value - divider * temp;
                temp.digits.pop_back();
            }
            return std::make_pair(result, value);
        }

        basicLong operator / (const basicLong& divider) const {
            return div(divider).first;
        }

        basicLong operator % (const basicLong& divider) const {
            return div(divider).second;
        }


        // сдвигает число на long_base^k влево
        basicLong operator << (int k) const {
            basicLong res = *this;
            while (k--) {
                res.digits.push_front(0);
            }
            return res;
        }

        basicLong firstHalf() const {
            int half = digits.size() >> 1;
            basicLong res;
            res.digits.resize(half);
            for (int i = 0; i < half; i++) {
                res.digits[i] = digits[i];
            }
            return res;
        }
        basicLong secondHalf() const {
            int half = digits.size() >> 1;
            basicLong res;
            res.digits.resize(half);
            for (int i = 0; i < half; i++) {
                res.digits[i] = digits[i + half];
            }
            return res;
        }
    };

    basicLong value;
    bool isNegative;

    // update
    elong& update() {
        if (value.digits.empty()) {
            value.digits.push_back(0);
        }
        return *this;
    }

    elong(const basicLong& value, bool isNegative) {
        this->value = value;
        this->isNegative = isNegative;
    }
    elong(basicLong&& value, bool isNegative) {
        this->value = std::move(value);
        this->isNegative = isNegative;
    }

    // turbo ++ and --

    void inc() {
        value.digits.front()++;
        if (value.digits.front() == long_base) {
            value.digits.front() = 0;

            int i;
            for (i = 1; i < value.digits.size() && value.digits[i] + 1 == long_base; i++) {
                value.digits[i] = 0;
            }
            if (i == value.digits.size()) {
                value.digits.push_back(1);
            }
            else {
                value.digits[i]++;
            }
        }
    }
    void dec() {
        value.digits.front()--;
        if (value.digits.front() == -1) {
            value.digits.front() = 0;

            long long i;
            for (i = 1; i < value.digits.size() && value.digits[i] == 0; i++) {
                value.digits[i] = long_base - 1;
            }

            if (i < value.digits.size()) {
                value.digits[i]--;
                if (value.digits.back() == 0) {
                    value.digits.pop_back();
                }
            }
            else {
                isNegative = !isNegative;
                value.digits[0] = 1;
            }
        }

        if (value.digits.size() == 1 && value.digits.front() == 0) {
            isNegative = false;
        }
    }

    void signedInc() {
        if (isNegative) {
            dec();
        }
        else {
            inc();
        }
    }
    void signedDec() {
        if (isNegative) {
            inc();
        }
        else {
            dec();
        }
    }

public:

    // default constructor
    elong() {
        value.digits.push_back(0);
        isNegative = false;
    }

    // convert constructors[

    elong(const std::string& String) {
        isNegative = String[0] == '-';
        value = isNegative ? basicLong(String.substr(1)) : basicLong(String);
    }
    elong(const char* String) {
        *this = elong(std::string(String));
    }

    template<typename T>
    elong(T number) {
        isNegative = number < 0;
        number *= isNegative ? -1 : 1;

        value.digits.push_back(number % long_base);
        number /= long_base;
        if (number > 0) {
            value.digits.push_back(number % long_base);
            number /= long_base;
            if (number > 0) {
                value.digits.push_back(number);
            }
        }
    }

    // ]

    //bool operators [

    bool operator == (const elong& Rhs) const {
        return isNegative != Rhs.isNegative ? false : value == Rhs.value;
    }
    bool operator != (const elong& Rhs) const {
        return !(*this == Rhs);
    }
    bool operator < (const elong& Rhs) const {
        if (isNegative) {
            return Rhs.isNegative ? -*this > -Rhs : true;
        }
        else {
            return Rhs.isNegative ? false : value < Rhs.value;
        }
    }
    bool operator > (const elong& Rhs) const {
        if (isNegative) {
            return Rhs.isNegative ? -*this < -Rhs : false;
        }
        else {
            return Rhs.isNegative ? true : value > Rhs.value;
        }
    }

    bool operator <= (const elong& Rhs) const {
        return !(*this > Rhs);
    }
    bool operator >= (const elong& Rhs) const {
        return !(*this < Rhs);
    }
    // ]

    elong operator - () const {
        return elong(value, !isNegative);
    }

    // operators [

    elong& operator ++() {
        signedInc();
        return *this;
    }
    elong operator ++(int) {
        auto temp = *this;
        signedInc();
        return temp;
    }
    elong& operator --() {
        signedDec();
        return *this;
    }
    elong operator --(int) {
        auto temp = *this;
        signedDec();
        return temp;
    }

    elong operator + (const elong& added) const {
        if (isNegative) {
            return added.isNegative ? -(-*this + (-added)) : added - (-*this);
        }
        else {
            return added.isNegative ? *this - (-added) : elong(std::move(value + added.value), false);
        }
    }
    elong& operator += (const elong& added) {
        return *this = *this + added;
    }

    elong operator - (const elong& subtrahend) const {
        if (subtrahend.isNegative) {
            return *this + (-subtrahend);
        }
        else if (isNegative) {
            return -(-*this + subtrahend);
        }
        else {
            return *this < subtrahend ? -(subtrahend - *this) : elong(std::move(value - subtrahend.value), false).update();
        }
    }
    elong& operator -= (const elong& subtrahend) {
        return *this = *this - subtrahend;
    }

    elong operator * (const elong& multiplied) const {
        return elong(std::move(value * multiplied.value), isNegative != multiplied.isNegative).update();
    }
    elong& operator *= (const elong& multiplied) {
        return *this = *this * multiplied;
    }

    elong operator / (const elong& divider) const {
        return elong(std::move(value / divider.value), isNegative != divider.isNegative).update();
    }
    elong& operator /= (const elong& divider) {
        return *this = *this / divider;
    }

    elong operator % (const elong& divider) const {
        return elong(std::move(value % divider.value), isNegative).update();
    }
    elong& operator %= (const elong& divider) {
        return *this = *this % divider;
    }

    // friend function :)
    friend std::ostream& operator << (std::ostream& out, const elong& var);

    elong mult_by_2() {
        basicLong result;
        result.digits.resize(value.digits.size() + 1);

        u64 c = 0, k;
        int i;
        for (i = 0; i < value.digits.size(); i++) {
            k = c;
            k += value.digits[i] * 2ULL;
            c = k / long_base;

            result.digits[i] = k - c * long_base;
        }
        if (c != 0) {
            result.digits.back() += c;
        }

        result.remove_leading_zeros();
        elong va;
        va.value = std::move(result);
        va.isNegative = isNegative;
        return va;
    }
};

template<typename T>
elong operator + (const T& a, const elong& b) {
    return elong(a) + b;
}
template<typename T>
elong operator - (const T& a, const elong& b) {
    return elong(a) - b;
}
template<typename T>
elong operator * (const T& a, const elong& b) {
    return elong(a) * b;
}
template<typename T>
elong operator / (const T& a, const elong& b) {
    return elong(a) / b;
}
template<typename T>
elong operator % (const T& a, const elong& b) {
    return elong(a) % b;
}


// output signed long integer
std::ostream& operator << (std::ostream& out, const elong& var) {
    if (var.isNegative) {
        out << '-';
    }
    out << var.value.digits.back();
    char k = out.fill('0');
    for (int i = var.value.digits.size() - 2; i >= 0; i--) {
        out << std::setw(long_length) << var.value.digits[i];
    }
    out.fill(k);
    return out;
}
// input signed long integer
std::istream& operator >> (std::istream& input, elong& var) {
    std::string str;
    input >> str;
    var = elong(str);
    return input;
}
