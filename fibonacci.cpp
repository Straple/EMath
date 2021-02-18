// fibonacci
// Developed by Mob

#include <vector>
#include "matrix_product.cpp"
using namespace std;


// Фибоначчи за log
template<typename T>
class fibonacci {
    typedef vector<vector<T>> matrix;

    matrix binPow(const matrix& a, const T& n) {
        if (n == 1) {
            return a;
        }
        else {
            auto z = binPow(a, n / 2);
            z *= z;

            return n % 2 == 0 ? z : z * a;
        }
    }
    matrix binPow(const matrix& a, const T& n, const T& mod) {
        if (n == 1) {
            return a;
        }
        else {
            auto z = binPow(a, n / 2);
            z = (z * z) % mod;

            return n % 2 == 0 ? z : (z * a) % mod;
        }
    }

    matrix getMatrix() const {
        return { {0, 1}, {1, 1} };
    }

public:
    T get(const T& n) {
        return binPow(getMatrix(), n)[0][0];
    }
    T get(const T& n, const T& mod) {
        return binPow(getMatrix(), n, mod)[0][0];
    }
};
