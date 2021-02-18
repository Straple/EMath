// gcd_machine
// Developed by Mob

#include<vector>
using namespace std;


class gcd_machine {
    vector<int> primes; // простые числа, по возрастанию
    vector<int> d; // d[x] - индекс минимального делителя числа x
    vector<int> rest; // rest[x] - x / p1^a1
    vector<int> term; // term[x] - p1^a1

    typedef uint64_t u64;

    // обновляет информацию для x
    inline void update(u64 x) {
        u64 p = primes[d[x]]; // наименьший простой делитель
        u64 y = x / p;

        rest[x] = (d[y] == d[x] ? rest[y] : y);

        term[x] = x / rest[x];
    }

    // O(n)
    void lineSieve(u64 n) {
        for (u64 x = 2; x <= n; x++) {
            if (d[x] == -1) { // новое простое
                d[x] = primes.size();
                primes.push_back(x);
            }
            for (u64 i = 0; i <= d[x] && x * primes[i] <= n; i++) {
                d[x * primes[i]] = i;
            }

            update(x); // обновить доп информацию
        }
    }

public:

    gcd_machine(u64 n) {
        d.resize(n + 1, -1);
        rest.resize(n + 1, 0);
        term.resize(n + 1, 0);
        lineSieve(n);
    }

    u64 gcd(u64 a, u64 b) {
        if (a > b) {
            swap(a, b);
        }
        if (a == 0) {
            return b;
        }
        else {
            u64 res = 1;
            while (a > 1 && b > 1) {
                if (primes[d[a]] < primes[d[b]]) {
                    a = rest[a];
                }
                else if (primes[d[a]] > primes[d[b]]) {
                    b = rest[b];
                }
                else {
                    u64 g = min(term[a], term[b]);
                    res *= g;
                    a /= g;
                    b /= g;
                }
            }
            return res;
        }
    }
};

