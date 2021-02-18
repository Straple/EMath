// prime sieve
// Developed by Mob

#include <vector>
using namespace std;

// O(n)
std::vector<int> lineSieve(long long n) {
    std::vector<int> primes, d(n + 1, -1);

    for (long long y = 2; y <= n; y++) {
        if (d[y] == -1) { // новое простое
            d[y] = primes.size();
            primes.push_back(y);
        }
        for (long long i = 0; i <= d[y] && y * primes[i] <= n; i++) {
            d[y * primes[i]] = i;
        }
    }

    return primes;
}

// дополнительная информация для линейного решета Эратосфена
struct sieve_info {
    std::vector<int>    primes; // primes[i] - простой делитель
    std::vector<int>         d; // d[x] - индекс минимального простого делителя числа x в массиве primes
    std::vector<int>       deg; // deg[x] - степень вхождения минимального простого делителя числа x
    std::vector<int>      rest; // rest[x] - x / p1^a1 - остаток числа x от p1^a1
    std::vector<int>      term; // term[x] - p1^a1
    std::vector<int>       phi; // количество взаимно простых с n
    std::vector<long long> div; // div[x] - количество делителей x
    std::vector<long long>  s1; // s1[x] - сумма делителей x

    sieve_info() {}
    sieve_info(long long n) {
        d.resize(n + 1, -1);

        deg.resize(n + 1, 0);
        rest.resize(n + 1, 0);
        term.resize(n + 1, 0);

        phi.resize(n + 1, 1);
        div.resize(n + 1, 1);

        s1.resize(n + 1, 0);
        s1[1] = 1;
    }

    // обновляет информацию для x
    void update(long long x) {
        long long p_x = primes[d[x]];
        long long y_x = x / p_x;

        deg[x] = d[y_x] == d[x] ? deg[y_x] + 1 : 1;

        rest[x] = d[y_x] == d[x] ? rest[y_x] : y_x;

        term[x] = x / rest[x];

        phi[x] = phi[rest[x]] * (term[x] / p_x) * (p_x - 1);

        div[x] = div[rest[x]] * (static_cast<long long>(deg[x]) + 1);

        s1[x] = s1[rest[x]] * ((p_x * term[x] - 1) / (p_x - 1));
    }
};

// O(n) + доп информация
sieve_info lineSieveInfo(long long n) {
    sieve_info info(n);

    for (long long y = 2; y <= n; y++) {
        if (info.d[y] == -1) { // новое простое
            info.d[y] = info.primes.size();
            info.primes.push_back(y);
        }
        for (long long i = 0; i <= info.d[y] && y * info.primes[i] <= n; i++) {
            info.d[y * info.primes[i]] = i;
        }

        info.update(y); // обновление доп информации
    }

    return info;
}
