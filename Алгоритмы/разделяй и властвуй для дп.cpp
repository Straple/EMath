#include <bits/stdc++.h>
#include <unordered_set>
#include <unordered_map>
#include <random>
using namespace std;

mt19937 rnd(time(NULL));

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

using s8 = int8_t;
using s16 = int16_t;
using s32 = int32_t;
using s64 = int64_t;

using ld = long double;

/*
* GLOBAL PROJECT SETTINGS
*/

#define HOME_MODE

//#define DEBUG
#define MULTITEST
#define STRESS_TEST // включает HOME_MODE, MULTITEST и тестирует решения

//#define TASK_FILE_NAME "supreme" // если есть, включает файловый ввод и вывод

/*
* SIMULATE SETTINGS
*/

#ifdef TASK_FILE_NAME
#define FILE_INPUT TASK_FILE_NAME".in"
#define FILE_OUTPUT TASK_FILE_NAME".out"
#endif

#ifdef STRESS_TEST
#define HOME_MODE // стресс-тест бывает только дома
#define MULTITEST // у нас есть несколько тестов для тестирования
#endif

#ifdef HOME_MODE
#undef FILE_INPUT
#undef FILE_OUTPUT
#define FILE_INPUT "input.txt"
#else
#undef DEBUG
#endif

/*
* SOLVE PROBLEM
*/

// даны n точек на прямой
// нужно не больше m отрезков, покрывающих [min x, max x],
// минимизировав при этом сумму квадратов их длин

// dp[i][j] - минимальная сумма квадратов длин j отрезков, покрывающих точки из x на префиксе i
vector<vector<int>> dp;

vector<int> x;

int cost(int i, int j) {
    return (x[i] - x[j]) * (x[i] - x[j]);
}

// разделяй и властвуй
// строим dp для i e[l, r]
// для этого будем рассматривать x'ы из [_l, _r]
// они являются оптимальными для данного отрезка
void slv(int l, int r, int _l, int _r, int k) {
    if (l > r) {
        return;
    }
    int t = (l + r) / 2, opt = _l;
    for (int i = _l; i < min(_r + 1, t); i++) {
        int val = dp[i][k - 1] + cost(i, t);
        if (val < dp[t][k]) {
            dp[t][k] = val;
            opt = i;
        }
    }

    slv(l, t - 1, _l, opt, k);
    slv(t + 1, r, opt, _r, k);
}

// O(nm * log n)
void solve(istream& cin, ostream& cout) {
    int n, m;
    cin >> n >> m;
    x.resize(n);
    for (int i = 0; i < n; i++) {
        cin >> x[i];
    }
    sort(x.begin(), x.end());

    dp = vector<vector<int>>(n + 1, vector<int>(m + 1, 1e9));
    dp[0][0] = 0;

    for (int k = 1; k <= m; k++) {
        slv(0, n - 1, 0, n - 1, k);
    }

    int res = 1e9;
    for (int j = 0; j <= m; j++) {
        res = min(res, dp[n - 1][j]);
    }
    cout << res << "\n";
}


void main_build() {

}

// O(n^2 * m)
void correct_solve(istream& cin, ostream& cout) {
    int n, m;
    cin >> n >> m;
    vector<int> x(n);
    for (int i = 0; i < n; i++) {
        cin >> x[i];
    }
    sort(x.begin(), x.end());

    // dp[i][j] - минимальная сумма квадратов длин j отрезков, покрывающих точки из x на префиксе i
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 1e9));
    dp[0][0] = 0;

    // O(n^2 * m)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int to = i + 1; to < n; to++) {
                // покрыть отрезком [xi, xto]
                dp[to][j + 1] = min(dp[to][j + 1], dp[i][j] + (x[to] - x[i]) * (x[to] - x[i]));
            }
        }
    }

    int res = 1e9;
    for (int j = 0; j <= m; j++) {
        res = min(res, dp[n - 1][j]);
    }
    cout << res << "\n";
}

#ifdef STRESS_TEST

// строит тест
string build_test() {
    stringstream cout;

    // вывести тест
    int n = rnd() % 100 + 1;
    int k = rnd() % n + 1;
    cout << n << " " << k << "\n";
    for (int i = 0; i < n; i++) {
        cout << (rnd() % 100) << " ";
    }

    return cout.str();
}

#endif


/*
* TEMPLATE MAIN
*/

int main() {
    //setlocale(LC_ALL, "Russian");
#ifdef FILE_INPUT
    ifstream cin(FILE_INPUT);
#endif
#ifdef FILE_OUTPUT
    ofstream cout(FILE_OUTPUT);
#endif
    ios::sync_with_stdio(false), cout.tie(nullptr), cin.tie(nullptr);

    main_build();

    s64 t = 1;
#ifdef MULTITEST
    cin >> t;
#endif

    while (t--) {
#ifdef STRESS_TEST
        string test_str_input = build_test(); // создаем тест

        stringstream input1(test_str_input);
        stringstream input2(test_str_input);

        // получаем ответы на тест
        stringstream wrong_out, correct_out;
        solve(input1, wrong_out);
        correct_solve(input2, correct_out);

        // сравниваем их
        if (wrong_out.str() != correct_out.str()) {
            cout << "[failed test]\n" << test_str_input << "\n";

            cout << "[wrong output]\n" << wrong_out.str() << "\n";

            cout << "[correct output]\n" << correct_out.str() << "\n";

            return 0;
        }
#else
        solve(cin, cout);
#endif
    }

#ifdef STRESS_TEST
    cout << "[all tests passed]\n";
#endif

    return 0;
}
