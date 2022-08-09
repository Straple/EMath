#include <bits/stdc++.h>
using namespace std;

using s64 = long long;
using u64 = unsigned long long;

//#define DEBUG

const long double eps = 1e-7;

int main() {
#ifdef DEBUG
    ifstream cin("input.txt");
#endif
    ios::sync_with_stdio(false), cout.tie(0), cin.tie(0);

    int n;
    cin >> n;

    vector<vector<long double>> A(n, vector<long double>(n + 1));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            cin >> A[i][j];
        }
    }

    vector<int> where(n, -1);
    for (int col = 0, row = 0; col < n && row < n; col++) {

        int sel = row;
        for (int i = row; i < n; i++) {
            if (abs(A[i][col]) > abs(A[sel][col])) {
                sel = i;
            }
        }

        if (abs(A[sel][col]) < eps) {
            continue;
        }

        for (int i = col; i <= n; i++) {
            swap(A[sel][i], A[row][i]);
        }
        where[col] = row;

        for (int i = 0; i < n; i++) {
            if (i != row) {
                long double c = A[i][col] / A[row][col];
                for (int j = col; j <= n; j++) {
                    A[i][j] -= A[row][j] * c;
                }
            }
        }

        row++;
    }

    vector<long double> ans(n, 0);

    for (int i = 0; i < n; i++) {
        if (where[i] != -1) {
            ans[i] = A[where[i]][n] / A[where[i]][i];
        }
    }
    for (int i = 0; i < n; i++) {

        long double sum = 0;
        for (int j = 0; j < n; j++) {
            sum += ans[j] * A[i][j];
        }

        if (abs(sum - A[i][n]) > eps) {
            cout << "impossible";
            return 0;
        }
    }

    for (int i = 0; i < n; i++) {
        if (where[i] == -1) {
            cout << "infinity";
            return 0;
        }
    }

    cout << "single\n";

    cout << fixed << setprecision(3);
    for (int i = 0; i < n; i++) {
        cout << ans[i] << " ";
    }



    return 0;
}