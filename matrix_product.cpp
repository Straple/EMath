// matrix product
// Developed by Mob

#include<vector>
using namespace std;

// O(n^3)
template<typename T>
vector<vector<T>> matrix_product(const vector<vector<T>>& a, const vector<vector<T>>& b) {
    vector<vector<T>> c(a.size(), vector<T>(b[0].size()));

    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < b[0].size(); j++) {
            for (int k = 0; k < b.size(); k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return c;
}
