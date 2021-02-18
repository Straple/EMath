// sparse table
// Developed by Mob

#include<vector>
using namespace std;

class sparse_table {

    vector<int> logs;

    vector<vector<int>> table;

    void computeLogs(int n) {
        logs.resize(n + 1);
        logs[1] = 0;
        for (s64 i = 2; i <= n; i++) {
            logs[i] = logs[i / 2] + 1;
        }
    }

    void buildSparseTable(int n, vector<int>& A) {
        for (int i = 0; i <= logs[n]; i++) {
            int curLen = 1 << i; // 2^i
            for (u32 j = 0; j + curLen <= n; j++) {
                if (curLen == 1) {
                    table[i][j] = A[j];
                }
                else {
                    table[i][j] = max(table[i - 1][j], table[i - 1][j + (curLen / 2)]);
                }
            }
        }
    }

public:

    sparse_table(vector<int>& A) {
        computeLogs(A.size());
        table.resize(logs[A.size()] + 1, vector<int>(A.size() + 1));
        buildSparseTable(A.size(), A);
    }

    int get(int l, int r) {
        int p = logs[r - l + 1];
        int pLen = 1 << p;
        return max(table[p][l], table[p][r - pLen + 1]);
    }
};
