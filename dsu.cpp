// dsu
// Developed by Mob

#include<vector>
using namespace std;

// Система непересекающихся множеств.
// Сжатие пути + ранговая эвристика на основе глубины деревьев
class dsu {
    std::vector<int> parent, rank;
public:

    dsu() {}
    dsu(int length) {
        parent.resize(length);
        rank.resize(length);
        for (int i = 0; i < length; i++) {
            makeSet(i);
        }
    }

    void makeSet(int v) {
        parent[v] = v;
        rank[v] = 0;
    }
    int findSet(int v) {
        return v == parent[v] ? v : (parent[v] = findSet(parent[v]));
    }
    void unionSets(int a, int b) {
        a = findSet(a);
        b = findSet(b);
        if (a != b) {
            if (rank[a] < rank[b]) {
                std::swap(a, b);
            }
            parent[b] = a;
            if (rank[a] == rank[b]) {
                rank[a]++;
            }
        }
    }
};
