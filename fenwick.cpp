// fenwick
// Developed by Mob

#include<vector>
using namespace std;

// Дерево Фенвика
template<typename T>
class fenwick {
    std::vector<T> t;
    int n;

public:
    fenwick() {}
    fenwick(int size) {
        n = size + 1;
        t.resize(n, 0);
    }
    fenwick(const std::vector<T>& array) {
        n = array.size() + 1;
        t.resize(n, 0);
        for (int i = 0; i < array.size(); i++) {
            update(i, array[i]);
        }
    }

    T get(int r) {
        r++;
        T res = 0;
        for (; r > 0; r -= r & -r) {
            res += t[r];
        }
        return res;
    }

    T get(int l, int r) {
        return get(r) - get(l - 1);
    }

    void update(int k, T val) {
        k++;
        for (; k < n; k += k & -k) {
            t[k] += val;
        }
    }
};
