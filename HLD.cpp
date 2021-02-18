// heavy light decomposition
// Developed by Mob

#include<vector>
using namespace std;

class hld {

    class seg_tree {
        int len;
        vector<int> T;

        int get(int v, int tl, int tr, int l, int r) {
            if (tl == l && tr == r) {
                return T[v];
            }
            else {
                int tm = (tl + tr) / 2;
                if (r <= tm) {
                    return get(v * 2 + 1, tl, tm, l, r);
                }
                else if (tm < l) {
                    return get(v * 2 + 2, tm + 1, tr, l, r);
                }
                else {
                    return max(get(v * 2 + 1, tl, tm, l, tm), get(v * 2 + 2, tm + 1, tr, tm + 1, r));
                }
            }
        }

        void update(int v, int tl, int tr, int pos, int add) {
            if (tl == tr) {
                T[v] += add;
            }
            else {
                int tm = (tl + tr) / 2;
                if (pos <= tm) {
                    update(v * 2 + 1, tl, tm, pos, add);
                }
                else {
                    update(v * 2 + 2, tm + 1, tr, pos, add);
                }
                T[v] = max(T[v * 2 + 1], T[v * 2 + 2]);
            }
        }

    public:

        seg_tree() {
            len = 0;
        }
        seg_tree(int n) {
            len = n;
            T.resize(4 * n);
        }

        int get(int l, int r) {
            return get(0, 0, len - 1, l, r);
        }

        void update(int pos, int add) {
            update(0, 0, len - 1, pos, add);
        }
    } Tree; // ДО


    vector<int> Size; // Size[v] = размер поддерева v
    vector<int> Parent; // Parent[v] = предок вершины v
    vector<int> Next; // Next[v] = сын для следущий вершины ДО, если Next[v] == -1, то она не входит в ДО
    vector<int> Depth; // Depth[v] = глубина вершины v
    vector<int> Num; // Num[v] = номер вершины в порядке "умного" обхода
    vector<int> Chain; // Chain[v] = номер цепочки, в которую входит вершина

    vector<int> Top; // Top[k] = верхняя вершина у цепочки k

    void dfs(int v, int prev, vector<vector<int>>& G) {
        Parent[v] = prev; // parent
        Depth[v] = Depth[prev] + 1; // depth

        // size, next
        Size[v] = 1;
        Next[v] = -1;
        for (auto& to : G[v]) {
            if (to != prev) {
                dfs(to, v, G); // dfs
                Size[v] += Size[to]; // size

                // next
                if (Next[v] == -1 || Size[to] > Size[Next[v]]) {
                    Next[v] = to; // берем самое тяжелое ребро
                }
            }
        }
    }

    int chain_cnt = 0;
    int num_cnt = 0;

    void build_hld(int v, int prev, vector<vector<int>>& G) {
        Chain[v] = chain_cnt; // chain
        Num[v] = num_cnt++; // num

        // top
        if (Top[chain_cnt] == -1) {
            Top[chain_cnt] = v;
        }

        // next
        if (Next[v] != -1) { // Эта вершина входит в ДО
            build_hld(Next[v], v, G);
        }

        for (auto& to : G[v]) {
            if (to != Next[v] && to != prev) {
                chain_cnt++; // new chain
                Top.push_back(-1);

                build_hld(to, v, G); // build hld
            }
        }
    }

public:

    hld(vector<vector<int>>& G) {
        int n = G.size();
        Size.resize(n);
        Parent.resize(n);
        Next.resize(n);
        Depth.resize(n);
        Num.resize(n);
        Chain.resize(n);

        Top.push_back(-1);

        Tree = seg_tree(n);

        dfs(0, 0, G);
        build_hld(0, 0, G);
    }

    void modify(int v, int add) {
        Tree.update(Num[v], add);
    }

    int get(int a, int b) {
        int res = 0;

        // пока не находимся в одной цепочке
        while (Chain[a] != Chain[b]) {
            // будем подниматься с самых нижних и релаксировать ответ

            if (Depth[Top[Chain[a]]] < Depth[Top[Chain[b]]]) {
                swap(a, b);
            }
            int begin = Top[Chain[a]];

            res = max(res, Tree.get(Num[begin], Num[a]));

            a = Parent[begin];
        }

        // теперь мы в одной цепочке. Возьмем в ней ответ

        if (Depth[a] > Depth[b]) {
            swap(a, b);
        }

        res = max(res, Tree.get(Num[a], Num[b]));
        return res;
    }
};
