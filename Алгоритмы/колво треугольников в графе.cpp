// дан граф. нужно посчитать колво треугольников (цикл длиной 3)
// O(n^1.5)
// n - колво вершин
// m - колво ребер
// g - граф
// edges - список ребра

int n, m;
cin >> n >> m;

vector<vector<int>> g(n);
vector<pair<int, int>> edges(m);

for (int i = 0; i < m; i++) {
    int u, v;
    cin >> u >> v;
    u--, v--;

    edges[i] = { u, v };

    g[v].push_back(u);
    g[u].push_back(v);
}

// убрать лишние ребра в графе
{
    vector<vector<int>> new_g(n);
    for (int i = 0; i < m; i++) {
        int u, v;
        tie(u, v) = edges[i];

        if (g[u].size() > g[v].size()) {
            new_g[v].push_back(u);
        }
        else if (g[u].size() == g[v].size()) {
            if (u > v) {
                swap(u, v);
            }
            new_g[u].push_back(v);
        }
        else {
            new_g[u].push_back(v);
        }
    }

    g = new_g;
}

vector<bool> have_v(n);
int ans = 0;

for (int v = 0; v < n; v++) {
    for (int u : g[v]) {
        have_v[u] = true;
    }

    // have_v[u] - смежна ли вершина с v

    for (int u : g[v]) {
        for (int w : g[u]) {
            // если есть ребро v->w
            ans += have_v[w];
        }
    }

    for (int u : g[v]) {
        have_v[u] = false;
    }
}

cout << ans << "\n";