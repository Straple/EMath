const int alp = 'z' - 'a' + 1;

struct trie {

    struct node {
        int next[alp];

        vector<int> set;

        int parent;
        char parent_char;

        int suff_link;

        node() {
            for (int i = 0; i < alp; i++) {
                next[i] = -1;
            }
            parent = suff_link = 0;
            parent_char = '?';
        }
    };

    vector<node> t;

    trie() {
        t.resize(1);
    }

    int get(int v, char c) {
        if (t[v].next[c - 'a'] == -1) {
            t.push_back(node());
            t.back().parent = v;
            t.back().parent_char = c;
            t[v].next[c - 'a'] = t.size() - 1;
        }
        return t[v].next[c - 'a'];
    }

    void ins(string s, int i) {
        int p = 0;
        for (char c : s) {
            p = get(p, c);
        }
        t[p].set.push_back(i);
    }

    void build() {
        queue<int> Q;
        Q.push(0);
        while (!Q.empty()) {
            int v = Q.front();
            Q.pop();

            auto& [next, set, parent, parent_char, suff_link] = t[v];

            // find suff_link
            {
                if (v == 0) {
                    suff_link = 0;
                }
                else {
                    int p = t[parent].suff_link;
                    while (p != 0 && t[p].next[parent_char - 'a'] == -1) {
                        p = t[p].suff_link;
                    }

                    int to = t[p].next[parent_char - 'a'];
                    if (to != -1 && to != v) {
                        suff_link = to;
                    }
                }
            }

            for (int i = 0; i < alp; i++) {
                if (next[i] != -1) {
                    Q.push(next[i]);
                }
            }
        }
    }

    int move(int v, char c) {
        while (v != 0 && t[v].next[c - 'a'] == -1) {
            v = t[v].suff_link;
        }
        v = t[v].next[c - 'a'];
        if (v == -1) {
            v = 0;
        }
        return v;
    }

    vector<int> get(int v) {
        vector<int> vs;
        while (v != 0) {
            vs += t[v].set;
            v = t[v].suff_link;
        }
        return vs;
    }
};

void solve(istream& cin, ostream& cout) {
    string t;
    int n;
    trie obj;
    vector<string> strs;
    {
        cin >> t;
        cin >> n;
        strs.resize(n);
        for (int i = 0; i < n; i++) {
            string s;
            cin >> s;
            strs[i] = s;
            obj.ins(s, i);
        }
        obj.build();
    }

    vector<vector<int>> vs(n);

    int v = 0;
    for (int i = 0; i < t.size(); i++) {
        v = obj.move(v, t[i]);

        vector<int> p = obj.get(v);
        for (int id : p) {
            vs[id].push_back(i - strs[id].size() + 1);
        }
    }

    for (int i = 0; i < n; i++) {
        cout << vs[i].size() << " ";
        for (int x : vs[i]) {
            cout << x + 1 << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}
