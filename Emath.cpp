// Emath 1.1 | Developed by Mr_Straym

#include<cstring>
#include<iomanip>
#include<algorithm>
#include<vector>
#include<stack>
#include<cmath>
#include<queue>
#include<random>
#include<ctime>
#include<set>
#include<deque>
#include<chrono>
#include<fstream>
#include<bitset>


// Data Structures: trie, edeque, segTree, hashTable, container, dsu, fenwick, lwf, PairingHeap, matrix, sparse_table
namespace dst {

    // в качестве значение возвращает ключ
    template<typename T>
    unsigned long long defaultHashFoo(const T& key) {
        return key;
    }

    // Динамическая хеш таблица на открытой адресации. Баланс скорости и памяти
    // Оптимизации Робин Гуд
    template<typename key_t, typename val_t, unsigned long long(*hashFunction)(const key_t&) = defaultHashFoo>
    class hashTable {

        // структура ячейки в хеш таблице
        struct cell {
            std::pair<key_t, val_t> p;
            bool free; // свободна

            // по умолчанию ячейка свободна. ЭТО ВАЖНО
            cell() {
                free = true;
            }
            cell(const key_t& key, const val_t& val) {
                free = false;
                p = std::make_pair(key, val);
            }
            cell(const std::pair<key_t, val_t>& p) {
                free = false;
                this->p = p;
            }
            cell(std::pair<key_t, val_t>&& p) {
                free = false;
                this->p = std::move(p);
            }
        };

        cell* A; // массив
        int length; // длина массива
        int count; // количество элементов

        // В конце массива A всегда будет находиться пустая ячейка. Для ускорения кода

        // хеширование ключа
        int hash(const key_t& key) const {
            return hashFunction(key) % length;
        }

        // берет размер массива с "запасом"
        int sizeInc() const {
            return 1.6 * length + 1;
        }

        // перестраивает структуру, увеличивая размер
        void rebuild() {
            hashTable temp(sizeInc()); // создаем новую хеш таблицу
            for (cell* it = A; it != A + length; it++) {
                if (!it->free) { // не пуста
                    temp.insert(it->p);
                }
            }
            delete[] A; // memory clear
            move(temp); // move
        }

        // копирует структуру, но не освобождает память
        void copy(const hashTable& source) {
            A = new cell[(length = source.length) + 1];
            count = source.count;
            for (int i = 0; i < length; i++) {
                A[i] = source.A[i];
            }
        }

        // обнуляет структуру. Изначально нужно иметь 1 элемент
        void zeroing() {
            A = new cell[2];
            length = 1;
            count = 0;
        }

        // возвращает индекс элемента с таким ключом или конец кластера если not found
        cell* findKey(const key_t& key) const {
            cell* it = A + hash(key);
            // пока ячейки не пусты и ключи не равны
            while (!it->free && it->p.first != key) {
                it++;
            }
            return it;
        }

        // move struct
        void move(hashTable& source) {
            // полностью копируем 
            A = source.A;
            length = source.length;
            count = source.count;
            // обнуляем source
            source.zeroing();
        }

    public:

        // конструктор по умолчанию
        hashTable() {
            zeroing();
        }
        // конструктор копирования
        hashTable(const hashTable& source) {
            copy(source);
        }
        // move constuctor
        hashTable(hashTable&& source) noexcept {
            move(source);
        }
        // деструктор
        ~hashTable() {
            delete[] A;
        }

        hashTable& operator = (const hashTable& source) {
            if (A != source.A) {
                delete[] A;
                copy(source);
            }
            return *this;
        }
        hashTable& operator = (hashTable&& source) {
            if (A != source.A) {
                delete[] A;
                move(source);
            }
            return *this;
        }

        // конструктор выделения
        hashTable(size_t newSize) {
            A = new cell[newSize + 1];
            length = newSize;
            count = 0;
        }

        // освобождает память
        void clear() {
            delete[] A;
            zeroing();
        }

        int size() const {
            return count;
        }
        bool empty() const {
            return count == 0;
        }

        // простой итератор
        class iterator {
            cell* it;
            cell* end; // ссылка на конец

            // идет вперед, пока не найдет занятую ячейку
            void findNext() {
                while (it != end && it->free) {
                    it++;
                }
            }

            iterator& inc() {
                it++;
                findNext();
                return *this;
            }

        public:
            iterator() {}
            iterator(cell* it, cell* end) {
                this->it = it;
                this->end = end;

                findNext();
            }

            iterator operator ++(int) {
                auto temp = *this;
                inc();
                return temp;
            }
            iterator& operator ++() {
                return inc();
            }

            std::pair<key_t, val_t>& operator *() {
                return it->p;
            }
            std::pair<key_t, val_t>* operator ->() {
                return &**this;
            }

            bool operator == (const iterator& Rhs) const {
                return it == Rhs.it;
            }
            bool operator != (const iterator& Rhs) const {
                return it != Rhs.it;
            }
        };

        iterator begin() {
            return iterator(A, A + length);
        }
        iterator end() {
            return iterator(A + length, A + length);
        }

        val_t& operator [](const key_t& key) {
            cell* it = findKey(key);
            if (it->free) { // если не нашли
                if (it == A + length) { // overflow
                    rebuild();
                    return this->operator[](key); // repeat
                }
                else { // add element
                    *it = cell(key, val_t());
                    count++;
                }
            }
            return it->p.second;
        }

    private:
        void insert(std::pair<key_t, val_t>& p) {
            cell* it = A + hash(p.first);
            while (!it->free) {
                // нашли богача
                if (hash(p.first) < hash(it->p.first)) {
                    // меняем значения и продолжаем
                    std::swap(p, it->p);
                }
                it++;
            }
            if (it == A + length) { // overflow
                rebuild();
                insert(p); // repeat
            }
            else {
                *it = cell(p);
                count++;
            }
        }
    public:
        // добавляет НОВЫЙ элемент
        // Оптимизация Робин Гуд
        void insert(const key_t& key, const val_t& value) {
            insert(std::make_pair(key, value));
        }
        void insert(key_t&& key, const val_t&& value) {
            insert(std::make_pair(std::move(key), std::move(value)));
        }

        // возвращает элемент с таким ключом или end если not found
        iterator find(const key_t& key) const {
            cell* it = findKey(key);
            return iterator(it->free ? A + length : it, A + length);
        }
    };

    // line with function(min/max)
    template<typename T, bool(*compare)(T, T)>
    class lwf {
        dst::edeque<T> d;
    public:
        int size() const {
            return d.size();
        }
        bool empty() const {
            return d.empty();
        }

        const T get() const {
            return d.front();
        }
        void push(const T& value) {
            while (!d.empty() && compare(value, d.back())) {
                d.pop_back();
            }
            d.push_back(value);
        }
        void erase(const T& key) {
            if (!d.empty() && d.front() == key) {
                d.pop_front();
            }
        }
    };
}

// Algorithm: Numbers, Graphs, Strings, Sort
namespace alg {

    // Graphs: dijkstra, mst, graphSets, graphCycles
    namespace gph {

        // find bridges
        struct bridge {

            vector<vector<int>> G;
            vector<bool> Vis;
            vector<int> tin, fup;
            int timer;


            map<pair<s64, s64>, s64> Map;
            vector<s64> Ans;

            void dfs(int v, int p = -1) {
                Vis[v] = true;
                tin[v] = fup[v] = timer++;

                for (auto to : G[v]) {
                    if (to != p) { // не идем в родителя

                        if (Vis[to]) {
                            fup[v] = min(fup[v], tin[to]);
                        }
                        else {
                            dfs(to, v);
                            fup[v] = min(fup[v], fup[to]);
                            if (fup[to] > tin[v])
                                Ans.push_back(Map[make_pair(v, to)]);
                        }
                    }
                }
            }

            void find_bridges() {
                timer = 0;

                Vis.resize(G.size(), false);
                tin.resize(G.size());
                fup.resize(G.size());

                for (int i = 0; i < G.size(); i++) {
                    if (!Vis[i]) {
                        dfs(i);
                    }
                }
            }
        };

        class lca {
            typedef std::vector<std::vector<s64>> graph;
            std::vector<s64> H, dfsList, first, tree;

            // dfs для инициализации LCA
            void dfs(s64 v, s64 h, s64 prev, const graph& G) {
                H[v] = h;
                dfsList.push_back(v);
                for (auto& to : G[v]) {
                    if (to != prev) {
                        dfs(to, h + 1, v, G);
                        dfsList.push_back(v);
                    }
                }
            }

            // строит дерево для LCA
            void buildTree(s64 v, s64 l, s64 r) {
                if (l == r) {
                    tree[v] = dfsList[l];
                }
                else {
                    s64 m = (l + r) / 2;
                    buildTree(v * 2, l, m);
                    buildTree(v * 2 + 1, m + 1, r);
                    tree[v] = H[tree[2 * v]] < H[tree[2 * v + 1]] ?
                        tree[2 * v] : tree[2 * v + 1];
                }
            }

            // строит LCA
            void build(s64 root, const graph& G) {
                s64 n = G.size();

                H.resize(n);
                dfsList.reserve(2 * n);

                dfs(root, 0, -1, G);

                s64 m = dfsList.size();
                tree.assign(dfsList.size() * 4 + 1, -1);
                buildTree(1, 0, m - 1);

                first.assign(n, -1);
                for (s64 i = 0; i < m; i++) {
                    s64 v = dfsList[i];
                    if (first[v] == -1) {
                        first[v] = i;
                    }
                }
            }

            // взятие минимума для LCA
            s64 treeMin(s64 v, s64 vl, s64 vr, s64 l, s64 r) {
                if (vl == l && vr == r) {
                    return tree[v];
                }
                else {
                    s64 vm = (vl + vr) / 2;
                    if (r <= vm) {
                        return treeMin(2 * v, vl, vm, l, r);
                    }
                    else if (l > vm) {
                        return treeMin(2 * v + 1, vm + 1, vr, l, r);
                    }
                    s64 cntL = treeMin(2 * v, vl, vm, l, vm);
                    s64 cntR = treeMin(2 * v + 1, vm + 1, vr, vm + 1, r);
                    return H[cntL] < H[cntR] ? cntL : cntR;
                }
            }

        public:

            // returns LCA(a, b)
            s64 get(s64 a, s64 b) {
                s64 left = first[a];
                s64 right = first[b];
                if (left > right) {
                    std::swap(left, right);
                }
                return treeMin(1, 0, dfsList.size() - 1, left, right);
            }

            // build LCA
            lca(s64 root, const graph& G) {
                build(root, G);
            }
        };

        // Алгоритм Дейкстры
        template<typename dist_t>
        class dijkstra {

            std::vector<dist_t> Dist;
            std::vector<int> Parent;
            int begin;

        public:

            // Алгоритм Дейкстры. O(m * log n)
            dijkstra(const std::vector<std::vector<std::pair<int, dist_t>>>& Graph, int begin, dist_t infinity) {
                this->begin = begin;

                Dist.resize(Graph.size(), infinity);
                Dist[begin] = 0;

                Parent.resize(Graph.size());

                dst::pairingHeap<std::pair<dist_t, int>> Q;
                Q.push(std::make_pair(0, begin));

                while (!Q.empty()) {
                    int v = Q.top().second;
                    dist_t cur_d = Q.top().first;
                    Q.pop();
                    if (cur_d > Dist[v]) {
                        continue;
                    }
                    for (int j = 0; j < Graph[v].size(); j++) {
                        int to = Graph[v][j].first;
                        dist_t len = Graph[v][j].second;
                        if (Dist[v] + len < Dist[to]) {
                            Dist[to] = Dist[v] + len;
                            Parent[to] = v;
                            Q.push(std::make_pair(Dist[to], to));
                        }
                    }
                }
            }

            // returns dist 
            dist_t getDist(size_t end) const {
                return Dist[end];
            }

            // returns path
            std::vector<int> getPath(size_t end) const {
                std::vector<int> res;
                while (end != begin) {
                    res.push_back(end);
                    end = Parent[end];
                }
                res.push_back(begin);
                std::reverse(res.begin(), res.end());
                return res;
            }
        };

        // Минимальное остовное дерево. Алгоритм Крускала с использованием dsu
        // Предполагается что граф связный
        // (dist, (u, v))
        template<typename dist_t>
        std::pair<dist_t, std::vector<std::pair<int, int>>> mst(std::vector<std::pair<dist_t, std::pair<int, int>>> Graph) {
            dist_t cost = 0;
            std::vector<std::pair<int, int>> res;

            sort(Graph.begin(), Graph.end());

            dst::dsu Dsu(Graph.size());

            for (int i = 0; i < Graph.size(); i++) {
                int u = Graph[i].second.first, v = Graph[i].second.second;
                dist_t len = Graph[i].first;

                if (Dsu.findSet(u) != Dsu.findSet(v)) {
                    cost += len;
                    res.push_back(Graph[i].second);
                    Dsu.unionSets(u, v);
                }
            }
            return std::make_pair(cost, res);
        }

        // Минимальное остовное дерево. Алгоритм Крускала с использованием dsu
        // Предполагается что граф связный
        // (to, dist)
        template<typename dist_t>
        std::pair<dist_t, std::vector<std::pair<int, int>>> mst(const std::vector<std::vector<std::pair<int, dist_t>>>& Graph) {
            // (dist_t, (u, v))
            std::vector<std::pair<dist_t, std::pair<int, int>>> g;
            for (int i = 0; i < Graph.size(); i++) {
                for (int j = 0; j < Graph[i].size(); j++) {
                    g.push_back(std::make_pair(Graph[i][j].second, std::make_pair(i, Graph[i][j].first)));
                }
            }

            return mst(g);
        }


        // возвращает множества графа
        std::vector<std::vector<int>> graphSets(const std::vector<std::vector<int>>& Graph) {
            dst::dsu Dsu(Graph.size());
            for (int i = 0; i < Graph.size(); i++) {
                for (int j = 0; j < Graph[i].size(); j++) {
                    Dsu.unionSets(i, Graph[i][j]);
                }
            }
            std::vector<std::vector<int>> Sets;
            dst::hashTable<int, int> Map;
            for (int i = 0; i < Graph.size(); i++) {
                int set = Dsu.findSet(i);
                if (Map.find(set) != Map.end()) { // если такое множество уже есть
                    Sets[Map[set]].push_back(i);
                }
                else { // новое множество
                    Map[set] = Sets.size();

                    Sets.push_back(std::vector<int>());
                    Sets.back().push_back(i);
                }
            }
            return Sets;
        }

        // возвращает множества графа
        template<typename dist_t>
        std::vector<std::vector<int>> graphSets(const std::vector<std::vector<std::pair<int, dist_t>>>& Graph) {
            dst::dsu Dsu(Graph.size());
            for (int i = 0; i < Graph.size(); i++) {
                for (int j = 0; j < Graph[i].size(); j++) {
                    Dsu.unionSets(i, Graph[i][j].first);
                }
            }
            std::vector<std::vector<int>> Sets;
            dst::hashTable<int, int> Map;
            for (int i = 0; i < Graph.size(); i++) {
                int set = Dsu.findSet(i);
                if (Map.find(set) != Map.end()) { // если такое множество уже есть
                    Sets[Map[set]].push_back(i);
                }
                else { // новое множество
                    Map[set] = Sets.size();

                    Sets.push_back(std::vector<int>());
                    Sets.back().push_back(i);
                }
            }
            return Sets;
        }


        // graphCycles

        enum class vertex_t {
            undefind, // не были
            visit, // были
            exit, // были, вышли
        };

        // добавляет новый цикл
        void addCycle(int begin, int end, const std::vector<int>& Parents, std::vector<std::vector<int>>& Cycles) {
            Cycles.push_back(std::vector<int>());
            for (int v = end; v != begin; v = Parents[v]) {
                Cycles.back().push_back(v);
            }
            Cycles.back().push_back(begin);
        }

        // (вершина, предыдущая вершина, родители, цвета, циклы, граф)
        void dfsFindCycles(int v, int prev, std::vector<int>& Parents, std::vector<vertex_t>& Type, std::vector<std::vector<int>>& Cycles, const std::vector<std::vector<int>>& G) {
            Type[v] = vertex_t::visit;

            // dfs
            {
                // проходим по всем ребрам
                for (int i = 0; i < G[v].size(); i++) {
                    // v -> to
                    int to = G[v][i];

                    // если мы не идем туда, где были ход назад
                    if (to != prev) {

                        // если еще не были там
                        if (Type[to] == vertex_t::undefind) {
                            Parents[to] = v;
                            dfsFindCycles(to, v, Parents, Type, Cycles, G);
                        }
                        // нашли цикл
                        else if (Type[to] == vertex_t::visit) {
                            addCycle(to, v, Parents, Cycles);
                        }
                    }
                }
            }

            Type[v] = vertex_t::exit;
        }

        // Находит все циклы в неориентированном графе
        std::vector<std::vector<int>> graphCycles(const std::vector<std::vector<int>>& G) {
            std::vector<int> Parents(G.size(), -1); // родители
            std::vector<vertex_t> Type(G.size(), vertex_t::undefind); // тип вершин
            std::vector<std::vector<int>> Cycles; // циклы

            for (int i = 0; i < G.size(); i++) {
                // если еще не были там
                if (Type[i] == vertex_t::undefind) {
                    dfsFindCycles(i, -1, Parents, Type, Cycles, G);
                }
            }

            return Cycles;
        }
    }

    // Strings: zFunction, strHash, strHashPref, findSubstr, cyclicShift, preFunction(prefix)
    namespace str {
        // z функция для строки за O(n)
        std::vector<int> zFunction(const std::string& str) {
            long long n = str.size();
            std::vector<int> z(n);
            for (long long i = 1, l = 0, r = 0; i < n; i++) {
                if (i <= r) {
                    z[i] = min(r - i + 1, static_cast<long long>(z[i - l]));
                }
                while (i + z[i] < n && str[z[i]] == str[i + z[i]]) {
                    z[i]++;
                }
                if (i + z[i] - 1 > r) {
                    l = i, r = i + z[i] - 1;
                }
            }
            return z;
        }

        std::vector<unsigned long long> hashPow;
        // Обновляет вектор множителей для хеша строк
        void updateHashPow(int strLen) {
            if (hashPow.size() < strLen) {
                size_t i = hashPow.size();
                hashPow.resize(strLen);
                if (i == 0) {
                    hashPow[0] = 1;
                    i++;
                }
                for (; i < strLen; i++) {
                    hashPow[i] = hashPow[i - 1] * static_cast<unsigned long long>(53);
                }
            }
        }

        // a..zA..Z0..9
        unsigned char SymbolCode[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 0, 0, 0, 0, 0, 0, 0, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        // Хеш строки
        unsigned long long strHash(const std::string& str) {
            updateHashPow(str.size());
            unsigned long long hash = 0;
            for (size_t i = 0; i < str.size(); i++) {
                hash += hashPow[i] * static_cast<unsigned long long>(SymbolCode[str[i]]);
            }
            return hash;
        }

        // Префиксный хеш строки
        std::vector<unsigned long long> strHashPref(const std::string& str) {
            updateHashPow(str.size());
            std::vector<unsigned long long> res(str.size());
            for (size_t i = 0; i < str.size(); i++) {
                res[i] = hashPow[i] * static_cast<unsigned long long>(SymbolCode[str[i]]);
                res[i] += (i ? res[i - 1] : 0);
            }
            return res;
        }

        // Находит все вхождения строки key в s
        std::vector<int> findSubstr(const std::string& s, const std::string& key) {
            auto sHashPref = strHashPref(s);
            auto keyHash = strHash(key);
            std::vector<int> res;
            for (size_t i = 0; i + key.size() - 1 < s.size(); i++) {
                unsigned long long curHash = sHashPref[i + key.size() - 1];
                curHash -= i ? sHashPref[i - 1] : 0;

                // приводим хэши к одной степени и сравниваем
                if (curHash == keyHash * hashPow[i]) {
                    res.push_back(i);
                }
            }
            return res;
        }

        // Возвращает наименьший циклический сдвиг
        int cyclicShift(std::string str) {
            str += str;
            int n = str.size(), i = 0, ans = 0;
            while (i < n / 2) {
                ans = i;
                int j = i + 1, k = i;
                while (j < n && str[k] <= str[j]) {
                    k = str[k] < str[j] ? i : k + 1;
                    j++;
                }
                while (i <= k) {
                    i += j - k;
                }
            }
            return ans + 1;
        }

        std::vector<int> preFunction(const std::string& str) {
            long long n = str.size();
            std::vector<int> result(n);
            for (long long i = 1; i < n; i++) {
                long long j = result[i - 1];
                while (j > 0 && str[i] != str[j]) {
                    j = result[j - 1];
                }
                if (str[i] == str[j]) {
                    j++;
                }
                result[i] = j;
            }
            return result;
        }
    }

    // Sort: bubbleSort, heapSort, mergeSort
    namespace sort {
        // пузырьковая сортировка. O(n^2)
        template<typename iter_t>
        void bubbleSort(const iter_t begin, const iter_t end) {
            iter_t temp = begin;
            temp++;
            while (temp != end) {
                iter_t it = temp, prev = temp;
                prev--;
                // пока не дошли до конца
                while (prev != begin && *it < *prev) {
                    std::swap(*it, *prev);
                    it = prev--;
                }
                if (*it < *prev) {
                    std::swap(*it, *prev);
                }
                temp++;
            }
        }

        // сортировка кучей. O(n * log n)
        template<typename T, typename iter_t>
        void heapSort(const iter_t begin, const iter_t end) {
            dst::pairingHeap<T> Heap;
            iter_t temp = begin;
            while (temp != end) {
                Heap.push(std::move(*temp));
                temp++;
            }
            temp = begin;
            while (temp != end) {
                *temp = std::move(Heap.top());
                Heap.pop();
                temp++;
            }
        }

        // begin     mid         end
        // sortArray1 + sortArray2
        template<typename T>
        void merge(T* begin, T* mid, T* end, T* A) {
            T* left = begin, * right = mid;
            int i = 0;
            // пока не вышли за пределы
            while (left != mid && right != end) {
                A[i++] = *left < *right ? std::move(*left++) : std::move(*right++);
            }
            while (left != mid) {
                A[i++] = std::move(*left++);
            }
            while (right != end) {
                A[i++] = std::move(*right++);
            }

            i = 0;
            while (begin != end) {
                *begin++ = std::move(A[i++]);
            }
        }

        template<typename T>
        void _mergeSort(T* begin, T* end, T* A) {
            if (begin + 1 != end) {
                T* mid = begin + ((end - begin) >> 1);
                // sort halves
                _mergeSort(begin, mid, A);
                _mergeSort(mid, end, A + (mid - begin));
                // merge halves
                merge(begin, mid, end, A);
            }
        }

        // сортировка слиянием. O(n * log n)
        template<typename T, typename iter_t>
        void mergeSort(const iter_t begin, const iter_t end) {
            size_t len = 0;
            iter_t temp = begin;
            { // find len
                while (temp != end) {
                    temp++;
                    len++;
                }
            }
            T* Res = new T[len], * A = new T[len];
            int i;
            { // move to Res
                i = 0;
                temp = begin;
                while (temp != end) {
                    Res[i++] = std::move(*temp);
                    temp++;
                }
            }
            // sort
            _mergeSort(Res, Res + len, A);
            { // move to container
                temp = begin;
                i = 0;
                while (temp != end) {
                    *temp = std::move(Res[i++]);
                    temp++;
                }
            }
            delete[] Res, A;
        }


        // !quickSort is under construction!
    }
}
