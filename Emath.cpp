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

// Developed by Mr_Straple

// Utils
// random, vector addition, vector output, clamp, range, abs
namespace utl {

    // random with random seed
    std::mt19937_64 random(std::random_device{}());

    // vector addition
    template<typename _Ty>
    std::vector<_Ty> operator + (const std::vector<_Ty>& a, const std::vector<_Ty>& b) {
        std::vector<_Ty> ret = a;
        ret.resize(a.size() + b.size());
        for (size_t i = 0; i < b.size(); i++) {
            ret[i + a.size()] = b[i];
        }
        return ret;
    }

    // Add to vector
    template<typename _Ty>
    std::vector<_Ty>& operator += (std::vector<_Ty>& a, const std::vector<_Ty>& b) {
        size_t len = a.size();
        a.resize(len + b.size());
        for (size_t i = 0; i < b.size(); i++) {
            a[i + len] = b[i];
        }
        return a;
    }

    template<typename _Ty>
    std::ostream& operator << (std::ostream& output, const std::vector<_Ty>& array) {
        for (int i = 0; i < array.size(); i++) {
            output << array[i] << " ";
        }
        return output << "\n";
    }

    template<typename _Ty1, typename _Ty2, typename _Ty3>
    _Ty2& clamp(const _Ty1& min, const _Ty2& value, const _Ty3& max) {
        return value > max ? static_cast<_Ty2>(max) :
            value < min ? static_cast<_Ty2>(min) :
            value;
    }

    template<typename _Ty>
    _Ty range(const _Ty& left, const _Ty& right) {
        return random() % (right - left + 1) + left;
    }

    template<typename _Ty>
    _Ty abs(const _Ty& value) {
        return value < 0 ? value * -1 : value;
    }
}
using namespace utl;

// Algorithm
// Numbers: sqrt, eratosthenesSieve, isPrime, factorize, epow, gcd, gcd(xy), lcm, fibonacci
// Graphs:  dijkstra
// Strings: zFunction, strHash, strHashPref, findSubstr, cyclicShift, preFunction(prefix)
namespace alg {

    // Метод Ньютона для поиска целочисленных корней
    template<typename _Ty>
    _Ty sqrt(const _Ty& n) {
        _Ty x = 1;
        bool decreased = false;
        while(true){
            _Ty nx = (x + n / x) / 2;
            if (x == nx || nx > x && decreased) {
                break;
            }
            decreased = nx < x;
            x = nx;
        }
        return x;
    }

    // Решето Эратосфена
    template<typename _Ty>
    std::vector<_Ty> eratosthenesSieve(_Ty n) {
        std::vector<bool> prime(n + 1, true);
        for (int i = 2; i * i <= n; i++) {
            if (prime[i]) {
                for (int j = i * i; j <= n; j += i) {
                    prime[j] = false;
                }
            }
        }

        std::vector<_Ty> result;
        for (int i = 2; i <= n; i++) {
            if (prime[i]) {
                result.push_back(i);
            }
        }
        return result;
    }

    // return true if n is prime
    template<typename _Ty>
    bool isPrime(_Ty n) {
        _Ty i;
        for (i = 2; i * i <= n && n % i != 0; i++) {}
        return i * i > n;
    }

    // factorizes the number
    // vector<a^b>
    template<typename _Ty>
    std::vector<std::pair<_Ty, int>> factorize(_Ty N) {
        std::vector<std::pair<_Ty, int>> result;
        for (_Ty i = 2; i * i <= N; i++) {
            if (N % i == 0) {
                result.push_back(std::make_pair(i, 0));
                while (N % i == 0) {
                    result.back().second++;
                    N /= i;
                }
            }
        }
        if (N != 1) {
            result.push_back(std::make_pair(N, 1));
        }
        return result;
    }

    // binary pow: a^n % mod
    template<typename T1, typename T2>
    T1 epow(const T1& a, const T2& n, const T1& mod) {
        if (n == 0) {
            return 1;
        }
        else {
            T1 z = epow(a, n / 2, mod);
            z = (z * z) % mod;

            return (n % 2 == 0) ? z : (z * a) % mod;
        }
    }

    // binary pow a^n
    template<typename T1, typename T2>
    T1 epow(const T1& a, const T2& n) {
        if (n == 0) {
            return 1;
        }
        else {
            T1 z = epow(a, n / 2);
            z = z * z;

            return (n % 2 == 0 ? z : z * a);
        }
    }

    // Алгоритм Евклида
    template<typename T>
    T gcd(T a, T b) {
        while (a > 0 && b > 0) {
            if (a > b) {
                a %= b;
            }
            else {
                b %= a;
            }
        }
        return std::max(a, b);
    }

    // Расширенный Алгоритм Евклида
    template<typename _Ty>
    _Ty gcd(_Ty a, _Ty b, _Ty& x, _Ty& y) {
        if (a == 0) {
            x = 0;
            y = 1;
            return b;
        }
        else {
            _Ty x1, y1;
            _Ty res = gcd(b % a, a, x1, y1);
            x = y1 - (b / a) * x1;
            y = x1;
            return res;
        }
    }

    // Линейные диофантовы уравнения с двумя переменными
    // a * x + b * y = c
    // a * xg + b * yg = g
    template<typename _Ty>
    bool solveLDE(_Ty a, _Ty b, _Ty c, _Ty& xg, _Ty& yg, _Ty& g) {
        g = gcd(abs(a), abs(b), xg, yg);
        if (c % g != 0) {
            return false;
        }
        xg *= c / g;
        yg *= c / g;
        xg *= a < 0 ? -1 : 1;
        yg *= b < 0 ? -1 : 1;
        return true;
    }

    //lcm(a, b) = a * b / gcd(a, b)
    template<typename T>
    T lcm(const T& a, const T& b) {
        return a * b / gcd(a, b);
    }

    // a^2 + b^2 = N
    // Представление числа в виде суммы двух квадратов
    template<typename _Ty>
    bool squarePres(const _Ty& N, _Ty& a, _Ty& b) {
        for (a = 2; a <= sqrt(N); a++) {
            b = sqrt(N - a * a);
            if (a * a + b * b == N) {
                return true;
            }
        }
        return false;
    }

    // Обратный элемент в кольце по модулю
    template<typename _Ty>
    bool reverseElement(_Ty a, _Ty m, _Ty &x){
        _Ty y;
        _Ty g = gcd(a, m, x, y);
        if (g != 1) {
            return false;
        }
        else {
            x = (x % m + m) % m;
            return true;
        }
    }

    // Фибоначчи за log
    template<typename T>
    class fibonacci {

        struct matrix {
            std::vector<std::vector<T>> a;

            matrix() {
                a = { {0, 1}, {1, 1} };
            }
            matrix(const matrix& other) {
                a = other.a;
            }
            matrix mult(const matrix& mult) const {
                matrix res;
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        res.a[i][j] = a[i][0] * mult.a[0][j] + a[i][1] * mult.a[1][j];
                    }
                }
                return res;
            }
            matrix mult(const matrix& mult, const T& mod) const {
                matrix res;
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        res.a[i][j] = ((a[i][0] * mult.a[0][j]) % mod + (a[i][1] * mult.a[1][j] % mod)) % mod;
                    }
                }
                return res;
            }
        };

        matrix binPow(const matrix& a, const T& n) {
            if (n == 1) {
                return a;
            }
            else {
                matrix z = binPow(a, n / 2);
                z = z.mult(z);

                return n % 2 == 0 ? z : z.mult(a);
            }
        }
        matrix binPow(const matrix& a, const T& n, const T& mod) {
            if (n == 1) {
                return a;
            }
            else {
                matrix z = binPow(a, n / 2, mod);
                z = z.mult(z, mod);

                return n % 2 == 0 ? z : z.mult(a, mod);
            }
        }

    public:

        T get(const T& n) {
            return binPow(matrix(), n).a[0][0];
        }
        T get(const T& n, const T& mod) {
            return binPow(matrix(), n, mod).a[0][0];
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
        template<typename inf_type>
        dijkstra(const std::vector<std::vector<std::pair<int, dist_t>>>& Graph, size_t newBegin, const inf_type& infinity) {
            begin = newBegin;

            Dist.resize(Graph.size(), static_cast<dist_t>(infinity));
            Dist[begin] = 0;

            Parent.resize(Graph.size());

            std::set<std::pair<dist_t, size_t>> Q;
            Q.insert(std::make_pair(Dist[begin], begin));

            while (!Q.empty()) {
                size_t v = Q.begin()->second;
                Q.erase(Q.begin());

                for (size_t j = 0; j < Graph[v].size(); j++) {
                    size_t to = Graph[v][j].first;
                    dist_t length = Graph[v][j].second, newDist = Dist[v] + length;

                    if (newDist < Dist[to]) {
                        Q.erase(std::make_pair(Dist[to], to));
                        Parent[to] = v;
                        Dist[to] = newDist;
                        Q.insert(std::make_pair(Dist[to], to));
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

    // z функция для строки за O(n)
    std::vector<int> zFunction(std::string str) {
        int n = str.size();
        std::vector<int> z(n);
        for (int i = 1, l = 0, r = 0; i < n; i++) {
            if (i <= r) {
                z[i] = std::min(r - i + 1, z[i - l]);
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
                hashPow[i] = hashPow[i - 1] * 53ULL;
            }
        }
    }

    // Хеш строки
    unsigned long long strHash(std::string str) {
        updateHashPow(str.size());
        unsigned long long hash = 0;
        for (size_t i = 0; i < str.size(); i++) {
            hash += static_cast<unsigned long long>(str[i] - 'a' + 1) * hashPow[i];
        }
        return hash;
    }

    // Префиксный хеш строки
    std::vector<unsigned long long> strHashPref(std::string str) {
        updateHashPow(str.size());
        std::vector<unsigned long long> res(str.size());
        for (size_t i = 0; i < str.size(); i++) {
            res[i] = static_cast<unsigned long long>(str[i] - 'a' + 1) * hashPow[i];
            res[i] += (i ? res[i - 1] : 0);
        }
        return res;
    }

    // Находит все вхождения строки key в s
    std::vector<int> findSubstr(std::string s, std::string key) {
        auto sHashPref = strHashPref(s);
        auto keyHash = strHash(key);
        std::vector<int> res;
        for (size_t i = 0; i + key.size() - 1 < s.size(); i++){
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

    std::vector<int> preFunction(std::string str) {
        int n = str.size();
        std::vector<int> result(n);
        for (int i = 1; i < n; i++) {
            int j = result[i - 1];
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

// Data Structures
// tos(tree of segments), edeque, container, hashTable, dsu, fenwick, lwf(line with function(min/max))
namespace dst {
    // tree of segments
    template<typename T>
    class tos {
        std::vector<T> Tree;
        T(*function)(T, T);
        T nothing;
        int length;

        // get(1, 0, n - 1, l, r)
        T get(int v, int tl, int tr, int l, int r) {
            if (l > r) {
                return nothing;
            }
            if (l == tl && r == tr) {
                return Tree[v];
            }
            int tm = (tl + tr) / 2;
            return (*function)(get(v * 2, tl, tm, l, std::min(r, tm)), get(v * 2 + 1, tm + 1, tr, std::max(l, tm + 1), r));
        }

        // build(1, 0, n - 1):  build Tree
        void build(int v, int tl, int tr, const std::vector<T>& Arr) {
            if (tl == tr) {
                Tree[v] = Arr[tl];
            }
            else {
                int tm = (tl + tr) / 2;
                build(v * 2, tl, tm, Arr);
                build(v * 2 + 1, tm + 1, tr, Arr);
                Tree[v] = (*function)(Tree[v * 2], Tree[v * 2 + 1]);
            }
        }

        // update(1, 0, n - 1, pos, new_val):
        void update(int v, int tl, int tr, int pos, const T& new_val) {
            if (tl == tr) {
                Tree[v] = new_val;
            }
            else {
                int tm = (tl + tr) / 2;
                if (pos <= tm) {
                    update(v * 2, tl, tm, pos, new_val);
                }
                else {
                    update(v * 2 + 1, tm + 1, tr, pos, new_val);
                }
                Tree[v] = (*function)(Tree[v * 2], Tree[v * 2 + 1]);
            }
        }

    public:

        tos() {}
        tos(T(*new_function)(T, T), T new_nothing, const std::vector<T>& Arr) {
            length = Arr.size();
            Tree.clear();
            Tree.resize(length * 4);
            function = new_function;
            nothing = new_nothing;
            build(1, 0, Arr.size() - 1, Arr);
        }

        // get on a segment [l, r]
        T get(int left, int right) {
            return get(1, 0, length - 1, left, right);
        }

        // update pos element value
        // Arr[pos] = new_val
        void update(int pos, const T& new_val) {
            update(1, 0, length - 1, pos, new_val);
        }
    };

    // deque on a circular array
    template<typename T>
    class edeque {
        T* A;
        int length_A, length, head, tail;

        bool overflow() const {
            return length + 1 >= length_A;
        }
        bool underflow() const {
            return length * 3 < length_A;
        }

        int sizeInc(int Size) const {
            return static_cast<int>(1.3 * Size + 1.0);
        }
        int sizeDec(int Size) const {
            return Size >> 1;
        }

        int inc(int index) const {
            index++;
            return index < length_A ? index : index - length_A;
        }
        int dec(int index) const {
            index--;
            return index >= 0 ? index : length_A - 1;
        }

        void rebuild(int new_size) {
            T* temp = new T[new_size];
            int index = 0;
            int length = size();
            while (index < length) {
                temp[index++] = A[head];
                head = inc(head);
            }
            delete[] A;
            A = temp;
            head = 0;
            tail = length;
            length_A = new_size;
        }

        void zeroing() {
            A = 0;
            length_A = length = head = tail = 0;
        }

        void create(int new_size) {
            length_A = sizeInc(new_size);
            length = new_size;
            A = new T[length_A];
            head = 0;
            tail = new_size;
        }

    public:

        edeque() {
            zeroing();
        }
        edeque(const edeque& new_value) {
            zeroing();
            this->operator=(new_value);
        }
        edeque(const std::vector<T>& new_value) {
            zeroing();
            this->operator=(new_value);
        }
        edeque(int new_size) {
            create(new_size);
        }
        edeque(int new_size, const T& fill_value) {
            create(new_size);
            for (int i = 0; i < new_size; i++) {
                A[i] = fill_value;
            }
        }
        ~edeque() {
            clear();
        }

        edeque& operator = (const edeque& new_value) {
            clear();
            A = new T[new_value.length_A];
            length = new_value.length;
            length_A = new_value.length_A;
            head = new_value.head;
            tail = new_value.tail;
            for (size_t i = 0; i < length_A; i++) {
                A[i] = new_value.A[i];
            }

            return *this;
        }
        edeque& operator = (const std::vector<T>& new_value) {
            clear();
            create(new_value.size());
            for (size_t i = 0; i < size(); i++) {
                this->operator[](i) = new_value[i];
            }
            return *this;
        }

        void clear() {
            delete[] A;
            zeroing();
        }

        void resize(int new_size) {
            edeque temp(new_size, T());
            for (size_t i = 0; i < std::min(size(), new_size); i++) {
                temp[i] = this->operator[](i);
            }
            *this = temp;
        }

        int size() const {
            return length;
        }
        bool empty() const {
            return length == 0;
        }

        void push_back(const T& value) {
            if (overflow()) {
                rebuild(sizeInc(length_A));
            }
            A[tail] = value;
            tail = inc(tail);
            length++;
        }
        void push_front(const T& value) {
            if (overflow()) {
                rebuild(sizeInc(length_A));
            }
            head = dec(head);
            A[head] = value;
            length++;
        }

        void pop_back() {
            if (underflow()) {
                rebuild(sizeDec(length_A));
            }
            tail = dec(tail);
            length--;
        }
        void pop_front() {
            if (underflow()) {
                rebuild(sizeDec(length_A));
            }
            head = inc(head);
            length--;
        }

        T& operator [](int index) const {
            index += head;
            return A[index < length_A ? index : index - length_A];
        }

        T& back() const {
            return A[dec(tail)];
        }
        T& front() const {
            return A[head];
        }

        std::vector<T> getArray() const {
            std::vector<T> res(size());
            for (size_t i = 0; i < size(); i++) {
                res[i] = this->operator[](i);
            }
            return res;
        }
    };

    /* Implicit cartesian tree
       O(1): empty, size
       O(log n): insert, erase, [index], back, front, push_back
       O(n): clear, get
       O(n * log n): copy constructor, fill constructor, sort
    */
    template<typename T>
    class container {

        struct node;
        typedef node* pnode;

        struct node {
            T value;
            pnode left, right;
            int size;
            long long prior;

            node() {
                left = right = 0;
            }
            node(const T& new_value, const pnode& new_left, const pnode& new_right, int new_size) {
                value = new_value;
                left = new_left;
                right = new_right;
                size = new_size;
                prior = utl::random();
            }
        };

        pnode root;

        void clear(pnode& branch) {
            if (branch->left) {
                clear(branch->left);
                delete branch->left;
            }
            if (branch->right) {
                clear(branch->right);
                delete branch->right;
            }
        }

        void get(const pnode& branch, std::vector<T>& Arr, size_t& index) const {
            if (branch->left) {
                get(branch->left, Arr, index);
            }
            Arr[index] = branch->value;
            index++;
            if (branch->right) {
                get(branch->right, Arr, index);
            }
        }

        // merge branches. O(h)
        void merge(pnode& branch, pnode left, pnode right) {
            if (!left || !right) {
                branch = left ? left : right;
            }
            else if (left->prior > right->prior) {
                merge(left->right, left->right, right);
                branch = left;
            }
            else {
                merge(right->left, left, right->left);
                branch = right;
            }
            resize(branch);
        }

        // split branches. O(h)
        void split(pnode branch, pnode& left, pnode& right, int key) {
            if (!branch) {
                left = right = 0;
                return;
            }
            else {

                if (key <= size(branch->left)) {
                    split(branch->left, left, branch->left, key);
                    right = branch;
                }
                else {
                    split(branch->right, branch->right, right, key - size(branch->left) - 1);
                    left = branch;
                }
                resize(branch);
            }
        }

        // return size in branch. O(1)
        int size(pnode branch) const {
            return branch ? branch->size : 0;
        }

        // update size in branch. O(1)
        void resize(pnode branch) {
            if (branch) {
                branch->size = size(branch->left) + size(branch->right) + 1;
            }
        }

        // find [index]. O(h)
        const pnode& find(const pnode& branch, int index) const {
            if (index < size(branch->left)) {
                return find(branch->left, index);
            }
            else {
                index -= size(branch->left);
                if (index == 0) {
                    return branch;
                }
                else {
                    return find(branch->right, index - 1);
                }
            }
        }

    public:

        // default constructor
        container() {
            root = 0;
        }

        // copy constructor. O(n * log n)
        container(const container<T>& source) {
            root = 0;

            *this = source;
        }

        // fill constructor. O(n * log n)
        container(int length, const T& value) {
            root = 0;

            while (length--) {
                push_back(value);
            }
        }

        // destructor. Clear all elements. O(n)
        ~container() {
            clear();
        }

        container<T>& operator = (const container<T>& source) {
            clear();
            for (size_t i = 0; i < source.size(); i++) {
                push_back(source[i]);
            }
            return *this;
        }

        // O(1)
        bool empty() const {
            return root == 0;
        }

        // Clear all elements. O(n)
        void clear() {
            if (root) {
                clear(root);
                delete root;
                root = 0;
            }
        }

        // return size in container. O(1)
        int size() const {
            return size(root);
        }

        // insert value in [index]. O(h)
        void insert(const T& value, int index) {
            pnode left, right;
            split(root, left, right, index);
            merge(left, left, new node(value, 0, 0, 0));
            merge(root, left, right);
        }

        // erase [index]. O(h)
        void erase(int index) {
            pnode left, right, mid;
            split(root, left, right, index);
            split(right, mid, right, 1);
            return merge(root, left, right);
        }

        T& operator [](int index) const {
            return find(root, index)->value;
        }

        void pop_back() {
            erase(size() - 1);
        }
        void pop_front() {
            erase(0);
        }

        T& back() const {
            return find(root, size() - 1)->value;
        }
        T& front() const {
            return find(root, 0)->value;
        }

        void push_back(const T& value) {
            insert(value, size());
        }
        void push_front(const T& value) {
            insert(value, 0);
        }

        // returns an array of elements. O(n)
        std::vector<T> get() const {
            std::vector<T> A(size());
            size_t index = 0;
            get(root, A, index);
            return A;
        }

        // sort elements. O(n * log n)
        void sort() {
            std::vector<T> A = get();
            std::sort(A.begin(), A.end());

            clear();
            for (size_t i = 0; i < A.size(); i++) {
                push_back(A[i]);
            }
        }
    };

    // Simple Hash Table
    template<typename T, size_t Length>
    class hashTable {
        typedef long long s64;

        s64 hash(s64 key) {
            return abs(key) % A.size();
        }

    public:
        // value, hash
        std::vector<std::vector<std::pair<T, s64>>> A;
        
        hashTable() {
            A.resize(Length);
        }
        template<typename _Ty2, size_t otherLength>
        hashTable(const hashTable<_Ty2, otherLength>& other) {
            A = other.A;
        }

        T& operator [](s64 key) {
            size_t it = hash(key);
            for (size_t i = 0; i < A[it].size(); i++) {
                if (A[it][i].second == key) {
                    return A[it][i].first;
                }
            }

            A[it].push_back(std::make_pair(T(), key));
            return A[it].back().first;
        }
        // erase element
        void erase(s64 key) {
            size_t it = hash(key);
            for (size_t i = 0; i < A[it].size(); i++) {
                if (A[it][i].second == key) {
                    A[it].erase(A[it].begin() + i);
                }
            }
        }
        // insert element
        void insert(const T& value, s64 key) {
            *this->operator[](key) = value;
        }
    };

    // Система непересекающихся множеств.
    // Сжатие пути + ранговая эвристика на основе размеров деревьев
    class dsu {
        std::vector<size_t> parent, size;

    public:

        dsu() {}
        dsu(size_t length) {
            resize(length);
        }

        void resize(size_t length) {
            parent.resize(length);
            size.resize(length);
        }

        void makeSet(size_t v) {
            parent[v] = v;
            size[v] = 1;
        }
        size_t findSet(size_t v) {
            if (v == parent[v]) {
                return v;
            }
            return parent[v] = findSet(parent[v]);
        }
        void unionSets(size_t a, size_t b) {
            a = findSet(a);
            b = findSet(b);
            if (a != b) {
                if (size[a] < size[b]) {
                    std::swap(a, b);
                }
                parent[b] = a;
                size[a] += size[b];
            }
        }
    };

    // Дерево Фенвика
    template<typename _Ty>
    class fenwick {
        _Ty* t;
        int n;
        _Ty(*function)(_Ty, _Ty);

    public:
        fenwick() {}
        fenwick(int size, _Ty(*newFunction)(_Ty, _Ty)) {
            n = size + 1;
            t = new _Ty[size + 1];
            for (int i = 0; i < n; i++) {
                t[i] = 0;
            }
            function = newFunction;
        }
        fenwick(const std::vector<_Ty>& array, _Ty(*newFunction)(_Ty, _Ty)) {
            this = fenwick(array.size(), newFunction);
            for (int i = 0; i < array.size(); i++) {
                update(i, array[i]);
            }
        }
        ~fenwick() {
            delete[] t;
        }

        _Ty get(int r) {
            r++;
            _Ty res = 0;
            for (; r > 0; r -= r & -r) {
                res += t[r];
            }
            return res;
        }

        _Ty get(int l, int r) {
            return get(r) - get(l - 1);
        }

        void update(int k, _Ty val) {
            k++;
            for (; k < n; k += k & -k) {
                t[k] += val;
            }
        }
    };

    // line with function(min/max)
    template<typename _Ty, bool(*compare)(_Ty, _Ty)>
    class lwf {
        dst::edeque<_Ty> d;
    public:
        _Ty get() {
            return d.front();
        }
        void push(_Ty value) {
            while (!d.empty() && compare(value , d.back())) {
                d.pop_back();
            }
            d.push_back(value);
        }
        void erase(_Ty key) {
            if (!d.empty() && d.front() == key) {
                d.pop_front();
            }
        }
    };
}

// Сomputational Geometry
// edouble, dot, line, circle, polygon
namespace mpg {

    long double eps = 1e-9, pi = acos(-1);

    //Superstructure over double. 
    //Compares a floating point number correctly
    class edouble {
        long double value;

    public:
        edouble() {
            value = 0.L;
        }
        template<typename T>
        edouble(const T& new_value) {
            value = static_cast<long double>(new_value);
        }

        bool operator == (const edouble& Rhs) const {
            return std::abs(value - Rhs.value) <= eps;
        }
        bool operator < (const edouble& Rhs) const {
            return value < Rhs.value - eps;
        }
        bool operator > (const edouble& Rhs) const {
            return value > Rhs.value + eps;
        }
        
        bool operator != (const edouble& Rhs) const {
            return !(*this == Rhs);
        }
        bool operator <= (const edouble& Rhs) const {
            return !(*this > Rhs);
        }
        bool operator >= (const edouble& Rhs) const {
            return !(*this < Rhs);
        }
    };

    template<typename T>
    bool operator == (const T& Lhs, const edouble& Rhs) {
        return edouble(Lhs) == Rhs;
    }
    template<typename T>
    bool operator < (const T& Lhs, const edouble& Rhs) {
        return edouble(Lhs) < Rhs;
    }
    template<typename T>
    bool operator > (const T& Lhs, const edouble& Rhs) {
        return edouble(Lhs) > Rhs;
    }

    template<typename T>
    bool operator != (const T& Lhs, const edouble& Rhs) {
        return edouble(Lhs) != Rhs;
    }
    template<typename T>
    bool operator <= (const T& Lhs, const edouble& Rhs) {
        return edouble(Lhs) <= Rhs;
    }
    template<typename T>
    bool operator >= (const T& Lhs, const edouble& Rhs) {
        return edouble(Lhs) >= Rhs;
    }

    // x, y
    struct dot {
        long double x, y;

        dot() {
            x = y = 0.L;
        }
        template<typename T>
        dot(const T& _x, const T& _y) {
            x = static_cast<long double>(_x);
            y = static_cast<long double>(_y);
        }

        // vector addition
        dot operator + (const dot& p) const {
            return dot(x + p.x, y + p.y);
        }
        dot& operator += (const dot& p) {
            return *this = *this + p;
        }

        // vector subtraction
        dot operator - (const dot& p) const {
            return dot(x - p.x, y - p.y);
        }
        dot& operator -= (const dot& p) {
            return *this = *this - p;
        }

        // multiplying a vector by a number k
        dot operator * (long double k) const {
            return dot{ x * k, y * k };
        }
        dot& operator *= (long double k) {
            return *this = *this * k;
        }

        // dividing a vector by a number k
        dot operator / (long double k) const {
            return dot{ x / k, y / k };
        }
        dot& operator /= (long double k) {
            return *this = *this / k;
        }

        // product of vectors
        long double operator % (const dot& p) const {
            return x * p.y - y * p.x;
        }
        // scalar product
        long double operator * (const dot& p) const {
            return x * p.x + y * p.y;
        }


        // vector len
        long double getLen() const {
            return sqrt(x * x + y * y);
        }
        // vector quare len
        long double getQuareLen() const {
            return x * x + y * y;
        }

        bool operator == (const dot& Rhs) const {
            return edouble(x) == Rhs.x && edouble(y) == Rhs.y;
        }
        bool operator != (const dot& Rhs) const {
            return !(*this == Rhs);
        }

        // самая левая, потом нижняя
        bool operator < (const dot& Rhs) const {
            return edouble(x) == Rhs.x ? y < Rhs.y : x < Rhs.x;
        }

        // vector normalize * mult
        dot normalize(long double mult = 1.0L) const {
            return *this / getLen() * mult;
        }
    };
    std::istream& operator >> (std::istream& input, dot& point) {
        return input >> point.x >> point.y;
    }
    std::ostream& operator << (std::ostream& output, const dot& point) {
        return output << point.x << " " << point.y;
    }

    // returns the angle between vectors
    long double getAngle(const dot& a, const dot& b) {
        return atan2(a % b, a * b);
    }
    // returns a non-negative angle between vectors
    long double getGoodAngle(const dot& a, const dot& b) {
        long double res = atan2(a % b, a * b);
        if (res < 0) {
            res += 2 * pi;
        }
        return res;
    }
    // returns a non - negative angle less than 180
    long double getVeryGoodAngle(const dot& a, const dot& b) {
        long double res = atan2(a % b, a * b);
        if (res < 0) {
            res += 2 * pi;
        }
        if (res > pi) {
            res = 2 * pi - res;
        }
        return res;
    }

    // A, B, C
    // A*x + B*y + c = 0
    struct line {
        long double a, b, c;

        line() {
            a = b = c = 0;
        }
        line(const dot& p0, const dot& p1) {
            a = p0.y - p1.y;
            b = p1.x - p0.x;
            c = -a * p0.x - b * p0.y;
        }
        template<typename T>
        line(const T& _a, const T& _b, const T& _c) {
            a = static_cast<long double>(_a);
            b = static_cast<long double>(_b);
            c = static_cast<long double>(_c);
        }

        // returns the perpendicular line through point
        line getPerp(const dot& point) const {
            long double A = -b, B = a;
            long double C = -A * point.x - B * point.y;
            return line(A, B, C);
        }

        // returns a parallel line at a distance
        // at negative distance will return a straight line from the other side
        line getParallel(long double dist) const {
            long double A = a, B = b;
            long double C = c + dist * sqrt(a * a + b * b);
            return line(A, B, C);
        }

        // returns the normalized vector line multiplied by mult
        dot getVector(long double mult = 1.0L) const {
            return dot(-b, a).normalize(mult);
        }


        // return the intersection point
        dot intersect(const line& Rhs) const {
            long double x, y;
            if (edouble(Rhs.b) != 0) {
                x = ((b * Rhs.c / Rhs.b - c) / (a - b * Rhs.a / Rhs.b));
                y = (-x * Rhs.a - Rhs.c) / Rhs.b;
            }
            else {
                x = -Rhs.c / Rhs.a;
                y = (-x * a - c) / b;
            }
            return dot(x, y);
        }

        // returns the intersection point of the perpendicular from the point
        dot perpIntersect(const dot& point) const {
            return intersect(getPerp(point));
        }


        // returns the distance between lineand point
        long double dist(const dot& point) const {
            return abs(a * point.x + b * point.y + c) / std::sqrt(a * a + b * b);
        }

        // returns the distance between two PARALLEL straight lines
        long double dist(const line& parallel_line) const {
            return abs(c - parallel_line.c) / sqrt(a * a + b * b);
        }


        // return true if lines is parallel
        bool isParallel(const line& Rhs) const {
            return edouble(a * Rhs.b - b * Rhs.a) == 0;
        }

        // return true if lines is perpendicular
        bool isPerp(const line& Rhs) {
            return edouble(a * Rhs.a + b * Rhs.b) == 0;
        }

        // return true if the point is on the line
        bool ison(const dot& point) const {
            return edouble(a * point.x + b * point.y + c) == 0;
        }
    };
    std::istream& operator >> (std::istream& input, line& Line) {
        return input >> Line.a >> Line.b >> Line.c;
    }
    std::ostream& operator << (std::ostream& output, const line& Line) {
        return output << Line.a << " " << Line.b << " " << Line.c;
    }

    // center, radius
    // (x - x0)^2 + (y - y0)^2 = r^2
    struct circle {
        dot center;
        long double radius;

        circle() {
            radius = 0;
        }
        circle(dot new_center, long double new_radius) {
            center = new_center;
            radius = new_radius;
        }
        // constructor of a circle by three points lying on it
        circle(const dot& p0, const dot& p1, const dot& p2) {
            dot point1 = p1 + (p0 - p1) * .5L;
            dot point2 = p2 + (p0 - p2) * .5L;

            line a(p0, p1), b(p0, p2);

            line l1(a.getPerp(point1)), l2(b.getPerp(point2));

            center = l1.intersect(l2);
            radius = (center - p0).getLen();
        }

        // returns a point on a circle
        // counterclockwise angle
        dot point(long double angle) {
            return center + dot(cos(angle), sin(angle)) * radius;
        }

        // returns the intersection points of two circles
        // if inf == true, then the circles have infinitely many intersection points
        std::vector<dot> intersect(const circle& Rhs, bool& inf) const {
            if (Rhs.center == center) {
                inf = Rhs.radius == radius;
                return {};
            }
            else {
                inf = false;

                dot vector(Rhs.center - center);
                std::vector<dot> A = circle(dot(0, 0), radius).intersect(
                    line(-2 * vector.x,
                        -2 * vector.y,
                        vector.x * vector.x + vector.y * vector.y + radius * radius - Rhs.radius * Rhs.radius));

                for (size_t i = 0; i < A.size(); i++) {
                    A[i] += center;
                }
                return A;
            }
        }

        // returns the intersection points of a line and a circle
        std::vector<dot> intersect(const line& Rhs) const {
            dot perp = Rhs.perpIntersect(center),
            delta = center - perp;

            long double quareRadius = radius * radius;
            long double quareDist = delta.getQuareLen();

            if (edouble(quareRadius) > quareDist) { // two points
                long double len = sqrt(quareRadius - quareDist);
                dot vector(Rhs.getVector(len));
                return { dot(perp + vector), dot(perp - vector) };
            }
            else if (edouble(quareRadius) == quareDist) { // one point
                return { perp };
            }
            else { // none point
                return {};
            }
        }

        // returns true if point lies on a circle
        bool ison(const dot& point) const {
            return edouble((center.x - point.x) * (center.x - point.x) +
                (center.y - point.y) * (center.y - point.y)) == radius * radius;
        }
    };
    std::istream& operator >> (std::istream& input, circle& Circle) {
        return input >> Circle.center >> Circle.radius;
    }
    std::ostream& operator << (std::ostream& output, const circle& Circle) {
        return output << Circle.center << " " << Circle.radius;
    }

    // сравнивает две точки для построения выпуклой оболочки
    bool compareConvexHull(const dot& Lhs, const dot& Rhs) {
        edouble vectorProduct = Lhs % Rhs;
        return vectorProduct == 0 ? edouble(Lhs.getQuareLen()) < Rhs.getQuareLen() : vectorProduct < 0;
    }

    // points
    class polygon {
        std::vector<dot> Points;

    public:
        polygon(){}
        polygon(const std::vector<dot>& newPoints) {
            Points = newPoints;
        }

        // Построение выпуклой оболочки за O(n * log n)
        // 39/41 tests. WA
        std::vector<dot> buildConvexHull() const{
            std::vector<dot> A = Points, result;
            sort(A.begin(), A.end());
            dot start = A[0];
            for (int i = 0; i < A.size(); i++) {
                A[i] -= start;
            }
            sort(A.begin() + 1, A.end(), compareConvexHull);

            result.push_back(A[0]);
            for (int i = 1; i < A.size(); i++) {
                while (result.size() > 1 && 
                    edouble((A[i] - result[result.size() - 2]) % (result.back() - result[result.size() - 2])) <= 0) {
                    result.pop_back();
                }
                result.push_back(A[i]);
            }

            for (int i = 0; i < result.size(); i++) {
                result[i] += start;
            }
            return result;
        }

        // Возвращает площадь многоугольника
        long double area() const{
            long double result = 0.L;
            for (int i = 0; i < Points.size(); i++)
            {
                dot p1 = i ? Points[i - 1] : Points.back(),
                    p2 = Points[i];
                result += (p1.x - p2.x) * (p1.y + p2.y);
            }
            return std::abs(result) * .5L;
        }

        // Возвращает периметр многоугольника
        long double perim() const {
            long double result = (Points[0] - Points.back()).getLen();
            for (int i = 1; i < Points.size(); i++) {
                result += (Points[i] - Points[i - 1]).getLen();
            }
            return result;
        }
    };
}

// Long Arithmetic
// elong
namespace emt {

    static const long long long_length = 18,
        long_base = alg::epow(10LL, long_length), 
        long_base_expansion = std::sqrt(long_base);

    // простое длинное число. НЕ ИСПОЛЬЗОВАТЬ!!!
    struct simpleLong {
        dst::edeque<long long> digits;

        simpleLong(dst::edeque<long long> newDigits) {
            digits = newDigits;
        }
        simpleLong(const std::string& str) {
            int i;
            for (i = str.size(); i >= long_length; i -= long_length) {
                digits.push_back(atoll(str.substr(i - long_length, long_length).c_str()));
            }
            if (i > 0) {
                digits.push_back(atoll(str.substr(0, i).c_str()));
            }
        }
        simpleLong() {}

        simpleLong& remove_leading_zeros() {
            while (!digits.empty() && digits.back() == 0) {
                digits.pop_back();
            }
            return *this;
        }

        enum class OP{
            less,
            equally,
            more,
        };

        template<typename T>
        OP comp(const T& a,const T& b) const{
            return a < b ? OP::less : OP::more;
        }

        // compare two numbers
        OP compare(const simpleLong& Rhs) const{
            if (digits.size() != Rhs.digits.size()) {
                return comp(digits.size(), Rhs.digits.size());
            }
            else {
                int i = digits.size() - 1;
                while (i >= 0 && digits[i] == Rhs.digits[i]) {
                    i--;
                }
                return i >= 0 ? comp(digits[i], Rhs.digits[i]) : OP::equally;
            }
        }
        
        bool operator < (const simpleLong& Rhs) const{
            return compare(Rhs) == OP::less;
        }
        bool operator == (const simpleLong& Rhs) const {
            return compare(Rhs) == OP::equally;
        }
        bool operator > (const simpleLong& Rhs) const {
            return compare(Rhs) == OP::more;
        }

        // return a + b
        simpleLong operator + (const simpleLong& added) const{
            simpleLong result = *this;
            bool k = 0;
            long long i = 0;
            for (i = 0; i < std::max(result.digits.size(), added.digits.size()) || k != 0; i++) {
                if (i == result.digits.size()) {
                    result.digits.push_back(0);
                }
                result.digits[i] += k + (i < added.digits.size() ? added.digits[i] : 0);
                k = result.digits[i] >= long_base;
                if (k != 0) {
                    result.digits[i] -= long_base;
                }
            }
            return result;
        }

        // return minuend - subtrahend
        simpleLong operator - (const simpleLong& subtrahend) const{
            simpleLong result = *this;

            bool k = 0;
            for (int i = 0; i < subtrahend.digits.size() || k != 0; i++) {
                result.digits[i] -= k + (i < subtrahend.digits.size() ? subtrahend.digits[i] : 0);
                k = result.digits[i] < 0;
                if (k != 0) {
                    result.digits[i] += long_base;
                }
            }
            return result.remove_leading_zeros();
        }

        // расширяет длину числа, уменьшая его модуль
        simpleLong expansion() const {
            simpleLong result;
            result.digits.resize(digits.size() * 2);
            for (int i = 0; i < digits.size(); i++) {
                result.digits[i * 2] = digits[i] % long_base_expansion;
                result.digits[i * 2 + 1] = digits[i] / long_base_expansion;
            }
            return result;
        }

        // сужает длину числа, увеличивая его модуль
        simpleLong reduction() const {
            simpleLong result;
            result.digits.resize(digits.size() / 2);
            for (int i = 0; i < result.digits.size(); i++) {
                result.digits[i] = digits[i * 2] + digits[i * 2 + 1] * long_base_expansion;
            }
            return result;
        }

        // return a * b
        simpleLong operator * (const simpleLong& mult) const {
            simpleLong result, mult1 = this->expansion(), mult2 = mult.expansion();
            result.digits.resize(2 * (digits.size() + mult.digits.size() + 1));

            long long c, k;
            int i, j;
            for (i = 0; i < mult1.digits.size(); i++) {
                c = 0;
                for (j = 0; j < mult2.digits.size() || c != 0; j++) {
                    k = result.digits[i + j] + (j < mult2.digits.size() ? mult1.digits[i] * 1LL * mult2.digits[j] : 0) + c;
                    result.digits[i + j] = k % long_base_expansion;
                    c = k / long_base_expansion;
                }
            }
            return result.reduction().remove_leading_zeros();
        }

        // return dividend / divider
        simpleLong operator / (const simpleLong& divider) const {
            simpleLong result, value, temp, t;

            int i = digits.size() - 1;
            while (i >= 0 && value < divider) {
                value.digits.push_front(digits[i]);
                i--;
            }
            value.digits.pop_front();
            i++;
            for (; i >= 0; i--) {
                value.remove_leading_zeros();

                value.digits.push_front(digits[i]);

                long long l = 0, r = long_base;
                while (l < r - 1) {
                    long long m = (l + r) / 2;

                    temp.digits.push_back(m);

                    t = divider * temp;

                    temp.digits.pop_back();

                    if (value < t) {
                        r = m;
                    }
                    else {
                        l = m;
                    }
                }

                result.digits.push_front(l);

                temp.digits.push_back(l);

                value = value - divider * temp;

                temp.digits.pop_back();
            }
            return result;
        }

        // return dividend % divider
        simpleLong operator % (const simpleLong& divider) const{
            simpleLong result, temp, t;

            int i = digits.size() - 1;
            while (i >= 0 && result < divider) {
                result.digits.push_front(digits[i]);
                i--;
            }
            result.digits.pop_front();
            i++;
            for (; i >= 0; i--) {
                result.remove_leading_zeros();
                result.digits.push_front(digits[i]);

                long long l = 0, r = long_base;
                while (l < r - 1) {
                    long long m = (l + r) / 2;

                    temp.digits.push_back(m);

                    t = divider * temp;

                    temp.digits.pop_back();

                    if (result < t) {
                        r = m;
                    }
                    else {
                        l = m;
                    }
                }

                temp.digits.push_back(l);

                result = result - divider * temp;

                temp.digits.pop_back();
            }
            return result;
        }
    };

    std::ostream& operator << (std::ostream& out, const simpleLong& var) {
        out << var.digits.back();
        char k = out.fill('0');
        for (int i = var.digits.size() - 2; i >= 0; i--) {
            out << std::setw(long_length) << var.digits[i];
        }
        out.fill(k);
        return out;
    }

    // signed long integer
    class elong {
    
        simpleLong value;
        bool isNegative;

        // update
        elong& update() {
            if (value.digits.empty()) {
                value.digits.push_back(0);
            }
            return *this;
        }

        elong(const simpleLong& newValue, bool newNegative) {
            value = newValue;
            isNegative = newNegative;
        }

        void inc() {
            value.digits.front()++;
            if (value.digits.front() == long_base) {
                value.digits.front() = 0;

                int i;
                for (i = 1; i < value.digits.size() && value.digits[i] + 1LL == long_base; i++) {
                    value.digits[i] = 0;
                }
                if (i == value.digits.size()) {
                    value.digits.push_back(1);
                }
                else {
                    value.digits[i]++;
                }
            }
        }
        void dec() {
            value.digits.front()--;
            if (value.digits.front() == -1) {
                value.digits.front() = 0;

                long long i;
                for (i = 1; i < value.digits.size() && value.digits[i] == 0; i++) {
                    value.digits[i] = long_base - 1;
                }

                if (i < value.digits.size()) {
                    value.digits[i]--;
                    if (value.digits.back() == 0) {
                        value.digits.pop_back();
                    }
                }
                else {
                    isNegative = !isNegative;
                    value.digits[0] = 1;
                }
            }

            if (value.digits.size() == 1 && value.digits.front() == 0) {
                isNegative = false;
            }
        }

        void signedInc() {
            if (isNegative) {
                dec();
            }
            else {
                inc();
            }
        }
        void signedDec() {
            if (isNegative) {
                inc();
            }
            else {
                dec();
            }
        }
        
    public:

        // default constructor
        elong() {
            isNegative = false;
            value.digits.push_back(0);
        };

        // convert constructors[

        elong(const std::string& String) {
            isNegative = String[0] == '-';
            value = isNegative ? simpleLong(String.substr(1)) : simpleLong(String);
        }
        elong(const char* String) {
            *this = elong(std::string(String));
        }

        template<typename _Ty>
        elong(_Ty number) {
            isNegative = number < 0;
            number = trunc(number);

            if (number == 0) {
                value.digits.push_back(0);
            }
            else {
                number = isNegative ? number * -1 : number;
                while (number > 0) {
                    value.digits.push_back(number - trunc(number / long_base) * long_base);
                    number = trunc(number / long_base);
                }
            }
        }
        // ]

        //bool operators [

        bool operator == (const elong& Rhs) const {
            return isNegative != Rhs.isNegative ? false : value == Rhs.value;
        }
        bool operator != (const elong& Rhs) const {
            return !(*this == Rhs);
        }
        bool operator < (const elong& Rhs) const {
            if (isNegative) {
                return Rhs.isNegative ? -*this > -Rhs : true;
            }
            else {
                return Rhs.isNegative ? false : value < Rhs.value;
            }
        }
        bool operator > (const elong& Rhs) const {
            if (isNegative) {
                return Rhs.isNegative ? -*this < -Rhs : false;
            }
            else {
                return Rhs.isNegative ? true : value > Rhs.value;
            }
        }

        bool operator <= (const elong& Rhs) const {
            return !(*this > Rhs);
        }
        bool operator >= (const elong& Rhs) const {
            return !(*this < Rhs);
        }
        // ]

        elong operator - () const {
            return elong(value, !isNegative);
        }

        // operators [

        elong& operator ++() {
            signedInc();
            return *this;
        }
        elong& operator ++(int) {
            signedInc();
            return *this;
        }
        elong& operator --() {
            signedDec();
            return *this;
        }
        elong& operator --(int) {
            signedDec();
            return *this;
        }

        elong operator + (const elong& added) const {
            if (isNegative) {
                return added.isNegative ? -(-*this + (-added)) : added - (-*this);
            }
            else {
                return added.isNegative ? *this - (-added) : elong(value + added.value, false);
            }
        }
        elong& operator += (const elong& added) {
            return *this = *this + added;
        }

        elong operator - (const elong& subtrahend) const {
            if (subtrahend.isNegative) {
                return *this + (-subtrahend);
            }
            else if (isNegative) {
                return -(-*this + subtrahend);
            }
            else {
                return *this < subtrahend ? -(subtrahend - *this) : elong(value - subtrahend.value, false).update();
            }
        }
        elong& operator -= (const elong& subtrahend) {
            return *this = *this - subtrahend;
        }

        elong operator * (const elong& multiplied) const {
            return elong(value * multiplied.value, isNegative != multiplied.isNegative).update();
        }
        elong& operator *= (const elong& multiplied) {
            return *this = *this * multiplied;
        }
        
        elong operator / (const elong& divider) const {
            return elong(value / divider.value, isNegative != divider.isNegative).update();
        }
        elong& operator /= (const elong& divider) {
            return *this = *this / divider;
        }

        elong operator % (const elong& divider) const {
            return divider == 2 ? elong(value.digits[0] & 1) : elong(value % divider.value, false).update();
        }
        elong& operator %= (const elong& divider) {
            return *this = *this % divider;
        }
        // ]

        elong sqrt() const {
            return alg::sqrt(*this);
        }

        // friend function :)
        friend std::ostream& operator << (std::ostream& out, const elong& var);
    };

    template<typename _Ty>
    elong operator + (const _Ty& a, const elong& b) {
        return elong(a) + b;
    }
    template<typename _Ty>
    elong operator - (const _Ty& a, const elong& b) {
        return elong(a) - b;
    }
    template<typename _Ty>
    elong operator * (const _Ty& a, const elong& b) {
        return elong(a) * b;
    }
    template<typename _Ty>
    elong operator / (const _Ty& a, const elong& b) {
        return elong(a) / b;
    }
    template<typename _Ty>
    elong operator % (const _Ty& a, const elong& b) {
        return elong(a) % b;
    }
    

    // output signed long integer
    std::ostream& operator << (std::ostream& out, const elong& var) {
        if (var.isNegative) {
            out << '-';
        }
        return out << var.value;
    }
    // input signed long integer
    std::istream& operator >> (std::istream& input, elong& var) {
        std::string str;
        input >> str;
        var = elong(str);
        return input;
    }
    
    elong tree_factor(int l, int r) {
        if (l == 4997 && r == 5000) {
            int k = 0;
            k++;
        }
        if (l < r - 1) {
            int m = (l + r) >> 1;
            return tree_factor(l, m) * tree_factor(m + 1, r);
        }
        else if (r - l == 1) {
            return 1LL * l * r;
        }
        else if (l == r) {
            return l;
        }
        else {
            return 1;
        }
    }

    elong factorial(int n) {
        return tree_factor(2, n);
    }
}

// Timepiece
// etime
namespace tmp{
    class etime {
        int time_begin;
        std::string text;

    public:
        etime(std::string newText) {
            start();
            text = newText;
        }

        void start() {
            time_begin = clock();
        }
        int operator *()const {
            return clock() - time_begin;
        }
        friend std::ostream& operator << (std::ostream& out, const etime& t);
    };
    std::ostream& operator << (std::ostream& out, const etime& t) {
        return out << "\n" + t.text + " " << *t << "\n\n";
    }
};



#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long ull;

int main() {
    ifstream cin("input.txt");
    

    return 0;
}