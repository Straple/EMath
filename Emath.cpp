﻿// Emath 1.1 | Developed by Mr_Straym

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

//std::ofstream debugLog("debug.txt");

// Utils: eclock, random, vector addition, clamp, range, roundTwo
namespace utl {

    using s8 = char;
    using u8 = unsigned char;
    using s16 = short;
    using u16 = unsigned short;
    using s32 = int;
    using u32 = unsigned int;
    using s64 = long long;
    using u64 = unsigned long long;

    // Макросы

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)
// возвращает значение между left и right
#define clamp(min, value, max) (value > max ? max : (value < min ? min : value))
// остаток от деления на степень 2
#define twoRemainder(value, degree) (value - ((value >> degree) << degree))

    class eclock {
        typedef std::chrono::steady_clock::time_point timePoint;

        timePoint begin;

        timePoint newNow() const {
            return std::chrono::high_resolution_clock::now();
        }

    public:
        eclock() {
            reset();
        }
        void reset() {
            begin = newNow();
        }

        // ms
        double count() const {
            std::chrono::duration<double> time = newNow() - begin;
            return time.count() * 1000;
        }
    };
    std::ostream& operator <<(std::ostream& output, const eclock& clock) {
        return output << clock.count() << "ms";
    }

    // vector addition
    template<typename T>
    std::vector<T> operator + (const std::vector<T>& a, const std::vector<T>& b) {
        std::vector<T> ret = a;
        ret.resize(a.size() + b.size());
        for (size_t i = 0; i < b.size(); i++) {
            ret[i + a.size()] = b[i];
        }
        return ret;
    }

    // Add to vector
    template<typename T>
    std::vector<T>& operator += (std::vector<T>& a, const std::vector<T>& b) {
        int len = a.size();
        a.resize(len + b.size());
        for (int i = 0; i < b.size(); i++) {
            a[i + len] = b[i];
        }
        return a;
    }

    template<typename T>
    T range(const T& left, const T& right) {
        static std::mt19937_64 random(42);
        return random() % (right - left + 1) + left;
    }

    // округление до степени 2
    size_t roundTwo(size_t v) {
        v--;

        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v |= v >> 32;

        return v + 1;
    }
}
using namespace utl;

// Data Structures: trie, edeque, segTree, hashTable, container, dsu, fenwick, lwf, PairingHeap
namespace dst {

    // шаблон
    struct trie {
        struct node;
        using pnode = node*;
        struct node {
            pnode Next[26];
            bool isEnd;
            node() {
                for (int i = 0; i < 26; i++) {
                    Next[i] = 0;
                }
                isEnd = false;
            }

            // лист?
            bool isLeaf() const {
                int i;
                // пока не вышли за пределы и следующего нет
                for (i = 0; i < 26 && Next[i] == 0; i++) {}
                return i == 26;
            }
        };

        pnode root = new node();

#define next()\
temp->Next[str[k] - 'a']


        // вернет 0 если такой строки нет
        pnode find(const std::string& str) {
            pnode temp = root;
            int k = 0;
            while (k < str.size() && next() != 0) {
                temp = next();
                k++;
            }
            // если мы не дошли до конца строки: 0 иначе temp
            return k < str.size() ? 0 : temp;
        }

        // добавляет строку и значение
        void insert(const std::string& str) {
            pnode temp = root;
            int k = 0;
            while (k < str.size() && next() != 0) {
                temp = next();
                k++;
            }
            while (k < str.size()) {
                temp = next() = new node();
                k++;
            }
            temp->isEnd = true;
        }
    };

    // deque on a circular array. 
    template<typename T>
    class edeque {
        // variables
        T* A; // массив
        int length_A; // размер А
        int length; // размер массива
        int head; // голова

        // проверяет нужно ли увеличить размер
        bool overflow() const {
            return length >= length_A;
        }
        // проверяет нужно ли уменьшить размер
        bool underflow() const {
            return (length << 1) < length_A;
        }


        // берет размер массива с "запасом"
        int sizeInc(int Size) const {
            return 1.4 * Size + 1.0;
        }


        // увеличивает индекс на 1
        int inc(int index) const {
            index++;
            return index < length_A ? index : 0;
        }
        // уменьшает индекс на 1
        int dec(int index) const {
            index--;
            return index >= 0 ? index : length_A - 1;
        }
        // дает индекс элемента в A
        int getIndex(int index) const {
            index += head;
            return index + (index < length_A ? 0 : -length_A);
        }


        // обнуляет структуру
        void zeroing() {
            A = 0;
            length_A = length = head = 0;
        }

        // Перестраивает массив. Урезаная версия resize для update/overflow/underflow
        void rebuild(int newSize) {
            // rebuild array
            {
                T* temp = new T[newSize];
                int i = 0;
                while (i < length) {
                    temp[i++] = std::move(A[head]);
                    head = inc(head);
                }
                delete[] A;
                A = temp;
            }
            head = 0;
            length_A = newSize;
        }

        // выделяет структуру с размером newSize
        void allocate(int newSize) {
            A = new T[length_A = sizeInc(length = newSize)];
            head = 0;
        }

        // полностью копирует структуру, но не очищает память
        void copy(const edeque& source) {
            A = new T[length_A = source.length_A];
            length = source.length;
            head = source.head;

            for (int i = 0; i < length_A; i++) {
                A[i] = source.A[i];
            }
        }

        // move structor
        void move(edeque& source) {
            // полностью скопировать
            A = source.A;
            head = source.head;
            length = source.length;
            length_A = source.length_A;
            // обнулить source
            source.zeroing();
        }

    public:

        // конструктор по умолчанию
        edeque() {
            zeroing();
        }
        // конструктор копирования
        edeque(const edeque& source) {
            copy(source);
        }
        // конструктор перемещения
        edeque(edeque&& source) noexcept {
            move(source);
        }
        // деструктор
        ~edeque() {
            delete[] A;
        }

        edeque& operator = (const edeque& source) {
            if (A != source.A) {
                delete[] A;
                copy(source);
            }
            return *this;
        }
        edeque& operator = (edeque&& source) noexcept {
            if (A != source.A) {
                delete[] A;
                move(source);
            }
            return *this;
        }

        // конструктор заполнения/выделения
        edeque(int newSize, const T fillValue = T()) {
            allocate(newSize);
            for (int i = 0; i < length; i++) {
                A[i] = fillValue;
            }
        }
        template<typename it>
        edeque(const it& begin, const it& end) {
            zeroing();
            it temp = begin;
            while (temp != end) {
                push_back(static_cast<T>(*temp++));
            }
        }

    private:
        template<typename array_t>
        void copyArray(const array_t& source) {
            auto it = source.begin();
            for (int i = 0; i < length; i++, it++) {
                A[i] = static_cast<T>(*it);
            }
        }

    public:

        edeque(std::initializer_list<T> list) {
            allocate(list.size());
            copyArray(list);
        }
        template<typename array_t>
        edeque(const array_t& source) {
            allocate(source.size());
            copyArray(source);
        }

    private:
        enum class OP {
            less,
            equally,
            more,
        };

        template<typename value_t>
        OP comp(const value_t& a, const value_t& b) const {
            return a < b ? OP::less : OP::more;
        }

        // compare two arrays
        OP compare(const edeque& Rhs) const {
            if (length != Rhs.length) {
                return comp(length, Rhs.length);
            }
            else {
                int i = 0;
                while (i < length && this->operator[](i) == Rhs[i]) { i++; }

                return i < length ? comp(this->operator[](i), Rhs[i]) : OP::equally;
            }
        }

    public:
        // compares

        bool operator < (const edeque& Rhs) const {
            return compare(Rhs) == OP::less;
        }
        bool operator == (const edeque& Rhs) const {
            return compare(Rhs) == OP::equally;
        }
        bool operator > (const edeque& Rhs) const {
            return compare(Rhs) == OP::more;
        }

        bool operator != (const edeque& Rhs) const {
            return !(*this == Rhs);
        }
        bool operator <= (const edeque& Rhs) const {
            return !(*this > Rhs);
        }
        bool operator >= (const edeque& Rhs) const {
            return !(*this < Rhs);
        }

        // освобождает память
        void clear() {
            delete[] A;
            zeroing();
        }

        // изменяет размер массива, добавляя недостающие элементы
        void resize(int newSize) {
            // array copy
            {
                T* temp = new T[sizeInc(newSize)];
                int i = 0, Min = min(newSize, length);
                while (i < Min) {
                    temp[i++] = std::move(A[head]);
                    head = inc(head);
                }
                while (i < newSize) {
                    temp[i++] = T();
                }
                delete[] A;
                A = temp;
            }
            length_A = sizeInc(length = newSize);
            head = 0;
        }

        size_t size() const {
            return length;
        }
        bool empty() const {
            return length == 0;
        }

    private:
        // Если нужно обновляет размер массива
        void updateOverflow() {
            if (overflow()) {
                rebuild(sizeInc(length));
            }
        }

    public:
        void push_back(const T& value) {
            updateOverflow();

            A[getIndex(length)] = value;
            length++;
        }
        void push_back(T&& value) {
            updateOverflow();

            A[getIndex(length)] = std::move(value);
            length++;
        }

        void push_front(const T& value) {
            updateOverflow();

            A[head = dec(head)] = value;
            length++;
        }
        void push_front(T&& value) {
            updateOverflow();

            A[head = dec(head)] = std::move(value);
            length++;
        }

    private:
        // Если нужно обновляет размер массива
        void updateUnderflow() {
            if (underflow()) {
                rebuild(sizeInc(length));
            }
        }

    public:
        void pop_back() {
            updateUnderflow();
            length--;
        }
        void pop_front() {
            updateUnderflow();
            head = inc(head);
            length--;
        }

        T& operator [](int index) {
            return A[getIndex(index)];
        }
        const T& operator[](int index) const {
            return A[getIndex(index)];
        }

        T& back() {
            return A[getIndex(length - 1)];
        }
        const T& back() const {
            return A[getIndex(length - 1)];
        }
        T& front() {
            return A[head];
        }
        const T& front() const {
            return A[head];
        }

        class iterator {
            int index; // индекс элемента в структуре
            edeque* source; // ресурс структуры
        public:
            // default consturctor
            iterator() {}
            iterator(int index, edeque* source) {
                this->index = index;
                this->source = source;
            }

            iterator& operator ++() {
                index++;
                return *this;
            }
            iterator operator ++(int) {
                auto temp = *this;
                index++;
                return temp;
            }
            iterator& operator --() {
                index--;
                return *this;
            }
            iterator operator --(int) {
                auto temp = *this;
                index--;
                return temp;
            }


            iterator operator + (int add) const {
                return iterator(index + add, source);
            }
            iterator& operator += (int add) {
                index += add;
                return *this;
            }

            iterator operator - (int sub) const {
                return iterator(index - sub, source);
            }
            iterator& operator -= (int sub) {
                index -= sub;
                return *this;
            }

            T& operator *() {
                return source->operator[](index);
            }
            const T& operator *() const {
                return source->operator[](index);
            }
            T* operator ->() {
                return &**this;
            }
            const T* operator ->() const {
                return &**this;
            }

            bool operator == (const iterator& Rhs) const {
                return index == Rhs.index && source == Rhs.source;
            }
            bool operator != (const iterator& Rhs) const {
                return !(*this == Rhs);
            }

            friend void edeque::insert(const iterator& where, const T& value);
            friend void edeque::erase(const iterator& where);
        };

        iterator begin() {
            return iterator(0, this);
        }
        iterator end() {
            return iterator(length, this);
        }

    private:
        // O(n / 2)
        void insert(const int index, const T& value) {
            // индекс находиться в правой половине
            if (index >= (length >> 1)) {
                push_back(value);
                for (int i = length - 1; i != index; i--) {
                    std::swap(this->operator[](i), this->operator[](i - 1));
                }
            }
            else {
                push_front(value);
                for (int i = 0; i != index; i++) {
                    std::swap(this->operator[](i), this->operator[](i + 1));
                }
            }
        }
        void insert(const int index, T&& value) {
            // индекс находиться в правой половине
            if (index >= (length >> 1)) {
                push_back(std::move(value));
                for (int i = length - 1; i != index; i--) {
                    std::swap(this->operator[](i), this->operator[](i - 1));
                }
            }
            else {
                push_front(std::move(value));
                for (int i = 0; i != index; i++) {
                    std::swap(this->operator[](i), this->operator[](i + 1));
                }
            }
        }

        // O(n / 2)
        void erase(int index) {
            // индекс находиться в правой половине
            if (index >= (length >> 1)) {
                while (index < length - 1) {
                    std::swap(this->operator[](index), this->operator[](index + 1));
                    index++;
                }
                pop_back();
            }
            else {
                while (index > 0) {
                    std::swap(this->operator[](index), this->operator[](index - 1));
                    index--;
                }
                pop_front();
            }
        }

    public:
        // O(n / 2)
        void insert(const iterator& where, const T& value) {
            insert(where.index, value);
        }
        // O(n / 2)
        void insert(const iterator& where, T&& value) {
            insert(where.index, std::move(value));
        }

        // O(n / 2)
        void erase(const iterator& where) {
            erase(where.index);
        }

        explicit operator std::vector<T>() const {
            std::vector<T> res(length);
            for (int i = 0; i < length; i++) {
                res[i] = this->operator[](i);
            }
            return res;
        }
    };

    /* matrix
    *
    *   col
    * +-----+
    * |#####| r
    * |#####| o
    * |#####| w
    * +-----+
    */
    template<typename T>
    class matrix {

        T* A; // память

        // обнуление
        void zeroing() {
            A = 0;
            rowLen = colLen = 0;
        }

        void copy(const matrix& source) {
            colLen = source.colLen;
            rowLen = source.rowLen;
            size_t n = colLen * rowLen;
            A = new T[n];
            for (int i = 0; i < n; i++) {
                A[i] = source.A[i];
            }
        }
        void move(matrix& source) {
            A = source.A;
            rowLen = source.rowLen;
            colLen = source.colLen;

            source.zeroing();
        }
    public:

        size_t rowLen; // колво строк
        size_t colLen; // колво столбцов

        matrix() {
            zeroing();
        }
        matrix(size_t rowLen, size_t colLen) {
            this->colLen = colLen;
            this->rowLen = rowLen;
            size_t n = colLen * rowLen;
            A = new T[n];
            for (int i = 0; i < n; i++) {
                A[i] = 0;
            }
        }
        matrix(std::initializer_list<std::initializer_list<T>> lists) {
            rowLen = lists.size();
            colLen = (*lists.begin()).size();
            A = new T[rowLen * colLen];
            T* temp = A;
            for (auto& it : lists) {
                for (auto& jt : it) {
                    *temp = jt;
                    temp++;
                }
            }
        }
        matrix(const matrix& source) {
            copy(source);
        }
        matrix(matrix&& source) noexcept {
            move(source);
        }
        ~matrix() {
            delete[] A;
        }

        void clear() {
            delete[] A;
            A = 0;
            colLen = rowLen = 0;
        }

        matrix& operator = (const matrix& source) {
            if (A != source.A) {
                delete[] A;
                copy(source);
            }
            return *this;
        }
        matrix& operator = (matrix&& source) noexcept {
            if (A != source.A) {
                delete[] A;
                move(source);
            }
            return *this;
        }

        T* operator [](size_t row) const {
            return A + row * colLen;
        }

        // O(n^3)
        matrix operator * (const matrix& mult) const {
            matrix c(rowLen, mult.colLen);

            if (colLen != mult.rowLen) {
                throw "error: bad matrix";
            }
            for (int i = 0; i < rowLen; i++) {
                for (int j = 0; j < mult.colLen; j++) {
                    for (int k = 0; k < mult.rowLen; k++) {
                        c[i][j] += (*this)[i][k] * mult[k][j];
                    }
                }
            }
            return c;
        }
        // O(n^2)
        matrix operator + (const matrix& add) const {
            if (colLen != add.colLen || rowLen != add.rowLen) {
                throw "error: matrix not equally";
            }
            matrix c(rowLen, colLen);
            for (int i = 0; i < rowLen; i++) {
                for (int j = 0; j < colLen; j++) {
                    c[i][j] = (*this)[i][j] + add[i][j];
                }
            }
            return c;
        }
        // O(n^2)
        matrix operator - (const matrix& sub) const {
            if (colLen != sub.colLen || rowLen != sub.rowLen) {
                throw "error: matrix not equally";
            }
            matrix c(rowLen, colLen);
            for (int i = 0; i < rowLen; i++) {
                for (int j = 0; j < colLen; j++) {
                    c[i][j] = (*this)[i][j] - sub[i][j];
                }
            }
            return c;
        }
        // O(n^2)
        matrix operator % (const T& mod) const {
            matrix c = *this;
            for (int i = 0; i < rowLen; i++) {
                for (int j = 0; j < colLen; j++) {
                    c[i][j] %= mod;
                }
            }
            return c;
        }

        matrix& operator *= (const matrix& mult) {
            return *this = *this * mult;
        }
        matrix& operator += (const matrix& add) {
            return *this = *this + add;
        }
        matrix& operator -= (const matrix& sub) {
            return *this = *this - sub;
        }
        matrix& operator %= (const T& mod) {
            return *this = *this % mod;
        }
    };

    // Дерево отрезков
    // взятия на отрезке и обновление элемента
    // calc - функция подсчета
    template<typename T, T(*calc)(T, T)>
    class segTree {
        std::vector<T> tree;
        int length;

        void build(int v, int tl, int tr, std::vector<T>& A) {
            if (tl == tr) {
                tree[v] = A[tl];
            }
            else {
                int tm = (tl + tr) / 2;
                build(v * 2, tl, tm, A);
                build(v * 2 + 1, tm + 1, tr, A);
                tree[v] = calc(tree[v * 2], tree[v * 2 + 1]);
            }
        }

        T get(int v, int tl, int tr, int l, int r) {
            if (tl == l && tr == r) {
                return tree[v];
            }
            else {
                int tm = (tl + tr) / 2;
                if (r <= tm) {
                    return get(v * 2, tl, tm, l, r);
                }
                else if (l > tm) {
                    return get(v * 2 + 1, tm + 1, tr, l, r);
                }
                return calc(get(2 * v, tl, tm, l, tm), get(2 * v + 1, tm + 1, tr, tm + 1, r));
            }
        }

        void update(int v, int tl, int tr, int pos, T val) {
            if (tl == tr) {
                tree[v] = val;
            }
            else {
                int tm = (tl + tr) / 2;
                if (pos <= tm) {
                    update(v * 2, tl, tm, pos, val);
                }
                else {
                    update(v * 2 + 1, tm + 1, tr, pos, val);
                }
                tree[v] = calc(tree[v * 2], tree[v * 2 + 1]);
            }
        }

        void resize(int length) {
            this->length = length;
            tree.resize(4 * length);
        }

    public:
        segTree() {}
        segTree(int length) {
            resize(length);
        }
        segTree(std::vector<T>& A) {
            resize(A.size());
            build(1, 0, A.size() - 1, A);
        }

        T get(int left, int right) {
            return get(1, 0, length - 1, left, right);
        }
        void update(int pos, T value) {
            update(1, 0, length - 1, pos, value);
        }
    };

    // Дерево отрезков с несоглассованными поддеревьями
    // взятие на отрезке и обновление на отрезке
    // modify - функция обновления на отрезке, calc - функция подсчета
    template<typename T, T(*modify)(T, T), T(*calc)(T, T)>
    class segTreeM {
        struct node {
            T value, add;
            node() {
                add = value = T();
            }
            node(T value, T add) {
                this->value = value;
                this->add = add;
            }
        };

        std::vector<node> tree;
        int length;

        void build(int v, int tl, int tr, const std::vector<T>& A) {
            if (tl == tr) {
                tree[v] = node(A[tl], tree[v].add);
            }
            else {
                int tm = (tl + tr) / 2;
                build(v * 2, tl, tm, A);
                build(v * 2 + 1, tm + 1, tr, A);
                tree[v] = node(calc(tree[v * 2].value, tree[v * 2 + 1].value), tree[v].add);
            }
        }

        // обновление на отрезке l...r value
        void update(int v, int tl, int tr, int l, int r, T& val) {
            if (tr < l || r < tl) {
                // отрезок не подходит
            }
            else if (l <= tl && tr <= r) {
                tree[v] = node(modify(val, tree[v].value), modify(val, tree[v].add));
            }
            else {
                int tm = (tl + tr) / 2;
                update(v * 2, tl, tm, l, r, val);
                update(v * 2 + 1, tm + 1, tr, l, r, val);
                tree[v].value = modify(calc(tree[v * 2].value, tree[v * 2 + 1].value), tree[v].add);
            }
        }

        // возвращает значение на отрезке
        T get(int v, int tl, int tr, int l, int r) {
            if (l == tl && r == tr) {
                return tree[v].value;
            }
            else { // l < r
                int tm = (tl + tr) / 2;
                if (r <= tm) {
                    return get(v * 2, tl, tm, l, r);
                }
                else if (l > tm) {
                    return get(v * 2 + 1, tm + 1, tr, l, r);
                }
                T cntL = get(2 * v, tl, tm, l, tm);
                T cntR = get(2 * v + 1, tm + 1, tr, tm + 1, r);

                return modify(calc(cntL, cntR), tree[v].add);
            }
        }

    public:

        segTreeM() {}
        segTreeM(std::vector<T>& A) {
            length = A.size();
            tree.resize(length * 4);
            build(1, 0, length - 1, A);
        }

        // get on a segment [l, r]
        T get(int l, int r) {
            return get(1, 0, length - 1, l, r);
        }

        // обновление на отрезке l...r значением value
        void update(int l, int r, T val) {
            update(1, 0, length - 1, l, r, val);
        }
    };

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

    // Implicit cartesian tree
    // O(1): empty, size
    // O(log n): insert, erase, [index], back, front, push_back
    // O(n): clear, get
    // O(n * log n): copy constructor, fill constructor, sort
    template<typename T>
    class container {

        struct node;
        typedef node* pnode;

        struct node {
            T value;
            pnode left = 0, right = 0;
            int size = 0;
            long long prior;

            node() {
                left = right = 0;
            }
            node(const T& value) {
                this->value = value;

                static std::mt19937_64 random(42);
                prior = random();
            }
        };

        pnode root = 0; // корень


        // clear branch
        void clearBranch(pnode& branch) {
            if (branch->left) { // если можно пойти в левый узел
                clearBranch(branch->left);
                delete branch->left;
            }
            if (branch->right) { // если можно пойти в правый узел
                clearBranch(branch->right);
                delete branch->right;
            }
        }

        // copy structor. Memory is clear
        void copy_memoryisClear(const container& source) {
            for (int i = 0; i < source.size(); i++) {
                push_back(source[i]);
            }
        }

        void move(container& source) {
            root = source.root;
            source.root = 0;
        }

    public:

        // default constructor
        container() {}
        // copy constructor. O(n * log n)
        container(const container& source) {
            copy_memoryisClear(source);
        }
        // move constuctor
        container(container&& source) noexcept {
            move(source);
        }
        // destructor. Clear all elements. O(n)
        ~container() {
            clear();
        }

        container& operator = (const container& source) {
            if (root != source.root) {
                clear();
                copy_memoryisClear(source);
            }
            return *this;
        }
        container& operator = (container&& source) noexcept {
            if (root != source.root) {
                clear();
                move(source);
            }
            return *this;
        }

        // fill constructor. O(n * log n)
        container(int length, const T& value) {
            while (length--) {
                push_back(value);
            }
        }

    private:
        template<typename array_t>
        void copyArray(const array_t& source) {
            for (auto it = source.begin(); it != source.end(); it++) {
                push_back(static_cast<T>(*it));
            }
        }
    public:

        container(std::initializer_list<T> list) {
            copyArray(list);
        }
        template<typename array_t>
        container(const array_t& A) {
            copyArray(A);
        }


        // return size in container. O(1)
        int size() const {
            return size(root);
        }
        // O(1)
        bool empty() const {
            return root == 0;
        }


        // Clear all elements. O(n)
        void clear() {
            if (root) {
                clearBranch(root);
                delete root;
                root = 0;
            }
        }

    private:
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

        // find [index]. O(h)
        const pnode& find(const pnode& branch, int index) const {
            if (index < size(branch->left)) {
                return find(branch->left, index);
            }
            else {
                index -= size(branch->left);
                return index == 0 ? branch : find(branch->right, index - 1);
            }
        }

    public:

        // insert value in [index]. O(h)
        void insert(const T& value, int index) {
            pnode left, right;
            split(root, left, right, index);
            merge(left, left, new node(value));
            merge(root, left, right);
        }

        // erase [index]. O(h)
        void erase(int index) {
            pnode left, right, mid;
            split(root, left, right, index);
            split(right, mid, right, 1);
            delete mid;
            merge(root, left, right);
        }


        const T& operator [](int index) const {
            return find(root, index)->value;
        }
        T& operator [](int index) {
            return find(root, index)->value;
        }
    };

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

    // куча с минимумом
    template<typename T>
    class pairingHeap {

        struct node;
        typedef node* pnode;

        // узел
        struct node {
            T key;
            pnode leftChild, nextSibling;

            node() {
                leftChild = nextSibling = 0;
            }
            node(const T& key, pnode leftChild, pnode nextSibling) {
                this->key = key;
                this->leftChild = leftChild;
                this->nextSibling = nextSibling;
            }
            node(T&& key, pnode leftChild, pnode nextSibling) {
                this->key = std::move(key);
                this->leftChild = leftChild;
                this->nextSibling = nextSibling;
            }

            // O(1)
            void addChild(pnode branch) {
                if (leftChild == 0) {
                    leftChild = branch;
                }
                else {
                    branch->nextSibling = leftChild;
                    leftChild = branch;
                }
            }
        };

        // слияние двух узлов. O(1)
        pnode merge(pnode A, pnode B) {
            if (A == 0) {
                return B;
            }
            else if (B == 0) {
                return A;
            }
            else if (A->key < B->key) {
                A->addChild(B);
                return A;
            }
            else {
                B->addChild(A);
                return B;
            }
        }

        // O(log n)
        pnode TwoPassMerge(pnode branch) {
            dst::edeque<pnode> stack;
            pnode B, newNode;
            while (!(branch == 0 || branch->nextSibling == 0)) {
                newNode = (B = branch->nextSibling)->nextSibling;

                branch->nextSibling = B->nextSibling = 0;
                stack.push_back(merge(branch, B));
                branch = newNode;
            }

            while (!stack.empty()) {
                branch = merge(stack.back(), branch);
                stack.pop_back();
            }
            return branch;
        }

        // O(log n)
        pnode TwoPassMerge_ERROR_STACKOVERFLOW(pnode branch) {
            if (branch == 0 || branch->nextSibling == 0) {
                return branch;
            }
            else {
                pnode B = branch->nextSibling;
                pnode newNode = B->nextSibling;

                branch->nextSibling = B->nextSibling = 0;

                return merge(merge(branch, B), TwoPassMerge(newNode));
            }
        }

        void clear(pnode branch) {
            if (branch->leftChild != 0) {
                clear(branch->leftChild);
                delete branch->leftChild;
            }
            if (branch->nextSibling) {
                clear(branch->nextSibling);
                delete branch->nextSibling;
            }
        }

        void copy(pnode branch, const pnode sourceBranch) {
            branch->key = sourceBranch->key;

            if (sourceBranch->leftChild != 0) {
                copy(branch->leftChild = new node(), sourceBranch->leftChild);
            }
            if (sourceBranch->nextSibling != 0) {
                copy(branch->nextSibling = new node(), sourceBranch->nextSibling);
            }
        }

        pnode root = 0; // корень

        void copy(const pairingHeap& source) {
            root = new node();
            copy(root, source.root);
        }

        void move(pairingHeap& source) {
            root = source.root;
            source.root = 0;
        }

    public:

        // default constructor
        pairingHeap() {}
        // copy constructor
        pairingHeap(const pairingHeap& source) {
            copy(source);
        }
        // move constructor
        pairingHeap(pairingHeap&& source) {
            move(source);
        }
        // destructor
        ~pairingHeap() {
            clear();
        }

        pairingHeap& operator = (const pairingHeap& source) {
            if (root != source.root) {
                clear();
                copy(source);
            }
            return *this;
        }
        pairingHeap& operator = (pairingHeap&& source) {
            if (root != source.root) {
                clear();
                move(source);
            }
            return *this;
        }

        void clear() {
            if (root) {
                clear(root);
                root = 0;
            }
        }

        bool empty() const {
            return root == 0;
        }

        // top. O(1)
        const T& top() const {
            return root->key;
        }

        // push key. O(1)
        void push(const T& key) {
            T addValue = key;
            push(std::move(addValue));
        }
        // push key. O(1)
        void push(T&& key) {
            root = merge(root, new node(std::move(key), 0, 0));
        }

        // амортизированно O(log n)
        void pop() {
            root = TwoPassMerge(root->leftChild);
        }

        // слияние двух очередей. O(1)
        void join(pairingHeap& other) {
            root = merge(root, other.root);
            other.root = 0;
        }
    };
}

// Algorithm: Numbers, Graphs, Strings, Sort
namespace alg {

    // Numbers: sqrt, namespace prime, factorize, epow, gcd, gcd(xy), lcm, fibonacci
    namespace nmb {
        // Метод Ньютона для поиска целочисленных корней
        template<typename T>
        T sqrt(const T& n) {
            T x = 1;
            bool decreased = false;
            while (true) {
                T nx = (x + n / x) / 2;
                if (x == nx || nx > x && decreased) {
                    break;
                }
                decreased = nx < x;
                x = nx;
            }
            return x;
        }

        // isPrime, lineSieve, sieve_info, lineSieveInfo, blockSieve, getPrimes
        namespace prime {
            // return true if n is prime. O(sqrt n)
            template<typename T>
            bool isPrime(const T& n) {
                T i, sqrtN = std::sqrt(n);
                // пока идем до корня и i не является делителем n
                for (i = 2; i <= sqrtN && n % i != 0; i++) {}
                return i > sqrtN;
            }

            // O(n)
            std::vector<int> lineSieve(long long n) {
                std::vector<int> primes, d(n + 1, -1);

                for (long long y = 2; y <= n; y++) {
                    if (d[y] == -1) { // новое простое
                        d[y] = primes.size();
                        primes.push_back(y);
                    }
                    for (long long i = 0; i <= d[y] && y * primes[i] <= n; i++) {
                        d[y * primes[i]] = i;
                    }
                }

                return primes;
            }

            // дополнительная информация для линейного решета Эратосфена
            struct sieve_info {
                std::vector<int>    primes; // primes[i] - простой делитель
                std::vector<int>         d; // d[x] - индекс минимального простого делителя числа x в массиве primes
                std::vector<int>       deg; // deg[x] - степень вхождения минимального простого делителя числа x
                std::vector<int>      rest; // rest[x] - x / p1^a1 - остаток числа x от p1^a1
                std::vector<int>      term; // term[x] - p1^a1
                std::vector<int>       phi; // количество взаимно простых с n
                std::vector<long long> div; // div[x] - количество делителей x
                std::vector<long long>  s1; // s1[x] - сумма делителей x

                sieve_info() {}
                sieve_info(long long n) {
                    d.resize(n + 1, -1);

                    deg.resize(n + 1, 0);
                    rest.resize(n + 1, 0);
                    term.resize(n + 1, 0);

                    phi.resize(n + 1, 1);
                    div.resize(n + 1, 1);

                    s1.resize(n + 1, 0);
                    s1[1] = 1;
                }

                // обновляет информацию для x
                void update(long long x) {
                    long long p_x = primes[d[x]];
                    long long y_x = x / p_x;

                    deg[x] = d[y_x] == d[x] ? deg[y_x] + 1 : 1;

                    rest[x] = d[y_x] == d[x] ? rest[y_x] : y_x;

                    term[x] = x / rest[x];

                    phi[x] = phi[rest[x]] * (term[x] / p_x) * (p_x - 1);

                    div[x] = div[rest[x]] * (static_cast<long long>(deg[x]) + 1);

                    s1[x] = s1[rest[x]] * ((p_x * term[x] - 1) / (p_x - 1));
                }
            };

            // O(n) + доп информация
            sieve_info lineSieveInfo(long long n) {
                sieve_info info(n);

                for (long long y = 2; y <= n; y++) {
                    if (info.d[y] == -1) { // новое простое
                        info.d[y] = info.primes.size();
                        info.primes.push_back(y);
                    }
                    for (long long i = 0; i <= info.d[y] && y * info.primes[i] <= n; i++) {
                        info.d[y * info.primes[i]] = i;
                    }

                    info.update(y); // обновление доп информации
                }

                return info;
            }

            // O(n * log log n). память sqrt(n). 
            // Находит все простые на отрезке [left, right]
            std::vector<int> blockSieve(long long left, long long right) {
                long long sqrtN = sqrt(right);
                long long s = max(sqrtN, static_cast<long long>(30));
                std::vector<int> primes;

                // просеивание до корня
                {
                    std::vector<bool> isPrime(sqrtN + 1, false);

                    for (long long i = 2; i <= sqrtN; i++) {
                        if (!isPrime[i]) { // простое
                            primes.push_back(i);
                            for (long long j = i * i; j <= sqrtN; j += i) {
                                isPrime[j] = true;
                            }
                        }
                    }
                }

                std::vector<int> result;

                for (long long k = left / s, maxk = right / s; k <= maxk; k++) {
                    std::vector<bool> isPrime(s + 1, false);

                    long long begin = k * s;
                    for (long long i = 0; i < primes.size(); ++i) {
                        long long beginIndex = (begin + primes[i] - 1) / primes[i];
                        long long j = max(beginIndex, static_cast<long long>(2)) * primes[i] - begin;
                        for (; j < s; j += primes[i]) {
                            isPrime[j] = true;
                        }
                    }
                    if (k == 0) { // ответ должен быть без 0 и 1
                        isPrime[0] = isPrime[1] = true;
                    }
                    for (long long i = 0; i < s && begin + i <= right; ++i) {
                        if (left <= begin + i && !isPrime[i]) { // простое
                            result.push_back(begin + i);
                        }
                    }
                }
                return result;
            }

            // O(n * log log n). primes [n^2, n^2 + n]
            std::vector<int> getPrimes(long long n) {
                std::vector<bool> isPrime(n + 1, true); // 0
                std::vector<bool> NisPrime(n + 1, true); // n * n 

                for (long long x = 2; x <= n; x++) {
                    if (isPrime[x]) {
                        for (long long y = x * x; y <= n; y += x) {
                            isPrime[y] = false;
                        }

                        long long y = (n * n / x) * x;
                        while (y < n * n) {
                            y += x;
                        }

                        while (y <= n * n + n) {
                            NisPrime[y - n * n] = false;
                            y += x;
                        }
                    }
                }

                std::vector<int> primes; // result
                for (long long i = 0; i <= n; i++) {
                    if (NisPrime[i]) {
                        primes.push_back(i + n * n);
                    }
                }
                return primes;
            }
        }

        // factorizes the number
        // vector<p^a>
        template<typename T>
        std::vector<std::pair<T, int>> factorize(T N) {
            std::vector<std::pair<T, int>> result;
            for (T i = 2; i * i <= N; i++) {
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
        template<typename value_type, typename power_type>
        value_type epow(const value_type& a, const power_type& n, const value_type& mod) {
            if (n == 0) {
                return 1;
            }
            else {
                value_type z = epow(a, n >> 1, mod);
                z = (z * z) % mod;

                return (n & 1) ? (z * a) % mod : z;
            }
        }

        // binary pow a^n
        template<typename value_type, typename power_type>
        value_type epow(const value_type& a, const power_type& n) {
            if (n == 0) {
                return 1;
            }
            else {
                value_type z = epow(a, n >> 1);
                z *= z;

                return (n & 1) ? z * a : z;
            }
        }

        // Алгоритм Евклида
        template<typename T>
        T gcd(T a, T b) {
            // a <= b ever
            if (a > b) {
                std::swap(a, b);
            }
            while (a != 0) {
                b %= a;
                // a > b
                std::swap(a, b);
                // a < b
            }
            return b;
        }

        // Расширенный Алгоритм Евклида
        template<typename T>
        T gcd(T a, T b, T& x, T& y) {
            if (a == 0) {
                x = 0;
                y = 1;
                return b;
            }
            else {
                T x1, y1;
                T res = gcd(b % a, a, x1, y1);
                x = y1 - (b / a) * x1;
                y = x1;
                return res;
            }
        }


        // Линейные диофантовы уравнения с двумя переменными
        // a * x + b * y = c
        // a * xg + b * yg = g
        template<typename T>
        bool solveLDE(T a, T b, T c, T& xg, T& yg, T& g) {
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
        template<typename T>
        bool squarePres(const T& N, T& a, T& b) {
            for (a = 2; a * a <= N; a++) {
                b = std::sqrt(N - a * a);
                if (a * a + b * b == N) {
                    return true;
                }
            }
            return false;
        }

        // Обратный элемент в кольце по модулю
        template<typename T>
        bool reverseElement(T a, T m, T& x) {
            T y;
            T g = gcd(a, m, x, y);
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
            dst::matrix<T> binPow(const dst::matrix<T>& a, const T& n) {
                if (n == 1) {
                    return a;
                }
                else {
                    auto z = binPow(a, n / 2);
                    z *= z;

                    return n % 2 == 0 ? z : z * a;
                }
            }
            dst::matrix<T> binPow(const dst::matrix<T>& a, const T& n, const T& mod) {
                if (n == 1) {
                    return a;
                }
                else {
                    auto z = binPow(a, n / 2);
                    z = (z * z) % mod;

                    return n % 2 == 0 ? z : (z * a) % mod;
                }
            }

            dst::matrix<T> getMatrix() const {
                return { {0, 1}, {1, 1} };
            }

        public:
            T get(const T& n) {
                return binPow(getMatrix(), n)[0][0];
            }
            T get(const T& n, const T& mod) {
                return binPow(getMatrix(), n, mod)[0][0];
            }
        };
    }

    // Graphs: dijkstra, mst, graphSets, graphCycles
    namespace gph {

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

// Long Arithmetic: elong
namespace emt {

    static const long long long_length = 18,
        long_base = alg::nmb::epow<long long>(10, long_length),
        long_base_expansion = std::sqrt(long_base);

    // signed long integer
    class elong {
        // unsigned basic long
        struct basicLong {
            dst::edeque<s64> digits; // цифры

            basicLong() {}
            basicLong(const dst::edeque<s64>& digits) {
                this->digits = digits;
            }
            basicLong(const std::string& str) {
                int i;
                for (i = str.size(); i >= long_length; i -= long_length) {
                    digits.push_back(atoll(str.substr(i - long_length, long_length).c_str()));
                }
                if (i > 0) {
                    digits.push_back(atoll(str.substr(0, i).c_str()));
                }
            }

            // удаляет ведущие нули
            basicLong& remove_leading_zeros() {
                while (!digits.empty() && digits.back() == 0) {
                    digits.pop_back();
                }
                return *this;
            }

        private:
            enum class OP {
                less,
                equally,
                more,
            };

#define compareValues(a, b) (a < b ? OP::less : OP::more)

            // сравнивает 2 числа
            OP compare(const basicLong& Rhs) const {
                // если у них разные длины
                if (digits.size() != Rhs.digits.size()) {
                    return compareValues(digits.size(), Rhs.digits.size());
                }
                else {
                    int i = digits.size() - 1;
                    while (i >= 0 && digits[i] == Rhs.digits[i]) {
                        i--;
                    }
                    return i >= 0 ? compareValues(digits[i], Rhs.digits[i]) : OP::equally;
                }
            }

        public:

            bool operator < (const basicLong& Rhs) const {
                return compare(Rhs) == OP::less;
            }
            bool operator == (const basicLong& Rhs) const {
                return compare(Rhs) == OP::equally;
            }
            bool operator > (const basicLong& Rhs) const {
                return compare(Rhs) == OP::more;
            }

            // операция сложения
            basicLong operator + (const basicLong& added) const {
                basicLong result = *this;
                bool k = 0;
                s64 i = 0;
                s64 len = max(result.digits.size(), added.digits.size());
                for (i = 0; i < len || k != 0; i++) {
                    if (i == result.digits.size()) {
                        result.digits.push_back(0);
                    }
                    result.digits[i] += k + (i < added.digits.size() ? added.digits[i] : 0);
                    k = result.digits[i] >= long_base;
                    result.digits[i] -= k * long_base;
                }
                return result;
            }

            // операция вычитания
            basicLong operator - (const basicLong& subtrahend) const {
                basicLong result = *this;
                bool k = 0;
                for (int i = 0; i < subtrahend.digits.size() || k != 0; i++) {
                    result.digits[i] -= k + (i < subtrahend.digits.size() ? subtrahend.digits[i] : 0);
                    k = result.digits[i] < 0;
                    result.digits[i] += k * long_base;
                }
                return result.remove_leading_zeros();
            }

            // расширяет длину числа, уменьшая его модуль
            basicLong expansion() const {
                basicLong result;
                result.digits.resize(digits.size() << 1);
                for (int i = 0; i < digits.size(); i++) {
                    result.digits[i << 1] = digits[i] % long_base_expansion;
                    result.digits[(i << 1) + 1] = digits[i] / long_base_expansion;
                }
                return result;
            }

            // сужает длину числа, увеличивая его модуль
            basicLong reduction() const {
                basicLong result;
                result.digits.resize(digits.size() >> 1);
                for (int i = 0; i < result.digits.size(); i++) {
                    result.digits[i] = digits[i << 1] + digits[(i << 1) + 1] * long_base_expansion;
                }
                return result;
            }

            // Простое умножение двух чисел за O(n^2)
            basicLong naiveMul(const basicLong& a, const basicLong& b) const {
                basicLong result, mult1 = a.expansion(), mult2 = b.expansion();
                result.digits.resize(mult1.digits.size() + mult2.digits.size() + 1);

                u64 c, k;
                int i, j;
                for (i = 0; i < mult1.digits.size(); i++) {
                    c = 0;
                    for (j = 0; j < mult2.digits.size() || c != 0; j++) {
                        // mul1[i] < 1e9 && mult2[j] < 1e9
                        // mult1[i] * mult2[j] < 1e18
                        k = result.digits[i + j] + (j < mult2.digits.size() ? mult1.digits[i] * mult2.digits[j] : 0) + c;
                        c = k / long_base_expansion;
                        result.digits[i + j] = k - c * long_base_expansion;
                    }
                }
                return result.reduction().remove_leading_zeros();
            }

            // O(n ^ 1.54)
            basicLong KaratsubaMul(const basicLong& x, const basicLong& y, int n) const {
                if (n <= 128) {
                    return naiveMul(x, y);
                }
                else {
                    int half = n >> 1;
                    // firstHalf - первая половина числа
                    basicLong Xl = x.firstHalf(), Xr = x.secondHalf(), Yl = y.firstHalf(), Yr = y.secondHalf();
                    basicLong sumX = Xl + Xr, sumY = Yl + Yr;

                    basicLong P1 = KaratsubaMul(Xl, Yl, half);
                    basicLong P2 = KaratsubaMul(Xr, Yr, half);
                    basicLong P3;
                    if (sumX.digits.size() != half || sumY.digits.size() != half) {
                        sumX.digits.resize(n);
                        sumY.digits.resize(n);
                        P3 = KaratsubaMul(sumX, sumY, n);
                    }
                    else {
                        P3 = KaratsubaMul(sumX, sumY, half);
                    }
                    // P1 * (base ^ n) + (P3 - P1 - P2) * (base ^ half) + P2
                    return (P1 + ((P3 - P1 - P2) << half) + (P2 << n)).remove_leading_zeros();
                }
            }

            basicLong operator * (const basicLong& mult) const {
                int k = max(this->digits.size(), mult.digits.size());
                if (k <= 128) {
                    return naiveMul(*this, mult);
                }
                else {
                    basicLong mult1 = *this, mult2 = mult;
                    int len = roundTwo(k);
                    mult1.digits.resize(len);
                    mult2.digits.resize(len);
                    return KaratsubaMul(mult1, mult2, len);
                }
            }

            // ( /, % )
            std::pair<basicLong, basicLong> div(const basicLong& divider) const {
                basicLong result, value, temp, t;

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
                        long long m = (l + r) >> 1;

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
                return std::make_pair(result, value);
            }

            basicLong operator / (const basicLong& divider) const {
                // если divider == 2
                if (divider.digits.size() == 1 && divider.digits.front() == 2) {
                    basicLong res = *this;
                    res.digits[0] >>= 1;
                    for (int i = 1; i < res.digits.size(); i++) {
                        if (res.digits[i] & 1) { // res.digits[i] % 2 
                            res.digits[i - 1] += long_base >> 1;
                        }
                        res.digits[i] >>= 1;
                    }
                    return res.remove_leading_zeros();
                }
                else {
                    return div(divider).first;
                }
            }

            basicLong operator % (const basicLong& divider) const {
                if (divider.digits.size() == 1 && divider.digits.front() == 2) {
                    basicLong res;
                    res.digits.push_back(digits.front() & 1);
                    return res;
                }
                else {
                    return div(divider).second;
                }
            }
            

            // сдвигает число на long_base^k влево
            basicLong operator << (int k) const {
                basicLong res = *this;
                while (k--) {
                    res.digits.push_front(0);
                }
                return res;
            }

            basicLong firstHalf() const {
                int half = digits.size() >> 1;
                basicLong res;
                res.digits.resize(half);
                for (int i = 0; i < half; i++) {
                    res.digits[i] = digits[i];
                }
                return res;
            }
            basicLong secondHalf() const {
                int half = digits.size() >> 1;
                basicLong res;
                res.digits.resize(half);
                for (int i = 0; i < half; i++) {
                    res.digits[i] = digits[i + half];
                }
                return res;
            }

            basicLong sqrt() const {
                basicLong left, right, two, one;
                right.digits.resize((digits.size() >> 1) + digits.size() & 1);
                for (int i = 0; i < right.digits.size(); i++) {
                    right.digits[i] = long_base - 1;
                }
                two.digits.push_back(2);
                one.digits.push_back(1);
                while (left + one < right) {
                    basicLong mid = (left + right) / two;
                    if (mid * mid > * this) {
                        right = mid;
                    }
                    else {
                        left = mid;
                    }
                }
                return left.remove_leading_zeros();
            }
        };

        basicLong value;
        bool isNegative;

        // update
        elong& update() {
            if (value.digits.empty()) {
                value.digits.push_back(0);
            }
            return *this;
        }

        elong(const basicLong& value, bool isNegative) {
            this->value = value;
            this->isNegative = isNegative;
        }
        elong(basicLong&& value, bool isNegative) {
            this->value = std::move(value);
            this->isNegative = isNegative;
        }

        void inc() {
            value.digits.front()++;
            if (value.digits.front() == long_base) {
                value.digits.front() = 0;

                int i;
                for (i = 1; i < value.digits.size() && value.digits[i] + 1 == long_base; i++) {
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
            value.digits.push_back(0);
            isNegative = false;
        }

        // convert constructors[

        elong(const std::string& String) {
            isNegative = String[0] == '-';
            value = isNegative ? basicLong(String.substr(1)) : basicLong(String);
        }
        elong(const char* String) {
            *this = elong(std::string(String));
        }

        template<typename T>
        elong(T number) {
            isNegative = number < 0;
            number *= isNegative ? -1 : 1;

            value.digits.push_back(number % long_base);
            number /= long_base;
            if (number > 0) {
                value.digits.push_back(number % long_base);
                number /= long_base;
                if (number > 0) {
                    value.digits.push_back(number);
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
        elong operator ++(int) {
            auto temp = *this;
            signedInc();
            return temp;
        }
        elong& operator --() {
            signedDec();
            return *this;
        }
        elong operator --(int) {
            auto temp = *this;
            signedDec();
            return temp;
        }

        elong operator + (const elong& added) const {
            if (isNegative) {
                return added.isNegative ? -(-*this + (-added)) : added - (-*this);
            }
            else {
                return added.isNegative ? *this - (-added) : elong(std::move(value + added.value), false);
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
                return *this < subtrahend ? -(subtrahend - *this) : elong(std::move(value - subtrahend.value), false).update();
            }
        }
        elong& operator -= (const elong& subtrahend) {
            return *this = *this - subtrahend;
        }

        elong operator * (const elong& multiplied) const {
            return elong(std::move(value * multiplied.value), isNegative != multiplied.isNegative).update();
        }
        elong& operator *= (const elong& multiplied) {
            return *this = *this * multiplied;
        }

        elong operator / (const elong& divider) const {
            return elong(std::move(value / divider.value), isNegative != divider.isNegative).update();
        }
        elong& operator /= (const elong& divider) {
            return *this = *this / divider;
        }

        elong operator % (const elong& divider) const {
            return elong(std::move(value % divider.value), isNegative).update();
        }
        elong& operator %= (const elong& divider) {
            return *this = *this % divider;
        }

        // ( /, % )
        std::pair<elong, elong> div(const elong& divider) const{
            auto p = std::move(value.div(divider.value));
            return std::make_pair(elong(std::move(p.first), isNegative != divider.isNegative).update(), elong(std::move(p.second), isNegative).update());
        }
        // ]

        // value * 2^k
        elong operator << (const size_t& k) const {
            return *this * alg::nmb::epow(elong(2), k);
        }
        elong& operator <<= (const size_t& k) {
            return *this = (*this << k);
        }

        // value / 2^k
        elong operator >> (const size_t& k) const {
            return *this / alg::nmb::epow(elong(2), k);
        }
        elong& operator >>= (const size_t& k) {
            return *this = (*this >> k);
        }

        elong sqrt() const {
            return elong(std::move(value.sqrt()), false).update();
        }

        // friend function :)
        friend std::ostream& operator << (std::ostream& out, const elong& var);

        // long > 0
        explicit operator bool() const {
            return !isNegative && value.digits[0] != 0;
        }
        explicit operator long long() const {
            return value.digits.back() * (isNegative ? -1 : 1);
        }
        explicit operator unsigned long long() const {
            if (value.digits.size() == 1) {
                return value.digits.back();
            }
            else {
                return static_cast<unsigned long long>(value.digits[1]) * long_base + value.digits[0];
            }
        }
    };

    template<typename T>
    elong operator + (const T& a, const elong& b) {
        return elong(a) + b;
    }
    template<typename T>
    elong operator - (const T& a, const elong& b) {
        return elong(a) - b;
    }
    template<typename T>
    elong operator * (const T& a, const elong& b) {
        return elong(a) * b;
    }
    template<typename T>
    elong operator / (const T& a, const elong& b) {
        return elong(a) / b;
    }
    template<typename T>
    elong operator % (const T& a, const elong& b) {
        return elong(a) % b;
    }

    // output signed long integer
    std::ostream& operator << (std::ostream& out, const elong& var) {
        if (var.isNegative) {
            out << '-';
        }
        out << var.value.digits.back();
        char k = out.fill('0');
        for (int i = var.value.digits.size() - 2; i >= 0; i--) {
            out << std::setw(long_length) << var.value.digits[i];
        }
        out.fill(k);
        return out;
    }
    // input signed long integer
    std::istream& operator >> (std::istream& input, elong& var) {
        std::string str;
        input >> str;
        var = elong(str);
        return input;
    }

    elong tree_factor(long long l, long long r) {
        if (l < r - 1) {
            long long m = (l + r) >> 1;
            return tree_factor(l, m) * tree_factor(m + 1, r);
        }
        else if (r - l == 1) {
            return l * r;
        }
        else if (l == r) {
            return l;
        }
        else {
            return 1;
        }
    }

    elong factorial(long long n) {
        return tree_factor(2, n);
    }
}

// Variables: u128, modVar
namespace var {

    const static emt::elong mod64 = "18446744073709551616";

    // unsigned int128. alpha version
    class u128 {

        // value = left + right * 2^64
        u64 left, right;

        void reduction(u64 A[4]) const {
            A[1] = left >> 32; 
            A[3] = right >> 32; 
            A[0] = left - (A[1] << 32); 
            A[2] = right - (A[3] << 32);
        }
        // naive multip
        u128 naiveMul(const u128& x, const u128& y) const {
            u64 result[] = { 0, 0, 0, 0 };
            static u64 mult1[4], mult2[4];

            x.reduction(mult1);
            y.reduction(mult2);

            u64 c, k;
            int i, j;
            // 10 умножений
            for (i = 0; i < 4; i++) {
                c = 0;
                for (j = 0; i + j < 4; j++) {
                    k = result[i + j] + mult1[i] * mult2[j] + c;
                    c = k >> 32;
                    result[i + j] = k - (c << 32);
                }
            }
            return u128(result[0] + (result[1] << 32), result[2] + (result[3] << 32));
        }

    public:

        u128(u64 l, u64 r) {
            left = l;
            right = r;
        } 
        u128() {
            left = right = 0;
        }
        // Подходят все встроенные типы чисел
        template<typename T>
        u128(const T& value) {
            left = value;
            right = value < 0 ? -1 : 0;
        }
        u128(const emt::elong& value) {
            auto var = value.div(mod64);
            left = var.second.operator size_t();
            right = var.first.operator size_t();
        }

        // возвращает (1 << bits), где bits принадлежит [0, 127]
        u128 getDegree2(size_t bits) const {
            return (bits < 64 ? u128(static_cast<u64>(1) << bits, 0) : u128(0, static_cast<u64>(1) << (bits - 64)));
        }


        bool operator == (const u128& Rhs) const {
            return left == Rhs.left && right == Rhs.right;
        }
        bool operator != (const u128& Rhs) const {
            return !(*this == Rhs);
        }
        bool operator <  (const u128& Rhs) const {
            return right == Rhs.right ? left < Rhs.left : right < Rhs.right;
        }
        bool operator >  (const u128& Rhs) const {
            return right == Rhs.right ? left > Rhs.left : right > Rhs.right;
        }

        bool operator <= (const u128& Rhs) const {
            return !(*this > Rhs);
        }
        bool operator >= (const u128& Rhs) const {
            return !(*this < Rhs);
        }


        u128& operator ++ () {
            left++;
            if (left == 0) {
                right++;
            }
            return *this;
        }
        u128  operator ++ (int) {
            u128 temp = *this;
            ++(*this);
            return temp;
        }

        u128& operator -- () {
            if (left == 0) {
                right--;
            }
            left--;
            return *this;
        }
        u128  operator -- (int) {
            u128 temp = *this;
            --(*this);
            return temp;
        }


        u128 operator + (const u128& add) const {
            u128 result(left + add.left, right + add.right);
            if (left < max(left, add.left)) {
                result.right++;
            }
            return result;
        }
        u128 operator - (const u128& sub) const {
            u128 result(left - sub.left, right - sub.right);
            if (left < result.left) {
                result.right--;
            }
            return result;
        }
        u128 operator * (const u128& mult) const {
            return naiveMul(*this, mult);
        }

        // 0.003ms - средний случай 
        // 0.1ms - худший случай
        // return (a / b, a % b)
        std::pair<u128, u128> division(u128 divValue, u128 divider) const {
            u128 result, temp;
            std::vector<u128> Divs; // max 128 el

            if (divValue >= divider) {
                // max 128 ветков
                while (divValue >= divider) {
                    divValue -= divider;

                    Divs.push_back(divider);

                    temp = divider + divider;
                    // если overflow
                    if (temp <= divider) {
                        break;
                    }
                    else {
                        divider = temp;
                    }
                }
                result += getDegree2(Divs.size()) - 1;
            }
            else {
                result = 0;
            }
            
            // раскрутка
            for (int i = Divs.size() - 1; i >= 0; i--) {
                if (divValue >= Divs[i]) {
                    divValue -= Divs[i];
                    result += getDegree2(i);
                }
            }
            return std::make_pair(result, divValue);
        }
        u128 operator / (const u128& div) const {
            return division(*this, div).first;
        }
        u128 operator % (const u128& div) const {
            return division(*this, div).second;
        }
        

        u128& operator += (const u128& add) {
            return *this = *this + add;
        }
        u128& operator -= (const u128& sub) {
            return *this = *this - sub;
        }
        u128& operator *= (const u128& mult) {
            return *this = *this * mult;
        }
        u128& operator /= (const u128& div) {
            return *this = *this / div;
        }
        u128& operator %= (const u128& div) {
            return *this = *this % div;
        }


        u128 operator & (const u128& k) const {
            return u128(left & k.left, right & k.right);
        }
        u128 operator | (const u128& k) const {
            return u128(left | k.left, right | k.right);
        }
        u128 operator ^ (const u128& k) const {
            return u128(left ^ k.left, right ^ k.right);
        }
        u128 operator ~ () const {
            return u128(~left, ~right);
        }


        u128& operator &= (const u128& k) {
            return *this = *this & k;
        }
        u128& operator |= (const u128& k) {
            return *this = *this | k;
        }
        u128& operator ^= (const u128& k) {
            return *this = *this ^ k;
        }

        

        // bits = [0, 127]
        // если биты вылезут за пределы, то они потеряются

        u128& operator <<= (u64 bits) {
            // bits / 64 == 1
            if ((bits >> 6) == 1) {
                // сдвинуть на 64 битов
                right = left;
                left = 0;
            }
            // (bits %= 64) != 0
            if ((bits = twoRemainder(bits, 6)) != 0) {
                right = (right << bits) | (left >> (64 - bits));
                left <<= bits;
            }
            return *this;
        }
        u128& operator >>= (u64 bits) {
            // bits / 64 == 1
            if (bits >> 6 == 1) {
                left = right;
                right = 0;
            }

            // (bits %= 64) != 0
            if ((bits = twoRemainder(bits, 6)) != 0) {
                left = (left >> bits) | (right << (64 - bits));
                right >>= bits;
            }
            return *this;
        }
        
        u128 operator << (u64 bits) const {
            return u128(*this) <<= bits;
        }
        u128 operator >> (u64 bits) const {
            return u128(*this) >>= bits;
        }

        u64 getLeft() const {
            return left;
        }
        u64 getRight() const {
            return right;
        }

        // меняет местами left и right
        u128 reverse() const {
            return u128(right, left);
        }

        template<typename T>
        explicit operator T() const {
            return (static_cast<T>(right) << 64) + left;
        }
        explicit operator emt::elong() const {
            return (static_cast<emt::elong>(right) * mod64) + left;
        }

        friend std::ostream& operator << (std::ostream& output, const u128& value);
    };
    // 0.4ms
    std::istream& operator >> (std::istream& input, u128& value) {
        emt::elong Long;
        input >> Long;
        value = u128(Long);
        return input;
    }
    // 0.5 - 1.5ms
    std::ostream& operator << (std::ostream& output, const u128& value) {
        if (value.right == 0) {
            output << value.left;
        }
        else {
            output << (value.operator emt::elong());
        }
        return output;
    }


    // value e [0, mod - 1]
    template<s64 mod>
    class modVar {

        s64 inverse(s64 n) const {
            s64 x, y;
            s64 d = gcd(n, mod, x, y);
            s64 inverse = (x % mod + mod) % mod;
            return inverse;
            //return epow(n, mod - 2);
        }

    public:

        s64 value;
        modVar() {
            value = 0;
        }
        template<typename T>
        modVar(const T& a) {
            value = static_cast<T>(a);
            value %= mod;
            if (value < 0) {
                value += mod;
            }
        }

        modVar operator + (const modVar& add) const {
            s64 res = value + add.value;
            if (res >= mod) {
                res -= mod;
            }
            return res;
        }
        modVar operator - (const modVar& sub) const {
            s64 res = value - sub.value;
            if (res < 0) {
                res += mod;
            }
            return res;
        }
        modVar operator * (const modVar& mult) const {
            return (value * mult.value) % mod;
        }
        modVar operator / (const modVar& div) const {
            return (value * inverse(div.value)) % mod;
        }

        modVar& operator += (const modVar& add) {
            return *this = *this + add;
        }
        modVar& operator -= (const modVar& sub) {
            return *this = *this - div;
        }
        modVar& operator *= (const modVar& mult) {
            return *this = *this * mult;
        }
        modVar& operator /= (const modVar& div) {
            return *this = *this / div;
        }

        bool operator == (const modVar& Rhs) const {
            return value == Rhs.value;
        }
        bool operator != (const modVar& Rhs) const {
            return !(*this == Rhs);
        }

        bool operator < (const modVar& Rhs) const {
            return value < Rhs.value;
        }
        bool operator > (const modVar& Rhs) const {
            return value > Rhs.value;
        }

        bool operator <= (const modVar& Rhs) const {
            return !(*this > Rhs.value);
        }
        bool operator >= (const modVar& Rhs) const {
            return !(*this < Rhs.value);
        }
    };
}

// Memory Manager: poolAllocator
namespace mem {
    // memory pool
    template<size_t bitLength>
    class poolAllocator {
        //memory pool
        char A[bitLength];

        struct freeSegment {
            void* begin; // начало 
            size_t length; // длина
            freeSegment(void* begin, size_t length) {
                this->begin = begin;
                this->length = length;
            }
            bool operator < (const freeSegment& Rhs) const {
                //      одинаковые длины     
                return length == Rhs.length ? begin < Rhs.begin : length < Rhs.length;
            }
        };

        struct segment {
            void* begin; // начало 
            mutable size_t length; // длина
            mutable bool free; // свободность
            segment(void* begin, size_t length, bool free) {
                this->begin = begin;
                this->length = length;
                this->free = free;
            }
            bool operator < (const segment& Rhs) const {
                return begin < Rhs.begin;
            }
        };

        std::set<freeSegment> freeSeg; // свободные отрезки памяти
        std::set<segment> Seg; // отрезки памяти
    public:

        poolAllocator() {
            freeSeg.insert(freeSegment(A, bitLength));
            Seg.insert(segment(A, bitLength, true));
        }

        // выделяет length элементов
        // O(log n). n - колво отрезков
        void* allocate(size_t length) {
            void* res;

            auto it = freeSeg.lower_bound(freeSegment(0, length));
            if (it != freeSeg.end()) { // нашли нужный отрезок
                segment newSegment(it->begin, length, false);// выделяем запрошенный отрезок
                size_t segLen = it->length; // длина прошлого отрезка

                // удаление отрезка
                freeSeg.erase(it);
                Seg.erase(newSegment);

                if (segLen != length) { // остался кусок отрезка
                    // закидываем остаточный свободный отрезок
                    freeSeg.insert(freeSegment(static_cast<char*>(newSegment.begin) + length, segLen - newSegment.length));
                    Seg.insert(segment(static_cast<char*>(newSegment.begin) + length, segLen - newSegment.length, true));
                }
                Seg.insert(newSegment);
                res = newSegment.begin;
            }
            else { // Pool overflow
                res = malloc(length);
            }
            return res;
        }

        // освобождает массив ptr
        // O(log n). n - колво отрезков
        void destroy(void* ptr) {
            // находим занятый отрезок с началом ptr 
            auto it = Seg.find(segment(ptr, 0, 0));
            it->free = true; // освобождаем отрезок

            // обьединяем свободных соседей

            auto right = it;
            right++;

            // есть правый сосед и он свободен
            if (right != Seg.end() && right->free) {
                it->length += right->length; // прибавляем отрезок
                // удаляем ненужный отрезок
                freeSeg.erase(freeSegment(right->begin, right->length));
                Seg.erase(right);
            }

            // есть левый сосед
            if (it != Seg.begin()) {
                auto left = it;
                left--;
                // он свободен
                if (left->free) {
                    // удаляем свободный отрезок
                    freeSeg.erase(freeSegment(left->begin, left->length));

                    // прибавляем отрезок
                    left->length += it->length;

                    // freeSeg(it) не нужно удалять, потому что его там нет
                    Seg.erase(it);
                    it = left;
                }
            }
            // добавляем свободный отрезок
            freeSeg.insert(freeSegment(it->begin, it->length));
        }

        // проверка принадлежности ссылки memory pool
        bool ison(void* ptr) const {
            return A <= ptr && ptr < A + bitLength;
        }
    };
};

//#define MEMORYPOOL
#ifdef MEMORYPOOL
mem::poolAllocator<1000000> memoryPool; // memory pool

void* operator new[](size_t length) {
    return length > 1e4 ? memoryPool.allocate(length) : malloc(length);
}

void operator delete[](void* ptr) {
    // если ссылка принадлежит memory pool
    if (memoryPool.ison(ptr)) {
        memoryPool.destroy(ptr);
    }
    else {
        free(ptr);
    }
}
#endif


#include<bits/stdc++.h>
using namespace std;
using namespace var;



int main() {
    ifstream cin("input.txt");

    dst::container<vector<int>> A;
    for (int i = 0; i < 1e6; i++) {
        A.insert({}, 0);
    }
    while (!A.empty()) {
        A.erase(0);
    }

    return 0;
}
