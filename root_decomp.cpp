// root_decomp
// Developed by Mob


#include<vector>
#include<iostream>
using namespace std;

template<typename T>
class root_decomp {

    int block_len; // длина блока

    vector<vector<T>> table; // таблица

    // строит table по вектору A
    // table должен быть пуст
    // O(n)
    void build(const vector<T>& A) {
        block_len = sqrt(A.size());

        // table resize sz(A) / block_len с округлением вверх
        table.resize((A.size() + block_len - 1) / block_len, vector<T>(block_len));

        for (int i = 0; i < table.size() - 1; i++) {
            for (int j = 0; j < block_len; j++) {
                table[i][j] = A[i * block_len + j];
            }
        }

        table.back().resize(A.size() - (table.size() - 1) * block_len);
        for (int j = 0; j < table.back().size(); j++) {
            table.back()[j] = A[(table.size() - 1) * block_len + j];
        }
    }

    // пробегает по table
    // возвращает номер блока
    // изменяет ind на номер числа в этом блоке
    // O(sqrt n)
    int run(int& ind) {
        int cpy_ind = ind;
        int block_i = 0;
        while (block_i < table.size() && ind >= table[block_i].size()) {
            ind -= table[block_i].size();
            block_i++;
        }
        return block_i;
    }


    // перестраивает структуру
    // O(n)
    void rebuild() {
        vector<T> A = get();
        table.clear();
        build(A);
    }


public:

    root_decomp() {
        block_len = 1;
    }

    // O(n)
    root_decomp(const vector<T>& A) {
        build(A);
    }

    // O(n)
    root_decomp(int len, const T& val) {
        block_len = sqrt(len);

        // table resize sz(A) / block_len с округлением вверх
        table.resize((len + block_len - 1) / block_len, vector<T>(block_len));

        for (int i = 0; i < table.size() - 1; i++) {
            for (int j = 0; j < block_len; j++) {
                table[i][j] = val;
            }
        }

        table.back().resize(len - (table.size() - 1) * block_len);
        for (int j = 0; j < table.back().size(); j++) {
            table.back()[j] = val;
        }
    }

    // O(sqrt n)
    int size() const {
        int sz = 0;
        for (int i = 0; i < table.size(); i++) {
            sz += table[i].size();
        }
        return sz;
    }
    // O(sqrt n)
    bool empty() const {
        return size() == 0;
    }


    // ind e [0, sz() - 1]
    // O(sqrt n)
    T& operator [](int ind) {
        int block_i = run(ind);
        return table[block_i][ind];
    }

    // ind e [0, sz() - 1]
    // O(sqrt n)
    const T& operator [](int ind) const {
        int block_i = run(ind);
        return table[block_i][ind];
    }


    // O(sqrt n)
    void insert(int ind, const T& val) {
        int block_i = run(ind);

        if (block_i == table.size()) {
            if (block_i > 0 && table[block_i - 1].size() < block_len * 2) {
                // можем положить в конец
                table[block_i - 1].push_back(val);
            }
            else {
                table.push_back({ val });

                if (table.size() > block_len * 2) { // колво блоков много
                    rebuild();
                }
            }
        }
        else {
            table[block_i].insert(table[block_i].begin() + ind, val);

            if (table[block_i].size() > block_len * 2) { 
                // разделить

                table.insert(table.begin() + block_i + 1, vector<T>(table[block_i].begin() + block_len, table[block_i].end()));

                table[block_i].resize(block_len);


                if (table.size() > block_len * 2) {
                    rebuild();
                }
            }
        }
    }

    // O(sqrt n)
    void erase(int ind) {
        int block_i = run(ind);

        table[block_i].erase(table[block_i].begin() + ind);
        if (table[block_i].empty()) {
            table.erase(table.begin() + block_i);

            if (table.size() < block_len / 2) {
                rebuild();
            }
        }
        else {
            if (table[block_i].size() < block_len / 2) {
                // попробовать объединить блоки

                // есть левый блок и объединение с ним не нарушает
                if (block_i > 0 && table[block_i - 1].size() + table[block_i].size() < 2 * block_len) {
                    table[block_i - 1].insert(table[block_i - 1].end(), table[block_i].begin(), table[block_i].end());
                    table.erase(table.begin() + block_i);
                }
                else if (block_i + 1 < table.size() && table[block_i].size() + table[block_i + 1].size() < 2 * block_len) {
                    table[block_i].insert(table[block_i].end(), table[block_i + 1].begin(), table[block_i + 1].end());
                    table.erase(table.begin() + block_i + 1);
                }
            }
        }
    }


    // возвращает вектор из элементов этой структуры
    // O(n)
    vector<T> get() {
        vector<T> A;
        A.reserve(size());
        for (int i = 0; i < table.size(); i++) {
            for (int j = 0; j < table[i].size(); j++) {
                A.push_back(table[i][j]);
            }
        }
        return A;
    }


    T& back() {
        return (*this)[size() - 1];
    }
    T& front() {
        return (*this)[0];
    }

    const T& back() const {
        return (*this)[size() - 1];
    }
    const T& front() const {
        return (*this)[0];
    }

    void pop_back() {
        erase(size() - 1);
    }
    void pop_front() {
        erase(0);
    }

    void push_back(const T& val) {
        insert(size(), val);
    }
    void push_front(const T& val) {
        insert(0, val);
    }
};
