// dsu
// Developed by Mob


#include<iostream>

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
    // destructor
    ~container() {
        delete root;
    }

    container& operator = (const container& source) {
        if (root != source.root) {
            delete root;
            copy_memoryisClear(source);
        }
        return *this;
    }
    container& operator = (container&& source) noexcept {
        if (root != source.root) {
            delete root;
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
