// pairing heap
// Developed by Mob

// куча с минимумом
template<typename T>
class pairing_heap {

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

    // ERROR_STACKOVERFLOW
    // O(log n)
    pnode TwoPassMerge(pnode branch) {
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

    void copy(const pairing_heap& source) {
        if (source.root != nullptr) {
            root = new node();
            copy(root, source.root);
        }
    }

    void move(pairing_heap& source) {
        root = source.root;
        source.root = 0;
    }

public:

    // default constructor
    pairing_heap() {}
    // copy constructor
    pairing_heap(const pairing_heap& source) {
        copy(source);
    }
    // move constructor
    pairing_heap(pairing_heap&& source) {
        move(source);
    }
    // destructor
    ~pairing_heap() {
        clear();
    }

    pairing_heap& operator = (const pairing_heap& source) {
        if (root != source.root) {
            clear();
            copy(source);
        }
        return *this;
    }
    pairing_heap& operator = (pairing_heap&& source) {
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
        root = merge(root, new node(key, 0, 0));
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
    void join(pairing_heap& other) {
        root = merge(root, other.root);
        other.root = 0;
    }
};
