// trie
// Developed by Mob

#include<string>

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
