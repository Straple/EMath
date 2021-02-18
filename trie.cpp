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

        // ����?
        bool isLeaf() const {
            int i;
            // ���� �� ����� �� ������� � ���������� ���
            for (i = 0; i < 26 && Next[i] == 0; i++) {}
            return i == 26;
        }
    };

    pnode root = new node();

#define next()\
temp->Next[str[k] - 'a']


    // ������ 0 ���� ����� ������ ���
    pnode find(const std::string& str) {
        pnode temp = root;
        int k = 0;
        while (k < str.size() && next() != 0) {
            temp = next();
            k++;
        }
        // ���� �� �� ����� �� ����� ������: 0 ����� temp
        return k < str.size() ? 0 : temp;
    }

    // ��������� ������ � ��������
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
