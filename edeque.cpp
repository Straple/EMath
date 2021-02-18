// edeque
// Developed by Mob


#include<algorithm>

// deque on a circular array.
template<typename T>
class edeque {
    // variables
    T* A; // ������
    int length_A; // ������ �
    int length; // ������ �������
    int head; // ������

    // ��������� ����� �� ��������� ������
    bool overflow() const {
        return length >= length_A;
    }
    // ��������� ����� �� ��������� ������
    bool underflow() const {
        return (length << 1) < length_A;
    }


    // ����� ������ ������� � "�������"
    int sizeInc(int Size) const {
        return 1.4 * Size + 1.0;
    }


    // ����������� ������ �� 1
    int inc(int index) const {
        index++;
        return index < length_A ? index : 0;
    }
    // ��������� ������ �� 1
    int dec(int index) const {
        index--;
        return index >= 0 ? index : length_A - 1;
    }
    // ���� ������ �������� � A
    int getIndex(int index) const {
        index += head;
        return index + (index < length_A ? 0 : -length_A);
    }


    // �������� ���������
    void zeroing() {
        A = 0;
        length_A = length = head = 0;
    }

    // ������������� ������. �������� ������ resize ��� update/overflow/underflow
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

    // �������� ��������� � �������� newSize
    void allocate(int newSize) {
        A = new T[length_A = sizeInc(length = newSize)];
        head = 0;
    }

    // ��������� �������� ���������, �� �� ������� ������
    void copy(const edeque& source) {
        A = new T[length_A = source.length_A];
        length = source.length;
        head = source.head;

        for (int i = 0; i < length_A; i++) {
            A[i] = source.A[i];
        }
    }

    // move structor
    void dst_move(edeque& source) {
        // ��������� �����������
        A = source.A;
        head = source.head;
        length = source.length;
        length_A = source.length_A;
        // �������� source
        source.zeroing();
    }

public:

    // ����������� �� ���������
    edeque() {
        zeroing();
    }
    // ����������� �����������
    edeque(const edeque& source) {
        copy(source);
    }
    // ����������� �����������
    edeque(edeque&& source) noexcept {
        dst_move(source);
    }
    // ����������
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
            dst_move(source);
        }
        return *this;
    }

    // ����������� ����������/���������
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

    template<typename array_t>
    edeque(const array_t& source) {
        allocate(source.size());

        auto it = source.begin();
        for (int i = 0; i < length; i++, it++) {
            A[i] = static_cast<T>(*it);
        }
    }

    // ����������� ������
    void clear() {
        delete[] A;
        zeroing();
    }

    // �������� ������ �������, �������� ����������� ��������
    void resize(int newSize) {
        // array copy
        {
            T* temp = new T[sizeInc(newSize)];
            int i = 0, Min = std::min(newSize, length);
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

    int size() const {
        return length;
    }
    bool empty() const {
        return length == 0;
    }

private:
    // ���� ����� ��������� ������ �������
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
    // ���� ����� ��������� ������ �������
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

    void insert(int index, T&& value) {
        // ������ ���������� � ������ ��������
        if (index >= (length >> 1)) {
            push_back(move(value));
            for (int i = length - 1; i != index; i--) {
                std::swap(this->operator[](i), this->operator[](i - 1));
            }
        }
        else {
            push_front(move(value));
            for (int i = 0; i != index; i++) {
                std::swap(this->operator[](i), this->operator[](i + 1));
            }
        }
    }

    void insert(int index, const T& value) {
        // ������ ���������� � ������ ��������
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

    void erase(int index) {
        // ������ ���������� � ������ ��������
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
};
