#include <vector>

template<typename T>
struct rank {
    static constexpr size_t value = 0;
};

static_assert(rank<int>::value == 0);

template<typename T>
struct rank<T*> {
    static constexpr size_t value = 1 + rank<T>::value;
};

static_assert(rank<int***>::value == 3);

template<typename T, size_t SIZE>
struct rank<T[SIZE]> {
    static constexpr size_t value = 1 + rank<T>::value;
};

static_assert(rank<int[10][2]>::value == 2);

template<typename T>
struct rank<T[]> {
    static constexpr size_t value = 1 + rank<T>::value;
};

static_assert(rank<int[]>::value == 1);


static_assert(rank<int* []>::value == 2);


static_assert(rank<std::vector<int>>::value == 0);

int main() {


    return 0;
}
