#include<cstring>
#include<iomanip>
#include<algorithm>
#include<vector>
#include<stack>
#include<cmath>
#include<queue>
#include<random>
#include<ctime>

// Developed by Mr_Straple

// Data Structures
// tos(tree of segments), edeque, container
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
        edeque(const edeque<T>& new_value) {
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

        edeque& operator = (const edeque<T>& new_value) {
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

        void clear() {
            delete[] A;
            zeroing();
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
            return index < length_A ? A[index] : A[index - length_A];
        }

        T& back() const {
            return A[dec(tail)];
        }
        T& front() const {
            return A[head];
        }
    };

    // For container
    std::mt19937_64 random(42);

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
                prior = random();
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
}

// Сomputational Geometry
// edouble, dot, line, circle
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
        edouble(const edouble& new_value) {
            value = new_value.value;
        }

        bool operator == (const edouble& Rhs) const {
            return std::abs(value - Rhs.value) < eps;
        }
        bool operator < (const edouble& Rhs) const {
            return value < Rhs.value - eps;
        }
        bool operator > (const edouble& Rhs) const {
            return value + eps > Rhs.value;
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

        long double operator + (const edouble& add) const {
            return value + add.value;
        }
        long double operator - (const edouble& subtrahend) const {
            return value - subtrahend.value;
        }
        long double operator * (const edouble& mult) const {
            return value * mult.value;
        }
        long double operator / (const edouble& div) const {
            return value / div.value;
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

    template<typename T>
    long double operator + (const T& value, const edouble& add) {
        return edouble(value) + add;
    }
    template<typename T>
    long double operator - (const T& value, const edouble& subtrahend) {
        return edouble(value) - subtrahend;
    }
    template<typename T>
    long double operator * (const T& value,const edouble& mult) {
        return edouble(value) * mult;
    }
    template<typename T>
    long double operator / (const T& value, const edouble& div) {
        return edouble(value) / div;
    }


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
            return std::abs(a * point.x + b * point.y + c) / std::sqrt(a * a + b * b);
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

    // dist^2 = r^2
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
}

// Long Arithmetic
// elong
namespace emt {

    static const long long long_base = 1000000000;

    // Long operators
    namespace lgp {
        void remove_leading_zeros(dst::edeque<int>& digits) {
            while (!digits.empty() && digits.back() == 0) {
                digits.pop_back();
            }
        }

        enum class OP{
            less,
            equally,
            more,
        };

        // compare two numbers
        OP compare(const dst::edeque<int>& a, const dst::edeque<int>& b) {
            if (a.size() != b.size()) {
                return a.size() < b.size() ? OP::less : OP::more;
            }
            else {
                int i = a.size() - 1;
                while (i >= 0 && a[i] == b[i]) {
                    i--;
                }
                if (i >= 0) {
                    return a[i] < b[i] ? OP::less : OP::more;
                }
                else {
                    return OP::equally;
                }
            }
        }
        
        bool less(const dst::edeque<int>& a, const dst::edeque<int>& b) {
            return compare(a, b) == OP::less;
        }
        bool equally(const dst::edeque<int>& a, const dst::edeque<int>& b) {
            return compare(a, b) == OP::equally;
        }
        bool more(const dst::edeque<int>& a, const dst::edeque<int>& b) {
            return compare(a, b) == OP::more;
        }

        // return a + b
        dst::edeque<int> addition(const dst::edeque<int>& a, const dst::edeque<int>& b) {
            dst::edeque<int> result = a;
            bool k = 0;
            int i = 0;
            for (i = 0; i < std::max(result.size(), b.size()) || k != 0; i++) {
                if (i == result.size()) {
                    result.push_back(0);
                }
                result[i] += k + (i < b.size() ? b[i] : 0);
                k = result[i] >= long_base;
                if (k != 0) {
                    result[i] -= long_base;
                }
            }
            return result;
        }

        // return minuend - subtrahend
        dst::edeque<int> subtraction(const dst::edeque<int>& minuend, const dst::edeque<int>& subtrahend) {
            dst::edeque<int> result = minuend;

            bool k = 0;
            int i;

            for (i = 0; i < subtrahend.size() || k != 0; i++) {
                result[i] -= k + (i < subtrahend.size() ? subtrahend[i] : 0);
                k = result[i] < 0;
                if (k != 0) {
                    result[i] += long_base;
                }
            }

            remove_leading_zeros(result);
            return result;
        }

        // return a * b
        dst::edeque<int> multiplication(const dst::edeque<int>& a, const dst::edeque<int>& b) {
            dst::edeque<int> result(a.size() + b.size() + 1, 0);

            long long c, k;
            int i, j;
            for (i = 0; i < a.size(); i++) {
                c = 0;
                for (j = 0; j < b.size() || c != 0; j++) {
                    k = result[i + j] + a[i] * 1LL * (j < b.size() ? b[j] : 0) + c;
                    result[i + j] = k % long_base;
                    c = k / long_base;
                }
            }

            remove_leading_zeros(result);
            return result;
        }

        // return dividend / divider
        dst::edeque<int> division(const dst::edeque<int>& dividend, const dst::edeque<int>& divider) {
            dst::edeque<int> result, value, temp, t;

            int i = dividend.size() - 1;
            while (i >= 0 && less(value, divider)) {
                value.push_front(dividend[i]);
                i--;
            }
            value.pop_front();
            i++;
            for (; i >= 0; i--) {
                remove_leading_zeros(value);

                value.push_front(dividend[i]);

                int l = 0, r = long_base;
                while (l < r - 1) {
                    int m = (l + r) / 2;

                    temp.push_back(m);

                    t = multiplication(divider, temp);

                    temp.pop_back();

                    if (less(value, t)) {
                        r = m;
                    }
                    else {
                        l = m;
                    }
                }

                result.push_front(l);

                temp.push_back(l);

                value = subtraction(value, multiplication(divider, temp));

                temp.pop_back();
            }
            return result;
        }

        // return dividend % divider
        dst::edeque<int> remainder_of_the_division(const dst::edeque<int>& dividend, const dst::edeque<int>& divider) {
            dst::edeque<int> result, temp, t;

            int i = dividend.size() - 1;
            while (i >= 0 && less(result, divider)) {
                result.push_front(dividend[i]);
                i--;
            }
            result.pop_front();
            i++;
            for (; i >= 0; i--) {
                remove_leading_zeros(result);
                result.push_front(dividend[i]);

                int l = 0, r = long_base;
                while (l < r - 1) {
                    int m = (l + r) / 2;

                    temp.push_back(m);

                    t = multiplication(divider, temp);

                    temp.pop_back();

                    if (less(result, t)) {
                        r = m;
                    }
                    else {
                        l = m;
                    }
                }

                temp.push_back(l);

                result = subtraction(result, multiplication(divider, temp));

                temp.pop_back();
            }
            return result;
        }

        void str_to_num(dst::edeque<int>& digits, const std::string& str) {
            int i;
            for (i = static_cast<int>(str.size()); i >= 9; i -= 9) {
                digits.push_back(atoi(str.substr(static_cast<long long>(i) - 9, 9).c_str()));
            }
            if (i > 0) {
                digits.push_back(atoi(str.substr(0, i).c_str()));
            }
        }

        void output(std::ostream& out, const dst::edeque<int>& digits) {
            out << digits.back();
            char k = out.fill('0');
            for (int i = digits.size() - 2; i >= 0; i--) {
                out << std::setw(9) << digits[i];
            }
            out.fill(k);
        }
    }

    // signed long integer
    class elong {
    
        dst::edeque<int> digits;
        bool is_negative;

        // update
        elong& update() {
            if (digits.empty()) {
                digits.push_back(0);
            }
            return *this;
        }

        elong(const dst::edeque<int>& _digits, bool _is_negative) {
            digits = _digits;
            is_negative = _is_negative;
        }

        void inc() {
            digits.front()++;
            if (digits.front() == long_base) {
                digits.front() = 0;

                int i;
                for (i = 1; i < digits.size() && static_cast<long long>(digits[i]) + 1 == long_base; i++) {
                    digits[i] = 0;
                }
                if (i == digits.size()) {
                    digits.push_back(1);
                }
                else {
                    digits[i]++;
                }
            }
        }
        void dec() {
            digits.front()--;
            if (digits.front() == -1) {
                digits.front() = 0;

                int i;
                for (i = 1; i < digits.size() && digits[i] == 0; i++) {
                    digits[i] = long_base - 1;
                }

                if (i < digits.size()) {
                    digits[i]--;
                    if (digits.back() == 0) {
                        digits.pop_back();
                    }
                }
                else {
                    is_negative = !is_negative;
                    digits[0] = 1;
                }
            }

            if (digits.size() == 1 && digits.front() == 0) {
                is_negative = false;
            }
        }

        void signedInc() {
            if (is_negative) {
                dec();
            }
            else {
                inc();
            }
        }
        void signedDec() {
            if (is_negative) {
                inc();
            }
            else {
                dec();
            }
        }
        
    public:

        // default constructor
        elong() {
            is_negative = false;
            digits.push_back(0);
        };

        // convert constructors[

        elong(const std::string& String) {
            is_negative = String[0] == '-';
            if (is_negative) {
                lgp::str_to_num(digits, String.substr(1));
            }
            else {
                lgp::str_to_num(digits, String);
            }
        }
        elong(const char* String) {
            *this = elong(std::string(String));
        }
        template<typename integer_type>
        elong(integer_type number) {
            is_negative = number < 0;

            if (number == 0) {
                digits.push_back(0);
            }
            else {
                long long value = static_cast<long long>(number);
                value = is_negative ? -value : value;
                while (value > 0) {
                    digits.push_back(value % long_base);
                    value /= long_base;
                }
            }
        }
        // ]

        //bool operators [

        bool operator == (const elong& Rhs) const {
            return is_negative != Rhs.is_negative ? false : lgp::equally(this->digits, Rhs.digits);
        }
        bool operator != (const elong& Rhs) const {
            return !(*this == Rhs);
        }
        bool operator < (const elong& Rhs) const {
            if (is_negative) {
                return Rhs.is_negative ? -*this > -Rhs : true;
            }
            else {
                return Rhs.is_negative ? false : lgp::less(this->digits, Rhs.digits);
            }
        }
        bool operator > (const elong& Rhs) const {
            if (is_negative) {
                return Rhs.is_negative ? -*this < -Rhs : false;
            }
            else {
                return Rhs.is_negative ? true : lgp::more(this->digits, Rhs.digits);
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
            return elong(digits, !is_negative);
        }

        // operators [

        const elong& operator ++() {
            signedInc();
            return *this;
        }
        const elong& operator ++(int) {
            signedInc();
            return *this;
        }
        const elong& operator --() {
            signedDec();
            return *this;
        }
        const elong& operator --(int) {
            signedDec();
            return *this;
        }

        elong operator + (const elong& added) const {
            if (is_negative) {
                return added.is_negative ? -(-*this + (-added)) : added - (-*this);
            }
            else {
                return added.is_negative ? *this - (-added) : elong(lgp::addition(digits, added.digits), false);
            }
        }
        elong& operator += (const elong& added) {
            return *this = *this + added;
        }

        elong operator - (const elong& subtrahend) const {
            if (subtrahend.is_negative) {
                return *this + (-subtrahend);
            }
            else if (is_negative) {
                return -(-*this + subtrahend);
            }
            else {
                return *this < subtrahend ? -(subtrahend - *this) : elong(lgp::subtraction(digits, subtrahend.digits), false).update();
            }
        }
        elong& operator -= (const elong& subtrahend) {
            return *this = *this - subtrahend;
        }

        elong operator * (const elong& multiplied) const {
            return elong(lgp::multiplication(digits, multiplied.digits), is_negative != multiplied.is_negative).update();
        }
        elong& operator *= (const elong& multiplied) {
            return *this = *this * multiplied;
        }
        
        elong operator / (const elong& divider) const {
            return elong(lgp::division(digits, divider.digits), is_negative != divider.is_negative).update();
        }
        elong& operator /= (const elong& divider) {
            return *this = *this / divider;
        }

        elong operator % (const elong& divider) const {
            return divider == 2 ? elong(digits[0] & 1) : elong(lgp::remainder_of_the_division(digits, divider.digits), false).update();
        }
        elong& operator %= (const elong& divider) {
            return *this = *this % divider;
        }
        // ]

        // friend function :)
        friend std::ostream& operator << (std::ostream& out, const elong& var);
    };

    template<typename num_variable>
    elong operator + (const num_variable& a, const elong& b) {
        return elong(a) + b;
    }
    template<typename num_variable>
    elong operator - (const num_variable& a, const elong& b) {
        return elong(a) - b;
    }
    template<typename num_variable>
    elong operator * (const num_variable& a, const elong& b) {
        return elong(a) * b;
    }
    template<typename num_variable>
    elong operator / (const num_variable& a, const elong& b) {
        return elong(a) / b;
    }
    template<typename num_variable>
    elong operator % (const num_variable& a, const elong& b) {
        return elong(a) % b;
    }
    

    // output signed long integer
    std::ostream& operator << (std::ostream& out, const elong& var) {
        if (var.is_negative) {
            out << '-';
        }
        lgp::output(out, var.digits);
        return out;
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

// Algorithm
namespace alg {
    // factorizes the number
    template<typename T>
    std::vector<T> factorize(const T& N) {
        std::vector<T> result;
        for (T i = 2; i * i <= N; i++) {
            while (N % i == 0) {
                result.push_back(i);
                N /= i;
            }
        }
        if (N != 1) {
            result.push_back(N);
        }
        return result;
    }

    // binary pow: a^n % mod
    template<typename T>
    T pow(const T& a, const T& n, const T& mod) {
        if (n == 0) {
            return 1;
        }
        else {
            T z = pow(a, n / 2, mod);
            z = (z * z) % mod;

            return (n % 2 == 0) ? z : (z * a) % mod;
        }
    }

    // binary pow a^n
    template<typename T>
    T pow(const T& a, const T& n) {
        if (n == 0) {
            return 1;
        }
        else {
            T z = pow(a, n / 2);
            z *= z;

            return (n % 2 == 0 ? z : z * a);
        }
    }

    // gcd(a, b)
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

    //lcm(a, b) = a * b / gcd(a, b)
    template<typename T>
    T lcm(const T& a, const T& b) {
        return a * b / gcd(a, b);
    }
}

// Timepiece
// etime
namespace tmp{
    class etime {
        int time_begin;

    public:
        etime() {
            start();
        }

        void start() {
            time_begin = clock();
        }
        int operator *()const {
            return clock() - time_begin;
        }
    };
    std::ostream& operator << (std::ostream& out, const etime& t) {
        return out << *t;
    }
};

// Utils
namespace utl {
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
        size_t len = a.size();
        a.resize(len + b.size());
        for (size_t i = 0; i < b.size(); i++) {
            a[i + len] = b[i];
        }
        return a;
    }

    template<typename T>
    T& clamp(const T& min, const T& value, const T& max) {
        return value > max ? max :
            value < min ? min :
            value;
    }
}
using namespace utl;
