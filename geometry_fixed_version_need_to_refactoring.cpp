// Сomputational Geometry: efloat, dot, line, circle, polygon, ConvexHull
// Developed by Mob

#include<cmath>
#include<iostream>
#include<vector>
#include<algorithm>
using namespace std;


typedef long double point_t;

// если true, то будут обрабатываться погрешности
#define PRECISION true

static const point_t eps = 1e-9, pi = acos(-1), inf = INFINITY;

// надстройка над point_t с правильным сравнением чисел с плавающей точкой
struct efloat {
    point_t value;

    efloat() {
        value = 0;
    }
    template<typename T>
    efloat(const T& value) {
        this->value = static_cast<point_t>(value);
    }

    bool operator == (const efloat& Rhs) const {
        return std::abs(value - Rhs.value) <= eps;
    }
    bool operator < (const efloat& Rhs) const {
        return value < Rhs.value - eps;
    }
    bool operator > (const efloat& Rhs) const {
        return value > Rhs.value + eps;
    }

    bool operator != (const efloat& Rhs) const {
        return !(*this == Rhs);
    }
    bool operator <= (const efloat& Rhs) const {
        return !(*this > Rhs);
    }
    bool operator >= (const efloat& Rhs) const {
        return !(*this < Rhs);
    }
};

#if PRECISION
#define EFLOAT efloat
#else
#define EFLOAT
#endif

// x, y
struct dot {
    point_t x, y;

    dot() {
        x = y = 0;
    }
    template<typename T1, typename T2>
    dot(const T1& x, const T2& y) {
        this->x = static_cast<point_t>(x);
        this->y = static_cast<point_t>(y);
    }

    dot operator + (const dot& p) const {
        return dot(x + p.x, y + p.y);
    }
    dot& operator += (const dot& p) {
        return *this = *this + p;
    }

    dot operator - (const dot& p) const {
        return dot(x - p.x, y - p.y);
    }
    dot& operator -= (const dot& p) {
        return *this = *this - p;
    }

    dot operator * (point_t k) const {
        return dot(x * k, y * k);
    }
    dot& operator *= (point_t k) {
        return *this = *this * k;
    }

    dot operator / (point_t k) const {
        return dot(x / k, y / k);
    }
    dot& operator /= (point_t k) {
        return *this = *this / k;
    }

    // векторное/косое произведение
    // Это площадь параллелограмма
    point_t operator % (const dot& p) const {
        return x * p.y - y * p.x;
    }
    // скалярное произведение
    point_t operator * (const dot& p) const {
        return x * p.x + y * p.y;
    }

    point_t getQuareLen() const {
        return x * x + y * y;
    }
    point_t getLen() const {
        return hypot(x, y);
    }

    bool operator == (const dot& Rhs) const {
        return EFLOAT(x) == Rhs.x && EFLOAT(y) == Rhs.y;
    }
    bool operator != (const dot& Rhs) const {
        return !(*this == Rhs);
    }

    // самая левая, потом самая нижняя
    bool operator < (const dot& Rhs) const {
        return EFLOAT(x) == Rhs.x ? EFLOAT(y) < Rhs.y : EFLOAT(x) < Rhs.x;
    }
    // самая правая, потом самая верхняя
    bool operator > (const dot& Rhs) const {
        return EFLOAT(x) == Rhs.x ? EFLOAT(y) > Rhs.y : EFLOAT(x) > Rhs.x;
    }

    dot normalize() const {
        if (*this == dot()) {
            return dot();
        }
        else {
            return *this / getLen();
        }
    }
    dot normalize(point_t mult) const {
        if (*this == dot()) {
            return dot();
        }
        else {
            return *this * (mult / getLen());
        }
    }
};

std::istream& operator >> (std::istream& input, dot& Dot) {
    return input >> Dot.x >> Dot.y;
}
std::ostream& operator << (std::ostream& output, const dot& Dot) {
    return output << Dot.x << " " << Dot.y;
}

// возвращает угол между векторами
point_t getAngle(const dot& a, const dot& b) {
    return atan2(a % b, a * b);
}
// возвращает неотрицательный угол между векторами
point_t getGoodAngle(const dot& a, const dot& b) {
    point_t res = getAngle(a, b);
    if (EFLOAT(res) < 0) {
        res += 2 * pi;
    }
    return res;
}
// возвращает неотрицательный угол меньше 180 между векторами
point_t getVeryGoodAngle(const dot& a, const dot& b) {
    point_t res = getGoodAngle(a, b);
    if (EFLOAT(res) > pi) {
        res = 2 * pi - res;
    }
    return res;
}

// a, b, c
class line {
    point_t a = 0, b = 0, c = 0; // ax * by + c = 0
    // a и b нормализованы!

    // во избежание ошибок a, b, c закрыты
    // можно написать сеттеры и геттеры

    void normalize() {
        point_t d = dot(a, b).getLen();
        if (EFLOAT(d) != 0) {
            a /= d;
            b /= d;
            c /= d;
        }
        else {
            //a = b = c = 0;
        }
        if (vector<efloat>{-a, -b, -c} < vector<efloat>{a, b, c}) {
            a = -a;
            b = -b;
            c = -c;
        }
    }

public:

    point_t get_a() const {
        return a;
    }
    point_t get_b() const {
        return b;
    }
    point_t get_c() const {
        return c;
    }

    line() {}
    // по двум точкам на прямой
    line(const dot& begin, const dot& end) {
        a = begin.y - end.y;
        b = end.x - begin.x;
        normalize();
        c = -a * begin.x - b * begin.y;
    }
    line(point_t A, point_t B, point_t C) {
        a = A;
        b = B;
        c = C;
        normalize();
    }

    // возвращает перпендикуляр из точки
    line getPerp(const dot& point) const {
        point_t A = -b, B = a;
        point_t C = -A * point.x - B * point.y;
        return line(A, B, C);
    }

    // возвращает параллельную прямую на расстоянии dist
    // если оно будет отрицательно, то вернет с другой стороны
    line getParallel(point_t dist) const {
        return line(a, b, c + dist /* * sqrt(a * a + b * b)*/);
    }

    // возвращает нормализованный вектор прямой умноженный на mult
    dot getVector(point_t mult = 1) const {
        return dot(-b, a) * mult /*.normalize(mult)*/;
    }

    // возвращает точку пересечения двух прямых
    dot intersect(const line& Rhs) const {
        point_t x, y;
        // ax + by + c = 0
        // by = -c - ax
        // y  = (-ax - c) / b

        if (EFLOAT(Rhs.b) != 0) {
            x = (b * Rhs.c / Rhs.b - c) / (a - b * Rhs.a / Rhs.b);
            y = (-x * Rhs.a - Rhs.c) / Rhs.b;
        }
        else {
            x = -Rhs.c / Rhs.a;
            y = (-x * a - c) / b;
        }
        return dot(x, y);
    }

    // возвращает точку пересечения перпендикуляра
    dot perpIntersect(const dot& point) const {
        return intersect(getPerp(point));
    }

    // отражает точки от прямой
    std::vector<dot> reflection(const std::vector<dot>& Dots) const {
        std::vector<dot> Result(Dots.size());
        for (int i = 0; i < Result.size(); i++) {
            Result[i] = Dots[i] + (perpIntersect(Dots[i]) - Dots[i]) * 2;
        }
        return Result;
    }

    // возвращает длину перпендикуляра
    point_t dist(const dot& point) const {
        return abs(a * point.x + b * point.y + c) /* / std::sqrt(a * a + b * b)*/;
    }
    point_t dist_NO_ABS(const dot& point) const {
        return a * point.x + b * point.y + c /* / std::sqrt(a * a + b * b)*/;
    }
    // возвращает расстояние между ПАРАЛЛЕЛЬНЫМИ прямыми
    point_t dist(const line& parallel) const {
        return abs(c - parallel.c) /* / sqrt(a * a + b * b)*/;
    }

    bool isParallel(const line& Rhs) const {
        return EFLOAT(a * Rhs.b - b * Rhs.a) == 0;
    }
    bool isPerp(const line& Rhs) const {
        return EFLOAT(a * Rhs.a + b * Rhs.b) == 0;
    }
    bool ison(const dot& point) const {
        return EFLOAT(a * point.x + b * point.y + c) == 0;
    }

    // прямые совпадают?
    bool operator == (const line& Rhs) const {
        return (EFLOAT(a) == Rhs.a && EFLOAT(b) == Rhs.b && EFLOAT(c) == Rhs.c) ||
            (EFLOAT(a) == -Rhs.a && EFLOAT(b) == -Rhs.b && EFLOAT(c) == -Rhs.c);
    }

    friend std::ostream& operator << (std::ostream& output, const line& Line);
};

std::istream& operator >> (std::istream& input, line& Line) {
    point_t A, B, C;
    input >> A >> B >> C;
    Line = line(A, B, C);
    return input;
}
std::ostream& operator << (std::ostream& output, const line& Line) {
    return output << Line.a << " " << Line.b << " " << Line.c;
}

// center, radius
struct circle {
    dot center;
    point_t radius;

    circle() {
        radius = 0;
    }
    circle(const dot& center, point_t radius) {
        this->center = center;
        this->radius = radius;
    }

    // returns a point on a circle
    // counterclockwise angle
    dot point(point_t angle) const {
        return center + dot(cos(angle), sin(angle)) * radius;
    }

    // 2*pi*R
    point_t getLength() const {
        return 2 * pi * radius;
    }
    // pi*R^2
    point_t getArea() const {
        return pi * radius * radius;
    }

    // пересечение окружности и прямой
    std::vector<dot> intersect(const line& Rhs) const {
        dot perp = Rhs.perpIntersect(center),
            delta = center - perp;

        point_t quareRadius = radius * radius;
        point_t quareDist = delta.getQuareLen();

        if (EFLOAT(quareRadius) > quareDist) { // two points
            point_t len = sqrt(quareRadius - quareDist);
            dot vector(Rhs.getVector(len));
            return { dot(perp + vector), dot(perp - vector) };
        }
        else if (EFLOAT(quareRadius) == quareDist) { // one point
            return { perp };
        }
        else { // none point
            return {};
        }
    }

    // пересечение двух окружностей
    // если вернет true то точек пересечения бесконечно много
    bool intersect(const circle& Rhs, std::vector<dot>& result) const {
        if (Rhs.center == center) {
            result.clear();
            return EFLOAT(radius) == Rhs.radius;
        }
        else {
            dot vector(Rhs.center - center);
            line l(-2 * vector.x, -2 * vector.y, vector.getQuareLen() + radius * radius - Rhs.radius * Rhs.radius);

            result = circle(dot(), radius).intersect(l);

            for (auto& point : result) {
                point += center;
            }

            return false;
        }
    }

    // проверка принадлежности точки окружности
    bool ison(const dot& point) const {
        return EFLOAT((center - point).getLen()) == radius;
    }
};

std::istream& operator >> (std::istream& input, circle& Circle) {
    return input >> Circle.center >> Circle.radius;
}
std::ostream& operator << (std::ostream& output, const circle& Circle) {
    return output << Circle.center << " " << Circle.radius;
}

// Многоугольник. points
struct polygon {
    std::vector<dot> Dots;

    polygon() {}
    polygon(const std::vector<dot>& Dots) {
        this->Dots = Dots;
    }

    // Возвращает площадь многоугольника
    point_t getArea() const {
        point_t result = 0;
        for (int i = 0; i < Dots.size(); i++) {
            dot p1 = i ? Dots[i - 1] : Dots.back(),
                p2 = Dots[i];
            result += (p1.x - p2.x) * (p1.y + p2.y);
        }
        return std::abs(result) * 0.5;
    }

    // Возвращает периметр многоугольника
    point_t getPerim() const {
        point_t result = (Dots[0] - Dots.back()).getLen();
        for (int i = 1; i < Dots.size(); i++) {
            result += (Dots[i] - Dots[i - 1]).getLen();
        }
        return result;
    }

    // Проверка на выпуклость многоугольника за O(n). 
    // WARNING !!No checked!!
    bool isConvexHull() const {
        int i = 2;
        while (i < Dots.size() && EFLOAT((Dots[i] - Dots[i - 2]) % (Dots[i - 1] - Dots[i - 2])) <= 0) {
            i++;
        }
        return i == Dots.size();
    }
};

// Выпуклая оболочка. polygon
struct convexHull {
    polygon poly; // точки выпуклой оболочки

    convexHull() {}
    // конструктор построения выпуклой оболочки
    convexHull(const polygon& newPolygon, bool isConvexHull = false) {
        poly = !isConvexHull ? buildConvexHull(newPolygon.Dots) : newPolygon;
    }

private:

    // Построение выпуклой оболочки за O(n * log n)
    std::vector<dot> buildConvexHull(std::vector<dot> Dots) {
        std::vector<dot> result;
        sort(Dots.begin(), Dots.end()); // Lhs.x == Rhs.x ? Lhs.y < Rhs.y : Lhs.x < Rhs.x

        dot start = Dots[0]; // начало
        for (auto& dot : Dots) {
            dot -= start;
        }

        sort(Dots.begin() + 1, Dots.end(), [](const dot& Lhs, const dot& Rhs) {
            point_t vectorProduct = Lhs % Rhs; // векторное/косое произведение = x * p.y - y * p.x
            if (EFLOAT(vectorProduct) == 0) { // лежат на одной прямой
                return EFLOAT(Lhs.getQuareLen()) < Rhs.getQuareLen();
            }
            else {
                return EFLOAT(vectorProduct) < 0;
            }
        });

        for (auto& dot : Dots) {
            while (result.size() > 1 &&
                EFLOAT((dot - result[result.size() - 2]) % (result.back() - result[result.size() - 2])) <= 0 /* <= или < тут уже философский вопрос*/) {
                result.pop_back(); // удаляем ненужные точки, которые точно не входят
            }
            result.push_back(dot);
        }

        for (auto& dot : result) {
            dot += start;
        }
        return result;
    }
};


//#include <bits/stdc++.h>
//#include <unordered_set>
#include<iostream>
#include<fstream>
#include<iomanip>
using namespace std;

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

using s8 = int8_t;
using s16 = int16_t;
using s32 = int32_t;
using s64 = int64_t;

using ld = long double;

/*
* GLOBAL PROJECT SETTINGS
*/

//#define HOME_MODE

#define DEBUG
//#define MULTITEST
//#define STRESS_TEST // включает HOME_MODE, MULTITEST и тестирует решения

//#define TASK_FILE_NAME "rag" // если есть, включает файловый ввод и вывод

/*
* SIMULATE SETTINGS
*/

#ifdef TASK_FILE_NAME
#define FILE_INPUT TASK_FILE_NAME".in"
#define FILE_OUTPUT TASK_FILE_NAME".out"
#endif

#ifdef STRESS_TEST
#define HOME_MODE // стресс-тест бывает только дома
#define MULTITEST // у нас есть несколько тестов для тестирования
#endif

#ifdef HOME_MODE
#undef FILE_INPUT
#undef FILE_OUTPUT
#define FILE_INPUT "input.txt"
#else
#undef DEBUG
#endif

/*
* SOLVE PROBLEM
*/

struct ray {
    dot begin, end;

    ray(dot begin, dot end) {
        this->begin = begin;
        this->end = end;
    }

    bool is_on(dot p) {
        return p == begin || (p - begin).normalize() == (end - begin).normalize();
    }

    bool is_inter(line l) {
        if (line(begin, end).isParallel(l)) {
            return false;
        }
        else {
            return is_on(l.intersect(line(begin, end)));
        }
    }

    ld dist(dot p) {
        dot inter = line(begin, end).perpIntersect(p);
        if (is_on(inter)) {
            return line(begin, end).dist(p);
        }
        else {
            return (p - begin).getLen();
        }
    }
};

struct segment {
    dot begin, end;

    segment(dot begin, dot end) {
        this->begin = begin;
        this->end = end;
    }

    bool is_on(dot p) {
        return EFLOAT((begin - p).getLen() + (end - p).getLen()) == (begin - end).getLen();
    }

    ld dist(dot p) {
        dot inter = line(begin, end).perpIntersect(p);

        if (is_on(inter)) {
            return line(begin, end).dist(p);
        }
        else {
            return min((p - begin).getLen(), (p - end).getLen());
        }
    }
    
};

ld dist_line_to_line(line a, line b) {
    if (a.isParallel(b)) {
        return a.dist(b);
    }
    else {
        return 0; // они пересекаются
    }
}

ld dist_line_to_ray(line l, ray r) {
    if (r.is_inter(l)) {
        return 0;
    }
    else {
        return l.dist(r.begin);
    }
}

ld dist_ray_to_ray(ray r1, ray r2) {
    line l1(r1.begin, r1.end), l2(r2.begin, r2.end);

    ld mn_dist = min(r1.dist(r2.begin), r2.dist(r1.begin));

    if (l1.isParallel(l2)) {
        return mn_dist;
    }
    else {
        dot inter = l1.intersect(l2);
        if (r1.is_on(inter) && r2.is_on(inter)) {
            return 0;
        }
        else {
            return mn_dist;
        }
    }
}

ld dist_line_to_segment(line l, segment s) {
    line line_s(s.begin, s.end);

    if (l.isParallel(line_s)) {
        return l.dist(line_s);
    }
    else {
        dot inter = l.intersect(line_s);

        if (s.is_on(inter)) {
            return 0;
        }
        else {
            return min(l.dist(s.begin), l.dist(s.end));
        }
    }
}

ld dist_ray_to_segment(ray r, segment s) {
    line line_r(r.begin, r.end), line_s(s.begin, s.end);
    
    if (r.is_inter(line_s)) {
        dot inter = line_r.intersect(line_s);

        if (s.is_on(inter)) {
            return 0;
        }
        else {
            return min(r.dist(s.begin), r.dist(s.end));
        }
    }
    else {
        if (line_r.isParallel(line_s)) {
            return min(r.dist(s.begin), r.dist(s.end));
        }
        else {
            dot inter = line_s.perpIntersect(r.begin);
            if (s.is_on(inter)) {
                return line_s.dist(r.begin);
            }
            else {
                return min(r.dist(s.begin), r.dist(s.end));
            }
        }
    }
}

ld dist_segment_to_segment(segment s1, segment s2) {
    line l1(s1.begin, s1.end), l2(s2.begin, s2.end);

    if (l1.isParallel(l2)) {
        return min({ s1.dist(s2.begin), s1.dist(s2.end), s2.dist(s1.begin), s2.dist(s1.end) });
    }
    else {

        dot inter = l1.intersect(l2);

        if (s1.is_on(inter) && s2.is_on(inter)) {
            return 0;
        }
        else {
            return min({ s1.dist(s2.begin), s1.dist(s2.end), s2.dist(s1.begin), s2.dist(s1.end) });
        }
    }
}

void solve(istream& cin, ostream& cout) {
    dot a, b, c, d;
    cin >> a >> b >> c >> d;


#define end << endl
    
    cout << fixed << setprecision(10);

#ifdef DEBUG

    cout << (a - c).getLen() << " от точки A до точки C\n";
    cout << segment(c, d).dist(a) <<  " от точки A до отрезка CD\n";
    cout << ray(c, d).dist(a) << " от точки A до луча CD\n";
    cout << line(c, d).dist(a) << " от точки A до прямой CD\n";
    cout << segment(a, b).dist(c) << " от отрезка AB до точки C\n";
    cout << dist_segment_to_segment(segment(a, b), segment(c, d)) << " от отрезка AB до отрезка CD\n";
    cout << dist_ray_to_segment(ray(c, d), segment(a, b)) << " от отрезка AB до луча CD\n";
    cout << dist_line_to_segment(line(c, d), segment(a, b))<< " от отрезка AB до прямой CD\n";
    cout << ray(a, b).dist(c) << " от луча AB до точки C\n";
    cout << dist_ray_to_segment(ray(a, b), segment(c, d)) << " от луча AB до отрезка CD\n";
    cout << dist_ray_to_ray(ray(a, b), ray(c, d)) << " от луча AB до луча CD\n";
    cout << dist_line_to_ray(line(c, d), ray(a, b)) << " от луча AB до прямой CD\n";
    cout << line(a, b).dist(c) << " от прямой AB до точки C\n";
    cout << dist_line_to_segment(line(a, b), segment(c, d)) << " от прямой AB до отрезка CD\n";
    cout << dist_line_to_ray(line(a, b), ray(c, d)) << " от прямой AB до луча CD\n";
    cout << dist_line_to_line(line(a, b), line(c, d)) << " от прямой AB до прямой CD\n";

#else
    cout << (a - c).getLen() end;
    cout << segment(c, d).dist(a) end;
    cout << ray(c, d).dist(a) end;
    cout << line(c, d).dist(a) end;
    cout << segment(a, b).dist(c) end;
    cout << dist_segment_to_segment(segment(a, b), segment(c, d)) end;
    cout << dist_ray_to_segment(ray(c, d), segment(a, b)) end;
    cout << dist_line_to_segment(line(c, d), segment(a, b)) end;
    cout << ray(a, b).dist(c) end;
    cout << dist_ray_to_segment(ray(a, b), segment(c, d)) end;
    cout << dist_ray_to_ray(ray(a, b), ray(c, d)) end;
    cout << dist_line_to_ray(line(c, d), ray(a, b)) end;
    cout << line(a, b).dist(c) end;
    cout << dist_line_to_segment(line(a, b), segment(c, d)) end;
    cout << dist_line_to_ray(line(a, b), ray(c, d)) end;
    cout << dist_line_to_line(line(a, b), line(c, d)) end;
#endif
}


void main_build() {
    
}

void correct_solve(istream& cin, ostream& cout) {
    
}

#ifdef STRESS_TEST

#include <random>
mt19937_64 rnd(42);

// строит тест
string build_test() {
    stringstream out;

    // вывести тест

    return out.str();
}

#endif


/*
* TEMPLATE MAIN
*/

int main() {
    setlocale(LC_ALL, "Russian");
#ifdef FILE_INPUT
    ifstream cin(FILE_INPUT);
#endif
#ifdef FILE_OUTPUT
    ofstream cout(FILE_OUTPUT);
#endif
    ios::sync_with_stdio(false), cout.tie(nullptr), cin.tie(nullptr);

    main_build();

    s64 t = 1;
#ifdef MULTITEST
    cin >> t;
#endif

    while (t--) {
#ifdef STRESS_TEST
        string test_str_input = build_test(); // создаем тест

        stringstream input1 = stringstream(test_str_input);
        stringstream input2 = stringstream(test_str_input);

        // получаем ответы на тест
        stringstream wrong_out, correct_out;
        solve(input1, wrong_out);
        correct_solve(input2, correct_out);

        // сравниваем их
        if (wrong_out.str() != correct_out.str()) {
            cout << "[failed test]\n" << test_str_input << "\n";
            
            cout << "[wrong output]\n" << wrong_out.str() << "\n";

            cout << "[correct output]\n" << correct_out.str() << "\n";

            return 0;
        }
#else
        solve(cin, cout);
#endif
    }

#ifdef STRESS_TEST
    cout << "[all tests passed]\n";
#endif

    return 0;
}
