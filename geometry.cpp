#include<cmath>
#include<iostream>
#include<vector>
#include<algorithm>


// Сomputational Geometry: efloat, dot, line, circle, polygon, ConvexHull
namespace mpg {
    typedef double point_t;

// если true, то будут обрабатываться погрешности
#define PRECISION true

    point_t eps = 1e-9, pi = acos(-1), inf = INFINITY;

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
            return sqrt(getQuareLen());
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

        dot normalize(point_t mult = 1) const {
            return *this * (mult / getLen());
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
        if (res < 0) {
            res += 2 * pi;
        }
        return res;
    }
    // возвращает неотрицательный угол меньше 180 между векторами
    point_t getVeryGoodAngle(const dot& a, const dot& b) {
        point_t res = getGoodAngle(a, b);
        if (res > pi) {
            res = 2 * pi - res;
        }
        return res;
    }

    // a, b, c
    class line {
        point_t a, b, c; // ax * by + c = 0
        // a и b нормализованы!

        // во избежание ошибок a, b, c закрыты
        // можно написать сеттеры и геттеры
    public:

        line() {
            a = b = c = 0;
        }
        // по двум точкам на прямой
        line(const dot& begin, const dot& end) {
            a = begin.y - end.y;
            b = end.x - begin.x;
            // normalize
            {
                point_t d = sqrt(a * a + b * b);
                if (EFLOAT(d) != 0) {
                    a /= d;
                    b /= d;
                }
            }
            c = -a * begin.x - b * begin.y;
        }
        line(point_t A, point_t B, point_t C) {
            a = A;
            b = B;
            c = C;
            // normalize
            {
                point_t d = sqrt(a * a + b * b);
                if(EFLOAT(d) != 0){
                    a /= d;
                    b /= d;
                    c /= d;
                }
            }
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
            return (EFLOAT(a) == Rhs.a  && EFLOAT(b) == Rhs.b  && EFLOAT(c) == Rhs.c) ||
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
            return EFLOAT((center - point).getQuareLen()) == radius * radius;
        }
    };
    std::istream& operator >> (std::istream& input, circle& Circle) {
        return input >> Circle.center >> Circle.radius;
    }
    std::ostream& operator << (std::ostream& output, const circle& Circle) {
        return output << Circle.center << " " << Circle.radius;
    }

    // луч
    struct ray {
        dot begin, // начало луча
            end; // конец луча

        ray(){}

        // проверка принадлежности точки лучу
        bool ison(const dot& point) const {
            return EFLOAT((end - begin) % (point - begin)) == 0 && // точка лежит на прямой
                EFLOAT((begin - point) * (begin - end)) >= 0;   // точка лежит на луче
        }
    };
    std::istream& operator >> (std::istream& input, ray& Ray) {
        return input >> Ray.begin >> Ray.end;
    }
    std::ostream& operator << (std::ostream& output, const ray& Ray) {
        return output << Ray.begin << " " << Ray.end;
    }

    // отрезок
    struct segment {
        dot begin, // начало отрезка
            end;   // конец отрезка

        segment(){}
        segment(const dot& begin, const dot& end) {
            this->begin = begin;
            this->end = end;
        }

        // проверка принадлежности точки отрезку
        bool ison(const dot& point) const {
            return EFLOAT((end - begin) % (point - begin)) == 0 && // точка лежит на прямой
                   EFLOAT((point - begin) * (point - end)) <= 0;   // точка лежит между begin и end
        }
        // проверка пересечения двух отрезков
        bool isIntersect(const segment& Rhs) const {
            // !!! warning no checked

            line l1(begin, end), l2(Rhs.begin, Rhs.end);
            if (l1 == l2) {
                return true;
            }
            else if (l1.isParallel(l2)) {
                return false;
            }
            else {
                dot inter = line(begin, end).intersect(line(Rhs.begin, Rhs.end));

                return ison(inter) && Rhs.ison(inter);
            }


            // 52/63 tests complited. 

            // если точка одного отрезка равна точке другого отрезка
            //if (begin == Rhs.begin || begin == Rhs.end || end == Rhs.begin || end == Rhs.end) {
            //    return true;
            //}
            //else {
            //    return EFLOAT(((begin - end) % (begin - Rhs.begin)) * ((begin - end) % (begin - Rhs.end))) < 0;//&&
            //        //EFLOAT(((Rhs.begin - Rhs.end) % (Rhs.begin - begin)) * ((Rhs.begin - Rhs.end) % (Rhs.begin - end))) < 0;
            //}
        }
    };
    std::istream& operator >> (std::istream& input, segment& seg) {
        return input >> seg.begin >> seg.end;
    }
    std::ostream& operator << (std::ostream& output, const segment& seg) {
        return output << seg.begin << " " << seg.end;
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
            sort(Dots.begin(), Dots.end());
            dot start = Dots[0];
            for (int i = 0; i < Dots.size(); i++) {
                Dots[i] -= start;
            }
            sort(Dots.begin() + 1, Dots.end(), [](const dot& Lhs, const dot& Rhs) {
                point_t vectorProduct = Lhs % Rhs;
                if (EFLOAT(vectorProduct) == 0) {
                    return EFLOAT(Lhs.getQuareLen()) < Rhs.getQuareLen();
                }
                else {
                    return EFLOAT(vectorProduct) < 0;
                }
            });

            result.push_back(Dots[0]);
            for (int i = 1; i < Dots.size(); i++) {
                while (result.size() > 1 &&
                    EFLOAT((Dots[i] - result[result.size() - 2]) % (result.back() - result[result.size() - 2])) <= 0) {
                    result.pop_back();
                }
                result.push_back(Dots[i]);
            }

            for (int i = 0; i < result.size(); i++) {
                result[i] += start;
            }
            return result;
        }
    };
}
