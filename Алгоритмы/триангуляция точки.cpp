// сколько точек находяться внутри выпуклой оболочки. O(n * log n)

ld eps = 1e-9;

struct dot {
    ld x, y;

    dot(ld X, ld Y) {
        x = X;
        y = Y;
    }
    dot() {
        x = y = 0;
    }

    dot operator + (const dot& p) const {
        return dot(x + p.x, y + p.y);
    }
    dot operator - (const dot& p) const {
        return dot(x - p.x, y - p.y);
    }
    dot operator * (ld k) const {
        return dot(x * k, y * k);
    }
    dot operator / (ld k) const {
        return dot(x / k, y / k);
    }

    // векторное/косое произведение
    // Это площадь параллелограмма
    ld operator % (const dot& p) const {
        return x * p.y - y * p.x;
    }
    // скалярное произведение
    ld operator * (const dot& p) const {
        return x * p.x + y * p.y;
    }

    ld sqrt_len() const {
        return x * x + y * y;
    }
    ld len() const {
        return hypot(x, y);
    }

    dot normalize() const {
        return (*this) / len();
    }
};

template<typename T>
int sign(const T& value) {
    if (value < -eps) {
        return -1;
    }
    else if (value > eps) {
        return 1;
    }
    else {
        return 0;
    }
}

bool is_in_triangle(const vector<dot>& pts, dot p) {
    int s1 = sign((pts[1] - pts[0]) % (p - pts[0]));
    int s2 = sign((pts[2] - pts[1]) % (p - pts[1]));
    int s3 = sign((pts[0] - pts[2]) % (p - pts[2]));

    return (s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0);
}

bool is_in_convex_hull(const vector<dot>& pts, dot p) {
    int n = pts.size();

    int tl = 1, tr = n - 1;
    while (tl < tr - 1) {
        int tm = (tl + tr) / 2;

        if ((pts[tm] - pts[0]) % (p - pts[0]) < 0) {
            tr = tm;
        }
        else {
            tl = tm;
        }
    }
    return is_in_triangle({ pts[0], pts[tl], pts[tr] }, p);
}

void solve(istream& cin, ostream& cout) {
    int n, m, k;
    cin >> n >> m >> k;

    vector<dot> pts(n);
    for (int i = 0; i < n; i++) {
        cin >> pts[i].x >> pts[i].y;
    }

    int cnt = 0;
    for (int j = 0; j < m; j++) {
        dot p;
        cin >> p.x >> p.y;
        cnt += is_in_convex_hull(pts, p);
    }

    cout << (cnt >= k ? "YES" : "NO") << "\n";
}
