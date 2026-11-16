// physics.hpp - фізична модель та чисельні методи (без залежності від SFML).
// Модель: затухаючий гармонічний осцилятор  m x'' + c x' + k x = 0.
#pragma once
#include <cmath>
#include <limits>

namespace phys {

struct Params {
    double m = 1.0;    // маса, кг
    double k = 25.0;   // жорсткість пружини, Н/м
    double c = 0.5;    // коефіцієнт в'язкого опору, кг/с
    double x0 = 1.0;   // початкове зміщення, м
    double v0 = 0.0;   // початкова швидкість, м/с
    double T = 10.0;   // тривалість моделювання, с
};

enum class Method { Euler = 0, Symplectic = 1, RK4 = 2 };
inline constexpr Method kMethods[3] = {Method::Euler, Method::Symplectic, Method::RK4};

inline const char* name(Method m) {
    switch (m) {
        case Method::Euler: return "euler";
        case Method::Symplectic: return "symp";
        default: return "rk4";
    }
}

struct State { double x, v; };

inline double energy(const State& s, const Params& p) {
    return 0.5 * p.m * s.v * s.v + 0.5 * p.k * s.x * s.x;
}

// ---------------------------------------------------------------- аналітика
// Точний розв'язок для коливального, критичного та надкритичного режимів.
inline State analytic(double t, const Params& p) {
    const double g = p.c / (2.0 * p.m), w2 = p.k / p.m, disc = g * g - w2;
    const double e = std::exp(-g * t);
    if (std::fabs(disc) < 1e-12) {                       // критичне затухання
        const double B = p.v0 + g * p.x0;
        return {e * (p.x0 + B * t), e * (p.v0 - g * B * t)};
    }
    if (disc < 0) {                                      // коливальний режим
        const double wd = std::sqrt(-disc), B = (p.v0 + g * p.x0) / wd;
        const double cs = std::cos(wd * t), sn = std::sin(wd * t);
        const double u = p.x0 * cs + B * sn;
        const double du = -p.x0 * wd * sn + B * wd * cs;
        return {e * u, e * (-g * u + du)};
    }
    const double s = std::sqrt(disc), l1 = -g + s, l2 = -g - s;   // надкритичний
    const double A = (p.v0 - l2 * p.x0) / (l1 - l2), B = p.x0 - A;
    return {A * std::exp(l1 * t) + B * std::exp(l2 * t),
            A * l1 * std::exp(l1 * t) + B * l2 * std::exp(l2 * t)};
}

// ------------------------------------------------------------ один крок
inline State step(Method mth, const State& s, const Params& p, double dt) {
    const double w2 = p.k / p.m, d2 = p.c / p.m;
    switch (mth) {
        case Method::Euler: {   // x_{n+1} = x_n + dt v_n ; v_{n+1} = v_n + dt a_n
            const double a = -w2 * s.x - d2 * s.v;
            return {s.x + dt * s.v, s.v + dt * a};
        }
        case Method::Symplectic: {   // спочатку v, потім x з НОВОЮ v
            const double v = s.v + dt * (-w2 * s.x - d2 * s.v);
            return {s.x + dt * v, v};
        }
        default: {   // класичний РК4
            const double h2 = 0.5 * dt, h6 = dt / 6.0;
            const double k1x = s.v,                 k1v = -w2 * s.x - d2 * s.v;
            const double k2x = s.v + h2 * k1v,      k2v = -w2 * (s.x + h2 * k1x) - d2 * (s.v + h2 * k1v);
            const double k3x = s.v + h2 * k2v,      k3v = -w2 * (s.x + h2 * k2x) - d2 * (s.v + h2 * k2v);
            const double k4x = s.v + dt * k3v,      k4v = -w2 * (s.x + dt * k3x) - d2 * (s.v + dt * k3v);
            return {s.x + h6 * (k1x + 2 * k2x + 2 * k3x + k4x),
                    s.v + h6 * (k1v + 2 * k2v + 2 * k3v + k4v)};
        }
    }
}

// n кроків; obs(i, state) викликається для i = 0..n. Повертає кінцевий стан.
template <class Obs>
inline State integrate(Method mth, const Params& p, double dt, long n, Obs&& obs) {
    State s{p.x0, p.v0};
    obs(0L, s);
    for (long i = 1; i <= n; ++i) {
        s = step(mth, s, p, dt);
        obs(i, s);
    }
    return s;
}

// ------------------------------------------------- стійкість (матриця кроку)
// Лінійна система => один крок: s_{n+1} = M s_n, M = 2x2. Метод стійкий <=> rho(M) <= 1.
struct Mat2 { double a, b, c, d; };
inline Mat2 mul(const Mat2& A, const Mat2& B) {
    return {A.a * B.a + A.b * B.c, A.a * B.b + A.b * B.d,
            A.c * B.a + A.d * B.c, A.c * B.b + A.d * B.d};
}
inline Mat2 stepMatrix(Method mth, const Params& p, double dt) {
    const double w2 = p.k / p.m, d2 = p.c / p.m;
    switch (mth) {
        case Method::Euler: return {1.0, dt, -w2 * dt, 1.0 - d2 * dt};
        case Method::Symplectic: return {1.0 - w2 * dt * dt, dt * (1.0 - d2 * dt), -w2 * dt, 1.0 - d2 * dt};
        default: {
            const Mat2 Z{0.0, dt, -w2 * dt, -d2 * dt};
            const Mat2 Z2 = mul(Z, Z), Z3 = mul(Z2, Z), Z4 = mul(Z3, Z);
            return {1.0 + Z.a + Z2.a / 2 + Z3.a / 6 + Z4.a / 24, Z.b + Z2.b / 2 + Z3.b / 6 + Z4.b / 24,
                    Z.c + Z2.c / 2 + Z3.c / 6 + Z4.c / 24, 1.0 + Z.d + Z2.d / 2 + Z3.d / 6 + Z4.d / 24};
        }
    }
}
// Спектральний радіус 2x2 матриці через слід і визначник.
inline double spectralRadius(const Mat2& M) {
    const double tr = M.a + M.d, det = M.a * M.d - M.b * M.c;
    const double disc = tr * tr / 4.0 - det;
    if (disc >= 0) {
        const double r = std::sqrt(disc);
        return std::fmax(std::fabs(tr / 2 + r), std::fabs(tr / 2 - r));
    }
    return std::sqrt(det);   // комплексна пара: |lambda|^2 = det
}
inline double rho(Method mth, const Params& p, double dt) {
    return spectralRadius(stepMatrix(mth, p, dt));
}

// Найбільший dt, для якого rho <= 1: сканування сітки + бісекція. 0 - якщо такого dt немає.
inline double dtMax(Method mth, const Params& p, double tol = 1e-9) {
    const double w0 = std::sqrt(p.k / p.m);
    const int N = 20000;
    const double lo0 = 1e-5, hi0 = 8.0 / w0;
    auto ok = [&](double d) { return rho(mth, p, d) <= 1.0 + tol; };
    if (!ok(lo0)) return 0.0;
    double prev = lo0;
    for (int i = 1; i < N; ++i) {
        const double d = lo0 + (hi0 - lo0) * i / (N - 1);
        if (!ok(d)) {
            double lo = prev, hi = d;
            for (int j = 0; j < 60; ++j) {
                const double mid = 0.5 * (lo + hi);
                (ok(mid) ? lo : hi) = mid;
            }
            return lo;
        }
        prev = d;
    }
    return std::numeric_limits<double>::infinity();
}

}  // namespace phys
