// experiments.cpp - серія експериментів для лабораторної роботи №5, варіант №34.
// Збірка:  g++ -O2 -std=c++17 experiments.cpp -o experiments
// Запуск:  ./experiments <каталог_для_CSV>
// Усі числа в тезах беруться з файлів, які створює ця програма.
#include "physics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace phys;

static std::string OUT = "data";
static volatile double g_sink = 0.0;   // щоб оптимізатор не викинув цикл при вимірюванні часу

static FILE* openCsv(const std::string& name, const char* header) {
    FILE* f = std::fopen((OUT + "/" + name).c_str(), "w");
    if (!f) { std::perror(name.c_str()); std::exit(1); }
    std::fprintf(f, "%s\n", header);
    return f;
}

// ------------------------------------------------------------- один запуск
struct Run {
    long n = 0;
    double abs_max = 0, abs_end = 0, rel = 0;
    double time = 0, time_std = 0;
    double E_end = 0, E_max = 0, rho = 0;
};

// Вимірювання часу: чистий цикл інтегрування (без збереження траєкторії).
// Кожен вимір = `inner` повторів, щоб тривалість значно перевищувала розрізнення годинника;
// `repeats` вимірів -> середнє та стандартне відхилення часу ОДНОГО прогону.
static void timeMethod(Method mth, const Params& p, double dt, long n, int repeats,
                       double& mean, double& sd) {
    const long inner = std::max(1L, 4000000L / n);
    auto once = [&]() {
        double acc = 0;
        for (long r = 0; r < inner; ++r) {
            Params q = p;
            double d = dt;
            asm volatile("" : "+m"(q.x0), "+m"(d));   // бар'єр: заборона виносу обчислень з циклу
            const State s = integrate(mth, q, d, n, [](long, const State&) {});
            acc += s.x;
        }
        g_sink = g_sink + acc;
    };
    once();   // прогрів
    std::vector<double> ts;
    for (int i = 0; i < repeats; ++i) {
        const auto t0 = std::chrono::steady_clock::now();
        once();
        const auto t1 = std::chrono::steady_clock::now();
        ts.push_back(std::chrono::duration<double>(t1 - t0).count() / inner);
    }
    double s = 0;
    for (double t : ts) s += t;
    mean = s / ts.size();
    double v = 0;
    for (double t : ts) v += (t - mean) * (t - mean);
    sd = ts.size() > 1 ? std::sqrt(v / (ts.size() - 1)) : 0.0;
}

static Run runCase(Method mth, const Params& p, double dt, int repeats) {
    Run r;
    r.n = std::lround(p.T / dt);
    double emax_x = 0, amax = 0, Emax = 0, E0 = energy({p.x0, p.v0}, p), Eend = 0;
    bool finite = true;
    integrate(mth, p, dt, r.n, [&](long i, const State& s) {
        const State e = analytic(i * dt, p);
        const double err = std::fabs(s.x - e.x);
        if (!std::isfinite(err)) finite = false;
        amax = std::max(amax, err);
        emax_x = std::max(emax_x, std::fabs(e.x));
        r.abs_end = err;
        const double E = energy(s, p);
        Emax = std::max(Emax, E);
        Eend = E;
    });
    r.abs_max = finite ? amax : INFINITY;
    r.rel = r.abs_max / emax_x;
    r.E_end = Eend / E0;
    r.E_max = Emax / E0;
    r.rho = rho(mth, p, dt);
    timeMethod(mth, p, dt, r.n, repeats, r.time, r.time_std);
    return r;
}

static double fitSlope(const std::vector<double>& x, const std::vector<double>& y) {
    const int n = (int)x.size();
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int i = 0; i < n; ++i) {
        const double a = std::log(x[i]), b = std::log(y[i]);
        sx += a; sy += b; sxx += a * a; sxy += a * b;
    }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

static const Params BASE{};   // m=1, k=25, c=0.5, x0=1, v0=0, T=10
static const int REPEATS = 30;

// ============================================ 0. перевірка еталона і реалізації
static void exp0() {
    std::printf("== EXP0 verification\n");
    FILE* f = openCsv("verify.csv", "check,case,value");
    // (a) аналітичний розв'язок проти RK4 з дуже малим кроком (незалежний еталон)
    for (double c : {0.5, 10.0, 20.0}) {
        Params p = BASE; p.c = c;
        const double dt = 1e-4;
        const long n = std::lround(p.T / dt);
        double dmax = 0;
        integrate(Method::RK4, p, dt, n, [&](long i, const State& s) {
            const State e = analytic(i * dt, p);
            dmax = std::max({dmax, std::fabs(s.x - e.x), std::fabs(s.v - e.v)});
        });
        std::printf("  analytic vs RK4(dt=1e-4), c=%g : max|dx,dv| = %.3e\n", c, dmax);
        std::fprintf(f, "analytic_vs_rk4_small_dt,c=%g,%.6e\n", c, dmax);
    }
    // (b) цикл step() проти матриці переходу M^n (перевірка формул стійкості)
    for (Method m : kMethods) {
        const double dt = 0.01; const long n = 300;
        State s{BASE.x0, BASE.v0};
        const Mat2 M = stepMatrix(m, BASE, dt);
        for (long i = 0; i < n; ++i) s = {M.a * s.x + M.b * s.v, M.c * s.x + M.d * s.v};
        const State q = integrate(m, BASE, dt, n, [](long, const State&) {});
        std::printf("  loop vs matrix %s: |dx| = %.3e\n", name(m), std::fabs(s.x - q.x));
        std::fprintf(f, "loop_vs_matrix,%s,%.6e\n", name(m), std::fabs(s.x - q.x));
    }
    std::fclose(f);
}

// ===================================================== 1. траєкторії різних методів
static void writeTraj(double dt) {
    char fn[64]; std::snprintf(fn, sizeof fn, "traj_dt%.3f.csv", dt);
    FILE* f = openCsv(fn, "t,analytic,euler,symp,rk4");
    const long n = std::lround(BASE.T / dt);
    std::vector<double> X[3];
    for (int m = 0; m < 3; ++m) {
        X[m].reserve(n + 1);
        integrate(kMethods[m], BASE, dt, n, [&](long, const State& s) { X[m].push_back(s.x); });
    }
    for (long i = 0; i <= n; ++i)
        std::fprintf(f, "%.10g,%.10e,%.10e,%.10e,%.10e\n", i * dt, analytic(i * dt, BASE).x, X[0][i], X[1][i], X[2][i]);
    std::fclose(f);
}

static void exp1() {
    std::printf("== EXP1 trajectories\n");
    for (double dt : {0.01, 0.05, 0.25, 0.5, 0.625}) writeTraj(dt);
    FILE* f = openCsv("analytic_fine.csv", "t,x");
    for (int i = 0; i <= 2000; ++i) {
        const double t = BASE.T * i / 2000.0;
        std::fprintf(f, "%.10g,%.12e\n", t, analytic(t, BASE).x);
    }
    std::fclose(f);
}

// =============================== 2. залежність похибки та часу від dt, порядок точності
static void exp2() {
    std::printf("== EXP2 dt sweep\n");
    const std::vector<double> dts = {0.25, 0.2, 0.1, 0.05, 0.025, 0.02, 0.01, 0.005, 0.0025, 0.001, 0.0005};
    FILE* f = openCsv("sweep.csv", "method,dt,n,abs_max,abs_end,rel,time_s,time_std_s,E_end_ratio,E_max_ratio,rho");
    std::vector<Run> rows[3];
    for (double dt : dts)
        for (int m = 0; m < 3; ++m) {
            Run r = runCase(kMethods[m], BASE, dt, REPEATS);
            rows[m].push_back(r);
            std::fprintf(f, "%s,%.6g,%ld,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.8f\n", name(kMethods[m]), dt, r.n,
                         r.abs_max, r.abs_end, r.rel, r.time, r.time_std, r.E_end, r.E_max, r.rho);
            std::printf("  dt=%-7g %-6s rel=%.3e  t=%9.4f us  rho=%.5f\n", dt, name(kMethods[m]), r.rel, r.time * 1e6, r.rho);
        }
    std::fclose(f);

    // глобальний нахил (емпіричний порядок точності) та локальні нахили
    FILE* fo = openCsv("orders.csv", "method,fit_range,slope");
    FILE* fl = openCsv("orders_local.csv", "method,dt_a,dt_b,slope");
    for (int m = 0; m < 3; ++m) {
        std::vector<double> X, Y;
        for (size_t i = 0; i < dts.size(); ++i) {
            const bool in = (m == 2) ? (dts[i] >= 0.0005 && dts[i] <= 0.05) : (dts[i] <= 0.01);
            if (in) { X.push_back(dts[i]); Y.push_back(rows[m][i].rel); }
        }
        const double s = fitSlope(X, Y);
        std::fprintf(fo, "%s,%s,%.4f\n", name(kMethods[m]), m == 2 ? "5e-4..5e-2" : "5e-4..1e-2", s);
        std::printf("  slope %s = %.4f\n", name(kMethods[m]), s);
        // dts спадають -> сусідні пари в порядку зростання dt
        for (int i = (int)dts.size() - 1; i > 0; --i) {
            const double a = dts[i], b = dts[i - 1];
            const double sl = std::log(rows[m][i - 1].rel / rows[m][i].rel) / std::log(b / a);
            std::fprintf(fl, "%s,%g,%g,%.4f\n", name(kMethods[m]), a, b, sl);
        }
    }
    std::fclose(fo); std::fclose(fl);
}

// ================================================================= 3. стійкість
static void exp3() {
    std::printf("== EXP3 stability\n");
    // (а) спектральний радіус як функція dt
    FILE* f = openCsv("rho_curve.csv", "dt,euler,symp,rk4");
    for (int i = 0; i < 700; ++i) {
        const double dt = 0.001 + (0.7 - 0.001) * i / 699.0;
        std::fprintf(f, "%.8g,%.8f,%.8f,%.8f\n", dt, rho(Method::Euler, BASE, dt), rho(Method::Symplectic, BASE, dt), rho(Method::RK4, BASE, dt));
    }
    std::fclose(f);

    // (б) теорія проти експерименту: 3000 кроків, зростання енергії => нестійкість
    FILE* fc = openCsv("classify.csv", "dt,method,rho,unstable_theory,unstable_empirical,max_abs_x_10s");
    int agree = 0, total = 0;
    for (double dt : {0.005, 0.015, 0.019, 0.021, 0.03, 0.1, 0.3, 0.38, 0.41, 0.5, 0.54, 0.57, 0.62})
        for (Method m : kMethods) {
            double E0 = energy({BASE.x0, BASE.v0}, BASE), Elast = 0;
            bool fin = true;
            const State s = integrate(m, BASE, dt, 3000, [&](long, const State& st) {
                Elast = energy(st, BASE);
                if (!std::isfinite(Elast)) fin = false;
            });
            (void)s;
            const bool emp = !fin || Elast > E0;
            const double r = rho(m, BASE, dt);
            const bool th = r > 1.0 + 1e-9;
            // амплітуда на інтервалі T = 10 с (для таблиці)
            double mx = 0;
            const long n10 = std::lround(BASE.T / dt);
            integrate(m, BASE, dt, n10, [&](long, const State& st) { mx = std::max(mx, std::fabs(st.x)); });
            agree += (emp == th); ++total;
            std::fprintf(fc, "%g,%s,%.8f,%d,%d,%.6e\n", dt, name(m), r, th, emp, mx);
        }
    std::fclose(fc);
    std::printf("  theory vs experiment agree: %d / %d\n", agree, total);

    // (в) критичний крок для базових параметрів
    FILE* fd = openCsv("dtmax_base.csv", "method,dt_max");
    for (Method m : kMethods) {
        const double d = dtMax(m, BASE);
        std::fprintf(fd, "%s,%.8f\n", name(m), d);
        std::printf("  dt_max(%s) = %.6f\n", name(m), d);
    }
    std::fclose(fd);

    // (г) dt_max залежно від коефіцієнта опору c
    FILE* fv = openCsv("dtmax_vs_c.csv", "c,euler,symp,rk4");
    for (double c : {0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 14.0, 20.0, 30.0}) {
        Params p = BASE; p.c = c;
        std::fprintf(fv, "%g,%.8f,%.8f,%.8f\n", c, dtMax(Method::Euler, p), dtMax(Method::Symplectic, p), dtMax(Method::RK4, p));
    }
    std::fclose(fv);
}

// ======================================================= 4. вплив фізичних параметрів
static void exp4() {
    std::printf("== EXP4 parameters (dt = 0.01)\n");
    const double dt = 0.01;
    FILE* f = openCsv("params.csv", "param,value,w0,w0dt,gamma,c_crit,regime,method,rel,abs_end,time_s,E_end_ratio,dt_max");
    struct Set { const char* name; std::vector<double> vals; };
    const std::vector<Set> sets = {{"m", {0.25, 0.5, 1.0, 2.0, 4.0}},
                                   {"k", {6.25, 25.0, 100.0, 400.0}},
                                   {"c", {0.0, 0.5, 2.0, 5.0, 10.0, 20.0}}};
    for (const auto& s : sets)
        for (double v : s.vals) {
            Params p = BASE;
            if (s.name[0] == 'm') p.m = v; else if (s.name[0] == 'k') p.k = v; else p.c = v;
            const double w0 = std::sqrt(p.k / p.m), ccr = 2 * std::sqrt(p.k * p.m);
            const char* reg = p.c < ccr - 1e-9 ? "oscillatory" : (std::fabs(p.c - ccr) < 1e-9 ? "critical" : "overdamped");
            std::printf("  %s=%g (w0*dt=%.3f, %s):", s.name, v, w0 * dt, reg);
            for (Method m : kMethods) {
                const Run r = runCase(m, p, dt, 5);
                std::fprintf(f, "%s,%g,%.6f,%.6f,%.6f,%.6f,%s,%s,%.6e,%.6e,%.6e,%.6e,%.8f\n", s.name, v, w0, w0 * dt,
                             p.c / (2 * p.m), ccr, reg, name(m), r.rel, r.abs_end, r.time, r.E_end, dtMax(m, p));
                std::printf("  %s %.2e", name(m), r.rel);
            }
            std::printf("\n");
        }
    std::fclose(f);
}

// ====================================== 5. енергія, фазовий портрет, монотонність
static void exp5() {
    std::printf("== EXP5 energy (c = 0) and phase portrait\n");
    Params p = BASE; p.c = 0.0; p.T = 100.0;
    const double E0 = energy({p.x0, p.v0}, p);
    FILE* fs = openCsv("undamped_summary.csv", "case,method,E_end_ratio,max_drift");

    // дрейф енергії, dt = 0.01 (усі три методи)
    {
        const double dt = 0.01; const long n = std::lround(p.T / dt);
        std::vector<double> D[3];
        for (int m = 0; m < 3; ++m) {
            double mx = 0, Eend = 0;
            integrate(kMethods[m], p, dt, n, [&](long i, const State& s) {
                const double E = energy(s, p), d = std::fabs(E / E0 - 1.0);
                if (i % 10 == 0) D[m].push_back(d);
                mx = std::max(mx, d); Eend = E;
            });
            std::fprintf(fs, "dt=0.01,%s,%.8e,%.8e\n", name(kMethods[m]), Eend / E0, mx);
            std::printf("  dt=0.01 %-6s E_end/E0=%.6g max|dE/E0|=%.4e\n", name(kMethods[m]), Eend / E0, mx);
        }
        FILE* f = openCsv("energy_dt0.010.csv", "t,euler,symp,rk4");
        for (size_t i = 0; i < D[0].size(); ++i) std::fprintf(f, "%.6g,%.6e,%.6e,%.6e\n", i * 10 * dt, D[0][i], D[1][i], D[2][i]);
        std::fclose(f);
    }
    // dt = 0.05: симплектичний та РК4
    {
        const double dt = 0.05; const long n = std::lround(p.T / dt);
        std::vector<double> D[2];
        for (int j = 0; j < 2; ++j) {
            const Method m = j == 0 ? Method::Symplectic : Method::RK4;
            double mx = 0, Eend = 0;
            integrate(m, p, dt, n, [&](long, const State& s) {
                const double E = energy(s, p), d = std::fabs(E / E0 - 1.0);
                D[j].push_back(d); mx = std::max(mx, d); Eend = E;
            });
            std::fprintf(fs, "dt=0.05,%s,%.8e,%.8e\n", name(m), Eend / E0, mx);
            std::printf("  dt=0.05 %-6s E_end/E0=%.6g max|dE/E0|=%.4e\n", name(m), Eend / E0, mx);
        }
        FILE* f = openCsv("energy_dt0.050.csv", "t,symp,rk4");
        for (size_t i = 0; i < D[0].size(); ++i) std::fprintf(f, "%.6g,%.6e,%.6e\n", i * dt, D[0][i], D[1][i]);
        std::fclose(f);
    }
    std::fclose(fs);

    // фазовий портрет (x, v/w0), dt = 0.05, T = 6
    {
        Params q = BASE; q.c = 0.0; q.T = 6.0;
        const double dt = 0.05, w0 = std::sqrt(q.k / q.m); const long n = std::lround(q.T / dt);
        FILE* f = openCsv("phase.csv", "method,t,x,v_over_w0");
        FILE* fr = openCsv("phase_radius.csv", "method,min_radius,max_radius");
        for (Method m : kMethods) {
            double rmin = 1e300, rmax = 0;
            integrate(m, q, dt, n, [&](long i, const State& s) {
                const double r = std::hypot(s.x, s.v / w0);
                rmin = std::min(rmin, r); rmax = std::max(rmax, r);
                std::fprintf(f, "%s,%.6g,%.8e,%.8e\n", name(m), i * dt, s.x, s.v / w0);
            });
            std::fprintf(fr, "%s,%.8f,%.8f\n", name(m), rmin, rmax);
        }
        std::fclose(f); std::fclose(fr);
    }
    // затухаюча система (базові параметри): частка кроків, на яких E зростає
    FILE* fm = openCsv("damped_monotonic.csv", "dt,method,frac_increase,E_end_ratio");
    for (double dt : {0.001, 0.01, 0.019, 0.025, 0.1, 0.25})
        for (Method m : kMethods) {
            const long n = std::lround(BASE.T / dt);
            const double E0b = energy({BASE.x0, BASE.v0}, BASE);
            double Eprev = E0b, Eend = E0b; long inc = 0;
            integrate(m, BASE, dt, n, [&](long i, const State& s) {
                const double E = energy(s, BASE);
                if (i > 0 && E - Eprev > 1e-12 * E0b) ++inc;
                Eprev = E; Eend = E;
            });
            std::fprintf(fm, "%g,%s,%.6f,%.8e\n", dt, name(m), (double)inc / n, Eend / E0b);
        }
    std::fclose(fm);
}

int main(int argc, char** argv) {
    if (argc > 1) OUT = argv[1];
    std::printf("Lab5 / variant 34: damped oscillator m=%g k=%g c=%g x0=%g v0=%g T=%g\n",
                BASE.m, BASE.k, BASE.c, BASE.x0, BASE.v0, BASE.T);
    exp0(); exp1(); exp2(); exp3(); exp4(); exp5();
    std::printf("DONE\n");
    return 0;
}
