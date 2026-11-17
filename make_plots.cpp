// make_plots.cpp - будує всі рисунки з CSV-файлів, створених experiments.cpp 
#define _USE_MATH_DEFINES
#include <cmath>
#include "svgplot.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using Row = std::map<std::string, std::string>;
static std::string IN = "data", OUT = "figs_svg";

static std::vector<Row> readCsv(const std::string& name) {
    std::ifstream f(IN + "/" + name);
    if (!f) { std::fprintf(stderr, "cannot open %s/%s\n", IN.c_str(), name.c_str()); std::exit(1); }
    std::string line; std::getline(f, line);
    std::vector<std::string> hdr; { std::stringstream ss(line); std::string c; while (std::getline(ss, c, ',')) hdr.push_back(c); }
    std::vector<Row> rows;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line); std::string c; Row r; size_t i = 0;
        while (std::getline(ss, c, ',') && i < hdr.size()) r[hdr[i++]] = c;
        rows.push_back(r);
    }
    return rows;
}
static std::vector<double> col(const std::vector<Row>& rows, const std::string& k, double scale = 1.0) {
    std::vector<double> v; for (auto& r : rows) v.push_back(std::atof(r.at(k).c_str()) * scale); return v;
}
static std::vector<Row> where(const std::vector<Row>& rows, const std::string& k, const std::string& v) {
    std::vector<Row> o; for (auto& r : rows) if (r.at(k) == v) o.push_back(r); return o;
}
static std::string uk(double v, int prec = 3) {   // число з десятковою комою
    char b[32]; std::snprintf(b, sizeof b, "%.*g", prec, v);
    std::string s = b; for (auto& c : s) if (c == '.') c = ','; return s;
}

static std::string fix2(double v) {   // рівно 2 знаки після коми
    char b[32]; std::snprintf(b, sizeof b, "%.2f", v);
    std::string s = b; for (auto& c : s) if (c == '.') c = ','; return s;
}

static const char* M[3] = {"euler", "symp", "rk4"};
static const char* NAME[3] = {"Явний Ейлер", "Симплектичний Ейлер", "Рунге–Кутта 4"};
static const char* COL[3] = {"#D55E00", "#009E73", "#0072B2"};
static const char MRK[3] = {'s', '^', 'o'};

static svg::Series line(const std::vector<double>& x, const std::vector<double>& y, int m, const std::string& label, bool marker) {
    svg::Series s; s.x = x; s.y = y; s.color = COL[m]; s.label = label; s.marker = marker ? MRK[m] : 0; s.lw = marker ? 1.5 : 1.3; return s;
}
static svg::Series ref(const std::vector<double>& x, const std::vector<double>& y, const std::string& label = "") {
    svg::Series s; s.x = x; s.y = y; s.color = "#000000"; s.dash = '-'; s.lw = 1.8; s.label = label; return s;
}
static void save(svg::Figure& f, const std::string& name) { f.write(OUT + "/" + name + ".svg"); std::printf("  %s.svg\n", name.c_str()); }

int main(int argc, char** argv) {
    if (argc > 1) IN = argv[1];
    if (argc > 2) OUT = argv[2];
    std::printf("make_plots: %s -> %s\n", IN.c_str(), OUT.c_str());

    // ---------------- Рис.1: траєкторії
    {
        svg::Figure fig(800, 270, 3);
        const auto fine = readCsv("analytic_fine.csv");
        const double dts[3] = {0.01, 0.05, 0.25};
        const char* fn[3] = {"traj_dt0.010.csv", "traj_dt0.050.csv", "traj_dt0.250.csv"};
        for (int i = 0; i < 3; ++i) {
            svg::Panel p; p.title = "dt = " + uk(dts[i]) + " с (ω₀dt = " + uk(5 * dts[i]) + ")"; p.xlabel = "t, с";
            if (i == 0) p.ylabel = "x, м";
            p.ymin = -1.6; p.ymax = (i == 0) ? 2.2 : 1.6; p.xmin = -0.3; p.xmax = 10.3;
            if (i == 0) p.legend = "tr";
            p.s.push_back(ref(col(fine, "t"), col(fine, "x"), i == 0 ? "Аналітичний" : ""));
            const auto tr = readCsv(fn[i]);
            for (int m = 0; m < 3; ++m) p.s.push_back(line(col(tr, "t"), col(tr, M[m]), m, i == 0 ? NAME[m] : "", false));
            fig.panels.push_back(p);
        }
        save(fig, "fig1_trajectories");
    }
    // ---------------- Рис.2 і 3: похибка, час, «ціна–точність»
    {
        const auto sw = readCsv("sweep.csv"); const auto ord = readCsv("orders.csv");
        svg::Figure fig(640, 300, 2);
        svg::Panel a, b;
        a.title = "а) Похибка від dt"; a.xlabel = "крок інтегрування dt, с"; a.ylabel = "відносна похибка (макс. за [0, T])";
        b.title = "б) Час обчислень від dt"; b.xlabel = "крок інтегрування dt, с"; b.ylabel = "середній час прогону, мкс";
        a.logx = a.logy = b.logx = b.logy = true; a.ymax = 1e12; a.invx = b.invx = true; a.legend = "tr"; b.legend = "tl";
        for (int m = 0; m < 3; ++m) {
            const auto r = where(sw, "method", M[m]);
            const double sl = std::atof(where(ord, "method", M[m]).at(0).at("slope").c_str());
            a.s.push_back(line(col(r, "dt"), col(r, "rel"), m, std::string(NAME[m]) + ", p = " + fix2(sl), true));
            b.s.push_back(line(col(r, "dt"), col(r, "time_s", 1e6), m, NAME[m], true));
        }
        svg::Series r1; r1.x = {5e-4, 5e-2}; r1.y = {1.5 * 5e-4, 1.5 * 5e-2}; r1.color = "#888888"; r1.dash = ':'; r1.lw = 1;
        svg::Series r4; r4.x = {5e-4, 5e-2}; r4.y = {160 * std::pow(5e-4, 4), 160 * std::pow(5e-2, 4)}; r4.color = "#888888"; r4.dash = ':'; r4.lw = 1;
        a.s.push_back(r1); a.s.push_back(r4);
        a.notes.push_back({3e-2, 2e-3, "∝ dt", "#666666", 10}); a.notes.push_back({1.2e-2, 3e-5, "∝ dt⁴", "#666666", 10});
        fig.panels = {a, b};
        save(fig, "fig2_error_time_vs_dt");

        svg::Figure f3(480, 340, 1);
        svg::Panel c; c.xlabel = "середній час прогону, мкс"; c.ylabel = "відносна похибка"; c.logx = c.logy = true; c.legend = "tr";
        for (int m = 0; m < 3; ++m) {
            std::vector<double> x, y;
            for (auto& r : where(sw, "method", M[m])) if (std::atof(r.at("rho").c_str()) <= 1.0) { x.push_back(std::atof(r.at("time_s").c_str()) * 1e6); y.push_back(std::atof(r.at("rel").c_str())); }
            c.s.push_back(line(x, y, m, NAME[m], true));
        }
        f3.panels = {c};
        save(f3, "fig3_error_vs_time");
    }
    // ---------------- Рис.4: стійкість, |x(t)| при різних dt
    {
        svg::Figure fig(800, 440, 3);
        const double dts[5] = {0.01, 0.05, 0.25, 0.5, 0.625};
        const char* fn[5] = {"traj_dt0.010.csv", "traj_dt0.050.csv", "traj_dt0.250.csv", "traj_dt0.500.csv", "traj_dt0.625.csv"};
        for (int i = 0; i < 5; ++i) {
            svg::Panel p; p.title = "dt = " + uk(dts[i]) + " с"; p.xlabel = "t, с"; p.logy = true; p.ymin = 1e-4; p.ymax = 1e8; p.floorY = 1e-12;
            p.xmin = -0.4; p.xmax = 10.4;
            if (i % 3 == 0) p.ylabel = "|x|, м (лог. шкала)";
            if (i == 0) p.legend = "tl";
            const auto tr = readCsv(fn[i]); auto t = col(tr, "t");
            auto absv = [](std::vector<double> v) { for (auto& x : v) x = std::fabs(x); return v; };
            p.s.push_back(ref(t, absv(col(tr, "analytic")), i == 0 ? "Аналітичний" : ""));
            for (int m = 0; m < 3; ++m) p.s.push_back(line(t, absv(col(tr, M[m])), m, i == 0 ? NAME[m] : "", false));
            fig.panels.push_back(p);
        }
        save(fig, "fig4_stability_dt");
    }
    // ---------------- Рис.5: ρ(dt) та dt_max(c)
    {
        svg::Figure fig(700, 300, 2);
        svg::Panel a, b;
        a.title = "а) ρ(dt), c = 0,5 кг/с"; a.xlabel = "dt, с"; a.ylabel = "спектральний радіус ρ матриці кроку"; a.ymin = 0.4; a.ymax = 1.6; a.legend = "bl";
        b.title = "б) dt_max залежно від c"; b.xlabel = "коефіцієнт опору c, кг/с"; b.ylabel = "максимальний стійкий крок dt_max, с"; b.logy = true; b.ymin = 5e-5; b.ymax = 1.2; b.floorY = 1e-4; b.legend = "cr";
        const auto rc = readCsv("rho_curve.csv"), dv = readCsv("dtmax_vs_c.csv");
        const char* SHORT[3] = {"Ейлер", "Симплектичний", "РК4"};
        for (int m = 0; m < 3; ++m) { a.s.push_back(line(col(rc, "dt"), col(rc, M[m]), m, SHORT[m], false)); b.s.push_back(line(col(dv, "c"), col(dv, M[m]), m, NAME[m], true)); }
        a.refs.push_back({false, 1.0, "#000000", '-'});
        b.refs.push_back({true, 10.0, "#888888", ':'});
        b.notes.push_back({10.4, 1.4e-4, "c_кр = 10", "#666666", 10});
        b.notes.push_back({1.5, 3e-4, "c = 0: Ейлер без стійкого кроку", COL[0], 9.5});
        fig.panels = {a, b};
        save(fig, "fig5_rho_dtmax");
    }
    // ---------------- Рис.6: вплив параметрів
    {
        svg::Figure fig(800, 280, 3);
        const auto pr = readCsv("params.csv");
        const char* names[3] = {"m", "k", "c"}; const char* lab[3] = {"маса m, кг", "жорсткість k, Н/м", "коеф. опору c, кг/с"};
        for (int i = 0; i < 3; ++i) {
            svg::Panel p; p.xlabel = lab[i]; p.logy = true; p.floorY = 1e-16;
            if (i == 0) { p.ylabel = "відносна похибка при dt = 0,01 с"; p.legend = "cr"; }
            for (int m = 0; m < 3; ++m) {
                const auto r = where(where(pr, "param", names[i]), "method", M[m]);
                std::vector<double> x; for (size_t k = 0; k < r.size(); ++k) x.push_back((double)k);
                p.s.push_back(line(x, col(r, "rel"), m, NAME[m], true));
                if (m == 0) for (size_t k = 0; k < r.size(); ++k) p.xticks.push_back({(double)k, uk(std::atof(r[k].at("value").c_str()))});
            }
            p.xmin = -0.3; p.xmax = (double)p.xticks.size() - 0.7;
            fig.panels.push_back(p);
        }
        save(fig, "fig6_params_error");
    }
    // ---------------- Рис.7: енергія (c = 0)
    {
        svg::Figure fig(700, 300, 2);
        svg::Panel a, b;
        a.title = "а) dt = 0,01 с"; b.title = "б) dt = 0,05 с (без Ейлера)";
        for (auto* p : {&a, &b}) { p->xlabel = "t, с"; p->ylabel = "|E(t)/E(0) − 1|"; p->logy = true; p->floorY = 1e-17; p->legend = "br"; }
        const auto e1 = readCsv("energy_dt0.010.csv"), e2 = readCsv("energy_dt0.050.csv");
        for (int m = 0; m < 3; ++m) a.s.push_back(line(col(e1, "t"), col(e1, M[m]), m, NAME[m], false));
        for (int m = 1; m < 3; ++m) b.s.push_back(line(col(e2, "t"), col(e2, M[m]), m, NAME[m], false));
        fig.panels = {a, b};
        save(fig, "fig7_energy");
    }
    // ---------------- Рис.8: фазові портрети (c = 0, dt = 0,05)
    {
        svg::Figure fig(800, 300, 3);
        const auto ph = readCsv("phase.csv");
        std::vector<double> cx, cy;
        for (int i = 0; i <= 400; ++i) { const double th = 2 * M_PI * i / 400; cx.push_back(std::cos(th)); cy.push_back(std::sin(th)); }
        for (int m = 0; m < 3; ++m) {
            svg::Panel p; p.title = NAME[m]; p.xlabel = "x, м"; p.equal = true;
            if (m == 0) { p.ylabel = "v / ω₀, м"; p.legend = "tl"; }
            p.s.push_back(ref(cx, cy, m == 0 ? "Аналітичний (коло)" : ""));
            const auto r = where(ph, "method", M[m]);
            p.s.push_back(line(col(r, "x"), col(r, "v_over_w0"), m, m == 0 ? NAME[m] : "", false));
            fig.panels.push_back(p);
        }
        save(fig, "fig8_phase");
    }
    std::printf("DONE\n");
    return 0;
}
