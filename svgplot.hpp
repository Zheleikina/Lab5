// svgplot.hpp - мінімальна бібліотека побудови графіків (вивід у SVG), без сторонніх залежностей.
// Підтримує: лінійні/логарифмічні осі, кілька панелей, маркери, пунктир, легенду, довільні підписи,
// довільні мітки на осі X (категорійні), однаковий масштаб осей (фазові портрети).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace svg {

struct Series {
    std::vector<double> x, y;
    std::string color = "#000000", label;
    char marker = 0;            // 'o', 's', '^' або 0
    char dash = 0;              // '-' пунктир, ':' крапки, 0 суцільна
    double lw = 1.6;
};
struct Note { double x, y; std::string text, color = "#555555"; double size = 11; };
struct Ref { bool vertical; double v; std::string color = "#000000"; char dash = '-'; };

struct Panel {
    std::string title, xlabel, ylabel;
    bool logx = false, logy = false, invx = false, equal = false;
    double xmin = NAN, xmax = NAN, ymin = NAN, ymax = NAN;
    double floorY = 0;          // для logy: значення нижче floorY підтягуються до floorY
    std::vector<Series> s;
    std::vector<Note> notes;
    std::vector<Ref> refs;
    std::string legend;         // "tl","tr","bl","br","cr" або ""
    std::vector<std::pair<double, std::string>> xticks;   // довільні мітки X
};

inline std::string esc(const std::string& t) {
    std::string o;
    for (char c : t) {
        if (c == '&') o += "&amp;"; else if (c == '<') o += "&lt;"; else if (c == '>') o += "&gt;"; else o += c;
    }
    return o;
}
// число з десятковою комою та справжнім мінусом
inline std::string numLabel(double v) {
    char b[64];
    std::snprintf(b, sizeof b, "%.6g", std::fabs(v) < 1e-12 ? 0.0 : v);
    std::string s = b;
    std::string o;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '.') o += ','; else if (s[i] == '-') o += "\xE2\x88\x92"; else o += s[i];
    }
    return o;
}

class Figure {
public:
    int W, H, cols;
    std::vector<Panel> panels;
    Figure(int w, int h, int c) : W(w), H(h), cols(c) {}

    void write(const std::string& path) const {
        FILE* f = std::fopen(path.c_str(), "w");
        if (!f) { std::perror(path.c_str()); return; }
        std::fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\" "
                        "font-family=\"DejaVu Sans, Arial, sans-serif\" font-size=\"11\">\n", W, H, W, H);
        std::fprintf(f, "<rect width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n", W, H);
        const int rows = ((int)panels.size() + cols - 1) / cols;
        const double cw = (double)W / cols, ch = (double)H / rows;
        for (size_t i = 0; i < panels.size(); ++i)
            drawPanel(f, panels[i], (int)i, (i % cols) * cw, (i / cols) * ch, cw, ch);
        std::fprintf(f, "</svg>\n");
        std::fclose(f);
    }

private:
    struct Axis {
        bool log; double lo, hi;
        double t(double v) const {   // 0..1
            if (log) return (std::log10(std::max(v, 1e-300)) - std::log10(lo)) / (std::log10(hi) - std::log10(lo));
            return (v - lo) / (hi - lo);
        }
    };

    static void autoRange(bool lg, double mn, double mx, double userLo, double userHi, double& lo, double& hi) {
        if (!(mx >= mn)) { mn = lg ? 1 : 0; mx = lg ? 10 : 1; }
        if (lg) { lo = mn / 1.3; hi = mx * 1.3; }
        else {
            double pad = (mx - mn) * 0.04; if (pad == 0) pad = 0.5;
            lo = mn - pad; hi = mx + pad;
        }
        if (!std::isnan(userLo)) lo = userLo;
        if (!std::isnan(userHi)) hi = userHi;
    }

    static std::vector<double> linTicks(double lo, double hi) {
        const double raw = (hi - lo) / 5.0, mag = std::pow(10.0, std::floor(std::log10(raw)));
        const double n = raw / mag, step = (n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10) * mag;
        std::vector<double> t;
        for (double v = std::ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) t.push_back(v);
        return t;
    }
    static std::vector<int> logDecades(double lo, double hi) {
        const int a = (int)std::ceil(std::log10(lo) - 1e-9), b = (int)std::floor(std::log10(hi) + 1e-9);
        const int step = std::max(1, (int)std::ceil((b - a + 1) / 7.0));
        std::vector<int> d;
        for (int k = a; k <= b; k += step) d.push_back(k);
        return d;
    }
    static std::string logLabel(int k) {
        char b[160];
        std::snprintf(b, sizeof b, "10<tspan dy=\"-5\" font-size=\"8.5\">%s</tspan>", numLabel(k).c_str());
        return b;
    }
    static const char* dashAttr(char d) { return d == '-' ? " stroke-dasharray=\"6 4\"" : d == ':' ? " stroke-dasharray=\"1.5 3\"" : ""; }

    void drawPanel(FILE* f, const Panel& p, int id, double ox, double oy, double cw, double ch) const {
        // межі даних
        double xmn = 1e300, xmx = -1e300, ymn = 1e300, ymx = -1e300;
        for (const auto& s : p.s)
            for (size_t i = 0; i < s.x.size() && i < s.y.size(); ++i) {
                const double x = s.x[i], y = s.y[i];
                if (!std::isfinite(x) || !std::isfinite(y)) continue;
                if (p.logx && x <= 0) continue;
                if (p.logy && y <= 0 && p.floorY <= 0) continue;
                const double yy = (p.logy && y < p.floorY) ? p.floorY : y;
                xmn = std::min(xmn, x); xmx = std::max(xmx, x); ymn = std::min(ymn, yy); ymx = std::max(ymx, yy);
            }
        Axis ax{p.logx, 0, 1}, ay{p.logy, 0, 1};
        autoRange(p.logx, xmn, xmx, p.xmin, p.xmax, ax.lo, ax.hi);
        autoRange(p.logy, ymn, ymx, p.ymin, p.ymax, ay.lo, ay.hi);

        // внутрішня область
        const double ml = 64, mr = 14, mt = 28, mb = 44;
        double ix = ox + ml, iy = oy + mt, iw = cw - ml - mr, ih = ch - mt - mb;
        if (p.equal) {
            const double sx = iw / (ax.hi - ax.lo), sy = ih / (ay.hi - ay.lo), sc = std::min(sx, sy);
            const double nw = sc * (ax.hi - ax.lo), nh = sc * (ay.hi - ay.lo);
            ix += (iw - nw) / 2; iy += (ih - nh) / 2; iw = nw; ih = nh;
        }
        auto X = [&](double v) { double t = ax.t(v); if (p.invx) t = 1 - t; return ix + std::clamp(t, -50.0, 50.0) * iw; };
        auto Y = [&](double v) { return iy + (1 - std::clamp(ay.t(v), -50.0, 50.0)) * ih; };

        std::fprintf(f, "<defs><clipPath id=\"c%d\"><rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\"/></clipPath></defs>\n",
                     id, ix, iy, iw, ih);
        // заголовок
        std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-weight=\"bold\" font-size=\"12.5\">%s</text>\n",
                     ix + iw / 2, oy + 17, esc(p.title).c_str());

        // сітка та мітки
        std::vector<std::pair<double, std::string>> xt, yt;
        if (!p.xticks.empty()) for (auto& t : p.xticks) xt.push_back({t.first, esc(t.second)});
        else if (p.logx) for (int k : logDecades(ax.lo, ax.hi)) xt.push_back({std::pow(10.0, k), logLabel(k)});
        else for (double v : linTicks(ax.lo, ax.hi)) xt.push_back({v, numLabel(v)});
        if (p.logy) for (int k : logDecades(ay.lo, ay.hi)) yt.push_back({std::pow(10.0, k), logLabel(k)});
        else for (double v : linTicks(ay.lo, ay.hi)) yt.push_back({v, numLabel(v)});
        for (auto& t : xt) {
            const double px = X(t.first);
            if (px < ix - 0.5 || px > ix + iw + 0.5) continue;
            std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#dddddd\" stroke-width=\"0.7\"/>\n", px, iy, px, iy + ih);
            std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" fill=\"#333\">%s</text>\n", px, iy + ih + 15, t.second.c_str());
        }
        for (auto& t : yt) {
            const double py = Y(t.first);
            if (py < iy - 0.5 || py > iy + ih + 0.5) continue;
            std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#dddddd\" stroke-width=\"0.7\"/>\n", ix, py, ix + iw, py);
            std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" fill=\"#333\">%s</text>\n", ix - 6, py + 4, t.second.c_str());
        }
        // осі (лише ліва та нижня)
        std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#222\"/>\n", ix, iy, ix, iy + ih);
        std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#222\"/>\n", ix, iy + ih, ix + iw, iy + ih);
        if (!p.xlabel.empty())
            std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\">%s</text>\n", ix + iw / 2, iy + ih + 34, esc(p.xlabel).c_str());
        if (!p.ylabel.empty())
            std::fprintf(f, "<text transform=\"translate(%.1f %.1f) rotate(-90)\" text-anchor=\"middle\">%s</text>\n", ox + 14, iy + ih / 2, esc(p.ylabel).c_str());

        // дані (з відсіканням по області графіка)
        std::fprintf(f, "<g clip-path=\"url(#c%d)\">\n", id);
        for (const auto& r : p.refs) {
            if (r.vertical) { const double px = X(r.v); std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"1\"%s/>\n", px, iy, px, iy + ih, r.color.c_str(), dashAttr(r.dash)); }
            else { const double py = Y(r.v); std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"1\"%s/>\n", ix, py, ix + iw, py, r.color.c_str(), dashAttr(r.dash)); }
        }
        for (const auto& s : p.s) {
            std::string pts; bool any = false;
            auto flush = [&]() {
                if (any) std::fprintf(f, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"%.2f\" stroke-linejoin=\"round\"%s points=\"%s\"/>\n", s.color.c_str(), s.lw, dashAttr(s.dash), pts.c_str());
                pts.clear(); any = false;
            };
            for (size_t i = 0; i < s.x.size() && i < s.y.size(); ++i) {
                double x = s.x[i], y = s.y[i];
                if (!std::isfinite(x) || !std::isfinite(y) || (p.logx && x <= 0)) { flush(); continue; }
                if (p.logy) { const double fl = p.floorY > 0 ? p.floorY : ay.lo / 10; if (y < fl) y = fl; }
                char b[48]; std::snprintf(b, sizeof b, "%.2f,%.2f ", X(x), Y(y)); pts += b; any = true;
            }
            flush();
            if (s.marker)
                for (size_t i = 0; i < s.x.size() && i < s.y.size(); ++i) {
                    double x = s.x[i], y = s.y[i];
                    if (!std::isfinite(x) || !std::isfinite(y) || (p.logx && x <= 0)) continue;
                    if (p.logy) { const double fl = p.floorY > 0 ? p.floorY : ay.lo / 10; if (y < fl) y = fl; }
                    marker(f, s.marker, X(x), Y(y), s.color);
                }
        }
        std::fprintf(f, "</g>\n");
        for (const auto& n : p.notes)
            std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" fill=\"%s\" font-size=\"%.1f\">%s</text>\n", X(n.x), Y(n.y), n.color.c_str(), n.size, esc(n.text).c_str());

        // легенда
        std::vector<const Series*> ls;
        for (const auto& s : p.s) if (!s.label.empty()) ls.push_back(&s);
        if (!p.legend.empty() && !ls.empty()) {
            const double lh = 16; const double bw = 22;
            double maxc = 0;
            for (auto* s : ls) {   // довжина в символах UTF-8 (не в байтах)
                double n = 0; for (unsigned char ch : s->label) if ((ch & 0xC0) != 0x80) n += 1;
                maxc = std::max(maxc, n);
            }
            const double w = bw + 8 + maxc * 5.6, h = lh * ls.size();
            double lx = ix + 8, ly = iy + 8;
            if (p.legend[1] == 'r') lx = ix + iw - w - 8;
            if (p.legend[0] == 'b') ly = iy + ih - h - 8;
            if (p.legend == "cr") { lx = ix + iw - w - 8; ly = iy + (ih - h) / 2; }
            for (size_t i = 0; i < ls.size(); ++i) {
                const double yy = ly + lh * i + 8;
                std::fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"%.2f\"%s/>\n", lx, yy, lx + bw, yy, ls[i]->color.c_str(), ls[i]->lw, dashAttr(ls[i]->dash));
                if (ls[i]->marker) marker(f, ls[i]->marker, lx + bw / 2, yy, ls[i]->color);
                std::fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" fill=\"#222\" font-size=\"10.5\">%s</text>\n", lx + bw + 6, yy + 3.6, esc(ls[i]->label).c_str());
            }
        }
    }

    static void marker(FILE* f, char m, double x, double y, const std::string& c) {
        if (m == 'o') std::fprintf(f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"3.3\" fill=\"%s\"/>\n", x, y, c.c_str());
        else if (m == 's') std::fprintf(f, "<rect x=\"%.2f\" y=\"%.2f\" width=\"6.4\" height=\"6.4\" fill=\"%s\"/>\n", x - 3.2, y - 3.2, c.c_str());
        else if (m == '^') std::fprintf(f, "<polygon points=\"%.2f,%.2f %.2f,%.2f %.2f,%.2f\" fill=\"%s\"/>\n", x, y - 4, x - 3.8, y + 3, x + 3.8, y + 3, c.c_str());
    }
};

}  // namespace svg
