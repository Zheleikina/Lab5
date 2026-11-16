// main.cpp - SFML-візуалізація (SFML 3). Фізика береться з physics.hpp, графіка лише малює стан.
// УВАГА: у середовищі, де готувалися тези, SFML відсутня, тому цей файл тут НЕ компілювався.
// Керування:  Up/Down - жорсткість k (+/-25)   Right/Left - опір c (+/-0.25)
//             ] / [   - крок dt (x2 / /2)       Space - пауза   R - перезапуск
// Три кулі (Ейлер, симплектичний Ейлер, РК4) і біла «аналітична» куля рухаються синхронно;
// внизу — графік x(t) за останні 10 с (sf::VertexArray + LineStrip, як траєкторія в методичці).
#include <SFML/Graphics.hpp>
#include <algorithm>
#include <optional>
#include <vector>
#include <string>
#include "physics.hpp"

int main() {
    using namespace phys;
    sf::RenderWindow window(sf::VideoMode({900, 800}), "SFML Lab5: numerical methods");
    window.setFramerateLimit(60);

    Params p;                       // m=1, k=25, c=0.5, x0=1
    double dt = 0.02, simT = 0.0, acc = 0.0;
    State s[3];
    std::vector<double> histT, histX[4];   // історія для графіка (по одній точці на кадр)
    auto reset = [&]() { for (auto& q : s) q = {p.x0, p.v0}; simT = 0.0; acc = 0.0; histT.clear(); for (auto& h : histX) h.clear(); };
    reset();
    bool paused = false;
    sf::Clock clock;
    const sf::Color col[3] = {sf::Color(213, 94, 0), sf::Color(0, 158, 115), sf::Color(0, 114, 178)};
    const float cx = 450.f, scale = 200.f, radius = 22.f;

    while (window.isOpen()) {
        // 1. Обробка введення
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) window.close();
            if (const auto* key = event->getIf<sf::Event::KeyPressed>()) {
                using K = sf::Keyboard::Key;
                // зміна k або c змінює модель -> перезапуск, щоб графік і точний розв'язок лишалися узгодженими
                if (key->code == K::Up) { p.k = std::clamp(p.k + 25.0, 1.0, 1000.0); reset(); }
                if (key->code == K::Down) { p.k = std::clamp(p.k - 25.0, 1.0, 1000.0); reset(); }
                if (key->code == K::Right) { p.c = std::clamp(p.c + 0.25, 0.0, 30.0); reset(); }
                if (key->code == K::Left) { p.c = std::clamp(p.c - 0.25, 0.0, 30.0); reset(); }
                if (key->code == K::RBracket) dt = std::min(dt * 2, 0.8);
                if (key->code == K::LBracket) dt = std::max(dt / 2, 0.0005);
                if (key->code == K::Space) paused = !paused;
                if (key->code == K::R) reset();
            }
        }
        // 2. Оновлення моделі: ФІКСОВАНИЙ крок dt, накопичувач реального часу
        acc += std::min(clock.restart().asSeconds(), 0.1f);
        while (!paused && acc >= dt) {
            for (int i = 0; i < 3; ++i) s[i] = step(kMethods[i], s[i], p, dt);
            simT += dt; acc -= dt;
        }
        const double xe = analytic(simT, p).x;
        if (!paused) {
            histT.push_back(simT);
            for (int i = 0; i < 3; ++i) histX[i].push_back(s[i].x);
            histX[3].push_back(xe);
            if (histT.size() > 2000) { histT.erase(histT.begin()); for (auto& h : histX) h.erase(h.begin()); }
        }
        window.setTitle("t=" + std::to_string(simT) + "  dt=" + std::to_string(dt) + "  k=" + std::to_string(p.k) +
                        "  c=" + std::to_string(p.c) + "  x_exact=" + std::to_string(xe));
        // 3. Відображення
        window.clear(sf::Color::Black);
        for (int i = 0; i < 4; ++i) {
            const float y = 120.f + 120.f * i;
            const double x = (i < 3) ? s[i].x : xe;
            const float px = cx + (float)std::clamp(x, -2.2, 2.2) * scale;   // обмеження, щоб розбіжність не «втекла» з вікна
            sf::RectangleShape rail({700.f, 2.f}); rail.setPosition({100.f, y}); rail.setFillColor(sf::Color(80, 80, 80));
            sf::CircleShape ball(radius); ball.setOrigin({radius, radius}); ball.setPosition({px, y});
            ball.setFillColor(i < 3 ? col[i] : sf::Color::White);
            window.draw(rail); window.draw(ball);
        }
        // графік x(t): область 700x220 px, вікно часу 10 с, масштаб 90 px/м
        const float gx = 100.f, gy = 560.f, gw = 700.f, gh = 220.f;
        sf::RectangleShape frame({gw, gh}); frame.setPosition({gx, gy});
        frame.setFillColor(sf::Color::Transparent); frame.setOutlineColor(sf::Color(120, 120, 120)); frame.setOutlineThickness(1.f);
        sf::RectangleShape zero({gw, 1.f}); zero.setPosition({gx, gy + gh / 2}); zero.setFillColor(sf::Color(70, 70, 70));
        window.draw(frame); window.draw(zero);
        const double t0 = std::max(0.0, simT - 10.0);
        for (int i = 3; i >= 0; --i) {          // спочатку «точна» крива, поверх неї методи
            sf::VertexArray line(sf::PrimitiveType::LineStrip);
            const sf::Color c = (i < 3) ? col[i] : sf::Color::White;
            for (size_t k = 0; k < histT.size(); ++k) {
                const float px = gx + (float)((histT[k] - t0) / 10.0) * gw;
                const float py = gy + gh / 2 - (float)std::clamp(histX[i][k], -1.2, 1.2) * 90.f;
                line.append(sf::Vertex{sf::Vector2f(px, py), c});
            }
            window.draw(line);
        }
        window.display();
    }
    return 0;
}
