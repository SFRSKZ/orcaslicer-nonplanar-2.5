// OrcaSlicer Nonplanar: G-code clean-up after the nonplanar output (port of np_postfix.py).
#include "NonplanarPostfix.hpp"

#include "../LocalesUtils.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

namespace Slic3r {

namespace {

constexpr double CELL = 0.2; // grid cell (mm)
constexpr double PI_  = 3.14159265358979323846;

// One parsed G-code line: the command word and its letter parameters.
struct Words
{
    std::string g;
    std::array<double, 26> v{};
    uint32_t has = 0;

    bool   is(char c) const { return (has >> (c - 'A')) & 1u; }
    double get(char c, double def) const { return is(c) ? v[c - 'A'] : def; }
    // true when every given parameter is one of `allowed`
    bool only(const char *allowed) const
    {
        for (int k = 0; k < 26; ++k)
            if (((has >> k) & 1u) && std::strchr(allowed, 'A' + k) == nullptr)
                return false;
        return true;
    }
    bool any_of(const char *letters) const
    {
        for (const char *p = letters; *p; ++p)
            if (is(*p)) return true;
        return false;
    }
};

bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r'; }

// Number in the form [-+]?\d*\.?\d+ starting at s[i]; returns its length (0 when there is none).
size_t number_len(const std::string &s, size_t i)
{
    size_t j = i;
    if (j < s.size() && (s[j] == '-' || s[j] == '+')) ++j;
    size_t digits_before = 0;
    while (j < s.size() && std::isdigit((unsigned char)s[j])) { ++j; ++digits_before; }
    size_t k = j;
    if (k < s.size() && s[k] == '.') {
        ++k;
        size_t digits_after = 0;
        while (k < s.size() && std::isdigit((unsigned char)s[k])) { ++k; ++digits_after; }
        if (digits_after > 0) return k - i;
    }
    return digits_before > 0 ? j - i : 0;
}

Words parse(const std::string &line)
{
    Words w;
    const size_t semi = line.find(';');
    const size_t end  = semi == std::string::npos ? line.size() : semi;
    size_t i = 0;
    while (i < end && is_space(line[i])) ++i;
    size_t j = i;
    while (j < end && !is_space(line[j])) ++j;
    if (i == j) return w;
    w.g.assign(line, i, j - i);
    const std::string rest = line.substr(j, end - j);
    for (size_t k = 0; k < rest.size(); ++k) {
        const char c = rest[k];
        if (c < 'A' || c > 'Z') continue;
        const size_t n = number_len(rest, k + 1);
        if (n == 0) continue;
        w.v[c - 'A'] = std::strtod(rest.substr(k + 1, n).c_str(), nullptr);
        w.has |= 1u << (c - 'A');
        k += n;
    }
    return w;
}

std::string trim(const std::string &s)
{
    size_t a = 0, b = s.size();
    while (a < b && (is_space(s[a]) || s[a] == '\n')) ++a;
    while (b > a && (is_space(s[b - 1]) || s[b - 1] == '\n')) --b;
    return s.substr(a, b - a);
}

bool starts_with(const std::string &s, const char *p) { return s.rfind(p, 0) == 0; }

bool is_move(const std::string &g) { return g == "G0" || g == "G1" || g == "G2" || g == "G3"; }
bool is_passive(const std::string &g) { return g == "M204" || g == "M73" || g == "M106" || g == "M107" || g == "M1002"; }

double arc_len(double x0, double y0, double x1, double y1, double i, double j, bool cw)
{
    const double cx = x0 + i, cy = y0 + j, r = std::hypot(i, j);
    const double a0 = std::atan2(y0 - cy, x0 - cx), a1 = std::atan2(y1 - cy, x1 - cx);
    double da = cw ? (a0 - a1) : (a1 - a0);
    while (da <= 0.) da += 2. * PI_;
    return r * da;
}

// The line without its E parameter; empty when nothing but "G1 [F..]" is left.
std::string strip_e(const std::string &line)
{
    const size_t semi = line.find(';');
    std::string code    = semi == std::string::npos ? line : line.substr(0, semi);
    std::string comment = semi == std::string::npos ? std::string() : line.substr(semi + 1);
    std::string body;
    for (size_t k = 0; k < code.size(); ++k) {
        if (code[k] == 'E') {
            const size_t n = number_len(code, k + 1);
            if (n > 0) {
                while (!body.empty() && is_space(body.back())) body.pop_back();
                k += n;
                continue;
            }
        }
        body += code[k];
    }
    while (!body.empty() && is_space(body.back())) body.pop_back();
    const Words w = parse(body);
    if (w.g == "G1" && w.only("F") && trim(body).find_first_of("XYZIJ") == std::string::npos)
        return std::string();
    return comment.empty() && semi == std::string::npos ? body : body + " ;" + comment;
}

// The line as a non-extruding move: an arc becomes a straight travel to its end.
std::string no_e(const std::string &line)
{
    const Words w = parse(line);
    if (w.g == "G2" || w.g == "G3") {
        std::string out = "G1";
        char buf[48];
        for (char c : { 'X', 'Y', 'Z' })
            if (w.is(c)) { std::snprintf(buf, sizeof(buf), " %c%g", c, w.get(c, 0.)); out += buf; }
        return out;
    }
    return strip_e(line);
}

void compact(std::vector<std::string> &lines, const std::vector<bool> &alive)
{
    size_t k = 0;
    for (size_t i = 0; i < lines.size(); ++i)
        if (alive[i]) {
            if (k != i) lines[k] = std::move(lines[i]);
            ++k;
        }
    lines.resize(k);
}

} // namespace

NonplanarPostfix::NonplanarPostfix(double min_len, double clearance, int initial_tool)
    : m_min_len(min_len), m_clear(clearance), m_tool(initial_tool)
{
    if (initial_tool >= 0 && initial_tool < 64)
        m_seen_tools |= uint64_t(1) << initial_tool;
}

uint64_t NonplanarPostfix::cell_key(int cx, int cy)
{
    return (uint64_t(uint32_t(cx)) << 32) | uint64_t(uint32_t(cy));
}

std::string NonplanarPostfix::process_layer(const std::string &gcode)
{
    if (gcode.empty())
        return gcode;
    CNumericLocalesSetter locales_setter;

    std::vector<std::string> lines;
    size_t start = 0;
    while (start < gcode.size()) {
        size_t nl = gcode.find('\n', start);
        if (nl == std::string::npos) nl = gcode.size();
        lines.emplace_back(gcode, start, nl - start);
        start = nl + 1;
    }
    const bool trailing_nl = gcode.back() == '\n';

    stage_lifts(lines);
    if (m_min_len > 0.)
        stage_dots(lines);
    if (m_clear > 0.)
        stage_hops(lines);

    std::string out;
    out.reserve(gcode.size() + 256);
    for (size_t i = 0; i < lines.size(); ++i) {
        out += lines[i];
        if (i + 1 < lines.size() || trailing_nl) out += '\n';
    }
    return out;
}

// ------------------------------------------------------------------ 0. lifts that go nowhere
void NonplanarPostfix::stage_lifts(std::vector<std::string> &lines)
{
    std::vector<bool> alive(lines.size(), true);
    Pos &p = m_s0;
    for (size_t i = 0; i < lines.size();) {
        const Words w = parse(lines[i]);
        if (w.g == "G1" && w.is('Z') && w.only("ZF") && w.get('Z', 0.) > p.z + 0.05 && i + 2 < lines.size()) {
            const Words w2 = parse(lines[i + 1]);
            const Words w3 = parse(lines[i + 2]);
            if (w2.g == "G1" && !w2.is('E') && !w2.is('Z') && w3.g == "G1" && w3.is('Z') && w3.only("ZF")
                && std::hypot(w2.get('X', p.x) - p.x, w2.get('Y', p.y) - p.y) < 0.05
                && std::abs(w3.get('Z', 0.) - p.z) < 0.05) {
                alive[i] = alive[i + 1] = false;
                p.x = w2.get('X', p.x);
                p.y = w2.get('Y', p.y);
                ++m_lifts;
                i += 2; // the lower move is kept and handled next
                continue;
            }
        }
        if (is_move(w.g)) {
            p.x = w.get('X', p.x);
            p.y = w.get('Y', p.y);
            p.z = w.get('Z', p.z);
        }
        ++i;
    }
    compact(lines, alive);
}

// ------------------------------------------------------------------ 1. tiny dots
void NonplanarPostfix::stage_dots(std::vector<std::string> &lines)
{
    struct Piece { std::vector<size_t> lines; double L = 0.; int layer = 0; std::string feat; bool tower = false; };
    std::vector<Piece> pieces;
    Piece cur;
    bool  open = false;
    auto end_piece = [&]() {
        if (open) pieces.push_back(std::move(cur));
        cur = Piece();
        open = false;
    };

    Pos &p = m_s1;
    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string s = trim(lines[i]);
        if (s.empty() || s[0] == ';') {
            if (starts_with(s, "; WIPE_TOWER_START")) { m_s1_tower = true; end_piece(); }
            else if (starts_with(s, "; WIPE_TOWER_END")) { m_s1_tower = false; end_piece(); }
            else if (starts_with(s, ";LAYER_CHANGE") || starts_with(s, "; CHANGE_LAYER")) { ++m_s1_layer; end_piece(); }
            else if (starts_with(s, "; FEATURE:")) m_s1_feat = trim(s.substr(10));
            continue;
        }
        const Words w = parse(lines[i]);
        if (is_move(w.g)) {
            const double nx = w.get('X', p.x), ny = w.get('Y', p.y), e = w.get('E', 0.);
            double d = std::hypot(nx - p.x, ny - p.y);
            if ((w.g == "G2" || w.g == "G3") && e > 0.)
                d = arc_len(p.x, p.y, nx, ny, w.get('I', 0.), w.get('J', 0.), w.g == "G2");
            if (e > 0. && d > 1e-4) {
                if (!open) { open = true; cur.layer = m_s1_layer; cur.feat = m_s1_feat; cur.tower = m_s1_tower; }
                cur.lines.push_back(i);
                cur.L += d;
            } else if (e == 0. && d <= 1e-3 && std::abs(w.get('Z', p.z) - p.z) < 0.05) {
                // feed-rate line, zero-length split or tiny settle inside a path
            } else {
                end_piece(); // travel, retract, unretract or Z move
            }
            p.x = nx; p.y = ny; p.z = w.get('Z', p.z);
        } else if (!is_passive(w.g)) {
            end_piece();
        }
    }
    end_piece();

    std::vector<bool> alive(lines.size(), true);
    for (const Piece &pc : pieces) {
        if (pc.L < m_min_len && !pc.tower && pc.layer > 1
            && pc.feat.find("Prime") == std::string::npos && pc.feat.find("Custom") == std::string::npos) {
            for (size_t k : pc.lines) {
                lines[k] = no_e(lines[k]);
                if (lines[k].empty()) alive[k] = false;
            }
            ++m_dots;
        }
    }

    // an unretract whose piece now extrudes nothing: drop it and its retract, so nothing is primed in place
    const size_t n = lines.size();
    for (size_t i = 0; i < n; ++i) {
        if (!alive[i]) continue;
        const Words w = parse(lines[i]);
        if (!(w.g == "G1" && w.get('E', 0.) > 0. && !w.any_of("XYZ")))
            continue;
        const double unret = w.get('E', 0.);
        double back = 0.;
        bool   extr = false;
        std::vector<size_t> ret;
        size_t j = i + 1;
        for (; j < n; ++j) {
            if (!alive[j]) continue;
            const std::string t = trim(lines[j]);
            if (starts_with(t, "; WIPE_TOWER") || starts_with(t, ";LAYER_CHANGE") || starts_with(t, "; CHANGE_LAYER")) break;
            const Words ww = parse(lines[j]);
            if (ww.g == "G1" || ww.g == "G2" || ww.g == "G3") {
                const double e = ww.get('E', 0.);
                if (e > 0.) { extr = true; break; }
                if (e < 0.) {
                    back += -e;
                    ret.push_back(j);
                    if (back >= unret - 1e-6) break;
                }
            } else if (!ww.g.empty() && !is_passive(ww.g)) {
                break;
            }
        }
        if (!extr && !ret.empty() && std::abs(back - unret) < 1e-4) {
            ret.insert(ret.begin(), i);
            for (size_t k : ret) {
                lines[k] = strip_e(lines[k]);
                if (lines[k].empty()) alive[k] = false;
            }
            i = j;
        }
    }
    compact(lines, alive);
}

// ------------------------------------------------------------------ 2. hop over another filament
void NonplanarPostfix::mark(int tool, double x0, double y0, double x1, double y1, double z0, double z1)
{
    const int n = std::max(1, int(std::hypot(x1 - x0, y1 - y0) / (CELL * 0.5)));
    for (int k = 0; k <= n; ++k) {
        const double t = double(k) / n;
        const double px = x0 + t * (x1 - x0), py = y0 + t * (y1 - y0);
        const float  pz = float(z0 + t * (z1 - z0));
        const int cx = int(std::floor(px / CELL)), cy = int(std::floor(py / CELL));
        for (int di = -1; di <= 1; ++di)
            for (int dj = -1; dj <= 1; ++dj) {
                auto &cell = m_top[cell_key(cx + di, cy + dj)];
                auto it = std::find_if(cell.begin(), cell.end(), [tool](const std::pair<int, float> &e) { return e.first == tool; });
                if (it == cell.end()) cell.emplace_back(tool, pz);
                else if (pz > it->second) it->second = pz;
            }
    }
}

bool NonplanarPostfix::danger(int tool, double x0, double y0, double x1, double y1, double z, double &hi) const
{
    const double L = std::hypot(x1 - x0, y1 - y0);
    if (L < 1.0) return false; // short hops inside a path do not drag across anything
    const int n = std::max(1, int(L / (CELL * 0.5)));
    bool found = false;
    for (int k = 0; k <= n; ++k) {
        const double t = double(k) / n;
        if (std::min(t, 1. - t) * L < 0.3) continue; // ends sit next to material anyway
        auto it = m_top.find(cell_key(int(std::floor((x0 + t * (x1 - x0)) / CELL)), int(std::floor((y0 + t * (y1 - y0)) / CELL))));
        if (it == m_top.end()) continue;
        for (const auto &[tl, zz] : it->second)
            if (tl != tool && z - zz < m_clear && (!found || zz > hi)) { hi = zz; found = true; }
    }
    return found;
}

void NonplanarPostfix::stage_hops(std::vector<std::string> &lines)
{
    std::vector<std::string> res;
    res.reserve(lines.size() + 16);
    Pos &p = m_s2;
    for (std::string &line : lines) {
        const std::string s = trim(line);
        if (starts_with(s, "; WIPE_TOWER_START")) m_s2_tower = true;
        else if (starts_with(s, "; WIPE_TOWER_END")) m_s2_tower = false;
        else if (starts_with(s, "; FEATURE:")) m_s2_feat = trim(s.substr(10));
        const Words w = parse(line);
        if (w.g.size() >= 2 && w.g[0] == 'T' && std::all_of(w.g.begin() + 1, w.g.end(), [](char c) { return std::isdigit((unsigned char)c); })) {
            const int t = std::atoi(w.g.c_str() + 1);
            if (t < 64) { m_tool = t; m_seen_tools |= uint64_t(1) << t; }
        }
        if (is_move(w.g)) {
            const double nx = w.get('X', p.x), ny = w.get('Y', p.y), nz = w.get('Z', p.z), e = w.get('E', 0.);
            const bool moved = std::abs(nx - p.x) > 1e-4 || std::abs(ny - p.y) > 1e-4;
            const bool several_tools = (m_seen_tools & (m_seen_tools - 1)) != 0;
            double hi = 0.;
            if (moved && e > 0. && m_tool >= 0 && !m_s2_tower && m_s2_feat.find("Prime") == std::string::npos) {
                mark(m_tool, p.x, p.y, nx, ny, p.z, nz);
            } else if (moved && e == 0. && !w.is('Z') && w.g == "G1" && m_tool >= 0 && !m_s2_tower && several_tools
                       && danger(m_tool, p.x, p.y, nx, ny, p.z, hi)) {
                char buf[96];
                std::snprintf(buf, sizeof(buf), "G1 Z%.3f F3000 ; np_postfix: hop over other filament", hi + m_clear);
                res.emplace_back(buf);
                res.push_back(std::move(line));
                std::snprintf(buf, sizeof(buf), "G1 Z%.3f F3000 ; np_postfix: back down", p.z);
                res.emplace_back(buf);
                ++m_hops;
                p.x = nx; p.y = ny;
                continue;
            }
            p.x = nx; p.y = ny; p.z = nz;
        }
        res.push_back(std::move(line));
    }
    lines = std::move(res);
}

std::string NonplanarPostfix::summary() const
{
    char buf[160];
    std::snprintf(buf, sizeof(buf), "; np_postfix: removed %zu empty lifts, dropped %zu dots, %zu hops\n", m_lifts, m_dots, m_hops);
    return buf;
}

} // namespace Slic3r
