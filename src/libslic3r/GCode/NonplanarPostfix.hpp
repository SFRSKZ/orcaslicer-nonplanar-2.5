// OrcaSlicer Nonplanar: G-code clean-up after the nonplanar output, a port of np_postfix.py.
// Runs on each layer's G-code right before it is written (so the preview and the time
// estimate see the cleaned G-code). Needs relative E.
//
//  0. Empty lifts: lift + travel of < 0.05 mm + lower (no retract) between two pieces of one
//     path is removed; each landing would leave a dot.
//  1. Tiny dots: a piece of extrusion between two travels / retracts shorter than min_len mm
//     is not printed (its lines lose their E). An unretract that then has nothing to print
//     before its retract goes too, with that retract.
//  2. Hop over another filament: a travel passing over material of another filament with less
//     than clearance mm between the nozzle and that material is lifted to clearance above it.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Slic3r {

class NonplanarPostfix
{
public:
    // min_len <= 0 turns step 1 off, clearance <= 0 turns step 2 off.
    NonplanarPostfix(double min_len, double clearance, int initial_tool);

    std::string process_layer(const std::string &gcode);

    // One G-code comment line with the counters, e.g. for the end of the file.
    std::string summary() const;

    size_t lifts() const { return m_lifts; }
    size_t dots()  const { return m_dots; }
    size_t hops()  const { return m_hops; }

private:
    struct Pos { double x = 0., y = 0., z = 0.; };

    void stage_lifts(std::vector<std::string> &lines);
    void stage_dots(std::vector<std::string> &lines);
    void stage_hops(std::vector<std::string> &lines);

    void   mark(int tool, double x0, double y0, double x1, double y1, double z0, double z1);
    bool   danger(int tool, double x0, double y0, double x1, double y1, double z, double &hi) const;
    static uint64_t cell_key(int cx, int cy);

    double m_min_len;
    double m_clear;

    // per-stage machine state, carried over from layer to layer
    Pos  m_s0;
    Pos  m_s1;
    bool m_s1_tower = false;
    int  m_s1_layer = 0;
    std::string m_s1_feat;
    Pos  m_s2;
    bool m_s2_tower = false;
    std::string m_s2_feat;
    int  m_tool = -1;
    uint64_t m_seen_tools = 0;

    // 0.2 mm grid: cell -> (tool, highest z printed there)
    std::unordered_map<uint64_t, std::vector<std::pair<int, float>>> m_top;

    size_t m_lifts = 0, m_dots = 0, m_hops = 0;
};

} // namespace Slic3r
