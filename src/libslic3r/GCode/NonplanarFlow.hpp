// OrcaSlicer Nonplanar: extrusion of a conformal path scaled to the real gap under it.
//
// A conformal shell is extruded with its nominal height, but where it lies on the
// steps of the planar layers below, the gap between the nozzle and the material
// under it can be much thinner (0.05 mm under a 0.12 mm shell): the same flow then
// over-extrudes there.  This map keeps the top of the material deposited so far on
// a 0.1 mm grid, under each bead's centre width only.  Each cell also remembers which
// "shell" (a planar layer, or a conformal surface at one offset) put the top there,
// so a conformal path reads the material below its own shell, not its neighbour lines.
#pragma once

#include <cmath>
#include <cstdint>
#include <unordered_map>

#include "../Point.hpp"

namespace Slic3r {

class NonplanarFlowMap
{
public:
    void clear() { m_cells.clear(); }

    // Segment a-b (G-code XY, bead top za / zb) of a bead of the given width, printed as part of `shell`.
    void deposit(const Vec2d &a, double za, const Vec2d &b, double zb, double width, uint64_t shell);

    // Top of the material under p that belongs to another shell; false when nothing is known.
    bool below(const Vec2d &p, uint64_t shell, double &z) const;

    // A planar layer is told apart by its height.
    static uint64_t planar_shell(double z) { return (uint64_t(1) << 62) | uint64_t(std::llround(z * 1000.)); }
    static uint64_t conformal_shell(int surface_id, float distance_to_top)
    {
        return (uint64_t(uint32_t(surface_id)) << 24) ^ uint64_t(uint32_t(std::lround(double(distance_to_top) * 1000.)) & 0xffffffu);
    }

private:
    struct Cell
    {
        float    top { -1e9f };      // top of the uppermost shell here, from its bead closest to the cell
        float    top_dist { 1e9f };  // distance of that bead's centre line
        uint64_t top_shell { ~uint64_t(0) };
        float    other { -1e9f };    // top of the material under top_shell
    };
    static constexpr double CELL = 0.1;
    static uint64_t key(int ix, int iy) { return (uint64_t(uint32_t(ix)) << 32) | uint64_t(uint32_t(iy)); }

    std::unordered_map<uint64_t, Cell> m_cells;
};

} // namespace Slic3r
