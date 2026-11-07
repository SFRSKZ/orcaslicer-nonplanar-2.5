// OrcaSlicer Nonplanar: extrusion of a conformal path scaled to the real gap under it.
#include "NonplanarFlow.hpp"

#include <algorithm>
#include <cmath>

namespace Slic3r {

void NonplanarFlowMap::deposit(const Vec2d &a, double za, const Vec2d &b, double zb, double width, uint64_t shell)
{
    const double r = 0.5 * width;
    if (r <= 0.)
        return;
    const Vec2d  d    = b - a;
    const double len2 = d.squaredNorm();
    const int ix0 = int(std::floor((std::min(a.x(), b.x()) - r) / CELL)), ix1 = int(std::floor((std::max(a.x(), b.x()) + r) / CELL));
    const int iy0 = int(std::floor((std::min(a.y(), b.y()) - r) / CELL)), iy1 = int(std::floor((std::max(a.y(), b.y()) + r) / CELL));
    for (int ix = ix0; ix <= ix1; ++ix)
        for (int iy = iy0; iy <= iy1; ++iy) {
            const Vec2d  c((ix + 0.5) * CELL, (iy + 0.5) * CELL);
            const double t    = len2 > 0. ? std::clamp((c - a).dot(d) / len2, 0., 1.) : 0.;
            const double dist = (a + d * t - c).norm();
            if (dist > r)
                continue;
            // the bead top where it passes the cell (not its highest end: on a slope that is off by the slope)
            const float z = float(za + (zb - za) * t);
            Cell &cell = m_cells[key(ix, iy)];
            if (cell.top_shell == shell) {
                // the neighbour line of the same shell is uphill or downhill: keep the closest one
                if (float(dist) < cell.top_dist) { cell.top = z; cell.top_dist = float(dist); }
            } else if (z >= cell.top - 0.02f) {
                // a new shell on top: the old top is now the material under it
                cell.other     = cell.top;
                cell.top       = z;
                cell.top_dist  = float(dist);
                cell.top_shell = shell;
            } else {
                // lower material printed later (deferred): it lies under the top shell
                cell.other = std::max(cell.other, z);
            }
        }
}

bool NonplanarFlowMap::below(const Vec2d &p, uint64_t shell, double &z) const
{
    auto it = m_cells.find(key(int(std::floor(p.x() / CELL)), int(std::floor(p.y() / CELL))));
    if (it == m_cells.end())
        return false;
    const Cell &cell = it->second;
    const float v = cell.top_shell == shell ? cell.other : cell.top;
    if (v < -1e8f)
        return false;
    z = v;
    return true;
}

} // namespace Slic3r
