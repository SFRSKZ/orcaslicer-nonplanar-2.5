// Nonplanar ironing: height of the flat nozzle tip so that no part of the tip
// disk goes below the surface ("drop cutter" with a flat end mill).
//
// Header-only and independent of the rest of libslic3r.  `Facet` must provide
// `vertex[3]` with float members x, y, z, as Slic3r::NonplanarFacet does.
// NonplanarIroning.cpp calls it with the facets of one NonplanarSurface
// (mesh + rim) that lie under the tip disk.
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace Slic3r {
namespace NonplanarTip {

struct Result {
    bool   found     = false;
    double z         = -std::numeric_limits<double>::max(); // highest surface point under the disk
    double slope_deg = 0.;                                    // slope of the facet that gave it
};

namespace detail {

struct P2 { double x, y; };

// Plane of a triangle: z = a*x + b*y + c.  False for (near) vertical facets,
// which have no height over an area; their top edge is shared with a
// neighbour that does.
template<class Facet>
inline bool plane_of(const Facet &f, double &a, double &b, double &c)
{
    const double x0 = f.vertex[0].x, y0 = f.vertex[0].y, z0 = f.vertex[0].z;
    const double x1 = f.vertex[1].x, y1 = f.vertex[1].y, z1 = f.vertex[1].z;
    const double x2 = f.vertex[2].x, y2 = f.vertex[2].y, z2 = f.vertex[2].z;
    const double ux = x1 - x0, uy = y1 - y0, uz = z1 - z0;
    const double vx = x2 - x0, vy = y2 - y0, vz = z2 - z0;
    const double nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    const double len = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (len <= 0. || std::abs(nz) < 1e-6 * len)
        return false;
    a = -nx / nz; b = -ny / nz; c = z0 - a * x0 - b * y0;
    return true;
}

inline bool inside_triangle(const P2 t[3], double x, double y, double eps)
{
    auto cross = [](const P2 &p, const P2 &q, double x, double y) {
        return (q.x - p.x) * (y - p.y) - (q.y - p.y) * (x - p.x);
    };
    const double d0 = cross(t[0], t[1], x, y), d1 = cross(t[1], t[2], x, y), d2 = cross(t[2], t[0], x, y);
    const bool neg = d0 < -eps || d1 < -eps || d2 < -eps;
    const bool pos = d0 >  eps || d1 >  eps || d2 >  eps;
    return !(neg && pos);
}

} // namespace detail

// Highest point of one facet within the disk (cx, cy, r).  Exact: a plane's
// maximum over the convex set triangle ∩ disk lies at one of its extreme
// points -- a triangle vertex inside the disk, an edge/circle crossing, or the
// circle point in the plane's uphill direction when that point is inside the
// triangle.
template<class Facet>
inline bool facet_max_in_disk(const Facet &f, double cx, double cy, double r, double &zmax, double &slope_deg)
{
    using namespace detail;
    double a, b, c;
    if (!plane_of(f, a, b, c))
        return false;
    const P2 t[3] = { { f.vertex[0].x, f.vertex[0].y }, { f.vertex[1].x, f.vertex[1].y }, { f.vertex[2].x, f.vertex[2].y } };
    const double r2 = r * r;
    bool   any  = false;
    double best = -std::numeric_limits<double>::max();
    auto take = [&](double x, double y) { best = std::max(best, a * x + b * y + c); any = true; };

    for (const P2 &v : t)
        if ((v.x - cx) * (v.x - cx) + (v.y - cy) * (v.y - cy) <= r2)
            take(v.x, v.y);
    for (int i = 0; i < 3; ++i) {
        const P2 &p = t[i], &q = t[(i + 1) % 3];
        const double dx = q.x - p.x, dy = q.y - p.y;
        const double A = dx * dx + dy * dy;
        if (A <= 0.)
            continue;
        const double B = 2. * (dx * (p.x - cx) + dy * (p.y - cy));
        const double C = (p.x - cx) * (p.x - cx) + (p.y - cy) * (p.y - cy) - r2;
        const double disc = B * B - 4. * A * C;
        if (disc < 0.)
            continue;
        const double s = std::sqrt(disc);
        for (double u : { (-B - s) / (2. * A), (-B + s) / (2. * A) })
            if (u >= 0. && u <= 1.)
                take(p.x + u * dx, p.y + u * dy);
    }
    const double g = std::sqrt(a * a + b * b);
    const double ux = g > 1e-12 ? a / g : 0., uy = g > 1e-12 ? b / g : 0.;
    if (inside_triangle(t, cx + r * ux, cy + r * uy, 1e-12))
        take(cx + r * ux, cy + r * uy);
    // Disk entirely inside the triangle and the plane flat: the centre.
    if (!any && inside_triangle(t, cx, cy, 1e-12))
        take(cx, cy);
    if (!any)
        return false;
    zmax      = best;
    slope_deg = std::atan(g) * 180. / M_PI;
    return true;
}

// Tip height over the disk (cx, cy, tip_radius) for the given facets.  Facets
// whose XY bounding box misses the disk are skipped cheaply, so the caller may
// pass the candidate list of a coarse grid.  The result is the surface height
// the tip may touch; the caller adds the press offset (<= 0 presses down).
template<class FacetRange>
inline Result tip_height(const FacetRange &facets, double cx, double cy, double tip_radius)
{
    Result res;
    for (const auto &f : facets) {
        const double xmin = std::min({ double(f.vertex[0].x), double(f.vertex[1].x), double(f.vertex[2].x) });
        const double xmax = std::max({ double(f.vertex[0].x), double(f.vertex[1].x), double(f.vertex[2].x) });
        const double ymin = std::min({ double(f.vertex[0].y), double(f.vertex[1].y), double(f.vertex[2].y) });
        const double ymax = std::max({ double(f.vertex[0].y), double(f.vertex[1].y), double(f.vertex[2].y) });
        if (xmax < cx - tip_radius || xmin > cx + tip_radius || ymax < cy - tip_radius || ymin > cy + tip_radius)
            continue;
        double z, slope;
        if (facet_max_in_disk(f, cx, cy, tip_radius, z, slope) && z > res.z) {
            res.found = true; res.z = z; res.slope_deg = slope;
        }
    }
    return res;
}

} // namespace NonplanarTip
} // namespace Slic3r
