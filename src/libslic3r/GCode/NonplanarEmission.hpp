// OrcaSlicer Nonplanar: ordering helpers for the emission of conformal shells (moved out of GCode.cpp).
#ifndef slic3r_GCode_NonplanarEmission_hpp_
#define slic3r_GCode_NonplanarEmission_hpp_

#include "../ExtrusionEntity.hpp"
#include "../ExtrusionEntityCollection.hpp"
#include "../NonplanarSchedule.hpp"
#include "../Polyline.hpp"
#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

namespace Slic3r {

// ---- Nonplanar top shells (ported from PrusaSlicer Nonplanar) ----
//
// Nonplanar surfaces are moved into a shared nominal layer.  The ordinary
// order is two-dimensional, so it may print a flat, higher part of the same
// nominal layer before a conformal perimeter which needs to pass beside it.
// Preserve normal ordering for purely planar work, but emit each conformal
// shell from its lower supporting shell towards its visible shell before the
// planar remainder.
struct NonplanarEmissionKey {
    bool  is_nonplanar = false;
    float shell_offset = -std::numeric_limits<float>::infinity();
};

inline void collect_nonplanar_emission_key(const ExtrusionEntity &entity, NonplanarEmissionKey &key)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath*>(&entity)) {
        if (path->nonplanar_source) {
            key.is_nonplanar = true;
            key.shell_offset = std::max(key.shell_offset, path->distance_to_top);
        }
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop*>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            collect_nonplanar_emission_key(path, key);
    } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath*>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            collect_nonplanar_emission_key(path, key);
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection*>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            collect_nonplanar_emission_key(*child, key);
    }
}

inline NonplanarEmissionKey nonplanar_emission_key(const ExtrusionEntity &entity)
{
    NonplanarEmissionKey key;
    collect_nonplanar_emission_key(entity, key);
    return key;
}

inline bool nonplanar_before_planar(const ExtrusionEntity &left, const ExtrusionEntity &right)
{
    const NonplanarEmissionKey a = nonplanar_emission_key(left);
    const NonplanarEmissionKey b = nonplanar_emission_key(right);
    if (a.is_nonplanar != b.is_nonplanar)
        return a.is_nonplanar;
    // Larger distance-to-top is the lower (supporting) conformal shell.
    if (a.is_nonplanar && std::abs(a.shell_offset - b.shell_offset) > EPSILON)
        return a.shell_offset > b.shell_offset;
    return false;
}

inline bool same_nonplanar_emission_group(const ExtrusionEntity &left, const ExtrusionEntity &right)
{
    const NonplanarEmissionKey a = nonplanar_emission_key(left);
    const NonplanarEmissionKey b = nonplanar_emission_key(right);
    return a.is_nonplanar == b.is_nonplanar &&
        (!a.is_nonplanar || std::abs(a.shell_offset - b.shell_offset) <= EPSILON);
}

inline bool has_nonplanar_source(const ExtrusionEntity &entity)
{
    return nonplanar_emission_key(entity).is_nonplanar;
}

// A conformal path carries absolute projected heights in Point3::z; planar
// material deferred to a later layer carries its own layer height there.
// (Z anti-aliasing offsets are excluded: ZAA and nonplanar layers are exclusive.)
inline bool path_has_absolute_z(const ExtrusionPath &path)
{
    return path.nonplanar_source ||
        (! path.z_contoured && ! path.polyline.points.empty() && path.polyline.points.front().z() > 0);
}

// Islands containing conformal paths get a deterministic order: supporting
// shells first, then the visible shell, then planar islands.  This is an
// island ordering, not a geometry filter.
struct NonplanarIslandKey {
    bool   has_nonplanar = false;
    float  shell_offset  = -std::numeric_limits<float>::infinity();
    size_t source_layer  = std::numeric_limits<size_t>::max();
};

inline void collect_nonplanar_island_key(const ExtrusionEntity &entity, NonplanarIslandKey &key)
{
    if (const auto *path = dynamic_cast<const ExtrusionPath*>(&entity)) {
        if (path->nonplanar_source) {
            key.has_nonplanar = true;
            key.shell_offset = std::max(key.shell_offset, path->distance_to_top);
            key.source_layer = std::min(key.source_layer, path->source_layer_id);
        }
    } else if (const auto *loop = dynamic_cast<const ExtrusionLoop*>(&entity)) {
        for (const ExtrusionPath &path : loop->paths)
            collect_nonplanar_island_key(path, key);
    } else if (const auto *multipath = dynamic_cast<const ExtrusionMultiPath*>(&entity)) {
        for (const ExtrusionPath &path : multipath->paths)
            collect_nonplanar_island_key(path, key);
    } else if (const auto *collection = dynamic_cast<const ExtrusionEntityCollection*>(&entity)) {
        for (const ExtrusionEntity *child : collection->entities)
            collect_nonplanar_island_key(*child, key);
    }
}

inline bool nonplanar_island_before(const NonplanarIslandKey &left, const NonplanarIslandKey &right)
{
    if (left.has_nonplanar != right.has_nonplanar)
        return left.has_nonplanar;
    if (!left.has_nonplanar)
        return false;
    // Larger distance-to-top is the lower supporting shell.  Source layer is
    // only a deterministic tie breaker; it does not change the geometry.
    if (std::abs(left.shell_offset - right.shell_offset) > EPSILON)
        return left.shell_offset > right.shell_offset;
    return left.source_layer < right.source_layer;
}

// XY arc length of a path (scaled), as measured by the deposition planner.
inline double polyline_length_scaled(const Polyline3 &polyline)
{
    double length = 0.;
    for (size_t i = 1; i < polyline.points.size(); ++ i)
        length += (polyline.points[i].to_point() - polyline.points[i-1].to_point()).cast<double>().norm();
    return length;
}

// Sub-polyline between XY arc lengths a and b; Z is interpolated.
inline Polyline3 sub_polyline(const Polyline3 &source, double a, double b)
{
    Polyline3 out;
    const Points3 &src = source.points;
    double walked = 0.;
    auto at = [](const Point3 &p, const Point3 &q, double t) {
        const Vec2d xy = p.to_point().cast<double>() + t * (q.to_point() - p.to_point()).cast<double>();
        return Point3(coord_t(std::round(xy.x())), coord_t(std::round(xy.y())),
                      coord_t(std::round(double(p.z()) + t * double(q.z() - p.z()))));
    };
    for (size_t i = 1; i < src.size(); ++ i) {
        const double len = (src[i].to_point() - src[i-1].to_point()).cast<double>().norm();
        const double s0 = walked, s1 = walked + len;
        walked = s1;
        if (s1 < a || len <= 0.)
            continue;
        if (s0 > b)
            break;
        if (out.points.empty())
            out.points.push_back(at(src[i-1], src[i], std::clamp((a - s0) / len, 0., 1.)));
        if (s1 <= b) {
            if (out.points.back() != src[i])
                out.points.push_back(src[i]);
        } else {
            out.points.push_back(at(src[i-1], src[i], (b - s0) / len));
            break;
        }
    }
    return out;
}

// Slot intervals of a path piece in its emitted direction: (slot, a, b) with
// arc lengths of the emitted polyline.  Slivers below 0.01 mm are dropped: a
// piece of almost no length would still cause a travel and descent.
inline std::vector<std::tuple<int, double, double>> plan_pieces(const NonplanarDeferral::PathPlan &plan, bool reversed, double length)
{
    std::vector<std::tuple<int, double, double>> out;
    const double scale = plan.length > 0. ? length / plan.length : 1.;
    for (size_t k = 0; k < plan.breaks.size(); ++ k) {
        const double a = plan.breaks[k].first * scale;
        const double b = (k + 1 < plan.breaks.size() ? plan.breaks[k + 1].first : plan.length) * scale;
        if (b - a <= scaled<double>(0.01))
            continue;
        out.emplace_back(plan.breaks[k].second, reversed ? length - b : a, reversed ? length - a : b);
    }
    if (reversed)
        std::reverse(out.begin(), out.end());
    return out;
}

} // namespace Slic3r

#endif // slic3r_GCode_NonplanarEmission_hpp_
