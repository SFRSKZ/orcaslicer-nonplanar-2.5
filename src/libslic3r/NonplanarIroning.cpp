// Nonplanar ironing: a very low flow pass over the visible conformal shell
// (distance_to_top == 0) that follows the surface.
//
// The flat nozzle tip must never go below the surface.  Centred on a slope its
// uphill edge would dig in by r * tan(slope) (0.12 mm at 12 deg with the A1 mini
// tip), so every point of a line is lifted to the highest surface point under
// the tip disk (NonplanarTipHeight.hpp, exact).  On a slope only the uphill edge
// touches; the lines run along the surface's contour so that edge keeps a
// constant contact, and parts steeper than nonplanar_ironing_max_angle are left
// out.  nonplanar_ironing_slope_press lowers the tip on a slope by a share of
// the gap under its middle (at most 0.05 mm) so the slopes are pressed too.
// The ironing area is the top shell inset by the tip radius, so the disk
// never reaches past the surface (a raised letter, a wall).
//
// The paths are conformal (nonplanar_source) and carry their absolute heights;
// they are added after the projection step and printed with the other ironing
// of the island, i.e. after the shell they iron.  The deposition planner sees
// them like any other conformal path.

#include "Layer.hpp"
#include "ClipperUtils.hpp"
#include "ExtrusionEntity.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Fill/FillBase.hpp"
#include "Flow.hpp"
#include "NonplanarSchedule.hpp"
#include "NonplanarSurface.hpp"
#include "NonplanarTipHeight.hpp"
#include "Print.hpp"

#include <boost/log/trivial.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>

namespace Slic3r {

namespace {

// Facets of one surface (mesh + rim) on a coarse grid, for the disk queries.
class TipFacetGrid {
public:
    explicit TipFacetGrid(const NonplanarSurface &surface)
    {
        for (const auto *facets : { &surface.mesh, &surface.rim })
            for (const auto &entry : *facets)
                add(&entry.second);
    }

    // Highest surface point under the disk (exact over the facets that reach it).
    NonplanarTip::Result tip_height(double cx, double cy, double r) const
    {
        NonplanarTip::Result res;
        for (int ix = cell(cx - r); ix <= cell(cx + r); ++ ix)
            for (int iy = cell(cy - r); iy <= cell(cy + r); ++ iy) {
                const auto it = m_cells.find(key(ix, iy));
                if (it == m_cells.end())
                    continue;
                for (const NonplanarFacet *f : it->second) {
                    double z, slope;
                    if (NonplanarTip::facet_max_in_disk(*f, cx, cy, r, z, slope) && z > res.z) {
                        res.found = true; res.z = z; res.slope_deg = slope;
                    }
                }
            }
        return res;
    }

    // Surface height at a single point (for the verification only).
    bool height_at(double x, double y, double &z) const
    {
        bool found = false;
        const auto it = m_cells.find(key(cell(x), cell(y)));
        if (it == m_cells.end())
            return false;
        for (const NonplanarFacet *f : it->second) {
            double zz, slope;
            if (NonplanarTip::facet_max_in_disk(*f, x, y, 1e-6, zz, slope) && (! found || zz > z)) {
                z = zz; found = true;
            }
        }
        return found;
    }

private:
    void add(const NonplanarFacet *f)
    {
        const auto &v = f->vertex;
        const double xmin = std::min({ double(v[0].x), double(v[1].x), double(v[2].x) });
        const double xmax = std::max({ double(v[0].x), double(v[1].x), double(v[2].x) });
        const double ymin = std::min({ double(v[0].y), double(v[1].y), double(v[2].y) });
        const double ymax = std::max({ double(v[0].y), double(v[1].y), double(v[2].y) });
        for (int ix = cell(xmin); ix <= cell(xmax); ++ ix)
            for (int iy = cell(ymin); iy <= cell(ymax); ++ iy)
                m_cells[key(ix, iy)].push_back(f);
    }
    static int cell(double v) { return int(std::floor(v / 2.0)); }
    static uint64_t key(int x, int y) { return (uint64_t(uint32_t(x)) << 32) | uint32_t(y); }
    std::unordered_map<uint64_t, std::vector<const NonplanarFacet*>> m_cells;
};

// Line direction along the surface's contour: perpendicular to the area-weighted
// mean slope direction of the shallow facets.  Returns false for a (nearly)
// flat or symmetric surface, where the normal ironing angle is used.
bool contour_angle(const NonplanarSurface &surface, double &angle)
{
    double sx = 0., sy = 0., area = 0.;
    for (const auto &entry : surface.mesh) {
        const NonplanarFacet &f = entry.second;
        const auto &v = f.vertex;
        const double ux = v[1].x - v[0].x, uy = v[1].y - v[0].y, wx = v[2].x - v[0].x, wy = v[2].y - v[0].y;
        const double a = 0.5 * std::abs(ux * wy - uy * wx);
        sx += a * f.normal.x; sy += a * f.normal.y; area += a;
    }
    if (area <= 0.)
        return false;
    sx /= area; sy /= area;
    // Mean slope below about 1 degree: no preferred direction.
    if (std::hypot(sx, sy) < 0.017)
        return false;
    angle = std::atan2(sy, sx) + 0.5 * M_PI;
    return true;
}

bool ironing_verify() { static const bool v = std::getenv("NONPLANAR_IRONING_VERIFY") != nullptr; return v; }

} // namespace

void Layer::make_nonplanar_ironing()
{
    const PrintObjectConfig &object_config = this->object()->config();
    if (! object_config.use_nonplanar_layers.value || ! object_config.nonplanar_ironing.value)
        return;

    // Radius of the flat face at the nozzle tip; the collision model keeps its own tip.
    const double tip_r     = 0.5 * object_config.nonplanar_ironing_tip_diameter.value;
    const double spacing   = object_config.nonplanar_ironing_spacing.value;
    const double max_angle = object_config.nonplanar_ironing_max_angle.value;
    const double z_offset  = object_config.nonplanar_ironing_z_offset.value;
    // Share of the gap under the middle of the tip (r * tan(slope) on a plane) that is pressed away.
    // The gap is measured, not taken from the facet slope: it changes smoothly along a line, while
    // the slope of the facet that holds the tip jumps from facet to facet.
    const double press     = 0.01 * object_config.nonplanar_ironing_slope_press.value;
    const double press_max = 0.05;  // mm
    const double step      = 0.1;   // sampling along a line (mm)
    const double min_len   = 1.0;   // shorter pieces are dropped (mm)

    // Material of the layers above this one, of every region (nonplanar strips included: they
    // lie in their home layers, above the material they cover).  Only needed around a visible
    // shell that is ironed, so it is built there and nowhere else: a union of every layer above
    // for every layer grew with the square of the layer count (Stanford Bunny, ~590 layers:
    // half an hour in this union alone).  Clipping to the shell's box keeps the result inside it.
    auto covered_above_near = [this](const BoundingBox &bbox) {
        Polygons covered;
        for (const Layer *l = this->upper_layer; l != nullptr; l = l->upper_layer)
            for (const LayerRegion *r : l->regions())
                for (const Surface &s : r->slices.surfaces)
                    if (get_extents(s.expolygon.contour).overlap(bbox))
                        polygons_append(covered, ClipperUtils::clip_clipper_polygons_with_subject_bbox(s.expolygon, bbox));
        if (! covered.empty())
            covered = union_(covered);
        return covered;
    };

    for (LayerRegion *layerm : m_regions) {
        const NonplanarSurfaces &surfaces = layerm->nonplanar_surfaces();
        if (surfaces.empty())
            continue;
        // Visible shell of each target surface in this (home) layer.
        std::map<int, ExPolygons> top_by_surface;
        for (const Surface &s : layerm->slices.surfaces)
            if (s.surface_type == stTopNonplanar && s.distance_to_top >= 0.f && s.distance_to_top < 1e-4f && s.source_surface_id >= 0)
                top_by_surface[s.source_surface_id].push_back(s.expolygon);
        if (top_by_surface.empty())
            continue;
        BoundingBox shell_bbox;
        for (const auto &[surface_id, area] : top_by_surface)
            for (const ExPolygon &e : area)
                shell_bbox.merge(get_extents(e));
        shell_bbox.offset(scale_(1.));
        const Polygons covered_above = covered_above_near(shell_bbox);

        const PrintRegionConfig &region_config = layerm->region().config();
        const int    extruder   = region_config.top_surface_filament_id.value;
        const double nozzle_dmr = this->object()->print()->config().nozzle_diameter.get_at(std::max(0, extruder - 1));
        // Same flow model as the planar ironing (Layer::make_ironing).
        const double height            = object_config.layer_height.value * 0.01 * object_config.nonplanar_ironing_flow.value;
        const double extrusion_height  = height * spacing / nozzle_dmr;
        const float  extrusion_width   = Flow::rounded_rectangle_extrusion_width_from_spacing(float(nozzle_dmr), float(extrusion_height));
        const double flow_mm3_per_mm   = nozzle_dmr * extrusion_height;

        // Material of another filament near this shell (e.g. the inlaid letter of a two-colour
        // keycap, flush with the top).  The tip carries this filament: running next to such a
        // part it would rub over it and smear this colour onto it, so the ironing keeps the
        // flat face of the tip away from it.  Layers down to the bottom of the shell's surfaces
        // are looked at, because a flush part ends in the layer of the surface below the tip.
        Polygons other_filament;
        {
            double span = 0.;
            for (const NonplanarSurface &s : surfaces)
                span = std::max(span, double(s.stats.max.z - s.stats.min.z));
            const double z_low = this->print_z - span - 2. * object_config.layer_height.value;
            // From the top down: material of another filament counts only where this filament's
            // material does not cover it from above (an inlaid letter covers the other filament's
            // layers under it; those must not keep the letter itself from being ironed).
            Polygons own_above;
            const auto &layers = this->object()->layers();
            for (auto it = layers.rbegin(); it != layers.rend(); ++ it) {
                const Layer *l = *it;
                if (l->print_z > this->print_z + EPSILON)
                    continue;
                if (l->print_z < z_low - EPSILON)
                    break;
                Polygons other_here, own_here;
                for (const LayerRegion *r : l->regions())
                    for (const Surface &s : r->slices.surfaces)
                        polygons_append(r->region().config().top_surface_filament_id.value != extruder ? other_here : own_here, to_polygons(s.expolygon));
                if (! other_here.empty())
                    polygons_append(other_filament, own_above.empty() ? other_here : diff(other_here, own_above));
                if (! own_here.empty())
                    own_above = union_(own_above, own_here);
            }
            // The flat face of the tip stays off it (plus 0.1 mm); past the face the nozzle rises.
            if (! other_filament.empty())
                other_filament = offset(union_(other_filament), float(scale_(tip_r + 0.1)));
        }

        for (auto &[surface_id, area] : top_by_surface) {
            const NonplanarSurface *surface = nullptr;
            for (const NonplanarSurface &s : surfaces)
                if (s.stable_id == surface_id) { surface = &s; break; }
            if (surface == nullptr)
                continue;
            // Keep the whole tip disk on this surface.
            ExPolygons iron_area = offset_ex(union_ex(area), -float(scale_(tip_r + 0.05)));
            if (! other_filament.empty())
                iron_area = diff_ex(iron_area, other_filament);
            // A shell that material printed later covers from above (the floor of the pocket
            // under an inlay) is not seen: ironing it only costs time.
            if (! iron_area.empty() && ! covered_above.empty())
                iron_area = opening_ex(diff_ex(iron_area, covered_above), float(scale_(0.1)));
            // Small areas are not worth a pass (and, ironing last, a tool change).
            if (const double min_area = object_config.nonplanar_ironing_min_area.value; min_area > 0.)
                iron_area.erase(std::remove_if(iron_area.begin(), iron_area.end(),
                    [min_area](const ExPolygon &e) { return e.area() * SCALING_FACTOR * SCALING_FACTOR < min_area; }), iron_area.end());
            if (iron_area.empty())
                continue;

            double angle = 0.;
            if (! contour_angle(*surface, angle))
                angle = region_config.ironing_angle.value * M_PI / 180.;

            std::unique_ptr<Fill> f(Fill::new_from_type(ipRectilinear));
            f->set_bounding_box(this->object()->bounding_box());
            f->layer_id        = this->id();
            f->z               = this->print_z;
            f->overlap         = 0;
            f->spacing         = spacing;
            f->angle           = float(angle);
            f->fixed_angle     = true;
            f->link_max_length = coord_t(scale_(3. * spacing));
            FillParams fill_params;
            fill_params.density   = 1.;
            fill_params.monotonic = true;

            const TipFacetGrid grid(*surface);
            size_t points = 0, skipped = 0;
            double min_clearance = std::numeric_limits<double>::max();
            ExtrusionPaths paths;

            auto flush = [&](Polyline3 &piece) {
                double len = 0.;
                for (size_t i = 1; i < piece.points.size(); ++ i)
                    len += std::hypot(unscale<double>(piece.points[i].x() - piece.points[i - 1].x()), unscale<double>(piece.points[i].y() - piece.points[i - 1].y()));
                if (piece.points.size() >= 2 && len >= min_len) {
                    ExtrusionPath path(erIroning, flow_mm3_per_mm, extrusion_width, float(extrusion_height));
                    path.polyline          = std::move(piece);
                    path.nonplanar_source  = true;
                    path.distance_to_top   = 0.f;
                    path.source_layer_id   = this->id();
                    path.source_surface_id = surface_id;
                    path.source_path_id    = next_nonplanar_path_id();
                    paths.emplace_back(std::move(path));
                }
                piece = Polyline3();
            };

            // Optional second pass across the first one.
            std::vector<double> pass_angles{ angle };
            if (object_config.nonplanar_ironing_cross.value)
                pass_angles.push_back(angle + 0.5 * M_PI);
            for (const double pass_angle : pass_angles)
            for (const ExPolygon &expoly : iron_area) {
                f->angle = float(pass_angle);
                Surface surface_fill(stTop, expoly);
                Polylines polylines;
                try {
                    polylines = f->fill_surface(&surface_fill, fill_params);
                } catch (InfillFailedException &) {
                }
                for (const Polyline &pl : polylines) {
                    Polyline3 piece;
                    for (size_t i = 1; i < pl.points.size(); ++ i) {
                        const Vec2d a(unscale<double>(pl.points[i - 1].x()), unscale<double>(pl.points[i - 1].y()));
                        const Vec2d b(unscale<double>(pl.points[i].x()), unscale<double>(pl.points[i].y()));
                        const double len = (b - a).norm();
                        const int    n   = std::max(1, int(std::ceil(len / step)));
                        for (int k = (i == 1 ? 0 : 1); k <= n; ++ k) {
                            const Vec2d p = a + (b - a) * (double(k) / n);
                            const NonplanarTip::Result tip = grid.tip_height(p.x(), p.y(), tip_r);
                            ++ points;
                            if (! tip.found || tip.slope_deg > max_angle + 1e-6) {
                                ++ skipped;
                                flush(piece);
                                continue;
                            }
                            double centre_z, gap = 0.;
                            if (press > 0. && grid.height_at(p.x(), p.y(), centre_z))
                                gap = std::max(0., tip.z - centre_z);
                            const double z = tip.z + z_offset - std::min(press_max, press * gap);
                            if (ironing_verify()) {
                                // Independent check: dense samples of the disk.
                                double highest = -std::numeric_limits<double>::max();
                                for (int u = -20; u <= 20; ++ u)
                                    for (int v = -20; v <= 20; ++ v) {
                                        const double dx = tip_r * u / 20., dy = tip_r * v / 20.;
                                        if (dx * dx + dy * dy > tip_r * tip_r) continue;
                                        double hz;
                                        if (grid.height_at(p.x() + dx, p.y() + dy, hz)) highest = std::max(highest, hz);
                                    }
                                if (highest > -1e30)
                                    min_clearance = std::min(min_clearance, z - highest);
                            }
                            piece.points.emplace_back(coord_t(scale_(p.x())), coord_t(scale_(p.y())), coord_t(scale_(z)));
                        }
                    }
                    flush(piece);
                }
            }

            if (ironing_verify()) {
                static std::mutex log_mutex;
                std::lock_guard<std::mutex> lock(log_mutex);
                std::cerr << "[IRON] layer " << this->id() << " z " << this->print_z << " surface " << surface_id
                          << " passes " << pass_angles.size() << " paths " << paths.size() << " points " << points << " skipped " << skipped
                          << " min_clearance " << (min_clearance < 1e30 ? min_clearance : 0.) << std::endl;
            }
            if (! paths.empty()) {
                ExtrusionEntityCollection *eec = new ExtrusionEntityCollection();
                // Monotonic like the planar ironing; the lines are already ordered.
                eec->no_sort = true;
                for (ExtrusionPath &path : paths)
                    eec->entities.push_back(new ExtrusionPath(std::move(path)));
                layerm->fills.entities.push_back(eec);
            }
        }
    }
}

} // namespace Slic3r
