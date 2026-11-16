#ifndef slic3r_NonplanarSurface_hpp_
#define slic3r_NonplanarSurface_hpp_

#include "libslic3r.h"
#include <memory>
#include "NonplanarFacet.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "ExPolygon.hpp"
#include "Geometry.hpp"
#include "ClipperUtils.hpp"

namespace Slic3r {

typedef struct {
  float x;
  float y;
  float z;
} mesh_vertex;

typedef struct {
  mesh_vertex    max;
  mesh_vertex    min;
} mesh_stats;

class NonplanarSurface;

typedef std::vector<NonplanarSurface> NonplanarSurfaces;

class NonplanarSurface
{
    public:
    // What a horizontal projection is used for.  Rim facets (the steeper fillet
    // between the target surface and a wall) must be part of the area a conformal
    // strip may cover -- that is what lets the strip reach the wall with the model's
    // own outline -- but they must not change the surface's extents, area or home
    // layer, which stay defined by the shallow facets alone.
    enum class Projection {
        Area,     // shallow facets only: stats, area / height checks
        Target,   // shallow + rim: strip generation, collision check
    };
    std::map<int, NonplanarFacet> mesh;
    // Facets kept only as projection targets (slope above the nonplanar angle).
    std::map<int, NonplanarFacet> rim;
    mesh_stats stats;
    // A piece of a surface rebuilt with a lower angle keeps the extents of the surface it came from.
    bool       has_origin = false;
    mesh_stats origin_stats {};
    int stable_id = -1;
    NonplanarSurface() {};
    ~NonplanarSurface() {};
    NonplanarSurface(std::map<int, NonplanarFacet> &_mesh);
    bool operator==(const NonplanarSurface& other) const;
    void calculate_stats();
    void translate(float x, float y, float z);
    void scale(float factor);
    void scale(float versor[3]);
    void rotate_z(float angle);
    void debug_output();
    NonplanarSurfaces group_surfaces();
    void mark_neighbor_surfaces(int id);
    bool check_max_printing_height(float height);
    void check_printable_surfaces(float max_angle);
    bool check_surface_area();
    ExPolygons horizontal_projection(Projection kind = Projection::Target) const;
    // Height of the target surface at (x, y), or of its nearest boundary point when
    // (x, y) lies outside the surface but within max_dist (unscaled mm).
    bool height_near(double x, double y, double max_dist, double &z) const;

private:
    // Grid over the facets for height_near(), built on first use.  It belongs to one
    // object and one mesh size; a copy or a changed mesh builds its own.
    struct HeightIndex {
        const NonplanarSurface              *owner { nullptr };
        size_t                               facets { 0 };
        double                               x0 { 0 }, y0 { 0 }, cell { 1. };
        int                                  nx { 0 }, ny { 0 };
        std::vector<const NonplanarFacet*>   order;   // map order
        std::vector<std::vector<int>>        cells;   // indices into order, ascending
    };
    mutable std::shared_ptr<HeightIndex> m_height_index;
    const HeightIndex& height_index() const;
};

// User-facing nonplanar errors share one format: what happened, where (height
// and layer as shown by the layer slider), what to do.  The technical detail
// goes to the log and, for internal checks, in short into the message.
// " at height 8.40 mm (layer 43)"; layer_id is 0-based.
std::string nonplanar_location(double z, size_t layer_id);
// Height written by nonplanar_location() into a slicing error message.
bool nonplanar_error_height(const std::string &message, double &z);
// Edge extension (setting "Nonplanar edge fix", on by default; NONPLANAR_EDGE_RING=0/1
// forces it off/on for tests): how far (mm) beyond the target surface a nonplanar shell
// may take over the thin planar remainder next to the wall, so the shell's outer line
// runs on the model's own outline.  Returns 0 when `enabled` is false.  0 = off.
double      nonplanar_edge_reach(bool enabled);
// Angle to suggest after a failure at `angle` (12 -> 10).
double      nonplanar_suggested_angle(double angle);
// A surface step a nonplanar line cannot follow (steeper than 12 degrees).
[[noreturn]] void throw_nonplanar_too_steep(double z, size_t layer_id, double angle, size_t object_id, const std::string &detail);
// An internal consistency check; not caused by a setting.
[[noreturn]] void throw_nonplanar_internal(const std::string &where, const std::string &detail, size_t object_id);
};

#endif
