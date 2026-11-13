#include "NonplanarSurface.hpp"
#include "Exception.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <boost/log/trivial.hpp>

namespace Slic3r {

NonplanarSurface::NonplanarSurface(std::map<int, NonplanarFacet> &_mesh)
{
    this->mesh = _mesh;
    this->calculate_stats();
}

bool 
NonplanarSurface::operator==(const NonplanarSurface& other) const {
   return (stats.min.x == other.stats.min.x &&
           stats.min.y == other.stats.min.y &&
           stats.min.z == other.stats.min.z &&
           stats.max.x == other.stats.max.x &&
           stats.max.y == other.stats.max.y &&
           stats.max.z == other.stats.max.z);
}

void
NonplanarSurface::calculate_stats()
{
    //calculate min and max values
    this->stats.min.x = 10000000;
    this->stats.min.y = 10000000;
    this->stats.min.z = 10000000;
    this->stats.max.x = -10000000;
    this->stats.max.y = -10000000;
    this->stats.max.z = -10000000;
    for(auto& facet : this->mesh) {
        // Edge-ring prototype: rim facets only exist so the extra ring can be
        // projected onto the steeper edge.  They must not change the surface's
        // extents (home-layer choice, max-height and area checks, shell offsets).
        if (facet.second.rim)
            continue;
        this->stats.min.x = std::min(this->stats.min.x, facet.second.stats.min.x);
        this->stats.min.y = std::min(this->stats.min.y, facet.second.stats.min.y);
        this->stats.min.z = std::min(this->stats.min.z, facet.second.stats.min.z);
        this->stats.max.x = std::max(this->stats.max.x, facet.second.stats.max.x);
        this->stats.max.y = std::max(this->stats.max.y, facet.second.stats.max.y);
        this->stats.max.z = std::max(this->stats.max.z, facet.second.stats.max.z);
    }
}

void
NonplanarSurface::translate(float x, float y, float z)
{
    m_height_index.reset();
    //translate all facets
    for(auto& facet : this->mesh) {
        facet.second.translate(x, y, z);
    }

    //translate min and max values
    this->stats.min.x += x;
    this->stats.min.y += y;
    this->stats.min.z += z;
    this->stats.max.x += x;
    this->stats.max.y += y;
    this->stats.max.z += z;
}

void
NonplanarSurface::scale(float factor)
{
    float versor[3];
    versor[0] = factor;
    versor[1] = factor;
    versor[2] = factor;
    this->scale(versor);
}

void
NonplanarSurface::scale(float versor[3])
{
    m_height_index.reset();
    //scale all facets
    for(auto& facet : this->mesh) {
        facet.second.scale(versor);
    }

    //scale min and max values
    this->stats.min.x *= versor[0];
    this->stats.min.y *= versor[1];
    this->stats.min.z *= versor[2];
    this->stats.max.x *= versor[0];
    this->stats.max.y *= versor[1];
    this->stats.max.z *= versor[2];
}

void
NonplanarSurface::rotate_z(float angle) {
  m_height_index.reset();
  double radian_angle = (angle / 180.0) * 3.1415927;
  double c = cos(radian_angle);
  double s = sin(radian_angle);

  for(auto& facet : this->mesh) {
    for(int j = 0; j < 3; j++) {
        double xold = facet.second.vertex[j].x;
        double yold = facet.second.vertex[j].y;
        facet.second.vertex[j].x = c * xold - s * yold;
        facet.second.vertex[j].y = s * xold + c * yold;
    }
    facet.second.calculate_stats();
  }

  this->calculate_stats();
}

void
NonplanarSurface::debug_output()
{
    std::cout << "Facets(" << this->mesh.size() << "): (min:X:" << this->stats.min.x << " Y:" << this->stats.min.y << " Z:" << this->stats.min.z <<
                           " max:X:" << this->stats.max.x << " Y:" << this->stats.max.y << " Z:" << this->stats.max.z << ")" << 
                           "Height " << this->stats.max.z - this->stats.min.z << std::endl;
    for(auto& facet : this->mesh) {
        std::cout << "triangle: (" << facet.first << ")(" << facet.second.marked << ") ";
        std::cout << " (" << (180*std::acos(facet.second.normal.z))/3.14159265 << "°)";

        std::cout << " | V0:";
        std::cout << " X:"<< facet.second.vertex[0].x;
        std::cout << " Y:"<< facet.second.vertex[0].y;
        std::cout << " Z:"<< facet.second.vertex[0].z;

        std::cout << " | V1:";
        std::cout << " X:"<< facet.second.vertex[1].x;
        std::cout << " Y:"<< facet.second.vertex[1].y;
        std::cout << " Z:"<< facet.second.vertex[1].z;

        std::cout << " | V2:";
        std::cout << " X:"<< facet.second.vertex[2].x;
        std::cout << " Y:"<< facet.second.vertex[2].y;
        std::cout << " Z:"<< facet.second.vertex[2].z;

        std::cout << " | Normal:";
        std::cout << " X:"<< facet.second.normal.x;
        std::cout << " Y:"<< facet.second.normal.y;
        std::cout << " Z:"<< facet.second.normal.z;

        //TODO check if neighbors exist
        // stl_neighbors* neighbors = mesh.stl.neighbors_start + facet.first;
        std::cout << " | Neighbors:";
        std::cout << " 0:"<< facet.second.neighbor[0];
        std::cout << " 1:"<< facet.second.neighbor[1];
        std::cout << " 2:"<< facet.second.neighbor[2];
        std::cout << std::endl;
    }
}

NonplanarSurfaces
NonplanarSurface::group_surfaces()
{
    // Splits off one connected group per round.  This used to recurse once per
    // group and overflowed the stack on meshes with thousands of small groups;
    // the result keeps the order of the recursive version (last group first,
    // this one last).
    NonplanarSurfaces nonplanar_surfaces;
    NonplanarSurface  rest;
    NonplanarSurface *current = this;
    while (! current->mesh.empty()) {
        current->mark_neighbor_surfaces(current->mesh.begin()->first);
        NonplanarSurface newSurface;
        for (std::map<int, NonplanarFacet>::iterator it = current->mesh.begin(); it != current->mesh.end();) {
            if ((*it).second.marked == false) {
                newSurface.mesh[(*it).first] = (*it).second;
                it = current->mesh.erase(it);
            } else
                ++it;
        }
        current->calculate_stats();
        nonplanar_surfaces.push_back(*current);
        if (newSurface.mesh.empty())
            break;
        rest    = std::move(newSurface);
        current = &rest;
    }
    std::reverse(nonplanar_surfaces.begin(), nonplanar_surfaces.end());
    return nonplanar_surfaces;
}

void
NonplanarSurface::mark_neighbor_surfaces(int id)
{
    // Flood fill with an explicit stack: the recursive version used one stack
    // frame per facet and overflowed the stack on finely meshed surfaces
    // (tens of thousands of facets), which crashed the slicer without a message.
    std::vector<int> todo { id };
    while (! todo.empty()) {
        const int current = todo.back();
        todo.pop_back();
        auto it = this->mesh.find(current);
        if (it == this->mesh.end() || it->second.marked)
            continue;
        it->second.marked = true;
        for (int j = 0; j < 3; j++)
            todo.push_back(it->second.neighbor[j]);
    }
}

bool
NonplanarSurface::check_max_printing_height(float height) 
{
    if ((this->stats.max.z - this->stats.min.z) > height ) {
        BOOST_LOG_TRIVIAL(trace) << "Surface removed: printheight too heigh (" << (this->stats.max.z - this->stats.min.z) << " mm)";
        return true;
    } else {
        return false;
    }
}

bool
NonplanarSurface::check_surface_area()
{
    //calculate surface area of nonplanar surface.
    float area = 0.0f;
    for(auto& facet : this->mesh) {
        // See calculate_stats(): rim facets are projection targets only.
        if (facet.second.rim)
            continue;
        area += facet.second.calculate_surface_area();
    }
    if (area < 20.0f) {
        BOOST_LOG_TRIVIAL(trace) << "Surface removed: area too small (" << area << " mm²)";
        return true;
    } else {
        return false;
    }
}

void
NonplanarSurface::check_printable_surfaces(float max_angle)
{
    //TODO do something
}

/* this will return scaled ExPolygons */
ExPolygons
NonplanarSurface::horizontal_projection(Projection kind) const
{
    Polygons pp;
    pp.reserve(this->mesh.size() + this->rim.size());
    auto append = [&pp](const std::map<int, NonplanarFacet> &facets) {
        for (auto &facet : facets) {
            Polygon p;
            p.points.resize(3);
            p.points[0] = Point(scale_(facet.second.vertex[0].x), scale_(facet.second.vertex[0].y));
            p.points[1] = Point(scale_(facet.second.vertex[1].x), scale_(facet.second.vertex[1].y));
            p.points[2] = Point(scale_(facet.second.vertex[2].x), scale_(facet.second.vertex[2].y));
            p.make_counter_clockwise();  // do this after scaling, as winding order might change while doing that
            pp.push_back(p);
        }
    };
    append(this->mesh);
    if (kind == Projection::Target)
        append(this->rim);

    // the offset factor was tuned using groovemount.stl
    return union_ex(offset(pp, scale_(0.01)));
}

const NonplanarSurface::HeightIndex&
NonplanarSurface::height_index() const
{
    if (m_height_index && m_height_index->owner == this && m_height_index->facets == this->mesh.size())
        return *m_height_index;
    auto idx = std::make_shared<HeightIndex>();
    idx->owner  = this;
    idx->facets = this->mesh.size();
    double x0 = 1e30, y0 = 1e30, x1 = -1e30, y1 = -1e30;
    for (const auto &entry : this->mesh) {
        const NonplanarFacet &f = entry.second;
        idx->order.push_back(&f);
        x0 = std::min(x0, double(f.stats.min.x)); y0 = std::min(y0, double(f.stats.min.y));
        x1 = std::max(x1, double(f.stats.max.x)); y1 = std::max(y1, double(f.stats.max.y));
    }
    if (! idx->order.empty()) {
        // About 1 mm cells, at most 1024 per side.
        idx->cell = std::max({ 1., (x1 - x0) / 1024., (y1 - y0) / 1024. });
        idx->x0 = x0; idx->y0 = y0;
        idx->nx = int((x1 - x0) / idx->cell) + 1;
        idx->ny = int((y1 - y0) / idx->cell) + 1;
        idx->cells.assign(size_t(idx->nx) * size_t(idx->ny), {});
        for (int i = 0; i < int(idx->order.size()); ++ i) {
            const NonplanarFacet &f = *idx->order[i];
            const int ix0 = std::clamp(int((f.stats.min.x - x0) / idx->cell), 0, idx->nx - 1), ix1 = std::clamp(int((f.stats.max.x - x0) / idx->cell), 0, idx->nx - 1);
            const int iy0 = std::clamp(int((f.stats.min.y - y0) / idx->cell), 0, idx->ny - 1), iy1 = std::clamp(int((f.stats.max.y - y0) / idx->cell), 0, idx->ny - 1);
            for (int iy = iy0; iy <= iy1; ++ iy)
                for (int ix = ix0; ix <= ix1; ++ ix)
                    idx->cells[size_t(iy) * size_t(idx->nx) + size_t(ix)].push_back(i);
        }
    }
    m_height_index = idx;
    return *m_height_index;
}

bool
NonplanarSurface::height_near(double x, double y, double max_dist, double &z) const
{
    double best = max_dist * max_dist;
    bool   found = false;
    // Candidates are the facets whose box (plus max_dist) contains (x, y), visited in
    // the order of the mesh map, so the result is the same as scanning all facets.
    const HeightIndex &idx = this->height_index();
    if (idx.order.empty())
        return false;
    std::vector<int> candidates;
    {
        const int ix0 = std::clamp(int(std::floor((x - max_dist - idx.x0) / idx.cell)), 0, idx.nx - 1);
        const int ix1 = std::clamp(int(std::floor((x + max_dist - idx.x0) / idx.cell)), 0, idx.nx - 1);
        const int iy0 = std::clamp(int(std::floor((y - max_dist - idx.y0) / idx.cell)), 0, idx.ny - 1);
        const int iy1 = std::clamp(int(std::floor((y + max_dist - idx.y0) / idx.cell)), 0, idx.ny - 1);
        for (int iy = iy0; iy <= iy1; ++ iy)
            for (int ix = ix0; ix <= ix1; ++ ix) {
                const std::vector<int> &c = idx.cells[size_t(iy) * size_t(idx.nx) + size_t(ix)];
                candidates.insert(candidates.end(), c.begin(), c.end());
            }
        std::sort(candidates.begin(), candidates.end());
        candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    }
    for (const int ci : candidates) {
        const NonplanarFacet &f = *idx.order[ci];
        if (x < f.stats.min.x - max_dist || x > f.stats.max.x + max_dist ||
            y < f.stats.min.y - max_dist || y > f.stats.max.y + max_dist)
            continue;
        const auto &a = f.vertex[0], &b = f.vertex[1], &c = f.vertex[2];
        const double den = (b.y-a.y)*(c.x-a.x) - (b.x-a.x)*(c.y-a.y);
        if (std::abs(den) > 1e-12) {
            const double w1 = ((y-a.y)*(c.x-a.x) - (x-a.x)*(c.y-a.y)) / den;
            const double w2 = ((b.y-a.y)*(x-a.x) - (b.x-a.x)*(y-a.y)) / den;
            const double w0 = 1.0 - w1 - w2;
            if (w0 >= -1e-7 && w1 >= -1e-7 && w2 >= -1e-7) {
                z = w0*a.z + w1*b.z + w2*c.z;
                return true;
            }
        }
        for (int i = 0; i < 3; ++ i) {
            const auto &p = f.vertex[i], &q = f.vertex[(i+1)%3];
            const double dx = q.x-p.x, dy = q.y-p.y, len2 = dx*dx + dy*dy;
            if (len2 < 1e-12) continue;
            const double t  = std::max(0.0, std::min(1.0, ((x-p.x)*dx + (y-p.y)*dy) / len2));
            const double ex = x-p.x-t*dx, ey = y-p.y-t*dy, d2 = ex*ex + ey*ey;
            if (d2 <= best) {
                best  = d2;
                z     = p.z + t*(q.z-p.z);
                found = true;
            }
        }
    }
    return found;
}

double nonplanar_edge_reach(bool enabled)
{
    // The environment is a test hook only: NONPLANAR_EDGE_RING=0 forces the fix off,
    // =1 forces it on regardless of the setting.  Both are read once; the setting
    // itself is evaluated on every call, so it follows the object's config live.
    static const int env_mode = []() {
        const char *on = std::getenv("NONPLANAR_EDGE_RING");
        if (on != nullptr && std::strcmp(on, "0") == 0)
            return 0;   // forced off
        if (on != nullptr && std::strcmp(on, "1") == 0)
            return 1;   // forced on
        return -1;      // follow the setting
    }();
    static const double env_reach = []() {
        const char  *v = std::getenv("NONPLANAR_EDGE_REACH");
        const double r = v != nullptr ? std::atof(v) : 0.0;
        return (r > 0.05 && r <= 3.0) ? r : 1.0;
    }();
    if (env_mode == 0)
        return 0.0;
    if (env_mode < 0 && !enabled)
        return 0.0;
    return env_reach;
}

std::string nonplanar_location(double z, size_t layer_id)
{
    char buf[80];
    std::snprintf(buf, sizeof(buf), " at height %.2f mm (layer %zu)", z, layer_id + 1);
    return buf;
}

bool nonplanar_error_height(const std::string &message, double &z)
{
    // Reads back the height written by nonplanar_location().
    const std::string key = " at height ";
    const size_t pos = message.find(key);
    if (pos == std::string::npos)
        return false;
    const char *begin = message.c_str() + pos + key.size();
    char *end = nullptr;
    z = std::strtod(begin, &end);
    return end != begin && std::strncmp(end, " mm", 3) == 0;
}

double nonplanar_suggested_angle(double angle)
{
    return angle > 10. + EPSILON ? 10. : std::max(1., std::floor(angle) - 2.);
}

void throw_nonplanar_too_steep(double z, size_t layer_id, double angle, size_t object_id, const std::string &detail)
{
    BOOST_LOG_TRIVIAL(error) << "Nonplanar surface too steep" << nonplanar_location(z, layer_id) << ": " << detail;
    char buf[400];
    std::snprintf(buf, sizeof(buf), "Nonplanar slicing failed%s: a nonplanar line would have to climb more than 12° there "
        "(a sharp edge or step in the surface). Lower \"Maximum nonplanar angle\" to %g° or less and slice again.",
        nonplanar_location(z, layer_id).c_str(), nonplanar_suggested_angle(angle));
    throw Slic3r::SlicingError(buf, object_id);
}

void throw_nonplanar_internal(const std::string &where, const std::string &detail, size_t object_id)
{
    BOOST_LOG_TRIVIAL(error) << "Nonplanar internal check failed" << where << ": " << detail;
    throw Slic3r::SlicingError("Nonplanar slicing stopped" + where + ": an internal check failed (" + detail + "). "
        "This is a bug, not a setting problem: please report it with the 3MF file. "
        "To slice now, turn off nonplanar layers for this object.", object_id);
}

}
