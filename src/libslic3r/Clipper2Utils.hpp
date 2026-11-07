#ifndef slic3r_Clipper2Utils_hpp_
#define slic3r_Clipper2Utils_hpp_

// OrcaSlicer Nonplanar: on 2.4.2 these were Clipper2 versions of the boolean ops, because
// Clipper 1 could hang a hole on the wrong outer contour.  Since 2.5.0 ClipperUtils runs on
// Clipper2 itself and fixes misnested holes, so the plain functions are used.
#include "ClipperUtils.hpp"

namespace Slic3r {
inline ExPolygons intersection_ex_2(const Polygons &subject, const Polygons &clip) { return intersection_ex(subject, clip); }
inline ExPolygons diff_ex_2(const Polygons &subject, const Polygons &clip) { return diff_ex(subject, clip); }
}

#endif
