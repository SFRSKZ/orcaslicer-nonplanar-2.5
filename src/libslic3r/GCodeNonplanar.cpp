// OrcaSlicer Nonplanar: deferred-material and final-ironing members of GCode, kept out of GCode.cpp
// so that merging a new OrcaSlicer release touches less nonplanar code.
#include "ArcFitter.hpp"
#include <tuple>
#include "BoundingBox.hpp"
#include "Circle.hpp"
#include "Config.hpp"
#include "GCode/GCodeProcessor.hpp"
#include "Flow.hpp"
#include "GCode/ThumbnailData.hpp"
#include "GCode/ToolOrdering.hpp"
#include "GCode/SpiralVase.hpp"
#include "GCode/PressureEqualizer.hpp"
#include "GCode/SmallAreaInfillFlowCompensator.hpp"
#include "GCode/CoolingBuffer.hpp"
#include "GCode/AdaptivePAProcessor.hpp"
#include "CustomGCode.hpp"
#include "GCode/TimelapsePosPicker.hpp"
#include "ExPolygon.hpp"
#include "GCode/WipePathHelpers.hpp"
#include "GCodeReader.hpp"
#include "GCodeWriter.hpp"
#include "Point.hpp"
#include "Line.hpp"
#include "ObjectID.hpp"
#include "Layer.hpp"
#include "Polygon.hpp"
#include "Polyline.hpp"
#include "PrintBase.hpp"
#include "ConnectedBodies.hpp"
#include "PrintConfig.hpp"
#include "enum_bitmask.hpp"
#include "libslic3r.h"
#include "I18N.hpp"
#include "IMEXHelpers.hpp"
#include "GCode.hpp"
#include <cstdio>
#include "Feature/FuzzySkin/FuzzySkin.hpp"
#include "Exception.hpp"
#include "LifecycleEvents.hpp"
#include "ExtrusionEntity.hpp"
#include "EdgeGrid.hpp"
#include "BeltTransform.hpp"
#include "Geometry.hpp"
#include "Geometry/ConvexHull.hpp"
#include "GCode/PrintExtents.hpp"
#include "GCode/Thumbnails.hpp"
#include "GCode/WipeTower.hpp"
#include "GCode/WipeTower2.hpp"
#include "ShortestPath.hpp"
#include "GCode/OrderingStrategies.hpp"
#include "Print.hpp"
#include "Utils.hpp"
#include "ClipperUtils.hpp"
#include "LocalesUtils.hpp"
#include "libslic3r/format.hpp"
#include "Time.hpp"
#include "GCode/ExtrusionProcessor.hpp"
#include <Eigen/Geometry>
#include <Shiny/ShinyMacros.h>
#include <algorithm>
#include <cctype>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdarg>
#include <cstdlib>
#include <chrono>
#include <iomanip>
#include <cstring>
#include <exception>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <math.h>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <set>
#include <optional>
#include <stdlib.h>
#include <string>
#include <system_error>
#include <unordered_set>
#include <type_traits>
#include <utility>
#include <string_view>
#include <sstream>
#include <regex>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/find.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/beast/core/detail/base64.hpp>
#include <boost/nowide/iostream.hpp>
#include <boost/nowide/cstdio.hpp>
#include <boost/nowide/cstdlib.hpp>
#include "SVG.hpp"
#include <tbb/parallel_for.h>
#include <vector>
#include "calib.hpp"
#include "libslic3r_version.h"
#include "GCode/BeltKinematics.hpp"
#include <Shiny/Shiny.h>
#include <stdio.h>
#include <assert.h>
#include "AABBTreeLines.hpp"
#include "Extruder.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "FilamentMixer.hpp"
#include "Format/STEP.hpp"
#include "Model.hpp"
#include "MultiNozzleUtils.hpp"
#include "Slicing.hpp"
#include "GCode/NonplanarEmission.hpp"

namespace Slic3r {

void GCode::queue_advanced_conformal(const Print &print)
{
    for (const auto &[object, plan] : m_nonplanar_deferral) {
        if (plan.advanced_paths == 0)
            continue;
        for (const Layer *layer : object->layers()) {
            const int base = 2 * int(layer->id());
            for (const LayerRegion *layerm : layer->regions()) {
                const PrintRegion &region = print.get_print_region(layerm->region().print_region_id());
                auto queue_paths = [&](const ExtrusionEntitiesPtr &entities, const std::string &description, FlowRole role) {
                    std::vector<const ExtrusionPath*> paths;
                    std::function<void(const ExtrusionEntity*)> collect = [&](const ExtrusionEntity *e) {
                        if (auto *path = dynamic_cast<const ExtrusionPath*>(e)) paths.push_back(path);
                        else if (auto *loop = dynamic_cast<const ExtrusionLoop*>(e)) for (const ExtrusionPath &q : loop->paths) paths.push_back(&q);
                        else if (auto *multi = dynamic_cast<const ExtrusionMultiPath*>(e)) for (const ExtrusionPath &q : multi->paths) paths.push_back(&q);
                        else if (auto *c = dynamic_cast<const ExtrusionEntityCollection*>(e)) for (const ExtrusionEntity *child : c->entities) collect(child);
                    };
                    for (const ExtrusionEntity *e : entities)
                        collect(e);
                    for (const ExtrusionPath *path : paths) {
                        if (! path->nonplanar_source)
                            continue;
                        bool reversed = false;
                        const auto *path_plan = plan.find(layer->id(), path->source_path_id, path->polyline, reversed);
                        if (path_plan == nullptr)
                            continue;
                        for (const auto &[slot, a, b] : plan_pieces(*path_plan, reversed, polyline_length_scaled(path->polyline))) {
                            if (slot >= base)
                                continue;
                            Polyline3 piece = sub_polyline(path->polyline, a, b);
                            if (piece.points.size() < 2)
                                continue;
                            double zmin = std::numeric_limits<double>::max();
                            for (const Point3 &q : piece.points)
                                zmin = std::min(zmin, unscale<double>(q.z()));
                            for (size_t instance = 0; instance < object->instances().size(); ++ instance)
                                m_deferred_extrusions[{ object, instance }].push_back({ slot, base, zmin, layer->id(), layer->print_z,
                                    std::make_shared<ExtrusionPath>(piece, *path), description, -1., &region,
                                    unsigned(std::max(0, int(region.extruder(role)) - 1)) });
                        }
                    }
                };
                queue_paths(layerm->perimeters.entities, "perimeter", frPerimeter);
                queue_paths(layerm->fills.entities, "infill", frSolidInfill);
            }
        }
    }
}

bool GCode::defer_planar_material(const ExtrusionEntity &entity, const std::string &description, double speed, std::string &gcode)
{
    std::vector<const ExtrusionPath*> paths;
    if (auto *path = dynamic_cast<const ExtrusionPath*>(&entity))
        paths.push_back(path);
    else if (auto *loop = dynamic_cast<const ExtrusionLoop*>(&entity))
        for (const ExtrusionPath &p : loop->paths) paths.push_back(&p);
    else if (auto *multi = dynamic_cast<const ExtrusionMultiPath*>(&entity))
        for (const ExtrusionPath &p : multi->paths) paths.push_back(&p);

    struct Planned { const ExtrusionPath *path; const NonplanarDeferral::PathPlan *plan; bool reversed; };
    std::vector<Planned> planned;
    bool any = false, cut = false;
    int  slot = -1;
    for (const ExtrusionPath *path : paths) {
        bool reversed = false;
        const auto *plan = m_deferral_plan->find(m_deferral_layer_id, path->source_path_id, path->polyline, reversed);
        planned.push_back({ path, plan, reversed });
        if (plan != nullptr) {
            any = true;
            cut |= plan->breaks.size() > 1;
            slot = std::max(slot, plan->breaks.front().second);
        }
    }
    const int base = 2 * int(m_deferral_layer_id);
    if (! any || (! cut && slot == base))
        return false;
    // A whole path emitted before its home layer was queued at export start.
    if (! cut && slot < base)
        return true;

    const double layer_z = this->layer()->print_z;
    auto queue = [&](std::shared_ptr<ExtrusionEntity> unit, int unit_slot) {
        double zmin = std::numeric_limits<double>::max();
        // A planar unit keeps its own layer height (absolute, like a conformal path).
        auto at_height = [&](ExtrusionPath &path) {
            if (! path.nonplanar_source)
                for (Point3 &p : path.polyline.points)
                    p.z() = scaled<coord_t>(layer_z);
            for (const Point3 &p : path.polyline.points)
                zmin = std::min(zmin, unscale<double>(p.z()));
        };
        if (auto *loop = dynamic_cast<ExtrusionLoop*>(unit.get())) {
            // The seam is placed now, in the loop's own layer, as
            // extrude_loop() would; the later emission keeps it.
            float seam_overhang = std::numeric_limits<float>::lowest();
            if (! m_config.spiral_mode && description == "perimeter")
                m_seam_placer.place_seam(m_layer, *loop, this->last_pos(), seam_overhang);
            else
                loop->split_at(this->last_pos(), false);
            for (ExtrusionPath &path : loop->paths) at_height(path);
        } else if (auto *multi = dynamic_cast<ExtrusionMultiPath*>(unit.get())) {
            for (ExtrusionPath &path : multi->paths) at_height(path);
        } else if (auto *path = dynamic_cast<ExtrusionPath*>(unit.get())) {
            at_height(*path);
        }
        m_deferred_extrusions[m_deferral_key].push_back({ unit_slot, base, zmin, m_deferral_layer_id, layer_z, std::move(unit),
            description, speed, m_current_region, m_writer.filament() ? m_writer.filament()->id() : 0u });
    };

    if (! cut) {
        // The whole unit (shell path, wall loop, path) moves; never cut.
        queue(std::shared_ptr<ExtrusionEntity>(entity.clone()), slot);
        return true;
    }

    // A refined wall loop: emitted as pieces by slot.
    std::vector<ExtrusionPath> keep;
    for (const Planned &item : planned) {
        const ExtrusionPath &path = *item.path;
        if (item.plan == nullptr) {
            keep.push_back(path);
            continue;
        }
        for (const auto &[piece_slot, a, b] : plan_pieces(*item.plan, item.reversed, polyline_length_scaled(path.polyline))) {
            Polyline3 piece = sub_polyline(path.polyline, a, b);
            if (piece.points.size() < 2)
                continue;
            if (piece_slot < base)
                continue;                   // queued at the start of the export
            if (piece_slot == base)
                keep.emplace_back(std::move(piece), path);
            else
                queue(std::make_shared<ExtrusionPath>(std::move(piece), path), piece_slot);
        }
    }
    for (ExtrusionPath &path : keep)
        gcode += this->extrude_path(path, description, speed);
    return true;
}

std::string GCode::emit_deferred_material(int slot, unsigned int extruder_id)
{
    std::string gcode;
    auto &queue = m_deferred_extrusions[m_deferral_key];
    std::vector<DeferredExtrusion> now;
    std::vector<DeferredExtrusion> later;
    for (DeferredExtrusion &d : queue)
        (d.slot == slot && d.extruder == extruder_id ? now : later).push_back(std::move(d));
    queue = std::move(later);
    // Within a slot the plan orders by regular emission time only.
    std::stable_sort(now.begin(), now.end(), [](const DeferredExtrusion &a, const DeferredExtrusion &b) { return a.base < b.base; });
    m_force_safe_travel = true;
    m_force_safe_travel_once = false;
    for (DeferredExtrusion &d : now) {
        if (d.region != nullptr)
            m_config.apply(d.region->config());
        m_deferred_region = d.region;
        char buf[160];
        sprintf(buf, ";PATH_DEFERRED from_layer=%zu from_layer_z=%.3f min_z=%.3f emitted_layer=%zu slot=%d\n", d.source_layer, d.source_layer_z, d.z, m_deferral_layer_id, slot);
        gcode += buf;
        const NonplanarDeferral *plan = m_deferral_plan;
        m_deferral_plan = nullptr;
        m_loop_seam_placed = true;
        m_emitting_deferred = true;
        gcode += this->extrude_entity(*d.entity, d.description, d.speed);
        m_emitting_deferred = false;
        gcode += ";PATH_DEFERRED_END\n";
        m_loop_seam_placed = false;
        m_deferral_plan = plan;
    }
    // The next travel leaves the deferred height and is lifted as well.
    m_force_safe_travel_once = true;
    return gcode;
}

std::string GCode::emit_deferred_after_layer(const PrintObject &print_object, unsigned int extruder_id)
{
    std::string gcode;
    if (m_deferral_plan == nullptr)
        return gcode;
    // Deferred lower material of this object goes after everything of this
    // layer, ascending in height.  After the last layer of the object, the
    // remaining slots follow in ascending order.
    const int  slot = 2 * int(m_deferral_layer_id) + 1;
    const bool last = m_deferral_layer_id == print_object.layers().back()->id();
    // The plan orders material by height only, not by filament: material of a filament
    // can be planned after a layer that filament does not print (e.g. the top of a
    // two-colour inlay deferred past the layer where the inlay ends).  Such material is
    // written at this filament's next turn on this object; on its last turn (no later
    // layer uses it) everything still queued for it is written now, as late as possible.
    bool used_later = false;
    if (m_deferral_tool_ordering != nullptr && m_deferral_layer_id < print_object.layers().size()) {
        const double z = print_object.layers()[m_deferral_layer_id]->print_z;
        for (const LayerTools &lt : m_deferral_tool_ordering->layer_tools())
            if (lt.print_z > z + EPSILON && std::find(lt.extruders.begin(), lt.extruders.end(), extruder_id) != lt.extruders.end()) {
                used_later = true;
                break;
            }
    } else
        used_later = true;
    auto it = m_deferred_extrusions.find(m_deferral_key);
    while (it != m_deferred_extrusions.end()) {
        int next = std::numeric_limits<int>::max();
        for (const DeferredExtrusion &d : it->second)
            if (d.extruder == extruder_id && (d.slot <= slot || last || ! used_later))
                next = std::min(next, d.slot);
        if (next == std::numeric_limits<int>::max())
            break;
        gcode += this->emit_deferred_material(next, extruder_id);
    }
    m_deferral_plan = nullptr;
    return gcode;
}

std::string GCode::emit_final_ironing(double print_z)
{
    std::string gcode;
    if (m_final_ironing.empty())
        return gcode;
    // The filament that is loaded irons first; the others follow, and the layer ends with the
    // filament it would have ended with, so the next layer's tool order (and the prime tower)
    // stays as planned.
    const unsigned int back_to = m_writer.filament() ? m_writer.filament()->id() : m_final_ironing.front().extruder;
    std::vector<unsigned int> order { back_to };
    for (const FinalIroning &f : m_final_ironing)
        if (std::find(order.begin(), order.end(), f.extruder) == order.end())
            order.push_back(f.extruder);
    gcode += "; NONPLANAR_IRONING_LAST_START\n";
    const Vec2d origin = m_origin;
    m_force_safe_travel      = true;
    m_force_safe_travel_once = false;
    m_emitting_final_ironing = true;
    for (unsigned int extruder : order) {
        bool any = false;
        for (const FinalIroning &f : m_final_ironing)
            any |= f.extruder == extruder;
        if (! any)
            continue;
        if (m_writer.need_toolchange(extruder)) {
            gcode += this->set_extruder(extruder, print_z);
            // The filament change ends retracted by its long toolchange retraction. Prime it here,
            // still at the change position away from the part, and retract normally: otherwise
            // the whole length would be pushed out at the start of the ironing, onto the surface.
            gcode += m_writer.unretract();
            gcode += m_writer.retract();
        }
        for (const FinalIroning &f : m_final_ironing) {
            if (f.extruder != extruder)
                continue;
            if (f.region != nullptr)
                m_config.apply(f.region->config());
            this->set_origin(f.origin);
            gcode += this->extrude_path(*f.path, f.description, f.speed);
        }
    }
    if (m_writer.need_toolchange(back_to)) {
        gcode += this->set_extruder(back_to, print_z);
        gcode += m_writer.unretract();
        gcode += m_writer.retract();
    }
    m_emitting_final_ironing = false;
    m_force_safe_travel      = false;
    m_force_safe_travel_once = true;
    this->set_origin(origin);
    gcode += "; NONPLANAR_IRONING_LAST_END\n";
    m_final_ironing.clear();
    return gcode;
}

} // namespace Slic3r
