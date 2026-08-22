/**
 * @file cut_builder.cpp
 * @brief Builds topological cuts for pretension sections.
 */

#include "cut_builder.h"
#include "../element/element.h"
#include "../solid/c3d4.h"
#include "../solid/c3d6.h"
#include "../solid/c3d8.h"
#include "../solid/c3d8r.h"
#include "../solid/c3d10.h"
#include "../solid/c3d13.h"
#include "../solid/c3d15.h"
#include "../solid/c3d20.h"
#include "../solid/c3d20r.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace fem::model {

namespace {

pretension::InterfacePair get_or_create_interface_pair(
    ModelData& model_data,
    pretension::PretensionSection& section,
    const Vec3& position,
    bool /*reverse_sides*/ = false) {
    constexpr Precision key_scale = Precision(1e9);
    const std::array<long long, 3> key = {
        std::llround(position(0) * key_scale),
        std::llround(position(1) * key_scale),
        std::llround(position(2) * key_scale)};
    const auto found = section.interface_node_cache.find(key);
    if (found != section.interface_node_cache.end()) {
        return found->second;
    }

    pretension::InterfacePair pair{
        model_data.append_node(position),
        model_data.append_node(position)};
    // Cache every geometric interface point in the canonical orientation
    // negative = side A, positive = side B. Callers cutting with a reversed
    // local axis may swap their local node assignment, but must not change the
    // globally stored constraint orientation.
    section.interface_node_cache.emplace(key, pair);
    section.interface_pairs.push_back(pair);
    return pair;
}

constexpr std::array<std::pair<Index, Index>, 12> c3d8_edges = {{
    {0, 1}, {3, 2}, {4, 5}, {7, 6},
    {0, 3}, {1, 2}, {4, 7}, {5, 6},
    {0, 4}, {1, 5}, {2, 6}, {3, 7}
}};

constexpr std::array<std::array<Index, 4>, 6> c3d8_faces = {{
    {{0, 1, 2, 3}}, {{4, 5, 6, 7}},
    {{0, 1, 5, 4}}, {{1, 2, 6, 5}},
    {{2, 3, 7, 6}}, {{3, 0, 4, 7}}
}};

using FaceKey = std::array<ID, 4>;
using TetFaceKey = std::array<ID, 3>;

constexpr std::array<std::pair<Index, Index>, 6> c3d4_edges = {{
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
}};

struct C3D4Quality {
    Precision volume = 0;
    Precision minimum_edge = 0;
    Precision maximum_edge = 0;
    Precision normalized_volume = 0;
    Precision edge_ratio = 0;
};

enum class NodePlaneSide {
    Negative,
    OnPlane,
    Positive
};

enum class TetCutType {
    None,
    VertexAligned,
    EdgeAligned,
    FaceAligned,
    OneThree,
    TwoTwo
};

struct TetCutPlan {
    ID element_id = -1;
    std::array<ID, 4> node_ids{};
    std::array<Precision, 4> distances{};
    std::array<NodePlaneSide, 4> sides{};
    TetCutType type = TetCutType::None;
    Precision characteristic_length = 0;
};

struct C3D4CutAnalysis {
    std::vector<TetCutPlan> plans;
    std::set<ID> snap_node_ids;
};

C3D4Quality c3d4_quality(
    ModelData& model_data,
    const std::array<ID, 4>& nodes) {
    std::array<Vec3, 4> positions{};
    for (Index i = 0; i < positions.size(); ++i) {
        positions[i] = model_data.positions->row_vec3(
            static_cast<Index>(nodes[i]));
    }

    C3D4Quality quality;
    quality.volume = std::abs((positions[1] - positions[0]).dot(
        (positions[2] - positions[0]).cross(positions[3] - positions[0]))) /
        Precision(6);
    quality.minimum_edge = std::numeric_limits<Precision>::max();
    for (const auto [a, b] : c3d4_edges) {
        const Precision edge = (positions[a] - positions[b]).norm();
        quality.minimum_edge = std::min(quality.minimum_edge, edge);
        quality.maximum_edge = std::max(quality.maximum_edge, edge);
    }
    if (quality.maximum_edge > std::numeric_limits<Precision>::epsilon()) {
        quality.normalized_volume = quality.volume /
            (quality.maximum_edge * quality.maximum_edge * quality.maximum_edge);
        quality.edge_ratio = quality.minimum_edge / quality.maximum_edge;
    }
    return quality;
}

template<std::size_t N>
Precision minimum_normalized_tet_volume(
    ModelData& model_data,
    const std::array<std::array<ID, 4>, N>& tets) {
    Precision result = std::numeric_limits<Precision>::max();
    for (const auto& tet : tets) {
        result = std::min(
            result, c3d4_quality(model_data, tet).normalized_volume);
    }
    return result;
}

template<std::size_t N>
Precision total_tet_volume(
    ModelData& model_data,
    const std::array<std::array<ID, 4>, N>& tets) {
    Precision result = 0;
    for (const auto& tet : tets) {
        result += c3d4_quality(model_data, tet).volume;
    }
    return result;
}

TetCutType classify_tet_cut(
    const std::array<NodePlaneSide, 4>& sides) {
    Index negative = 0;
    Index on_plane = 0;
    Index positive = 0;
    for (const NodePlaneSide side : sides) {
        negative += side == NodePlaneSide::Negative;
        on_plane += side == NodePlaneSide::OnPlane;
        positive += side == NodePlaneSide::Positive;
    }

    if (on_plane == 3) return TetCutType::FaceAligned;
    if (negative == 0 || positive == 0) return TetCutType::None;
    if (on_plane == 2) return TetCutType::EdgeAligned;
    if (on_plane == 1) return TetCutType::VertexAligned;
    if (negative == 2 && positive == 2) return TetCutType::TwoTwo;
    return TetCutType::OneThree;
}

C3D4CutAnalysis analyze_c3d4_cuts(
    ModelData& model_data,
    const Vec3& plane_point,
    const Vec3& axis,
    Precision snap_ratio) {
    C3D4CutAnalysis analysis;
    std::map<ID, unsigned> near_plane_incidence;

    // First collect snap decisions by original node ID. No topology or
    // coordinates are changed during this pass.
    for (const ElementPtr& element : model_data.elements) {
        auto* c3d4 = element == nullptr ? nullptr : dynamic_cast<C3D4*>(element.get());
        if (c3d4 == nullptr) continue;

        std::array<Precision, 4> distances{};
        bool has_negative = false;
        bool has_positive = false;
        for (Index i = 0; i < 4; ++i) {
            const Vec3 position = model_data.positions->row_vec3(
                static_cast<Index>(c3d4->node_ids[i]));
            distances[i] = axis.dot(position - plane_point);
            has_negative = has_negative || distances[i] < 0;
            has_positive = has_positive || distances[i] > 0;
        }
        const Precision characteristic_length =
            c3d4_quality(model_data, c3d4->node_ids).maximum_edge;
        const Precision snap_tolerance =
            snap_ratio * characteristic_length;
        for (Index i = 0; i < 4; ++i) {
            if (std::abs(distances[i]) <= snap_tolerance) {
                unsigned& incidence = near_plane_incidence[c3d4->node_ids[i]];
                if (has_negative) incidence |= 1u;
                if (has_positive) incidence |= 2u;
            }
        }
    }
    for (const auto& [node_id, incidence] : near_plane_incidence) {
        if (incidence == 3u) analysis.snap_node_ids.insert(node_id);
    }

    // Reclassify every original tetrahedron with the global snap decisions.
    for (const ElementPtr& element : model_data.elements) {
        auto* c3d4 = element == nullptr ? nullptr : dynamic_cast<C3D4*>(element.get());
        if (c3d4 == nullptr) continue;

        TetCutPlan plan;
        plan.element_id = c3d4->elem_id;
        plan.node_ids = c3d4->node_ids;
        plan.characteristic_length =
            c3d4_quality(model_data, c3d4->node_ids).maximum_edge;
        for (Index i = 0; i < 4; ++i) {
            const ID node_id = c3d4->node_ids[i];
            const Vec3 position = model_data.positions->row_vec3(
                static_cast<Index>(node_id));
            plan.distances[i] = axis.dot(position - plane_point);
            if (analysis.snap_node_ids.count(node_id) > 0) {
                plan.sides[i] = NodePlaneSide::OnPlane;
            } else if (plan.distances[i] < 0) {
                plan.sides[i] = NodePlaneSide::Negative;
            } else if (plan.distances[i] > 0) {
                plan.sides[i] = NodePlaneSide::Positive;
            } else {
                plan.sides[i] = NodePlaneSide::OnPlane;
            }
        }
        plan.type = classify_tet_cut(plan.sides);
        if (plan.type != TetCutType::None) {
            analysis.plans.push_back(plan);
        }
    }
    return analysis;
}

Precision tetra_volume(const std::array<Vec3, 4>& points) {
    return std::abs((points[1] - points[0]).dot(
        (points[2] - points[0]).cross(points[3] - points[0]))) /
        Precision(6);
}

Precision tetra_normalized_volume(const std::array<Vec3, 4>& points) {
    Precision maximum_edge = 0;
    for (const auto [a, b] : c3d4_edges) {
        maximum_edge = std::max(maximum_edge, (points[a] - points[b]).norm());
    }
    if (maximum_edge <= std::numeric_limits<Precision>::epsilon()) return 0;
    return tetra_volume(points) /
        (maximum_edge * maximum_edge * maximum_edge);
}

void validate_vertex_aligned_plans(
    ModelData& model_data,
    const pretension::PretensionSection& section,
    const C3D4CutAnalysis& analysis,
    const Vec3& plane_point,
    const Vec3& axis) {
    Index checked = 0;
    Index invalid = 0;
    Precision maximum_relative_volume_error = 0;
    Precision minimum_normalized_volume = std::numeric_limits<Precision>::max();

    for (const TetCutPlan& plan : analysis.plans) {
        if (plan.type != TetCutType::VertexAligned) continue;

        Index on_plane = 0;
        Index single = 0;
        std::array<Index, 2> opposite{};
        Index opposite_count = 0;
        Index positive_count = 0;
        for (Index i = 0; i < 4; ++i) {
            if (plan.sides[i] == NodePlaneSide::OnPlane) on_plane = i;
            else if (plan.sides[i] == NodePlaneSide::Positive) ++positive_count;
        }
        const NodePlaneSide single_side = positive_count == 1
            ? NodePlaneSide::Positive : NodePlaneSide::Negative;
        for (Index i = 0; i < 4; ++i) {
            if (plan.sides[i] == single_side) single = i;
            else if (plan.sides[i] != NodePlaneSide::OnPlane) {
                opposite[opposite_count++] = i;
            }
        }

        std::array<Vec3, 4> original{};
        for (Index i = 0; i < 4; ++i) {
            original[i] = model_data.positions->row_vec3(
                static_cast<Index>(plan.node_ids[i]));
        }
        const Vec3 snapped = original[on_plane] -
            axis * axis.dot(original[on_plane] - plane_point);
        original[on_plane] = snapped;

        std::array<Vec3, 2> intersections{};
        for (Index i = 0; i < 2; ++i) {
            const Index other = opposite[i];
            const Precision parameter = plan.distances[single] /
                (plan.distances[single] - plan.distances[other]);
            intersections[i] = original[single] +
                parameter * (original[other] - original[single]);
        }

        const std::array<std::array<Vec3, 4>, 3> children = {{
            {{original[single], snapped, intersections[0], intersections[1]}},
            {{original[opposite[0]], original[opposite[1]], snapped, intersections[0]}},
            {{original[opposite[1]], snapped, intersections[0], intersections[1]}}
        }};
        const Precision original_volume = tetra_volume(original);
        Precision child_volume = 0;
        bool plan_valid = original_volume > std::numeric_limits<Precision>::epsilon();
        for (const auto& child : children) {
            const Precision volume = tetra_volume(child);
            child_volume += volume;
            minimum_normalized_volume = std::min(
                minimum_normalized_volume, tetra_normalized_volume(child));
            plan_valid = plan_valid &&
                volume > std::numeric_limits<Precision>::epsilon();
        }
        const Precision relative_error = original_volume > 0
            ? std::abs(child_volume - original_volume) / original_volume
            : std::numeric_limits<Precision>::max();
        maximum_relative_volume_error = std::max(
            maximum_relative_volume_error, relative_error);
        plan_valid = plan_valid && relative_error <= Precision(1e-10);
        invalid += !plan_valid;
        ++checked;
    }

    if (checked == 0) return;
    logging::info(true,
                  "PRETENSIONSECTION '", section.name,
                  "': validated ", checked,
                  " vertex-aligned C3D4 plan(s): invalid = ", invalid,
                  ", maximum relative volume error = ",
                  maximum_relative_volume_error,
                  ", minimum normalized child volume = ",
                  minimum_normalized_volume);
}

void validate_edge_aligned_plans(
    ModelData& model_data,
    const pretension::PretensionSection& section,
    const C3D4CutAnalysis& analysis,
    const Vec3& plane_point,
    const Vec3& axis) {
    Index checked = 0;
    Index invalid = 0;
    Precision maximum_relative_volume_error = 0;
    Precision minimum_normalized_volume = std::numeric_limits<Precision>::max();

    for (const TetCutPlan& plan : analysis.plans) {
        if (plan.type != TetCutType::EdgeAligned) continue;
        std::array<Index, 2> on_plane{};
        Index on_plane_count = 0;
        Index negative = 0;
        Index positive = 0;
        std::array<Vec3, 4> original{};
        for (Index i = 0; i < 4; ++i) {
            original[i] = model_data.positions->row_vec3(
                static_cast<Index>(plan.node_ids[i]));
            if (plan.sides[i] == NodePlaneSide::OnPlane) {
                on_plane[on_plane_count++] = i;
                original[i] -= axis * axis.dot(original[i] - plane_point);
            } else if (plan.sides[i] == NodePlaneSide::Negative) negative = i;
            else positive = i;
        }
        const Precision parameter = plan.distances[negative] /
            (plan.distances[negative] - plan.distances[positive]);
        const Vec3 intersection = original[negative] +
            parameter * (original[positive] - original[negative]);
        const std::array<std::array<Vec3, 4>, 2> children = {{
            {{original[negative], original[on_plane[0]],
              original[on_plane[1]], intersection}},
            {{original[positive], original[on_plane[0]],
              original[on_plane[1]], intersection}}
        }};
        const Precision original_volume = tetra_volume(original);
        Precision child_volume = 0;
        bool plan_valid = original_volume > std::numeric_limits<Precision>::epsilon();
        for (const auto& child : children) {
            const Precision volume = tetra_volume(child);
            child_volume += volume;
            minimum_normalized_volume = std::min(
                minimum_normalized_volume, tetra_normalized_volume(child));
            plan_valid = plan_valid &&
                volume > std::numeric_limits<Precision>::epsilon();
        }
        const Precision relative_error = original_volume > 0
            ? std::abs(child_volume - original_volume) / original_volume
            : std::numeric_limits<Precision>::max();
        maximum_relative_volume_error = std::max(
            maximum_relative_volume_error, relative_error);
        plan_valid = plan_valid && relative_error <= Precision(1e-10);
        invalid += !plan_valid;
        ++checked;
    }
    if (checked == 0) return;
    logging::info(true,
                  "PRETENSIONSECTION '", section.name,
                  "': validated ", checked,
                  " edge-aligned C3D4 plan(s): invalid = ", invalid,
                  ", maximum relative volume error = ",
                  maximum_relative_volume_error,
                  ", minimum normalized child volume = ",
                  minimum_normalized_volume);
}

void log_c3d4_cut_analysis(
    const pretension::PretensionSection& section,
    const C3D4CutAnalysis& analysis) {
    std::array<Index, 6> counts{};
    for (const TetCutPlan& plan : analysis.plans) {
        ++counts[static_cast<std::size_t>(plan.type)];
    }
    logging::info(true,
                  "PRETENSIONSECTION '", section.name,
                  "': C3D4 cut analysis with snap ratio ", section.snap_ratio,
                  ": snap nodes = ", analysis.snap_node_ids.size(),
                  ", vertex-aligned = ", counts[static_cast<std::size_t>(TetCutType::VertexAligned)],
                  ", edge-aligned = ", counts[static_cast<std::size_t>(TetCutType::EdgeAligned)],
                  ", face-aligned = ", counts[static_cast<std::size_t>(TetCutType::FaceAligned)],
                  ", 1/3 = ", counts[static_cast<std::size_t>(TetCutType::OneThree)],
                  ", 2/2 = ", counts[static_cast<std::size_t>(TetCutType::TwoTwo)]);
}

FaceKey make_face_key(const std::array<ID, 8>& nodes,
                     const std::array<Index, 4>& face) {
    FaceKey key = {
        nodes[face[0]], nodes[face[1]], nodes[face[2]], nodes[face[3]]};
    std::sort(key.begin(), key.end());
    return key;
}

constexpr std::array<std::array<Index, 3>, 4> c3d4_faces = {{
    {{0, 1, 2}}, {{0, 3, 1}}, {{1, 3, 2}}, {{2, 3, 0}}
}};

TetFaceKey make_face_key(const std::array<ID, 4>& nodes,
                         const std::array<Index, 3>& face) {
    TetFaceKey key = {
        nodes[face[0]], nodes[face[1]], nodes[face[2]]};
    std::sort(key.begin(), key.end());
    return key;
}

int c3d4_side_of_element(
    ModelData& model_data,
    const std::array<ID, 4>& nodes,
    const Vec3& plane_point,
    const Vec3& axis,
    Precision tolerance) {
    bool has_negative = false;
    bool has_positive = false;
    for (const ID node_id : nodes) {
        const Vec3 position = model_data.positions->row_vec3(
            static_cast<Index>(node_id));
        const Precision distance = axis.dot(position - plane_point);
        has_negative = has_negative || distance < -tolerance;
        has_positive = has_positive || distance > tolerance;
    }
    if (has_negative && !has_positive) return -1;
    if (has_positive && !has_negative) return 1;
    return 0;
}

int c3d8_side_of_element(
    ModelData& model_data,
    const std::array<ID, 8>& nodes,
    const Vec3& plane_point,
    const Vec3& axis,
    Precision tolerance) {
    bool has_negative = false;
    bool has_positive = false;
    for (const ID node_id : nodes) {
        const Vec3 position = model_data.positions->row_vec3(static_cast<Index>(node_id));
        const Precision distance = axis.dot(position - plane_point);
        has_negative = has_negative || distance < -tolerance;
        has_positive = has_positive || distance > tolerance;
    }
    if (has_negative && !has_positive) return -1;
    if (has_positive && !has_negative) return 1;
    return 0;
}

void split_face_aligned_c3d8(
    ModelData& model_data,
    pretension::PretensionSection& section,
    const Vec3& plane_point,
    const Vec3& axis,
    Precision tolerance) {
    struct FaceRecord {
        ID element_id = -1;
        int side = 0;
        FaceKey nodes{};
    };

    std::map<FaceKey, std::vector<FaceRecord>> faces;
    for (const ElementPtr& element : model_data.elements) {
        auto* c3d8 = element == nullptr ? nullptr : dynamic_cast<C3D8*>(element.get());
        if (c3d8 == nullptr) continue;

        const auto nodes = c3d8->node_ids;
        const int side = c3d8_side_of_element(
            model_data, nodes, plane_point, axis, tolerance);
        if (side == 0) continue;

        for (const auto& face : c3d8_faces) {
            bool on_plane = true;
            for (const Index local_node : face) {
                const Vec3 position = model_data.positions->row_vec3(
                    static_cast<Index>(nodes[local_node]));
                const Precision distance = axis.dot(position - plane_point);
                if (std::abs(distance) > tolerance) {
                    on_plane = false;
                    break;
                }
            }
            if (on_plane) {
                const FaceKey key = make_face_key(nodes, face);
                faces[key].push_back({c3d8->elem_id, side, key});
            }
        }
    }

    std::set<ID> interface_nodes;
    for (const auto& [key, records] : faces) {
        bool has_negative = false;
        bool has_positive = false;
        for (const auto& record : records) {
            has_negative = has_negative || record.side < 0;
            has_positive = has_positive || record.side > 0;
        }
        if (!has_negative || !has_positive) continue;

        for (const ID node_id : key) {
            const Vec3 position = model_data.positions->row_vec3(
                static_cast<Index>(node_id));
            const auto pair = get_or_create_interface_pair(model_data, section, position);
            interface_nodes.insert(node_id);
            section.side_a_nodes.push_back(pair.side_a);
            section.side_b_nodes.push_back(pair.side_b);
        }
    }

    if (interface_nodes.empty()) return;

    for (const ElementPtr& element : model_data.elements) {
        auto* c3d8 = element == nullptr ? nullptr : dynamic_cast<C3D8*>(element.get());
        if (c3d8 == nullptr) continue;
        if (c3d8_side_of_element(
                model_data, c3d8->node_ids, plane_point, axis, tolerance) <= 0) {
            continue;
        }

        for (const ID node_id : interface_nodes) {
            c3d8->replace_node(
                node_id,
                section.interface_node_cache.at({
                    std::llround(model_data.positions->row_vec3(static_cast<Index>(node_id))(0) * Precision(1e9)),
                    std::llround(model_data.positions->row_vec3(static_cast<Index>(node_id))(1) * Precision(1e9)),
                    std::llround(model_data.positions->row_vec3(static_cast<Index>(node_id))(2) * Precision(1e9))}).side_b);
        }
    }
}

void split_face_aligned_c3d4(
    ModelData& model_data,
    pretension::PretensionSection& section,
    const Vec3& plane_point,
    const Vec3& axis,
    Precision tolerance) {
    struct FaceRecord {
        int side = 0;
    };

    std::map<TetFaceKey, std::vector<FaceRecord>> faces;
    for (const ElementPtr& element : model_data.elements) {
        auto* c3d4 = element == nullptr ? nullptr : dynamic_cast<C3D4*>(element.get());
        if (c3d4 == nullptr) continue;

        const auto nodes = c3d4->node_ids;
        const int side = c3d4_side_of_element(
            model_data, nodes, plane_point, axis, tolerance);
        if (side == 0) continue;

        for (const auto& face : c3d4_faces) {
            bool on_plane = true;
            for (const Index local_node : face) {
                const Vec3 position = model_data.positions->row_vec3(
                    static_cast<Index>(nodes[local_node]));
                if (std::abs(axis.dot(position - plane_point)) > tolerance) {
                    on_plane = false;
                    break;
                }
            }
            if (on_plane) {
                faces[make_face_key(nodes, face)].push_back({side});
            }
        }
    }

    std::map<ID, pretension::InterfacePair> node_pairs;
    for (const auto& [face, records] : faces) {
        bool has_negative = false;
        bool has_positive = false;
        for (const auto& record : records) {
            has_negative = has_negative || record.side < 0;
            has_positive = has_positive || record.side > 0;
        }
        if (!has_negative || !has_positive) continue;

        for (const ID node_id : face) {
            if (node_pairs.find(node_id) != node_pairs.end()) continue;
            const Vec3 position = model_data.positions->row_vec3(
                static_cast<Index>(node_id));
            const auto pair = get_or_create_interface_pair(
                model_data, section, position);
            node_pairs.emplace(node_id, pair);
            section.side_a_nodes.push_back(pair.side_a);
            section.side_b_nodes.push_back(pair.side_b);
        }
    }

    if (node_pairs.empty()) return;

    for (const ElementPtr& element : model_data.elements) {
        auto* c3d4 = element == nullptr ? nullptr : dynamic_cast<C3D4*>(element.get());
        if (c3d4 == nullptr || c3d4_side_of_element(
                model_data, c3d4->node_ids, plane_point, axis, tolerance) <= 0) {
            continue;
        }
        for (const auto& [node_id, pair] : node_pairs) {
            c3d4->replace_node(node_id, pair.side_b);
        }
    }
}

void log_c3d4_interface_quality(
    ModelData& model_data,
    const pretension::PretensionSection& section) {
    std::set<ID> interface_nodes;
    for (const auto& pair : section.interface_pairs) {
        interface_nodes.insert(pair.side_a);
        interface_nodes.insert(pair.side_b);
    }

    Precision minimum_volume = std::numeric_limits<Precision>::max();
    Precision minimum_edge = std::numeric_limits<Precision>::max();
    Precision minimum_normalized_volume = std::numeric_limits<Precision>::max();
    Precision minimum_edge_ratio = std::numeric_limits<Precision>::max();
    ID minimum_volume_element = -1;
    ID minimum_edge_element = -1;
    Index checked_elements = 0;
    Index failed_elements = 0;
    constexpr Precision normalized_volume_warning = Precision(1e-4);
    constexpr Precision edge_ratio_warning = Precision(0.02);

    for (const ElementPtr& element : model_data.elements) {
        auto* c3d4 = element == nullptr ? nullptr : dynamic_cast<C3D4*>(element.get());
        if (c3d4 == nullptr) continue;

        bool touches_interface = false;
        for (const ID node_id : c3d4->node_ids) {
            touches_interface = touches_interface || interface_nodes.count(node_id) > 0;
        }
        if (!touches_interface) continue;

        const auto& nodes = c3d4->node_ids;
        const C3D4Quality quality = c3d4_quality(model_data, nodes);
        const Precision volume = quality.volume;
        minimum_volume = std::min(minimum_volume, volume);
        if (volume == minimum_volume) minimum_volume_element = c3d4->elem_id;
        if (quality.minimum_edge < minimum_edge) {
            minimum_edge = quality.minimum_edge;
            minimum_edge_element = c3d4->elem_id;
        }
        minimum_normalized_volume = std::min(
            minimum_normalized_volume, quality.normalized_volume);
        minimum_edge_ratio = std::min(minimum_edge_ratio, quality.edge_ratio);
        if (quality.normalized_volume < normalized_volume_warning ||
            quality.edge_ratio < edge_ratio_warning) {
            const auto origin = section.cut_quality_origins.find(c3d4->elem_id);
            const auto type = section.cut_quality_types.find(c3d4->elem_id);
            logging::warning(false,
                             "PRETENSIONSECTION '", section.name,
                             "': C3D4 quality gate failed for element ",
                             c3d4->elem_id, " (origin ",
                             origin == section.cut_quality_origins.end()
                                 ? c3d4->elem_id : origin->second,
                             ", cut type ",
                             type == section.cut_quality_types.end()
                                 ? std::string("unknown") : type->second,
                             "): normalized volume = ",
                             quality.normalized_volume,
                             ", edge ratio = ", quality.edge_ratio);
            ++failed_elements;
        }
        ++checked_elements;
    }

    if (checked_elements == 0) return;
    logging::info(true,
                  "PRETENSIONSECTION '", section.name,
                  "': C3D4 interface quality over ", checked_elements,
                  " elements: minimum volume = ", minimum_volume,
                  " (element ", minimum_volume_element,
                  "), minimum edge = ", minimum_edge,
                  " (element ", minimum_edge_element,
                  "), minimum normalized volume = ", minimum_normalized_volume,
                  ", minimum edge ratio = ", minimum_edge_ratio,
                  ", quality-gate failures = ", failed_elements);
}

void log_c3d10_interface_quality(
    ModelData& model_data,
    const pretension::PretensionSection& section) {
    std::set<ID> interface_nodes;
    for (const auto& pair : section.interface_pairs) {
        interface_nodes.insert(pair.side_a);
        interface_nodes.insert(pair.side_b);
    }
    Precision minimum_normalized_volume = std::numeric_limits<Precision>::max();
    Precision minimum_edge_ratio = std::numeric_limits<Precision>::max();
    ID minimum_element = -1;
    Index checked = 0;
    Index failed = 0;
    for (const ElementPtr& element : model_data.elements) {
        auto* c3d10 = element == nullptr ? nullptr : dynamic_cast<C3D10*>(element.get());
        if (c3d10 == nullptr) continue;
        bool touches_interface = false;
        for (const ID node_id : c3d10->node_ids) {
            touches_interface = touches_interface || interface_nodes.count(node_id) > 0;
        }
        if (!touches_interface) continue;
        const std::array<ID, 4> vertices = {
            c3d10->node_ids[0], c3d10->node_ids[1],
            c3d10->node_ids[2], c3d10->node_ids[3]};
        const C3D4Quality quality = c3d4_quality(model_data, vertices);
        if (quality.normalized_volume < minimum_normalized_volume) {
            minimum_normalized_volume = quality.normalized_volume;
            minimum_element = c3d10->elem_id;
        }
        minimum_edge_ratio = std::min(minimum_edge_ratio, quality.edge_ratio);
        if (quality.normalized_volume < Precision(1e-4) ||
            quality.edge_ratio < Precision(0.02)) {
            logging::warning(false,
                             "PRETENSIONSECTION '", section.name,
                             "': C3D10 corner quality gate failed for element ",
                             c3d10->elem_id, ": normalized volume = ",
                             quality.normalized_volume,
                             ", edge ratio = ", quality.edge_ratio);
            ++failed;
        }
        ++checked;
    }
    if (checked == 0) return;
    logging::info(true,
                  "PRETENSIONSECTION '", section.name,
                  "': C3D10 interface corner quality over ", checked,
                  " elements: minimum normalized volume = ",
                  minimum_normalized_volume, " (element ", minimum_element,
                  "), minimum edge ratio = ", minimum_edge_ratio,
                  ", quality-gate failures = ", failed);
}

void validate_interface_pair_connectivity(
    const ModelData& model_data,
    const pretension::PretensionSection& section) {
    std::map<ID, Index> use_count;
    std::set<ID> side_a_nodes;
    std::set<ID> side_b_nodes;
    for (const auto& pair : section.interface_pairs) {
        side_a_nodes.insert(pair.side_a);
        side_b_nodes.insert(pair.side_b);
    }
    Index mixed_elements = 0;
    for (const ElementPtr& element : model_data.elements) {
        if (element == nullptr) continue;
        bool uses_side_a = false;
        bool uses_side_b = false;
        for (Index local_node = 0;
             local_node < static_cast<Index>(element->n_nodes()); ++local_node) {
            const ID node_id = element->nodes()[local_node];
            ++use_count[node_id];
            uses_side_a = uses_side_a || side_a_nodes.count(node_id) > 0;
            uses_side_b = uses_side_b || side_b_nodes.count(node_id) > 0;
        }
        if (uses_side_a && uses_side_b) {
            logging::warning(false,
                             "PRETENSIONSECTION '", section.name,
                             "': element ", element->elem_id,
                             " connects both interface sides");
            ++mixed_elements;
        }
    }
    Index disconnected = 0;
    for (const auto& pair : section.interface_pairs) {
        const Index side_a_uses = use_count[pair.side_a];
        const Index side_b_uses = use_count[pair.side_b];
        if (side_a_uses == 0 || side_b_uses == 0) {
            logging::warning(false,
                             "PRETENSIONSECTION '", section.name,
                             "': disconnected interface pair ", pair.side_a,
                             "/", pair.side_b, " has element-use counts ",
                             side_a_uses, "/", side_b_uses);
            ++disconnected;
        }
    }
    logging::error(disconnected == 0,
                   "CutBuilder: ", disconnected,
                   " pretension interface pair(s) are not connected on both sides");
    logging::error(mixed_elements == 0,
                   "CutBuilder: ", mixed_elements,
                   " element(s) connect both pretension interface sides");

}

std::vector<std::pair<ID, ID>> find_crossing_edges(
    ModelData& model_data,
    const std::array<ID, 8>& node_ids,
    const Vec3& plane_point,
    const Vec3& axis,
    Precision tolerance) {
    std::vector<std::pair<ID, ID>> crossing_edges;
    for (const auto [local_a, local_b] : c3d8_edges) {
        const Vec3 position_a = model_data.positions->row_vec3(
            static_cast<Index>(node_ids[local_a]));
        const Vec3 position_b = model_data.positions->row_vec3(
            static_cast<Index>(node_ids[local_b]));
        const Precision distance_a = axis.dot(position_a - plane_point);
        const Precision distance_b = axis.dot(position_b - plane_point);

        if ((distance_a < -tolerance && distance_b > tolerance) ||
            (distance_a > tolerance && distance_b < -tolerance)) {
            crossing_edges.emplace_back(node_ids[local_a], node_ids[local_b]);
        }
    }
    return crossing_edges;
}

void add_element_to_matching_sets(
    ModelData& model_data,
    ID original_element_id,
    ID new_element_id) {
    for (auto& entry : model_data.elem_sets) {
        const auto& set = entry.second;
        if (set == nullptr) {
            continue;
        }

        const bool contains_original = std::find(
            set->begin(), set->end(), original_element_id) != set->end();
        if (contains_original) {
            set->add(new_element_id);
        }
    }

    // Section regions are compiled snapshots of the element sets.  Pretension
    // cuts happen after compilation, so extending only elem_sets would leave
    // generated children without the section of their parent element.
    for (const auto& section : model_data.sections) {
        if (section == nullptr || section->region_ == nullptr) continue;
        const auto& region = section->region_;
        const bool contains_original = std::find(
            region->begin(), region->end(), original_element_id) != region->end();
        if (contains_original) region->add(new_element_id);
    }
}

std::array<ID, 4> oriented_c3d4_nodes(
    ModelData& model_data,
    std::array<ID, 4> nodes);

void replace_c3d4_node(C3D4& element, ID old_node, ID new_node) {
    for (ID& node_id : element.node_ids) {
        if (node_id == old_node) node_id = new_node;
    }
}

bool apply_aligned_c3d4_plans(
    ModelData& model_data,
    pretension::PretensionSection& section,
    const C3D4CutAnalysis& analysis,
    const Vec3& plane_point,
    const Vec3& axis) {
    if (analysis.snap_node_ids.empty()) return false;

    std::map<ID, pretension::InterfacePair> snapped_pairs;
    for (const ID node_id : analysis.snap_node_ids) {
        Vec3 position = model_data.positions->row_vec3(static_cast<Index>(node_id));
        position -= axis * axis.dot(position - plane_point);
        for (Index component = 0; component < 3; ++component) {
            (*model_data.positions)(static_cast<Index>(node_id), component) =
                position(component);
        }
        const auto pair = get_or_create_interface_pair(model_data, section, position);
        snapped_pairs.emplace(node_id, pair);
        section.side_a_nodes.push_back(pair.side_a);
        section.side_b_nodes.push_back(pair.side_b);
    }

    std::set<ID> aligned_elements;
    for (const TetCutPlan& plan : analysis.plans) {
        if (plan.type == TetCutType::VertexAligned ||
            plan.type == TetCutType::EdgeAligned) {
            aligned_elements.insert(plan.element_id);
        }
    }

    // Reconnect every unsplit incident tetrahedron to exactly one side of the
    // duplicated snapped node. Mixed-side elements are handled below.
    for (const ElementPtr& element : model_data.elements) {
        auto* c3d4 = element == nullptr ? nullptr : dynamic_cast<C3D4*>(element.get());
        if (c3d4 == nullptr || aligned_elements.count(c3d4->elem_id) > 0) continue;
        bool has_negative = false;
        bool has_positive = false;
        for (const ID node_id : c3d4->node_ids) {
            if (analysis.snap_node_ids.count(node_id) > 0) continue;
            const Precision distance = axis.dot(
                model_data.positions->row_vec3(static_cast<Index>(node_id)) - plane_point);
            has_negative = has_negative || distance < 0;
            has_positive = has_positive || distance > 0;
        }
        if (has_negative == has_positive) continue;
        for (const auto& [old_node, pair] : snapped_pairs) {
            replace_c3d4_node(*c3d4, old_node,
                has_negative ? pair.side_a : pair.side_b);
        }
    }

    for (const TetCutPlan& plan : analysis.plans) {
        if (plan.type != TetCutType::VertexAligned) continue;
        auto* original_element = dynamic_cast<C3D4*>(
            model_data.elements[static_cast<std::size_t>(plan.element_id)].get());
        logging::error(original_element != nullptr,
                       "CutBuilder: missing vertex-aligned C3D4 element");

        Index on_plane = 0;
        Index single = 0;
        std::array<Index, 2> opposite{};
        Index opposite_count = 0;
        Index positive_count = 0;
        for (Index i = 0; i < 4; ++i) {
            if (plan.sides[i] == NodePlaneSide::OnPlane) on_plane = i;
            else if (plan.sides[i] == NodePlaneSide::Positive) ++positive_count;
        }
        const NodePlaneSide single_side = positive_count == 1
            ? NodePlaneSide::Positive : NodePlaneSide::Negative;
        for (Index i = 0; i < 4; ++i) {
            if (plan.sides[i] == single_side) single = i;
            else if (plan.sides[i] != NodePlaneSide::OnPlane) {
                opposite[opposite_count++] = i;
            }
        }

        const auto snapped_pair = snapped_pairs.at(plan.node_ids[on_plane]);
        std::array<pretension::InterfacePair, 2> intersection_pairs{};
        for (Index i = 0; i < 2; ++i) {
            const Vec3 single_position = model_data.positions->row_vec3(
                static_cast<Index>(plan.node_ids[single]));
            const Vec3 opposite_position = model_data.positions->row_vec3(
                static_cast<Index>(plan.node_ids[opposite[i]]));
            const Precision parameter = plan.distances[single] /
                (plan.distances[single] - plan.distances[opposite[i]]);
            intersection_pairs[i] = get_or_create_interface_pair(
                model_data, section,
                single_position + parameter * (opposite_position - single_position));
            section.side_a_nodes.push_back(intersection_pairs[i].side_a);
            section.side_b_nodes.push_back(intersection_pairs[i].side_b);
        }

        const bool single_positive = single_side == NodePlaneSide::Positive;
        const ID single_interface = single_positive
            ? snapped_pair.side_b : snapped_pair.side_a;
        const ID opposite_interface = single_positive
            ? snapped_pair.side_a : snapped_pair.side_b;
        const auto interface_for_single = [&](Index i) {
            return single_positive ? intersection_pairs[i].side_b
                                   : intersection_pairs[i].side_a;
        };
        const auto interface_for_opposite = [&](Index i) {
            return single_positive ? intersection_pairs[i].side_a
                                   : intersection_pairs[i].side_b;
        };

        std::array<std::array<ID, 4>, 3> children = {{
            {plan.node_ids[single], single_interface,
             interface_for_single(0), interface_for_single(1)},
            {plan.node_ids[opposite[0]], plan.node_ids[opposite[1]],
             opposite_interface, interface_for_opposite(0)},
            {plan.node_ids[opposite[1]], opposite_interface,
             interface_for_opposite(0), interface_for_opposite(1)}
        }};
        for (auto& child : children) child = oriented_c3d4_nodes(model_data, child);

        auto first = std::make_shared<C3D4>(plan.element_id, children[0]);
        first->_model_data = &model_data;
        model_data.elements[static_cast<std::size_t>(plan.element_id)] = first;
        section.cut_quality_origins[plan.element_id] = plan.element_id;
        section.cut_quality_types[plan.element_id] = "vertex-aligned";
        for (Index i = 1; i < children.size(); ++i) {
            const ID new_id = model_data.next_free_element_id();
            add_element_to_matching_sets(model_data, plan.element_id, new_id);
            model_data.insert_element(std::make_shared<C3D4>(new_id, children[i]));
            section.cut_quality_origins[new_id] = plan.element_id;
            section.cut_quality_types[new_id] = "vertex-aligned";
        }
    }


    for (const TetCutPlan& plan : analysis.plans) {
        if (plan.type != TetCutType::EdgeAligned) continue;
        std::array<Index, 2> on_plane{};
        Index on_plane_count = 0;
        Index negative = 0;
        Index positive = 0;
        for (Index i = 0; i < 4; ++i) {
            if (plan.sides[i] == NodePlaneSide::OnPlane) {
                on_plane[on_plane_count++] = i;
            } else if (plan.sides[i] == NodePlaneSide::Negative) negative = i;
            else positive = i;
        }
        const Vec3 negative_position = model_data.positions->row_vec3(
            static_cast<Index>(plan.node_ids[negative]));
        const Vec3 positive_position = model_data.positions->row_vec3(
            static_cast<Index>(plan.node_ids[positive]));
        const Precision parameter = plan.distances[negative] /
            (plan.distances[negative] - plan.distances[positive]);
        const auto crossing_pair = get_or_create_interface_pair(
            model_data, section,
            negative_position + parameter * (positive_position - negative_position));
        section.side_a_nodes.push_back(crossing_pair.side_a);
        section.side_b_nodes.push_back(crossing_pair.side_b);

        const auto edge_pair_0 = snapped_pairs.at(plan.node_ids[on_plane[0]]);
        const auto edge_pair_1 = snapped_pairs.at(plan.node_ids[on_plane[1]]);
        auto negative_tet = oriented_c3d4_nodes(model_data, {
            plan.node_ids[negative], edge_pair_0.side_a,
            edge_pair_1.side_a, crossing_pair.side_a});
        auto positive_tet = oriented_c3d4_nodes(model_data, {
            plan.node_ids[positive], edge_pair_0.side_b,
            edge_pair_1.side_b, crossing_pair.side_b});

        auto first = std::make_shared<C3D4>(plan.element_id, negative_tet);
        first->_model_data = &model_data;
        model_data.elements[static_cast<std::size_t>(plan.element_id)] = first;
        section.cut_quality_origins[plan.element_id] = plan.element_id;
        section.cut_quality_types[plan.element_id] = "edge-aligned";
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, plan.element_id, new_id);
        model_data.insert_element(std::make_shared<C3D4>(new_id, positive_tet));
        section.cut_quality_origins[new_id] = plan.element_id;
        section.cut_quality_types[new_id] = "edge-aligned";
    }
    return true;
}

void split_c3d8_x_plane(
    ModelData& model_data,
    pretension::PretensionSection& section,
    C3D8& element,
    const Vec3& plane_point,
    const Vec3& axis) {
    const ID element_id = element.elem_id;
    const bool reduced_integration = dynamic_cast<C3D8R*>(&element) != nullptr;
    auto old_nodes = element.node_ids;
    const bool x_positive = axis.isApprox(Vec3::UnitX(), Precision(1e-12));
    const bool x_negative = axis.isApprox(-Vec3::UnitX(), Precision(1e-12));
    const bool y_positive = axis.isApprox(Vec3::UnitY(), Precision(1e-12));
    const bool y_negative = axis.isApprox(-Vec3::UnitY(), Precision(1e-12));
    const bool z_positive = axis.isApprox(Vec3::UnitZ(), Precision(1e-12));
    const bool z_negative = axis.isApprox(-Vec3::UnitZ(), Precision(1e-12));
    const bool x_axis = x_positive || x_negative;
    const bool y_axis = y_positive || y_negative;
    const bool z_axis = z_positive || z_negative;
    logging::error(x_axis || y_axis || z_axis,
                   "CutBuilder: only global X, Y or Z axes are supported yet");
    const bool reverse_sides = x_negative || y_negative || z_negative;
    const Vec3 split_axis = x_axis ? Vec3::UnitX() :
                            (y_axis ? Vec3::UnitY() : Vec3::UnitZ());

    std::array<ID, 4> negative_nodes{};
    std::array<ID, 4> positive_nodes{};
    std::array<std::pair<ID, ID>, 4> crossing_edges{};

    const auto configure_groups = [&]() {
        if (x_axis) {
            negative_nodes = {old_nodes[0], old_nodes[3], old_nodes[4], old_nodes[7]};
            positive_nodes = {old_nodes[1], old_nodes[2], old_nodes[5], old_nodes[6]};
            crossing_edges = {{{old_nodes[0], old_nodes[1]},
                               {old_nodes[3], old_nodes[2]},
                               {old_nodes[4], old_nodes[5]},
                               {old_nodes[7], old_nodes[6]}}};
        } else if (y_axis) {
            negative_nodes = {old_nodes[0], old_nodes[1], old_nodes[5], old_nodes[4]};
            positive_nodes = {old_nodes[3], old_nodes[2], old_nodes[6], old_nodes[7]};
            crossing_edges = {{{old_nodes[0], old_nodes[3]},
                               {old_nodes[1], old_nodes[2]},
                               {old_nodes[5], old_nodes[6]},
                               {old_nodes[4], old_nodes[7]}}};
        } else {
            negative_nodes = {old_nodes[0], old_nodes[1], old_nodes[2], old_nodes[3]};
            positive_nodes = {old_nodes[4], old_nodes[5], old_nodes[6], old_nodes[7]};
            crossing_edges = {{{old_nodes[0], old_nodes[4]},
                               {old_nodes[1], old_nodes[5]},
                               {old_nodes[2], old_nodes[6]},
                               {old_nodes[3], old_nodes[7]}}};
        }
    };
    configure_groups();

    const auto group_is_on_expected_side = [&](const auto& nodes, bool negative) {
        for (const ID node_id : nodes) {
            const Vec3 position = model_data.positions->row_vec3(static_cast<Index>(node_id));
            const Precision distance = split_axis.dot(position - plane_point);
            if (negative ? !(distance < 0) : !(distance > 0)) {
                return false;
            }
        }
        return true;
    };

    // Gmsh may orient a valid C3D8 with the two axial node layers reversed.
    // Reorder only the local working copy so the split logic always receives
    // the canonical negative-side/positive-side connectivity.
    const bool canonical = group_is_on_expected_side(negative_nodes, true) &&
                           group_is_on_expected_side(positive_nodes, false);
    const bool reversed = group_is_on_expected_side(negative_nodes, false) &&
                          group_is_on_expected_side(positive_nodes, true);
    if (!canonical && reversed) {
        if (x_axis) {
            old_nodes = {old_nodes[1], old_nodes[0], old_nodes[3], old_nodes[2],
                         old_nodes[5], old_nodes[4], old_nodes[7], old_nodes[6]};
        } else if (y_axis) {
            old_nodes = {old_nodes[3], old_nodes[2], old_nodes[1], old_nodes[0],
                         old_nodes[7], old_nodes[6], old_nodes[5], old_nodes[4]};
        } else {
            // Swap the z-layers and reverse their in-plane winding.  The
            // winding reversal preserves the positive C3D8 Jacobian while
            // putting the lower layer into slots 0..3.
            old_nodes = {old_nodes[4], old_nodes[7], old_nodes[6], old_nodes[5],
                         old_nodes[0], old_nodes[3], old_nodes[2], old_nodes[1]};
        }
        configure_groups();
    }

    for (ID node_id : negative_nodes) {
        const Vec3 position = model_data.positions->row_vec3(static_cast<Index>(node_id));
        logging::error(split_axis.dot(position - plane_point) < 0,
                       "CutBuilder: C3D8 node ordering is not compatible with the cut plane");
    }
    for (ID node_id : positive_nodes) {
        const Vec3 position = model_data.positions->row_vec3(static_cast<Index>(node_id));
        logging::error(split_axis.dot(position - plane_point) > 0,
                       "CutBuilder: C3D8 node ordering is not compatible with the cut plane");
    }

    const auto detected_edges = find_crossing_edges(
        model_data, old_nodes, plane_point, split_axis,
        Precision(100) * std::numeric_limits<Precision>::epsilon());
    logging::error(detected_edges.size() == 4,
                   "CutBuilder: C3D8 cut must intersect exactly four edges");

    std::array<ID, 4> side_a_interface{};
    std::array<ID, 4> side_b_interface{};

    for (Index i = 0; i < crossing_edges.size(); ++i) {
        const auto [node_a, node_b] = crossing_edges[i];
        const Vec3 position_a = model_data.positions->row_vec3(static_cast<Index>(node_a));
        const Vec3 position_b = model_data.positions->row_vec3(static_cast<Index>(node_b));
        const Precision distance_a = split_axis.dot(position_a - plane_point);
        const Precision distance_b = split_axis.dot(position_b - plane_point);
        const Precision denominator = distance_a - distance_b;
        logging::error(std::abs(denominator) > std::numeric_limits<Precision>::epsilon(),
                       "CutBuilder: crossing edge is parallel to the cut plane");
        const Precision parameter = distance_a / denominator;
        const Vec3 intersection = position_a + parameter * (position_b - position_a);

        const auto pair = get_or_create_interface_pair(model_data, section, intersection);
        side_a_interface[i] = reverse_sides ? pair.side_b : pair.side_a;
        side_b_interface[i] = reverse_sides ? pair.side_a : pair.side_b;
    }

    std::array<ID, 8> left_nodes{};
    std::array<ID, 8> right_nodes{};
    if (x_axis) {
        left_nodes = {old_nodes[0], side_a_interface[0], side_a_interface[1], old_nodes[3],
                      old_nodes[4], side_a_interface[2], side_a_interface[3], old_nodes[7]};
        right_nodes = {side_b_interface[0], old_nodes[1], old_nodes[2], side_b_interface[1],
                       side_b_interface[2], old_nodes[5], old_nodes[6], side_b_interface[3]};
    } else if (y_axis) {
        left_nodes = {old_nodes[0], old_nodes[1], side_a_interface[1], side_a_interface[0],
                      old_nodes[4], old_nodes[5], side_a_interface[2], side_a_interface[3]};
        right_nodes = {side_b_interface[0], side_b_interface[1], old_nodes[2], old_nodes[3],
                       side_b_interface[3], side_b_interface[2], old_nodes[6], old_nodes[7]};
    } else {
        left_nodes = {old_nodes[0], old_nodes[1], old_nodes[2], old_nodes[3],
                      side_a_interface[0], side_a_interface[1], side_a_interface[2], side_a_interface[3]};
        right_nodes = {side_b_interface[0], side_b_interface[1], side_b_interface[2], side_b_interface[3],
                       old_nodes[4], old_nodes[5], old_nodes[6], old_nodes[7]};
    }

    const ID new_element_id = model_data.next_free_element_id();
    add_element_to_matching_sets(model_data, element_id, new_element_id);

    ElementPtr left_element;
    if (reduced_integration) {
        left_element = std::make_shared<C3D8R>(element_id, left_nodes);
    } else {
        left_element = std::make_shared<C3D8>(element_id, left_nodes);
    }
    left_element->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = left_element;

    ElementPtr right_element;
    if (reduced_integration) {
        right_element = std::make_shared<C3D8R>(new_element_id, right_nodes);
    } else {
        right_element = std::make_shared<C3D8>(new_element_id, right_nodes);
    }
    model_data.insert_element(std::move(right_element));
    logging::info(reduced_integration,
                  "PRETENSIONSECTION '", section.name,
                  "': preserved C3D8R reduced integration for split element ",
                  element_id, " and child ", new_element_id);

    const auto& section_side_a = reverse_sides ? side_b_interface : side_a_interface;
    const auto& section_side_b = reverse_sides ? side_a_interface : side_b_interface;
    section.side_a_nodes.insert(
        section.side_a_nodes.end(), section_side_a.begin(), section_side_a.end());
    section.side_b_nodes.insert(
        section.side_b_nodes.end(), section_side_b.begin(), section_side_b.end());
}

std::array<ID, 4> oriented_c3d4_nodes(
    ModelData& model_data,
    std::array<ID, 4> nodes) {
    const Vec3 p0 = model_data.positions->row_vec3(static_cast<Index>(nodes[0]));
    const Vec3 p1 = model_data.positions->row_vec3(static_cast<Index>(nodes[1]));
    const Vec3 p2 = model_data.positions->row_vec3(static_cast<Index>(nodes[2]));
    const Vec3 p3 = model_data.positions->row_vec3(static_cast<Index>(nodes[3]));
    const Precision volume6 = (p1 - p0).dot((p2 - p0).cross(p3 - p0));
    if (volume6 < 0) {
        std::swap(nodes[2], nodes[3]);
    }
    logging::error(std::abs(volume6) > std::numeric_limits<Precision>::epsilon(),
                   "CutBuilder: generated C3D4 has zero volume (nodes ",
                   nodes[0], ", ", nodes[1], ", ", nodes[2], ", ", nodes[3], ")");
    return nodes;
}

std::vector<ID> tetrahedralize_c3d8(
    ModelData& model_data,
    C3D8& element,
    std::map<FaceKey, std::pair<ID, ID>>& face_diagonal_cache) {
    const ID element_id = element.elem_id;
    const auto& n = element.node_ids;
    constexpr std::array<std::pair<Index, Index>, 4> body_diagonals = {{
        {0, 6}, {1, 7}, {2, 4}, {3, 5}
    }};
    constexpr std::array<std::array<int, 3>, 8> bits = {{
        {{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{0, 1, 0}},
        {{0, 0, 1}}, {{1, 0, 1}}, {{1, 1, 1}}, {{0, 1, 1}}
    }};
    const auto adjacent = [&](Index a, Index b) {
        return std::abs(bits[a][0] - bits[b][0]) +
               std::abs(bits[a][1] - bits[b][1]) +
               std::abs(bits[a][2] - bits[b][2]) == 1;
    };
    const auto diagonal = [](ID a, ID b) {
        return std::pair<ID, ID>{std::min(a, b), std::max(a, b)};
    };
    struct Candidate {
        std::array<std::array<ID, 4>, 6> tets{};
        std::map<FaceKey, std::pair<ID, ID>> diagonals;
        Precision score = 0;
        bool compatible = true;
    };
    std::array<Candidate, 4> candidates{};
    for (Index candidate_index = 0;
         candidate_index < candidates.size(); ++candidate_index) {
        const auto [a, b] = body_diagonals[candidate_index];
        std::array<Index, 6> ring{};
        Index ring_count = 0;
        for (Index node = 0; node < 8; ++node) {
            if (node != a && node != b) ring[ring_count++] = node;
        }
        std::array<Index, 6> cycle{};
        cycle[0] = ring[0];
        for (Index i = 1; i < cycle.size(); ++i) {
            for (const Index node : ring) {
                bool used = false;
                for (Index j = 0; j < i; ++j) used = used || cycle[j] == node;
                if (!used && adjacent(cycle[i - 1], node)) {
                    cycle[i] = node;
                    break;
                }
            }
        }
        Candidate& candidate = candidates[candidate_index];
        for (Index i = 0; i < cycle.size(); ++i) {
            candidate.tets[i] = oriented_c3d4_nodes(
                model_data, {n[a], n[b], n[cycle[i]],
                             n[cycle[(i + 1) % cycle.size()]]});
        }
        for (const auto& face : c3d8_faces) {
            const FaceKey key = make_face_key(n, face);
            const auto first = diagonal(n[face[0]], n[face[2]]);
            const auto second = diagonal(n[face[1]], n[face[3]]);
            bool uses_first = false;
            for (const auto& tet : candidate.tets) {
                uses_first = uses_first ||
                    (std::find(tet.begin(), tet.end(), first.first) != tet.end() &&
                     std::find(tet.begin(), tet.end(), first.second) != tet.end());
            }
            candidate.diagonals.emplace(key, uses_first ? first : second);
        }
        candidate.score = minimum_normalized_tet_volume(
            model_data, candidate.tets);
        for (const auto& [face, selected_diagonal] : candidate.diagonals) {
            const auto cached = face_diagonal_cache.find(face);
            if (cached != face_diagonal_cache.end() &&
                cached->second != selected_diagonal) {
                candidate.compatible = false;
            }
        }
    }

    const Candidate* selected = nullptr;
    for (const Candidate& candidate : candidates) {
        if (!candidate.compatible || candidate.score <= 0) continue;
        if (selected == nullptr || candidate.score > selected->score) {
            selected = &candidate;
        }
    }
    logging::error(selected != nullptr,
                   "CutBuilder: no face-compatible C3D8 tetrahedralization exists");
    for (const auto& [face, selected_diagonal] : selected->diagonals) {
        face_diagonal_cache.emplace(face, selected_diagonal);
    }
    const auto& tets = selected->tets;

    std::vector<ID> ids;
    ids.reserve(tets.size());
    auto first = std::make_shared<C3D4>(element_id, tets[0]);
    first->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = first;
    ids.push_back(element_id);

    for (Index i = 1; i < tets.size(); ++i) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D4>(new_id, tets[i]));
        ids.push_back(new_id);
    }
    return ids;
}

std::vector<ID> tetrahedralize_c3d6(
    ModelData& model_data,
    C3D6& element,
    std::map<FaceKey, std::pair<ID, ID>>& face_diagonal_cache) {
    const ID element_id = element.elem_id;
    const auto& n = element.node_ids;
    constexpr std::array<std::array<Index, 3>, 6> permutations = {{
        {{0, 1, 2}}, {{0, 2, 1}}, {{1, 0, 2}},
        {{1, 2, 0}}, {{2, 0, 1}}, {{2, 1, 0}}
    }};

    struct Candidate {
        std::array<std::array<ID, 4>, 3> tets{};
        std::map<FaceKey, std::pair<ID, ID>> diagonals;
        Precision score = 0;
        bool compatible = true;
    };
    const auto face_key = [&](Index a, Index b) {
        FaceKey key = {n[a], n[b], n[a + 3], n[b + 3]};
        std::sort(key.begin(), key.end());
        return key;
    };
    const auto diagonal = [&](ID a, ID b) {
        return std::pair<ID, ID>{std::min(a, b), std::max(a, b)};
    };

    std::array<Candidate, 6> candidates{};
    for (Index candidate_index = 0;
         candidate_index < candidates.size(); ++candidate_index) {
        const auto [a, b, c] = permutations[candidate_index];
        Candidate& candidate = candidates[candidate_index];
        candidate.tets = {{
            oriented_c3d4_nodes(model_data, {n[a], n[b], n[c], n[a + 3]}),
            oriented_c3d4_nodes(model_data, {n[b], n[c], n[a + 3], n[b + 3]}),
            oriented_c3d4_nodes(model_data, {n[c], n[a + 3], n[b + 3], n[c + 3]})
        }};
        candidate.diagonals.emplace(
            face_key(a, b), diagonal(n[b], n[a + 3]));
        candidate.diagonals.emplace(
            face_key(b, c), diagonal(n[c], n[b + 3]));
        candidate.diagonals.emplace(
            face_key(c, a), diagonal(n[c], n[a + 3]));
        candidate.score = minimum_normalized_tet_volume(
            model_data, candidate.tets);
        for (const auto& [face, selected_diagonal] : candidate.diagonals) {
            const auto cached = face_diagonal_cache.find(face);
            if (cached != face_diagonal_cache.end() &&
                cached->second != selected_diagonal) {
                candidate.compatible = false;
            }
        }
    }

    const Candidate* selected = nullptr;
    for (const Candidate& candidate : candidates) {
        if (!candidate.compatible || candidate.score <= 0) continue;
        if (selected == nullptr || candidate.score > selected->score) {
            selected = &candidate;
        }
    }
    logging::error(selected != nullptr,
                   "CutBuilder: no face-compatible C3D6 tetrahedralization exists");
    Index reused_face_diagonals = 0;
    for (const auto& [face, selected_diagonal] : selected->diagonals) {
        const auto [unused, inserted] = face_diagonal_cache.emplace(
            face, selected_diagonal);
        (void) unused;
        reused_face_diagonals += !inserted;
    }
    logging::info(reused_face_diagonals > 0,
                  "CutBuilder: C3D6 element ", element_id,
                  " reused ", reused_face_diagonals,
                  " cached shared-face diagonal(s)");
    const auto& tets = selected->tets;
    logging::error(total_tet_volume(model_data, tets) >
                       std::numeric_limits<Precision>::epsilon(),
                   "CutBuilder: C3D6 tetrahedralization has zero volume");

    std::vector<ID> ids;
    ids.reserve(tets.size());
    auto first = std::make_shared<C3D4>(element_id, tets[0]);
    first->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = first;
    ids.push_back(element_id);
    for (Index i = 1; i < tets.size(); ++i) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D4>(new_id, tets[i]));
        ids.push_back(new_id);
    }
    return ids;
}

bool is_c3d5(const C3D8& element) {
    return element.node_ids[4] == element.node_ids[5] &&
        element.node_ids[4] == element.node_ids[6] &&
        element.node_ids[4] == element.node_ids[7];
}

std::set<ID> regular_c3d8_tetrahedralization_closure(
    const ModelData& model_data,
    const std::vector<ID>& cut_elements) {
    std::map<ID, std::array<FaceKey, 6>> element_faces;
    std::map<FaceKey, std::vector<ID>> face_elements;
    for (const ElementPtr& element : model_data.elements) {
        const auto* c3d8 = element == nullptr
            ? nullptr : dynamic_cast<const C3D8*>(element.get());
        if (c3d8 == nullptr || is_c3d5(*c3d8)) continue;
        std::array<FaceKey, 6> faces{};
        for (Index face_index = 0; face_index < c3d8_faces.size(); ++face_index) {
            faces[face_index] = make_face_key(
                c3d8->node_ids, c3d8_faces[face_index]);
            face_elements[faces[face_index]].push_back(c3d8->elem_id);
        }
        element_faces.emplace(c3d8->elem_id, faces);
    }

    std::set<ID> closure;
    std::deque<ID> pending;
    for (const ID element_id : cut_elements) {
        if (element_faces.find(element_id) != element_faces.end() &&
            closure.insert(element_id).second) {
            pending.push_back(element_id);
        }
    }
    while (!pending.empty()) {
        const ID element_id = pending.front();
        pending.pop_front();
        for (const FaceKey& face : element_faces.at(element_id)) {
            for (const ID neighbor_id : face_elements.at(face)) {
                if (closure.insert(neighbor_id).second) {
                    pending.push_back(neighbor_id);
                }
            }
        }
    }
    return closure;
}

std::array<FaceKey, 6> c3d8_face_keys(const C3D8& element) {
    std::array<FaceKey, 6> faces{};
    for (Index face_index = 0; face_index < c3d8_faces.size(); ++face_index) {
        faces[face_index] = make_face_key(
            element.node_ids, c3d8_faces[face_index]);
    }
    return faces;
}

void validate_c3d8_tetrahedralization_closure(
    const ModelData& model_data,
    const std::map<ID, std::array<FaceKey, 6>>& parent_faces,
    const std::map<ID, std::vector<ID>>& parent_tetrahedra) {
    std::map<FaceKey, std::vector<ID>> shared_faces;
    for (const auto& [parent_id, faces] : parent_faces) {
        for (const FaceKey& face : faces) {
            shared_faces[face].push_back(parent_id);
        }
    }

    constexpr std::array<std::array<Index, 3>, 4> tet_faces = {{
        {{0, 1, 2}}, {{0, 1, 3}}, {{0, 2, 3}}, {{1, 2, 3}}
    }};
    Index checked_faces = 0;
    for (const auto& [quad_face, parents] : shared_faces) {
        if (parents.size() < 2) continue;
        const FaceKey& face_nodes = quad_face;
        logging::error(parents.size() == 2,
                       "CutBuilder: non-manifold C3D8 closure face shared by ",
                       parents.size(), " elements");
        std::array<std::set<TetFaceKey>, 2> triangulations;
        for (Index side = 0; side < parents.size(); ++side) {
            const auto children = parent_tetrahedra.find(parents[side]);
            logging::error(children != parent_tetrahedra.end(),
                           "CutBuilder: missing tetrahedra for C3D8 closure parent ",
                           parents[side]);
            for (const ID tetra_id : children->second) {
                const auto* tetra = dynamic_cast<const C3D4*>(
                    model_data.elements[static_cast<std::size_t>(tetra_id)].get());
                logging::error(tetra != nullptr,
                               "CutBuilder: expected C3D4 child in closure validation");
                for (const auto& local_face : tet_faces) {
                    TetFaceKey triangle = {
                        tetra->node_ids[local_face[0]],
                        tetra->node_ids[local_face[1]],
                        tetra->node_ids[local_face[2]]};
                    const bool on_quad = std::all_of(
                        triangle.begin(), triangle.end(), [&](ID node_id) {
                            return std::find(face_nodes.begin(), face_nodes.end(),
                                             node_id) != face_nodes.end();
                        });
                    if (on_quad) {
                        std::sort(triangle.begin(), triangle.end());
                        const bool inserted = triangulations[side].insert(triangle).second;
                        logging::error(inserted,
                                       "CutBuilder: duplicate triangle on C3D8 closure face");
                    }
                }
            }
            logging::error(triangulations[side].size() == 2,
                           "CutBuilder: C3D8 closure quad face did not produce exactly two triangles");
        }
        logging::error(triangulations[0] == triangulations[1],
                       "CutBuilder: incompatible triangulation on shared C3D8 closure face");
        ++checked_faces;
    }
    logging::info(checked_faces > 0,
                  "CutBuilder: validated ", checked_faces,
                  " shared C3D8 closure face(s) with matching triangle topology");
}

std::vector<ID> tetrahedralize_c3d5(
    ModelData& model_data,
    C3D8& element,
    std::map<FaceKey, std::pair<ID, ID>>& face_diagonal_cache) {
    logging::error(is_c3d5(element),
                   "CutBuilder: invalid degenerate C3D5 representation");
    const ID element_id = element.elem_id;
    const auto& n = element.node_ids;
    FaceKey base_face = {n[0], n[1], n[2], n[3]};
    std::sort(base_face.begin(), base_face.end());
    const auto diagonal = [](ID a, ID b) {
        return std::pair<ID, ID>{std::min(a, b), std::max(a, b)};
    };

    struct Candidate {
        std::array<std::array<ID, 4>, 2> tets{};
        std::pair<ID, ID> diagonal{};
        Precision score = 0;
    };
    std::array<Candidate, 2> candidates{};
    candidates[0].tets = {{
        {{n[0], n[1], n[2], n[4]}},
        {{n[0], n[2], n[3], n[4]}}
    }};
    candidates[0].diagonal = diagonal(n[0], n[2]);
    candidates[1].tets = {{
        {{n[0], n[1], n[3], n[4]}},
        {{n[1], n[2], n[3], n[4]}}
    }};
    candidates[1].diagonal = diagonal(n[1], n[3]);
    const auto cached = face_diagonal_cache.find(base_face);
    Candidate* selected = nullptr;
    for (Candidate& candidate : candidates) {
        for (auto& tet : candidate.tets) {
            tet = oriented_c3d4_nodes(model_data, tet);
        }
        candidate.score = minimum_normalized_tet_volume(
            model_data, candidate.tets);
        if (candidate.score <= 0 ||
            (cached != face_diagonal_cache.end() &&
             cached->second != candidate.diagonal)) {
            continue;
        }
        if (selected == nullptr || candidate.score > selected->score) {
            selected = &candidate;
        }
    }
    logging::error(selected != nullptr,
                   "CutBuilder: no valid C3D5 tetrahedralization exists");
    logging::error(total_tet_volume(model_data, selected->tets) >
                       std::numeric_limits<Precision>::epsilon(),
                   "CutBuilder: C3D5 tetrahedralization has zero volume");
    logging::info(cached != face_diagonal_cache.end(),
                  "CutBuilder: C3D5 element ", element_id,
                  " reused its cached base-face diagonal");
    face_diagonal_cache.emplace(base_face, selected->diagonal);

    std::vector<ID> ids{element_id};
    auto first = std::make_shared<C3D4>(element_id, selected->tets[0]);
    first->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = first;
    const ID second_id = model_data.next_free_element_id();
    add_element_to_matching_sets(model_data, element_id, second_id);
    model_data.insert_element(std::make_shared<C3D4>(
        second_id, selected->tets[1]));
    ids.push_back(second_id);
    return ids;
}

FaceKey c3d5_base_face_key(const C3D8& element) {
    FaceKey key = {
        element.node_ids[0], element.node_ids[1],
        element.node_ids[2], element.node_ids[3]};
    std::sort(key.begin(), key.end());
    return key;
}

std::set<ID> c3d5_c3d8_tetrahedralization_closure(
    const ModelData& model_data,
    const std::vector<ID>& cut_elements) {
    std::map<ID, std::vector<FaceKey>> element_faces;
    std::map<FaceKey, std::vector<ID>> face_elements;
    for (const ElementPtr& element : model_data.elements) {
        const auto* c3d8 = element == nullptr
            ? nullptr : dynamic_cast<const C3D8*>(element.get());
        if (c3d8 == nullptr) continue;
        std::vector<FaceKey> faces;
        if (is_c3d5(*c3d8)) {
            faces.push_back(c3d5_base_face_key(*c3d8));
        } else {
            for (const auto& local_face : c3d8_faces) {
                faces.push_back(make_face_key(c3d8->node_ids, local_face));
            }
        }
        element_faces.emplace(c3d8->elem_id, faces);
        for (const FaceKey& face : faces) {
            face_elements[face].push_back(c3d8->elem_id);
        }
    }

    std::set<ID> closure;
    std::deque<ID> pending;
    for (const ID element_id : cut_elements) {
        const auto found = element_faces.find(element_id);
        const auto* c3d8 = found == element_faces.end() ? nullptr
            : dynamic_cast<const C3D8*>(model_data.elements[
                static_cast<std::size_t>(element_id)].get());
        if (c3d8 != nullptr && is_c3d5(*c3d8) &&
            closure.insert(element_id).second) {
            pending.push_back(element_id);
        }
    }
    while (!pending.empty()) {
        const ID element_id = pending.front();
        pending.pop_front();
        for (const FaceKey& face : element_faces.at(element_id)) {
            for (const ID neighbor_id : face_elements.at(face)) {
                if (closure.insert(neighbor_id).second) {
                    pending.push_back(neighbor_id);
                }
            }
        }
    }
    return closure;
}

std::array<FaceKey, 3> c3d6_quad_face_keys(
    const std::array<ID, 6>& nodes) {
    constexpr std::array<std::pair<Index, Index>, 3> column_pairs = {{
        {0, 1}, {1, 2}, {2, 0}
    }};
    std::array<FaceKey, 3> keys{};
    for (Index face_index = 0; face_index < keys.size(); ++face_index) {
        const auto [a, b] = column_pairs[face_index];
        keys[face_index] = {
            nodes[a], nodes[b], nodes[a + 3], nodes[b + 3]};
        std::sort(keys[face_index].begin(), keys[face_index].end());
    }
    return keys;
}

std::set<ID> c3d6_tetrahedralization_closure(
    const ModelData& model_data,
    const std::vector<ID>& cut_elements) {
    std::map<ID, std::array<FaceKey, 3>> element_faces;
    std::map<FaceKey, std::vector<ID>> face_elements;
    for (const ElementPtr& element : model_data.elements) {
        const auto* c3d6 = element == nullptr
            ? nullptr : dynamic_cast<const C3D6*>(element.get());
        if (c3d6 == nullptr) continue;
        const auto faces = c3d6_quad_face_keys(c3d6->node_ids);
        element_faces.emplace(c3d6->elem_id, faces);
        for (const FaceKey& face : faces) {
            face_elements[face].push_back(c3d6->elem_id);
        }
    }

    std::set<ID> closure;
    std::deque<ID> pending;
    for (const ID element_id : cut_elements) {
        if (element_faces.find(element_id) != element_faces.end() &&
            closure.insert(element_id).second) {
            pending.push_back(element_id);
        }
    }
    while (!pending.empty()) {
        const ID element_id = pending.front();
        pending.pop_front();
        for (const FaceKey& face : element_faces.at(element_id)) {
            for (const ID neighbor_id : face_elements.at(face)) {
                if (closure.insert(neighbor_id).second) {
                    pending.push_back(neighbor_id);
                }
            }
        }
    }
    return closure;
}

std::set<ID> c3d20_linearization_closure(
    const ModelData& model_data,
    const std::vector<ID>& cut_elements) {
    std::map<ID, std::array<FaceKey, 6>> element_faces;
    std::map<FaceKey, std::vector<ID>> face_elements;
    for (const ElementPtr& element : model_data.elements) {
        if (element == nullptr) continue;
        std::array<ID, 8> corners{};
        if (const auto* c3d20 = dynamic_cast<const C3D20*>(element.get())) {
            std::copy_n(c3d20->node_ids.begin(), corners.size(), corners.begin());
        } else if (const auto* c3d20r = dynamic_cast<const C3D20R*>(element.get())) {
            std::copy_n(c3d20r->node_ids.begin(), corners.size(), corners.begin());
        } else if (const auto* c3d8 = dynamic_cast<const C3D8*>(element.get())) {
            const bool degenerate_pyramid =
                c3d8->node_ids[4] == c3d8->node_ids[5] &&
                c3d8->node_ids[4] == c3d8->node_ids[6] &&
                c3d8->node_ids[4] == c3d8->node_ids[7];
            if (degenerate_pyramid) continue;
            corners = c3d8->node_ids;
        } else {
            continue;
        }

        std::array<FaceKey, 6> faces{};
        for (Index face_index = 0; face_index < c3d8_faces.size(); ++face_index) {
            faces[face_index] = make_face_key(corners, c3d8_faces[face_index]);
            face_elements[faces[face_index]].push_back(element->elem_id);
        }
        element_faces.emplace(element->elem_id, faces);
    }

    std::set<ID> closure;
    std::deque<ID> pending;
    for (const ID element_id : cut_elements) {
        if (element_faces.find(element_id) != element_faces.end() &&
            closure.insert(element_id).second) {
            pending.push_back(element_id);
        }
    }
    while (!pending.empty()) {
        const ID element_id = pending.front();
        pending.pop_front();
        for (const FaceKey& face : element_faces.at(element_id)) {
            for (const ID neighbor_id : face_elements.at(face)) {
                if (closure.insert(neighbor_id).second) {
                    pending.push_back(neighbor_id);
                }
            }
        }
    }
    return closure;
}

std::array<ID, 10> make_c3d10_from_tet(
    ModelData& model_data,
    const std::array<ID, 4>& vertices,
    const std::array<ID, 4>& original_vertices,
    const std::array<ID, 6>& original_mids,
    std::map<std::pair<ID, ID>, ID>& midpoint_cache) {
    std::array<ID, 10> nodes{};
    std::copy(vertices.begin(), vertices.end(), nodes.begin());
    constexpr std::array<std::pair<Index, Index>, 6> edges = {{
        {0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}
    }};

    constexpr std::array<std::pair<Index, Index>, 6> original_edges = {{
        {0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}
    }};

    for (Index edge_index = 0; edge_index < edges.size(); ++edge_index) {
        const ID node_a = vertices[edges[edge_index].first];
        const ID node_b = vertices[edges[edge_index].second];
        const auto key = std::minmax(node_a, node_b);
        ID midpoint = -1;

        for (Index original_edge = 0; original_edge < original_edges.size(); ++original_edge) {
            const auto original_key = std::minmax(
                original_vertices[original_edges[original_edge].first],
                original_vertices[original_edges[original_edge].second]);
            if (key == original_key) {
                midpoint = original_mids[original_edge];
                break;
            }
        }

        if (midpoint < 0) {
            const auto cached = midpoint_cache.find(key);
            if (cached != midpoint_cache.end()) {
                midpoint = cached->second;
            } else {
                const Vec3 position_a = model_data.positions->row_vec3(
                    static_cast<Index>(node_a));
                const Vec3 position_b = model_data.positions->row_vec3(
                    static_cast<Index>(node_b));
                midpoint = model_data.append_node((position_a + position_b) * Precision(0.5));
                midpoint_cache.emplace(key, midpoint);
            }
        }
        nodes[4 + edge_index] = midpoint;
    }
    return nodes;
}

std::vector<ID> tetrahedralize_c3d15(
    ModelData& model_data,
    C3D15& element,
    std::map<FaceKey, std::pair<ID, ID>>& face_diagonal_cache) {
    const ID element_id = element.elem_id;
    const auto& n = element.node_ids;
    constexpr std::array<std::array<Index, 3>, 6> permutations = {{
        {{0, 1, 2}}, {{0, 2, 1}}, {{1, 0, 2}},
        {{1, 2, 0}}, {{2, 0, 1}}, {{2, 1, 0}}
    }};
    struct Candidate {
        std::array<std::array<ID, 4>, 3> vertices{};
        std::map<FaceKey, std::pair<ID, ID>> diagonals;
        Precision score = 0;
        bool compatible = true;
    };
    const auto face_key = [&](Index a, Index b) {
        FaceKey key = {n[a], n[b], n[a + 3], n[b + 3]};
        std::sort(key.begin(), key.end());
        return key;
    };
    const auto diagonal = [](ID a, ID b) {
        return std::pair<ID, ID>{std::min(a, b), std::max(a, b)};
    };
    std::array<Candidate, 6> candidates{};
    for (Index i = 0; i < candidates.size(); ++i) {
        const auto [a, b, c] = permutations[i];
        Candidate& candidate = candidates[i];
        candidate.vertices = {{
            oriented_c3d4_nodes(model_data, {n[a], n[b], n[c], n[a + 3]}),
            oriented_c3d4_nodes(model_data, {n[b], n[c], n[a + 3], n[b + 3]}),
            oriented_c3d4_nodes(model_data, {n[c], n[a + 3], n[b + 3], n[c + 3]})
        }};
        candidate.diagonals.emplace(face_key(a, b), diagonal(n[b], n[a + 3]));
        candidate.diagonals.emplace(face_key(b, c), diagonal(n[c], n[b + 3]));
        candidate.diagonals.emplace(face_key(c, a), diagonal(n[c], n[a + 3]));
        candidate.score = minimum_normalized_tet_volume(
            model_data, candidate.vertices);
        for (const auto& [face, selected_diagonal] : candidate.diagonals) {
            const auto cached = face_diagonal_cache.find(face);
            candidate.compatible = candidate.compatible &&
                (cached == face_diagonal_cache.end() ||
                 cached->second == selected_diagonal);
        }
    }
    const Candidate* selected = nullptr;
    for (const Candidate& candidate : candidates) {
        if (!candidate.compatible || candidate.score <= 0) continue;
        if (selected == nullptr || candidate.score > selected->score) selected = &candidate;
    }
    logging::error(selected != nullptr,
                   "CutBuilder: no face-compatible C3D15 tetrahedralization exists");
    for (const auto& [face, selected_diagonal] : selected->diagonals) {
        face_diagonal_cache.emplace(face, selected_diagonal);
    }

    std::vector<ID> ids{element_id};
    auto first = std::make_shared<C3D4>(element_id, selected->vertices[0]);
    first->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = first;
    for (Index i = 1; i < selected->vertices.size(); ++i) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D4>(
            new_id, selected->vertices[i]));
        ids.push_back(new_id);
    }
    return ids;
}

std::vector<ID> tetrahedralize_c3d13(
    ModelData& model_data,
    C3D13& element,
    std::map<FaceKey, std::pair<ID, ID>>& face_diagonal_cache) {
    const ID element_id = element.elem_id;
    const auto& n = element.node_ids;
    FaceKey base_face = {n[0], n[1], n[2], n[3]};
    std::sort(base_face.begin(), base_face.end());

    struct Candidate {
        std::array<std::array<ID, 4>, 2> vertices{};
        std::pair<ID, ID> diagonal{};
        Precision score = 0;
        bool compatible = true;
    };
    std::array<Candidate, 2> candidates{};
    candidates[0].vertices = {{
        oriented_c3d4_nodes(model_data, {n[0], n[1], n[2], n[4]}),
        oriented_c3d4_nodes(model_data, {n[0], n[2], n[3], n[4]})}};
    candidates[0].diagonal = std::minmax(n[0], n[2]);
    candidates[1].vertices = {{
        oriented_c3d4_nodes(model_data, {n[0], n[1], n[3], n[4]}),
        oriented_c3d4_nodes(model_data, {n[1], n[2], n[3], n[4]})}};
    candidates[1].diagonal = std::minmax(n[1], n[3]);

    const auto cached = face_diagonal_cache.find(base_face);
    const Candidate* selected = nullptr;
    for (Candidate& candidate : candidates) {
        candidate.score = minimum_normalized_tet_volume(
            model_data, candidate.vertices);
        candidate.compatible = cached == face_diagonal_cache.end() ||
            cached->second == candidate.diagonal;
        if (!candidate.compatible || candidate.score <= 0) continue;
        if (selected == nullptr || candidate.score > selected->score) {
            selected = &candidate;
        }
    }
    logging::error(selected != nullptr,
                   "CutBuilder: no face-compatible C3D13 tetrahedralization exists");
    face_diagonal_cache.emplace(base_face, selected->diagonal);

    std::vector<ID> ids{element_id};
    auto first = std::make_shared<C3D4>(element_id, selected->vertices[0]);
    first->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = first;
    const ID new_id = model_data.next_free_element_id();
    add_element_to_matching_sets(model_data, element_id, new_id);
    model_data.insert_element(std::make_shared<C3D4>(
        new_id, selected->vertices[1]));
    ids.push_back(new_id);
    return ids;
}

void split_c3d4_three_edge_cut(
    ModelData& model_data,
    pretension::PretensionSection& section,
    C3D4& element,
    const Vec3& plane_point,
    const Vec3& axis,
    bool reverse_sides) {
    const ID element_id = element.elem_id;
    const auto old_nodes = element.node_ids;
    std::array<Precision, 4> distances{};
    std::array<Index, 4> positive{};
    std::array<Index, 4> negative{};
    Index positive_count = 0;
    Index negative_count = 0;

    for (Index i = 0; i < 4; ++i) {
        const Vec3 position = model_data.positions->row_vec3(
            static_cast<Index>(old_nodes[i]));
        distances[i] = axis.dot(position - plane_point);
        if (distances[i] > 0) {
            positive[positive_count++] = i;
        } else if (distances[i] < 0) {
            negative[negative_count++] = i;
        }
    }

    logging::error(positive_count == 1 && negative_count == 3,
                   "CutBuilder: C3D4 currently supports one positive and three negative nodes");

    const bool positive_single = positive_count == 1;
    const Index single = positive_single ? positive[0] : negative[0];
    const auto& three = positive_single ? negative : positive;

    std::array<ID, 3> side_a_interface{};
    std::array<ID, 3> side_b_interface{};
    for (Index i = 0; i < 3; ++i) {
        const ID single_node = old_nodes[single];
        const ID three_node = old_nodes[three[i]];
        const Vec3 single_position = model_data.positions->row_vec3(
            static_cast<Index>(single_node));
        const Vec3 three_position = model_data.positions->row_vec3(
            static_cast<Index>(three_node));
        const Precision single_distance = distances[single];
        const Precision three_distance = distances[three[i]];
        const Precision parameter = single_distance /
            (single_distance - three_distance);
        const Vec3 intersection = single_position +
            parameter * (three_position - single_position);

        const auto pair = get_or_create_interface_pair(
            model_data, section, intersection, reverse_sides);
        side_a_interface[i] = reverse_sides ? pair.side_b : pair.side_a;
        side_b_interface[i] = reverse_sides ? pair.side_a : pair.side_b;
    }

    std::array<std::array<ID, 4>, 3> side_a_tets{};
    std::array<ID, 4> side_b_tet{};
    if (positive_single) {
        side_a_tets = {{
            {old_nodes[three[0]], old_nodes[three[1]], old_nodes[three[2]], side_a_interface[0]},
            {old_nodes[three[1]], old_nodes[three[2]], side_a_interface[0], side_a_interface[1]},
            {old_nodes[three[2]], side_a_interface[0], side_a_interface[1], side_a_interface[2]}
        }};
        side_b_tet = {old_nodes[single], side_b_interface[0],
                      side_b_interface[1], side_b_interface[2]};
    }

    for (auto& tet : side_a_tets) {
        tet = oriented_c3d4_nodes(model_data, tet);
    }
    side_b_tet = oriented_c3d4_nodes(model_data, side_b_tet);

    Precision minimum_child_quality = c3d4_quality(
        model_data, side_b_tet).normalized_volume;
    for (const auto& tet : side_a_tets) {
        minimum_child_quality = std::min(
            minimum_child_quality,
            c3d4_quality(model_data, tet).normalized_volume);
    }
    logging::info(minimum_child_quality < Precision(1e-4),
                  "CutBuilder: low-quality C3D4 1/3 split candidate from element ",
                  element_id, ": minimum normalized child volume = ",
                  minimum_child_quality);

    auto side_a_element = std::make_shared<C3D4>(element_id, side_a_tets[0]);
    side_a_element->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = side_a_element;
    section.cut_quality_origins[element_id] = element_id;
    section.cut_quality_types[element_id] = "1/3";

    for (Index i = 1; i < side_a_tets.size(); ++i) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D4>(new_id, side_a_tets[i]));
        section.cut_quality_origins[new_id] = element_id;
        section.cut_quality_types[new_id] = "1/3";
    }
    const ID side_b_id = model_data.next_free_element_id();
    add_element_to_matching_sets(model_data, element_id, side_b_id);
    model_data.insert_element(std::make_shared<C3D4>(side_b_id, side_b_tet));
    section.cut_quality_origins[side_b_id] = element_id;
    section.cut_quality_types[side_b_id] = "1/3";

    section.side_a_nodes.insert(section.side_a_nodes.end(),
                                side_a_interface.begin(), side_a_interface.end());
    section.side_b_nodes.insert(section.side_b_nodes.end(),
                                side_b_interface.begin(), side_b_interface.end());
}

void split_c3d4_four_edge_cut(
    ModelData& model_data,
    pretension::PretensionSection& section,
    C3D4& element,
    const Vec3& plane_point,
    const Vec3& axis,
    Index& alternative_diagonal_count) {
    const ID element_id = element.elem_id;
    const auto old_nodes = element.node_ids;
    std::array<Precision, 4> distances{};
    std::array<Index, 2> positive{};
    std::array<Index, 2> negative{};
    Index positive_count = 0;
    Index negative_count = 0;

    for (Index i = 0; i < 4; ++i) {
        const Vec3 position = model_data.positions->row_vec3(
            static_cast<Index>(old_nodes[i]));
        distances[i] = axis.dot(position - plane_point);
        if (distances[i] > 0) {
            positive[positive_count++] = i;
        } else if (distances[i] < 0) {
            negative[negative_count++] = i;
        }
    }

    logging::error(positive_count == 2 && negative_count == 2,
                   "CutBuilder: C3D4 four-edge cut requires two nodes on each side");

    std::array<std::array<ID, 2>, 2> side_a_interface{};
    std::array<std::array<ID, 2>, 2> side_b_interface{};

    for (Index n = 0; n < 2; ++n) {
        for (Index p = 0; p < 2; ++p) {
            const ID negative_node = old_nodes[negative[n]];
            const ID positive_node = old_nodes[positive[p]];
            const Vec3 negative_position = model_data.positions->row_vec3(
                static_cast<Index>(negative_node));
            const Vec3 positive_position = model_data.positions->row_vec3(
                static_cast<Index>(positive_node));
            const Precision parameter = distances[negative[n]] /
                (distances[negative[n]] - distances[positive[p]]);
            const Vec3 intersection = negative_position +
                parameter * (positive_position - negative_position);

            const auto pair = get_or_create_interface_pair(model_data, section, intersection);
            side_a_interface[n][p] = pair.side_a;
            side_b_interface[n][p] = pair.side_b;
        }
    }

    const auto make_side_a = [&](Index first, Index second) {
        return std::array<std::array<ID, 4>, 3>{{
            {old_nodes[negative[first]], old_nodes[negative[second]],
             side_a_interface[first][0], side_a_interface[first][1]},
            {old_nodes[negative[second]], side_a_interface[first][0],
             side_a_interface[first][1], side_a_interface[second][1]},
            {old_nodes[negative[second]], side_a_interface[first][0],
             side_a_interface[second][1], side_a_interface[second][0]}
        }};
    };
    const auto make_side_b = [&](Index first, Index second) {
        return std::array<std::array<ID, 4>, 3>{{
            {old_nodes[positive[first]], old_nodes[positive[second]],
             side_b_interface[0][first], side_b_interface[1][first]},
            {old_nodes[positive[second]], side_b_interface[0][first],
             side_b_interface[1][first], side_b_interface[1][second]},
            {old_nodes[positive[second]], side_b_interface[0][first],
             side_b_interface[1][second], side_b_interface[0][second]}
        }};
    };

    auto side_a_tets = make_side_a(0, 1);
    auto side_b_tets = make_side_b(0, 1);
    const auto alternative_side_a = make_side_a(1, 0);
    const auto alternative_side_b = make_side_b(1, 0);

    const Precision original_volume = c3d4_quality(model_data, old_nodes).volume;
    const auto candidate_is_valid = [&](const auto& side_a, const auto& side_b) {
        const Precision child_volume = total_tet_volume(model_data, side_a) +
            total_tet_volume(model_data, side_b);
        const Precision relative_error = original_volume > 0
            ? std::abs(child_volume - original_volume) / original_volume
            : std::numeric_limits<Precision>::max();
        return minimum_normalized_tet_volume(model_data, side_a) > 0 &&
            minimum_normalized_tet_volume(model_data, side_b) > 0 &&
            relative_error <= Precision(1e-10);
    };
    logging::error(candidate_is_valid(side_a_tets, side_b_tets),
                   "CutBuilder: invalid primary C3D4 2/2 split candidate");
    const bool alternative_valid = candidate_is_valid(
        alternative_side_a, alternative_side_b);
    const Precision primary_score = std::min(
        minimum_normalized_tet_volume(model_data, side_a_tets),
        minimum_normalized_tet_volume(model_data, side_b_tets));
    const Precision alternative_score = alternative_valid
        ? std::min(
            minimum_normalized_tet_volume(model_data, alternative_side_a),
            minimum_normalized_tet_volume(model_data, alternative_side_b))
        : Precision(0);
    if (alternative_score > primary_score) {
        side_a_tets = alternative_side_a;
        side_b_tets = alternative_side_b;
        ++alternative_diagonal_count;
    }

    for (auto& tet : side_a_tets) {
        tet = oriented_c3d4_nodes(model_data, tet);
    }
    for (auto& tet : side_b_tets) {
        tet = oriented_c3d4_nodes(model_data, tet);
    }

    auto side_a_element = std::make_shared<C3D4>(element_id, side_a_tets[0]);
    side_a_element->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = side_a_element;
    section.cut_quality_origins[element_id] = element_id;
    section.cut_quality_types[element_id] = "2/2";

    for (Index i = 1; i < side_a_tets.size(); ++i) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D4>(new_id, side_a_tets[i]));
        section.cut_quality_origins[new_id] = element_id;
        section.cut_quality_types[new_id] = "2/2";
    }
    for (const auto& tet : side_b_tets) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D4>(new_id, tet));
        section.cut_quality_origins[new_id] = element_id;
        section.cut_quality_types[new_id] = "2/2";
    }

    for (const auto& pair : side_a_interface) {
        section.side_a_nodes.push_back(pair[0]);
        section.side_a_nodes.push_back(pair[1]);
    }
    for (const auto& pair : side_b_interface) {
        section.side_b_nodes.push_back(pair[0]);
        section.side_b_nodes.push_back(pair[1]);
    }
}

void split_c3d4_element(
    ModelData& model_data,
    pretension::PretensionSection& section,
    ID element_id,
    const Vec3& plane_point,
    const Vec3& axis,
    Precision tolerance,
    Index& alternative_diagonal_count) {
    auto* c3d4 = dynamic_cast<C3D4*>(
        model_data.elements[static_cast<std::size_t>(element_id)].get());
    logging::error(c3d4 != nullptr,
                   "CutBuilder: expected a C3D4 after tetrahedralization");

    Index positive_count = 0;
    for (Index local_node = 0; local_node < 4; ++local_node) {
        const Vec3 position = model_data.positions->row_vec3(
            static_cast<Index>(c3d4->node_ids[local_node]));
        if (axis.dot(position - plane_point) > tolerance) {
            ++positive_count;
        }
    }
    if (positive_count == 0 || positive_count == 4) {
        return;
    }
    if (positive_count == 2) {
        split_c3d4_four_edge_cut(model_data, section, *c3d4,
                                 plane_point, axis, alternative_diagonal_count);
    } else if (positive_count == 3) {
        split_c3d4_three_edge_cut(model_data, section, *c3d4,
                                  plane_point, -axis, true);
    } else {
        split_c3d4_three_edge_cut(model_data, section, *c3d4,
                                  plane_point, axis, false);
    }
}

void split_c3d10_three_edge_cut(
    ModelData& model_data,
    pretension::PretensionSection& section,
    C3D10& element,
    const Vec3& plane_point,
    const Vec3& axis) {
    const ID element_id = element.elem_id;
    const auto old_nodes = element.node_ids;
    const std::array<ID, 4> original_vertices = {
        old_nodes[0], old_nodes[1], old_nodes[2], old_nodes[3]};
    const std::array<ID, 6> original_mids = {
        old_nodes[4], old_nodes[5], old_nodes[6],
        old_nodes[7], old_nodes[8], old_nodes[9]};
    std::array<Precision, 4> distances{};
    std::array<Index, 4> positive{};
    std::array<Index, 4> negative{};
    Index positive_count = 0;
    Index negative_count = 0;

    for (Index i = 0; i < 4; ++i) {
        const Vec3 position = model_data.positions->row_vec3(
            static_cast<Index>(original_vertices[i]));
        distances[i] = axis.dot(position - plane_point);
        if (distances[i] > 0) {
            positive[positive_count++] = i;
        } else if (distances[i] < 0) {
            negative[negative_count++] = i;
        }
    }

    logging::error(positive_count == 1 && negative_count == 3,
                   "CutBuilder: C3D10 currently supports one positive and three negative vertices");

    for (Index edge = 0; edge < original_mids.size(); ++edge) {
        const std::array<std::pair<Index, Index>, 6> edge_map = {{
            {0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}
        }};
        const auto [a, b] = edge_map[edge];
        const Vec3 pa = model_data.positions->row_vec3(static_cast<Index>(original_vertices[a]));
        const Vec3 pb = model_data.positions->row_vec3(static_cast<Index>(original_vertices[b]));
        const Vec3 pm = model_data.positions->row_vec3(static_cast<Index>(original_mids[edge]));
        logging::error(pm.isApprox((pa + pb) * Precision(0.5), Precision(1e-10)),
                       "CutBuilder: curved C3D10 edges are not supported yet");
    }

    const Index single = positive[0];
    const auto& three = negative;
    std::array<ID, 3> side_a_interface{};
    std::array<ID, 3> side_b_interface{};
    for (Index i = 0; i < 3; ++i) {
        const ID single_node = original_vertices[single];
        const ID three_node = original_vertices[three[i]];
        const Vec3 single_position = model_data.positions->row_vec3(static_cast<Index>(single_node));
        const Vec3 three_position = model_data.positions->row_vec3(static_cast<Index>(three_node));
        const Precision parameter = distances[single] /
            (distances[single] - distances[three[i]]);
        const Vec3 intersection = single_position +
            parameter * (three_position - single_position);
        const auto pair = get_or_create_interface_pair(model_data, section, intersection);
        side_a_interface[i] = pair.side_a;
        side_b_interface[i] = pair.side_b;
    }

    const std::array<std::array<ID, 4>, 3> side_a_vertices = {{
        {original_vertices[three[0]], original_vertices[three[1]], original_vertices[three[2]], side_a_interface[0]},
        {original_vertices[three[1]], original_vertices[three[2]], side_a_interface[0], side_a_interface[1]},
        {original_vertices[three[2]], side_a_interface[0], side_a_interface[1], side_a_interface[2]}
    }};
    const std::array<ID, 4> side_b_vertices = {
        original_vertices[single], side_b_interface[0],
        side_b_interface[1], side_b_interface[2]};
    std::map<std::pair<ID, ID>, ID> side_a_midpoints;
    std::map<std::pair<ID, ID>, ID> side_b_midpoints;
    std::array<std::array<ID, 10>, 3> side_a_nodes{};
    for (Index i = 0; i < side_a_vertices.size(); ++i) {
        const auto vertices = oriented_c3d4_nodes(model_data, side_a_vertices[i]);
        side_a_nodes[i] = make_c3d10_from_tet(
            model_data, vertices, original_vertices, original_mids, side_a_midpoints);
    }
    const auto oriented_side_b = oriented_c3d4_nodes(model_data, side_b_vertices);
    const auto side_b_nodes = make_c3d10_from_tet(
        model_data, oriented_side_b, original_vertices, original_mids, side_b_midpoints);

    auto first = std::make_shared<C3D10>(element_id, side_a_nodes[0]);
    first->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = first;
    for (Index i = 1; i < side_a_nodes.size(); ++i) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D10>(new_id, side_a_nodes[i]));
    }
    const ID side_b_id = model_data.next_free_element_id();
    add_element_to_matching_sets(model_data, element_id, side_b_id);
    model_data.insert_element(std::make_shared<C3D10>(side_b_id, side_b_nodes));

    section.side_a_nodes.insert(section.side_a_nodes.end(),
                                side_a_interface.begin(), side_a_interface.end());
    section.side_b_nodes.insert(section.side_b_nodes.end(),
                                side_b_interface.begin(), side_b_interface.end());
}

void split_c3d10_four_edge_cut(
    ModelData& model_data,
    pretension::PretensionSection& section,
    C3D10& element,
    const Vec3& plane_point,
    const Vec3& axis) {
    const ID element_id = element.elem_id;
    const auto old_nodes = element.node_ids;
    const std::array<ID, 4> original_vertices = {
        old_nodes[0], old_nodes[1], old_nodes[2], old_nodes[3]};
    const std::array<ID, 6> original_mids = {
        old_nodes[4], old_nodes[5], old_nodes[6],
        old_nodes[7], old_nodes[8], old_nodes[9]};
    std::array<Precision, 4> distances{};
    std::array<Index, 2> positive{};
    std::array<Index, 2> negative{};
    Index positive_count = 0;
    Index negative_count = 0;

    for (Index i = 0; i < 4; ++i) {
        const Vec3 position = model_data.positions->row_vec3(
            static_cast<Index>(original_vertices[i]));
        distances[i] = axis.dot(position - plane_point);
        if (distances[i] > 0) {
            positive[positive_count++] = i;
        } else if (distances[i] < 0) {
            negative[negative_count++] = i;
        }
    }
    logging::error(positive_count == 2 && negative_count == 2,
                   "CutBuilder: C3D10 four-edge cut requires two vertices on each side");

    const std::array<std::pair<Index, Index>, 6> edge_map = {{
        {0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}
    }};
    for (Index edge = 0; edge < edge_map.size(); ++edge) {
        const auto [a, b] = edge_map[edge];
        const Vec3 pa = model_data.positions->row_vec3(static_cast<Index>(original_vertices[a]));
        const Vec3 pb = model_data.positions->row_vec3(static_cast<Index>(original_vertices[b]));
        const Vec3 pm = model_data.positions->row_vec3(static_cast<Index>(original_mids[edge]));
        logging::error(pm.isApprox((pa + pb) * Precision(0.5), Precision(1e-10)),
                       "CutBuilder: curved C3D10 edges are not supported yet");
    }

    std::array<std::array<ID, 2>, 2> side_a_interface{};
    std::array<std::array<ID, 2>, 2> side_b_interface{};
    for (Index n = 0; n < 2; ++n) {
        for (Index p = 0; p < 2; ++p) {
            const ID negative_node = original_vertices[negative[n]];
            const ID positive_node = original_vertices[positive[p]];
            const Vec3 negative_position = model_data.positions->row_vec3(static_cast<Index>(negative_node));
            const Vec3 positive_position = model_data.positions->row_vec3(static_cast<Index>(positive_node));
            const Precision parameter = distances[negative[n]] /
                (distances[negative[n]] - distances[positive[p]]);
            const Vec3 intersection = negative_position +
                parameter * (positive_position - negative_position);
            const auto pair = get_or_create_interface_pair(model_data, section, intersection);
            side_a_interface[n][p] = pair.side_a;
            side_b_interface[n][p] = pair.side_b;
        }
    }

    const std::array<std::array<ID, 4>, 3> side_a_vertices = {{
        {original_vertices[negative[0]], original_vertices[negative[1]], side_a_interface[0][0], side_a_interface[0][1]},
        {original_vertices[negative[1]], side_a_interface[0][0], side_a_interface[0][1], side_a_interface[1][1]},
        {original_vertices[negative[1]], side_a_interface[0][0], side_a_interface[1][1], side_a_interface[1][0]}
    }};
    const std::array<std::array<ID, 4>, 3> side_b_vertices = {{
        {original_vertices[positive[0]], original_vertices[positive[1]], side_b_interface[0][0], side_b_interface[1][0]},
        {original_vertices[positive[1]], side_b_interface[0][0], side_b_interface[1][0], side_b_interface[1][1]},
        {original_vertices[positive[1]], side_b_interface[0][0], side_b_interface[1][1], side_b_interface[0][1]}
    }};

    std::map<std::pair<ID, ID>, ID> side_a_midpoints;
    std::map<std::pair<ID, ID>, ID> side_b_midpoints;
    std::array<std::array<ID, 10>, 3> side_a_nodes{};
    std::array<std::array<ID, 10>, 3> side_b_nodes{};
    for (Index i = 0; i < 3; ++i) {
        side_a_nodes[i] = make_c3d10_from_tet(
            model_data, oriented_c3d4_nodes(model_data, side_a_vertices[i]),
            original_vertices, original_mids, side_a_midpoints);
        side_b_nodes[i] = make_c3d10_from_tet(
            model_data, oriented_c3d4_nodes(model_data, side_b_vertices[i]),
            original_vertices, original_mids, side_b_midpoints);
    }

    auto first = std::make_shared<C3D10>(element_id, side_a_nodes[0]);
    first->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = first;
    for (Index i = 1; i < 3; ++i) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D10>(new_id, side_a_nodes[i]));
    }
    for (const auto& nodes : side_b_nodes) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D10>(new_id, nodes));
    }

    for (const auto& pair : side_a_interface) {
        section.side_a_nodes.push_back(pair[0]);
        section.side_a_nodes.push_back(pair[1]);
    }
    for (const auto& pair : side_b_interface) {
        section.side_b_nodes.push_back(pair[0]);
        section.side_b_nodes.push_back(pair[1]);
    }
}

void split_c3d10_element(
    ModelData& model_data,
    pretension::PretensionSection& section,
    ID element_id,
    const Vec3& plane_point,
    const Vec3& axis,
    Precision tolerance) {
    auto* c3d10 = dynamic_cast<C3D10*>(model_data.elements[
        static_cast<std::size_t>(element_id)].get());
    logging::error(c3d10 != nullptr,
                   "CutBuilder: expected C3D10 element");
    Index positive_count = 0;
    Index negative_count = 0;
    for (Index local_node = 0; local_node < 4; ++local_node) {
        const Vec3 position = model_data.positions->row_vec3(
            static_cast<Index>(c3d10->node_ids[local_node]));
        const Precision distance = axis.dot(position - plane_point);
        positive_count += distance > tolerance;
        negative_count += distance < -tolerance;
    }
    if (positive_count == 0 || negative_count == 0) return;
    if (positive_count == 2) {
        split_c3d10_four_edge_cut(
            model_data, section, *c3d10, plane_point, axis);
    } else if (positive_count == 3) {
        const std::size_t pair_start = section.interface_pairs.size();
        const std::size_t node_start = section.side_a_nodes.size();
        split_c3d10_three_edge_cut(
            model_data, section, *c3d10, plane_point, -axis);
        for (std::size_t i = pair_start; i < section.interface_pairs.size(); ++i) {
            std::swap(section.interface_pairs[i].side_a,
                      section.interface_pairs[i].side_b);
        }
        for (std::size_t i = node_start; i < section.side_a_nodes.size(); ++i) {
            std::swap(section.side_a_nodes[i], section.side_b_nodes[i]);
        }
    } else {
        split_c3d10_three_edge_cut(
            model_data, section, *c3d10, plane_point, axis);
    }
}

} // namespace

void CutBuilder::split(
    ModelData& model_data,
    pretension::PretensionSection& section) {
    logging::error(model_data.positions != nullptr,
                   "CutBuilder: POSITION field is not initialized");

    const Precision axis_norm = section.axis_direction.norm();
    logging::error(axis_norm > std::numeric_limits<Precision>::epsilon(),
                   "CutBuilder: pretension axis must not be zero");

    const Vec3 axis = section.axis_direction / axis_norm;
    const Vec3 plane_point = section.axis_origin + section.cut_coordinate * axis;
    const Precision tolerance = Precision(100) * std::numeric_limits<Precision>::epsilon();

    section.cut_element_ids.clear();
    section.side_a_nodes.clear();
    section.side_b_nodes.clear();
    section.interface_pairs.clear();
    section.interface_node_cache.clear();
    section.cut_quality_origins.clear();
    section.cut_quality_types.clear();

    // Build the future snapping/splitting plan from the untouched C3D4 mesh.
    // This is deliberately read-only until all aligned cut cases use the same
    // plan in one mutation pass.
    const C3D4CutAnalysis c3d4_analysis = analyze_c3d4_cuts(
        model_data, plane_point, axis, section.snap_ratio);
    log_c3d4_cut_analysis(section, c3d4_analysis);
    validate_vertex_aligned_plans(
        model_data, section, c3d4_analysis, plane_point, axis);
    validate_edge_aligned_plans(
        model_data, section, c3d4_analysis, plane_point, axis);
    const bool applied_aligned_snapping = apply_aligned_c3d4_plans(
        model_data, section, c3d4_analysis, plane_point, axis);
    if (applied_aligned_snapping) {
        logging::info(true,
                      "PRETENSIONSECTION '", section.name,
                      "': applied aligned C3D4 snapping");
    }

    // Handle a mesh-aligned cut before looking for elements crossed by the
    // plane. In this case the plane coincides with an existing element face
    // and no element has nodes on both sides of it.
    split_face_aligned_c3d4(
        model_data, section, plane_point, axis, tolerance);
    split_face_aligned_c3d8(
        model_data, section, plane_point, axis, tolerance);

    for (const ElementPtr& element : model_data.elements) {
        if (element == nullptr) {
            continue;
        }

        Precision minimum_distance = std::numeric_limits<Precision>::max();
        Precision maximum_distance = std::numeric_limits<Precision>::lowest();

        for (Index local_node = 0;
             local_node < static_cast<Index>(element->n_nodes());
             ++local_node) {
            const ID node_id = element->nodes()[local_node];
            const Vec3 position = model_data.positions->row_vec3(static_cast<Index>(node_id));
            const Precision distance = axis.dot(position - plane_point);

            minimum_distance = std::min(minimum_distance, distance);
            maximum_distance = std::max(maximum_distance, distance);
        }

        if (minimum_distance < -tolerance && maximum_distance > tolerance) {
            section.cut_element_ids.push_back(element->elem_id);
        }
    }

    const auto cut_elements = section.cut_element_ids;
    Index alternative_diagonal_count = 0;
    std::map<FaceKey, std::pair<ID, ID>> c3d6_face_diagonal_cache;
    const bool axis_aligned_cut =
        axis.isApprox(Vec3::UnitX(), Precision(1e-12)) ||
        axis.isApprox(-Vec3::UnitX(), Precision(1e-12)) ||
        axis.isApprox(Vec3::UnitY(), Precision(1e-12)) ||
        axis.isApprox(-Vec3::UnitY(), Precision(1e-12)) ||
        axis.isApprox(Vec3::UnitZ(), Precision(1e-12)) ||
        axis.isApprox(-Vec3::UnitZ(), Precision(1e-12));

    const std::set<ID> c3d20_closure = c3d20_linearization_closure(
        model_data, cut_elements);
    Index cut_c3d20_count = 0;
    Index linearized_c3d20_count = 0;
    for (const ID element_id : c3d20_closure) {
        ElementPtr& element = model_data.elements[static_cast<std::size_t>(element_id)];
        const auto* c3d20 = dynamic_cast<const C3D20*>(element.get());
        const auto* c3d20r = dynamic_cast<const C3D20R*>(element.get());
        if (c3d20 == nullptr && c3d20r == nullptr) continue;
        std::array<ID, 8> corner_nodes{};
        if (c3d20 != nullptr) {
            std::copy_n(c3d20->node_ids.begin(), corner_nodes.size(),
                        corner_nodes.begin());
        } else {
            std::copy_n(c3d20r->node_ids.begin(), corner_nodes.size(),
                        corner_nodes.begin());
        }
        const bool is_cut = std::find(cut_elements.begin(), cut_elements.end(),
                                      element_id) != cut_elements.end();
        cut_c3d20_count += is_cut;
        ++linearized_c3d20_count;
        logging::warning(false,
                         "CutBuilder: ", c3d20r != nullptr ? "C3D20R" : "C3D20",
                         is_cut ? " cut element " : " closure element ", element_id,
                         " is converted to C3D8; quadratic order",
                         c3d20r != nullptr ? " and reduced integration are" : " is",
                         " not preserved");
        auto linear = std::make_shared<C3D8>(element_id, corner_nodes);
        linear->_model_data = &model_data;
        element = std::move(linear);
    }
    logging::info(linearized_c3d20_count > 0,
                  "PRETENSIONSECTION '", section.name,
                  "': linearized ", linearized_c3d20_count,
                  " C3D20/C3D20R element(s) in a regular hex closure of ",
                  c3d20_closure.size(), " element(s) for ",
                  cut_c3d20_count, " cut seed(s)");

    std::set<ID> cut_regular_c3d8_elements;
    if (!axis_aligned_cut) {
        const std::set<ID> c3d8_closure =
            regular_c3d8_tetrahedralization_closure(model_data, cut_elements);
        for (const ID element_id : cut_elements) {
            if (c3d8_closure.find(element_id) != c3d8_closure.end()) {
                cut_regular_c3d8_elements.insert(element_id);
            }
        }
        std::map<ID, std::vector<ID>> c3d8_tetrahedra;
        std::map<ID, std::array<FaceKey, 6>> c3d8_parent_faces;
        for (const ID element_id : c3d8_closure) {
            ElementPtr& element =
                model_data.elements[static_cast<std::size_t>(element_id)];
            auto* c3d8 = dynamic_cast<C3D8*>(element.get());
            logging::error(c3d8 != nullptr && !is_c3d5(*c3d8),
                           "CutBuilder: invalid regular C3D8 closure element");
            logging::warning(dynamic_cast<C3D8R*>(c3d8) == nullptr,
                             "CutBuilder: C3D8R closure element ", element_id,
                             " is tetrahedralized to C3D4 and loses reduced integration");
            c3d8_parent_faces.emplace(element_id, c3d8_face_keys(*c3d8));
            c3d8_tetrahedra.emplace(
                element_id,
                tetrahedralize_c3d8(
                    model_data, *c3d8, c3d6_face_diagonal_cache));
        }
        validate_c3d8_tetrahedralization_closure(
            model_data, c3d8_parent_faces, c3d8_tetrahedra);
        logging::info(!c3d8_closure.empty(),
                      "PRETENSIONSECTION '", section.name,
                      "': tetrahedralized regular C3D8 closure with ",
                      c3d8_closure.size(), " element(s) for ",
                      cut_regular_c3d8_elements.size(), " cut seed(s)");

        const C3D4CutAnalysis post_tetra_analysis = analyze_c3d4_cuts(
            model_data, plane_point, axis,
            Precision(1000) * std::numeric_limits<Precision>::epsilon());
        Index post_tetra_aligned_count = 0;
        for (const TetCutPlan& plan : post_tetra_analysis.plans) {
            if (plan.type == TetCutType::VertexAligned ||
                plan.type == TetCutType::EdgeAligned ||
                plan.type == TetCutType::FaceAligned) {
                ++post_tetra_aligned_count;
            }
        }
        if (post_tetra_aligned_count > 0) {
            apply_aligned_c3d4_plans(
                model_data, section, post_tetra_analysis, plane_point, axis);
        }
        for (const TetCutPlan& plan : post_tetra_analysis.plans) {
            if (plan.type == TetCutType::OneThree || plan.type == TetCutType::TwoTwo) {
                split_c3d4_element(model_data, section, plan.element_id,
                                   plane_point, axis, tolerance,
                                   alternative_diagonal_count);
            }
        }
    }

    const std::set<ID> c3d6_closure = c3d6_tetrahedralization_closure(
        model_data, cut_elements);
    std::set<ID> cut_c3d6_elements;
    for (const ID element_id : cut_elements) {
        if (c3d6_closure.find(element_id) != c3d6_closure.end()) {
            cut_c3d6_elements.insert(element_id);
        }
    }
    std::map<ID, std::vector<ID>> c3d6_tetrahedra;
    for (const ID element_id : c3d6_closure) {
        ElementPtr& element = model_data.elements[static_cast<std::size_t>(element_id)];
        auto* c3d6 = dynamic_cast<C3D6*>(element.get());
        logging::error(c3d6 != nullptr,
                       "CutBuilder: invalid C3D6 closure element");
        c3d6_tetrahedra.emplace(
            element_id,
            tetrahedralize_c3d6(
                model_data, *c3d6, c3d6_face_diagonal_cache));
    }
    logging::info(!c3d6_closure.empty(),
                  "PRETENSIONSECTION '", section.name,
                  "': tetrahedralized C3D6 closure with ",
                  c3d6_closure.size(), " element(s) for ",
                  cut_c3d6_elements.size(), " cut seed(s)");
    for (const ID element_id : cut_c3d6_elements) {
        for (const ID tetra_id : c3d6_tetrahedra.at(element_id)) {
            split_c3d4_element(model_data, section, tetra_id,
                               plane_point, axis, tolerance,
                               alternative_diagonal_count);
        }
    }

    const std::set<ID> c3d5_closure = c3d5_c3d8_tetrahedralization_closure(
        model_data, cut_elements);
    std::set<ID> cut_c3d5_closure_elements;
    for (const ID element_id : cut_elements) {
        if (c3d5_closure.find(element_id) != c3d5_closure.end()) {
            cut_c3d5_closure_elements.insert(element_id);
        }
    }
    std::map<ID, std::vector<ID>> c3d5_tetrahedra;
    for (const ID element_id : c3d5_closure) {
        ElementPtr& element = model_data.elements[static_cast<std::size_t>(element_id)];
        auto* c3d8 = dynamic_cast<C3D8*>(element.get());
        logging::error(c3d8 != nullptr,
                       "CutBuilder: invalid C3D5/C3D8 closure element");
        logging::warning(is_c3d5(*c3d8),
                         "CutBuilder: regular C3D8 closure element ",
                         element_id,
                         dynamic_cast<C3D8R*>(c3d8) == nullptr
                             ? " is tetrahedralized to C3D4"
                             : " is tetrahedralized to C3D4 and loses reduced integration");
        c3d5_tetrahedra.emplace(
            element_id,
            is_c3d5(*c3d8)
                ? tetrahedralize_c3d5(
                    model_data, *c3d8, c3d6_face_diagonal_cache)
                : tetrahedralize_c3d8(
                    model_data, *c3d8, c3d6_face_diagonal_cache));
    }
    logging::info(!c3d5_closure.empty(),
                  "PRETENSIONSECTION '", section.name,
                  "': tetrahedralized C3D5/C3D8 closure with ",
                  c3d5_closure.size(), " element(s) for ",
                  cut_c3d5_closure_elements.size(), " cut element(s)");
    for (const ID element_id : cut_c3d5_closure_elements) {
        for (const ID tetra_id : c3d5_tetrahedra.at(element_id)) {
            split_c3d4_element(model_data, section, tetra_id,
                               plane_point, axis, tolerance,
                               alternative_diagonal_count);
        }
    }

    for (ID element_id : cut_elements) {
        if (cut_regular_c3d8_elements.find(element_id) !=
            cut_regular_c3d8_elements.end()) {
            continue;
        }
        if (cut_c3d6_elements.find(element_id) != cut_c3d6_elements.end()) {
            continue;
        }
        if (cut_c3d5_closure_elements.find(element_id) !=
            cut_c3d5_closure_elements.end()) {
            continue;
        }
        ElementPtr& element = model_data.elements[static_cast<std::size_t>(element_id)];
        auto* c3d6 = dynamic_cast<C3D6*>(element.get());
        if (c3d6 != nullptr) {
            logging::error(false,
                           "CutBuilder: C3D6 cut element missing from closure");
        }
        auto* c3d8 = dynamic_cast<C3D8*>(element.get());
        if (c3d8 != nullptr) {
            if (is_c3d5(*c3d8)) {
                logging::error(false,
                               "CutBuilder: C3D5 cut element missing from closure");
            }
            if (axis_aligned_cut) {
                split_c3d8_x_plane(model_data, section, *c3d8,
                                   plane_point, axis);
            } else {
                logging::warning(dynamic_cast<C3D8R*>(c3d8) == nullptr,
                                 "CutBuilder: non-axis-aligned C3D8R cut is "
                                 "tetrahedralized to C3D4 and cannot preserve "
                                 "reduced integration");
                const auto tetrahedra = tetrahedralize_c3d8(
                    model_data, *c3d8, c3d6_face_diagonal_cache);
                for (const ID tetra_id : tetrahedra) {
                    split_c3d4_element(model_data, section, tetra_id,
                                       plane_point, axis, tolerance,
                                       alternative_diagonal_count);
                }
            }
            continue;
        }

        auto* c3d4 = dynamic_cast<C3D4*>(element.get());
        if (c3d4 != nullptr) {
            split_c3d4_element(model_data, section, c3d4->elem_id,
                               plane_point, axis, tolerance,
                               alternative_diagonal_count);
            continue;
        }

        auto* c3d15 = dynamic_cast<C3D15*>(element.get());
        if (c3d15 != nullptr) {
            logging::warning(false,
                             "CutBuilder: C3D15 cut element ", element_id,
                             " is tetrahedralized to three C3D4 elements; "
                             "quadratic order is not preserved");
            const auto tetrahedra = tetrahedralize_c3d15(
                model_data, *c3d15, c3d6_face_diagonal_cache);
            for (const ID tetra_id : tetrahedra) {
                split_c3d4_element(model_data, section, tetra_id,
                                   plane_point, axis, tolerance,
                                   alternative_diagonal_count);
            }
            continue;
        }

        auto* c3d13 = dynamic_cast<C3D13*>(element.get());
        if (c3d13 != nullptr) {
            logging::warning(false,
                             "CutBuilder: C3D13 cut element ", element_id,
                             " is tetrahedralized to two C3D4 elements; "
                             "quadratic order is not preserved");
            const auto tetrahedra = tetrahedralize_c3d13(
                model_data, *c3d13, c3d6_face_diagonal_cache);
            for (const ID tetra_id : tetrahedra) {
                split_c3d4_element(model_data, section, tetra_id,
                                   plane_point, axis, tolerance,
                                   alternative_diagonal_count);
            }
            continue;
        }

        auto* c3d10 = dynamic_cast<C3D10*>(element.get());
        if (c3d10 != nullptr) {
            split_c3d10_element(model_data, section, element_id,
                                plane_point, axis, tolerance);
            continue;
        }

        logging::error(false,
                       "CutBuilder: only C3D4, C3D5, C3D6, C3D8, C3D10, C3D13, C3D15, C3D20 and C3D20R cuts are implemented yet");
    }

    log_c3d4_interface_quality(model_data, section);
    log_c3d10_interface_quality(model_data, section);
    validate_interface_pair_connectivity(model_data, section);

    logging::info(alternative_diagonal_count > 0,
                  "PRETENSIONSECTION '", section.name,
                  "': selected the alternative quality-optimized diagonal for ",
                  alternative_diagonal_count, " C3D4 2/2 cut(s)");

    logging::info(true,
                  "PRETENSIONSECTION '", section.name,
                  "': detected ", section.cut_element_ids.size(),
                  " cut element(s) and created ",
                  section.side_a_nodes.size(),
                  " side-A / ",
                  section.side_b_nodes.size(),
                  " side-B interface nodes");
}

} // namespace fem::model
