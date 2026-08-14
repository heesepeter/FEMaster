/**
 * @file cut_builder.cpp
 * @brief Builds topological cuts for pretension sections.
 */

#include "cut_builder.h"
#include "../element/element.h"
#include "../solid/c3d4.h"
#include "../solid/c3d8.h"
#include "../solid/c3d10.h"

#include <algorithm>
#include <array>
#include <cmath>
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
    bool reverse_sides = false) {
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
    if (reverse_sides) {
        std::swap(pair.side_a, pair.side_b);
    }
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
    ID minimum_volume_element = -1;
    ID minimum_edge_element = -1;
    Index checked_elements = 0;

    for (const ElementPtr& element : model_data.elements) {
        auto* c3d4 = element == nullptr ? nullptr : dynamic_cast<C3D4*>(element.get());
        if (c3d4 == nullptr) continue;

        bool touches_interface = false;
        for (const ID node_id : c3d4->node_ids) {
            touches_interface = touches_interface || interface_nodes.count(node_id) > 0;
        }
        if (!touches_interface) continue;

        const auto& nodes = c3d4->node_ids;
        const Vec3 p0 = model_data.positions->row_vec3(static_cast<Index>(nodes[0]));
        const Vec3 p1 = model_data.positions->row_vec3(static_cast<Index>(nodes[1]));
        const Vec3 p2 = model_data.positions->row_vec3(static_cast<Index>(nodes[2]));
        const Vec3 p3 = model_data.positions->row_vec3(static_cast<Index>(nodes[3]));
        const Precision volume = std::abs((p1 - p0).dot((p2 - p0).cross(p3 - p0))) /
            Precision(6);
        minimum_volume = std::min(minimum_volume, volume);
        if (volume == minimum_volume) minimum_volume_element = c3d4->elem_id;

        constexpr std::array<std::pair<Index, Index>, 6> edges = {{
            {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
        }};
        for (const auto [a, b] : edges) {
            const Precision edge = (model_data.positions->row_vec3(
                static_cast<Index>(nodes[a])) - model_data.positions->row_vec3(
                static_cast<Index>(nodes[b]))).norm();
            if (edge < minimum_edge) {
                minimum_edge = edge;
                minimum_edge_element = c3d4->elem_id;
            }
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
                  " (element ", minimum_edge_element, ")");
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
}

void split_c3d8_x_plane(
    ModelData& model_data,
    pretension::PretensionSection& section,
    C3D8& element,
    const Vec3& plane_point,
    const Vec3& axis) {
    const ID element_id = element.elem_id;
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

    auto left_element = std::make_shared<C3D8>(element_id, left_nodes);
    left_element->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = left_element;

    auto right_element = std::make_shared<C3D8>(new_element_id, right_nodes);
    model_data.insert_element(std::move(right_element));

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
                   "CutBuilder: generated C3D4 has zero volume");
    return nodes;
}

std::vector<ID> tetrahedralize_c3d8(
    ModelData& model_data,
    C3D8& element) {
    const ID element_id = element.elem_id;
    const auto& n = element.node_ids;
    const std::array<std::array<ID, 4>, 5> raw_tets = {{
        {n[0], n[1], n[3], n[4]},
        {n[1], n[2], n[3], n[6]},
        {n[1], n[3], n[4], n[6]},
        {n[1], n[4], n[5], n[6]},
        {n[3], n[4], n[6], n[7]}
    }};

    std::array<std::array<ID, 4>, 5> tets{};
    for (Index i = 0; i < tets.size(); ++i) {
        tets[i] = oriented_c3d4_nodes(model_data, raw_tets[i]);
    }

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

    auto side_a_element = std::make_shared<C3D4>(element_id, side_a_tets[0]);
    side_a_element->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = side_a_element;

    for (Index i = 1; i < side_a_tets.size(); ++i) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D4>(new_id, side_a_tets[i]));
    }
    const ID side_b_id = model_data.next_free_element_id();
    add_element_to_matching_sets(model_data, element_id, side_b_id);
    model_data.insert_element(std::make_shared<C3D4>(side_b_id, side_b_tet));

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
    const Vec3& axis) {
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

    std::array<std::array<ID, 4>, 3> side_a_tets = {{
        {old_nodes[negative[0]], old_nodes[negative[1]], side_a_interface[0][0], side_a_interface[0][1]},
        {old_nodes[negative[1]], side_a_interface[0][0], side_a_interface[0][1], side_a_interface[1][1]},
        {old_nodes[negative[1]], side_a_interface[0][0], side_a_interface[1][1], side_a_interface[1][0]}
    }};
    std::array<std::array<ID, 4>, 3> side_b_tets = {{
        {old_nodes[positive[0]], old_nodes[positive[1]], side_b_interface[0][0], side_b_interface[1][0]},
        {old_nodes[positive[1]], side_b_interface[0][0], side_b_interface[1][0], side_b_interface[1][1]},
        {old_nodes[positive[1]], side_b_interface[0][0], side_b_interface[1][1], side_b_interface[0][1]}
    }};

    for (auto& tet : side_a_tets) {
        tet = oriented_c3d4_nodes(model_data, tet);
    }
    for (auto& tet : side_b_tets) {
        tet = oriented_c3d4_nodes(model_data, tet);
    }

    auto side_a_element = std::make_shared<C3D4>(element_id, side_a_tets[0]);
    side_a_element->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = side_a_element;

    for (Index i = 1; i < side_a_tets.size(); ++i) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D4>(new_id, side_a_tets[i]));
    }
    for (const auto& tet : side_b_tets) {
        const ID new_id = model_data.next_free_element_id();
        add_element_to_matching_sets(model_data, element_id, new_id);
        model_data.insert_element(std::make_shared<C3D4>(new_id, tet));
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
    Precision tolerance) {
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
                                 plane_point, axis);
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
    for (ID element_id : cut_elements) {
        ElementPtr& element = model_data.elements[static_cast<std::size_t>(element_id)];
        auto* c3d8 = dynamic_cast<C3D8*>(element.get());
        if (c3d8 != nullptr) {
            const bool axis_aligned =
                axis.isApprox(Vec3::UnitX(), Precision(1e-12)) ||
                axis.isApprox(-Vec3::UnitX(), Precision(1e-12)) ||
                axis.isApprox(Vec3::UnitY(), Precision(1e-12)) ||
                axis.isApprox(-Vec3::UnitY(), Precision(1e-12)) ||
                axis.isApprox(Vec3::UnitZ(), Precision(1e-12)) ||
                axis.isApprox(-Vec3::UnitZ(), Precision(1e-12));
            if (axis_aligned) {
                split_c3d8_x_plane(model_data, section, *c3d8,
                                   plane_point, axis);
            } else {
                const auto tetrahedra = tetrahedralize_c3d8(model_data, *c3d8);
                for (const ID tetra_id : tetrahedra) {
                    split_c3d4_element(model_data, section, tetra_id,
                                       plane_point, axis, tolerance);
                }
            }
            continue;
        }

        auto* c3d4 = dynamic_cast<C3D4*>(element.get());
        if (c3d4 != nullptr) {
            split_c3d4_element(model_data, section, c3d4->elem_id,
                               plane_point, axis, tolerance);
            continue;
        }

        auto* c3d10 = dynamic_cast<C3D10*>(element.get());
        if (c3d10 != nullptr) {
            Index positive_count = 0;
            for (Index local_node = 0; local_node < 4; ++local_node) {
                const Vec3 position = model_data.positions->row_vec3(
                    static_cast<Index>(c3d10->node_ids[local_node]));
                if (axis.dot(position - plane_point) > tolerance) {
                    ++positive_count;
                }
            }
            if (positive_count == 2) {
                split_c3d10_four_edge_cut(model_data, section, *c3d10,
                                          plane_point, axis);
            } else if (positive_count == 3) {
                const std::size_t pair_start = section.interface_pairs.size();
                const std::size_t node_start = section.side_a_nodes.size();
                split_c3d10_three_edge_cut(model_data, section, *c3d10,
                                           plane_point, -axis);
                for (std::size_t i = pair_start; i < section.interface_pairs.size(); ++i) {
                    std::swap(section.interface_pairs[i].side_a,
                              section.interface_pairs[i].side_b);
                }
                for (std::size_t i = node_start; i < section.side_a_nodes.size(); ++i) {
                    std::swap(section.side_a_nodes[i], section.side_b_nodes[i]);
                }
            } else {
                split_c3d10_three_edge_cut(model_data, section, *c3d10,
                                           plane_point, axis);
            }
            continue;
        }

        logging::error(false,
                       "CutBuilder: only C3D8, C3D4 and supported C3D10 cuts are implemented yet");
    }

    log_c3d4_interface_quality(model_data, section);

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
