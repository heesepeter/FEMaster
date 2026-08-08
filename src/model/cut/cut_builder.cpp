/**
 * @file cut_builder.cpp
 * @brief Builds topological cuts for pretension sections.
 */

#include "cut_builder.h"
#include "../element/element.h"
#include "../solid/c3d8.h"

#include <algorithm>
#include <array>
#include <limits>
#include <utility>

namespace fem::model {

namespace {

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
    C3D8& element) {
    const ID element_id = element.elem_id;
    const auto old_nodes = element.node_ids;
    const Precision plane_x = section.axis_origin(0) + section.cut_coordinate;

    const std::array<ID, 4> negative_nodes = {
        old_nodes[0], old_nodes[3], old_nodes[4], old_nodes[7]};
    const std::array<ID, 4> positive_nodes = {
        old_nodes[1], old_nodes[2], old_nodes[5], old_nodes[6]};

    for (ID node_id : negative_nodes) {
        const Vec3 position = model_data.positions->row_vec3(static_cast<Index>(node_id));
        logging::error(position(0) < plane_x,
                       "CutBuilder: C3D8 node ordering is not compatible with the X-plane split");
    }
    for (ID node_id : positive_nodes) {
        const Vec3 position = model_data.positions->row_vec3(static_cast<Index>(node_id));
        logging::error(position(0) > plane_x,
                       "CutBuilder: C3D8 node ordering is not compatible with the X-plane split");
    }

    const std::array<std::pair<ID, ID>, 4> crossing_edges = {{
        {old_nodes[0], old_nodes[1]},
        {old_nodes[3], old_nodes[2]},
        {old_nodes[4], old_nodes[5]},
        {old_nodes[7], old_nodes[6]}
    }};

    std::array<ID, 4> side_a_interface{};
    std::array<ID, 4> side_b_interface{};

    for (Index i = 0; i < crossing_edges.size(); ++i) {
        const auto [node_a, node_b] = crossing_edges[i];
        const Vec3 position_a = model_data.positions->row_vec3(static_cast<Index>(node_a));
        const Vec3 position_b = model_data.positions->row_vec3(static_cast<Index>(node_b));
        const Vec3 intersection = (position_a + position_b) * Precision(0.5);

        side_a_interface[i] = model_data.append_node(intersection);
        side_b_interface[i] = model_data.append_node(intersection);

        section.interface_pairs.push_back({
            side_a_interface[i],
            side_b_interface[i]});
    }

    const std::array<ID, 8> left_nodes = {
        old_nodes[0], side_a_interface[0], side_a_interface[1], old_nodes[3],
        old_nodes[4], side_a_interface[2], side_a_interface[3], old_nodes[7]};

    const std::array<ID, 8> right_nodes = {
        side_b_interface[0], old_nodes[1], old_nodes[2], side_b_interface[1],
        side_b_interface[2], old_nodes[5], old_nodes[6], side_b_interface[3]};

    const ID new_element_id = model_data.next_free_element_id();
    add_element_to_matching_sets(model_data, element_id, new_element_id);

    auto left_element = std::make_shared<C3D8>(element_id, left_nodes);
    left_element->_model_data = &model_data;
    model_data.elements[static_cast<std::size_t>(element_id)] = left_element;

    auto right_element = std::make_shared<C3D8>(new_element_id, right_nodes);
    model_data.insert_element(std::move(right_element));

    section.side_a_nodes.insert(
        section.side_a_nodes.end(), side_a_interface.begin(), side_a_interface.end());
    section.side_b_nodes.insert(
        section.side_b_nodes.end(), side_b_interface.begin(), side_b_interface.end());
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

    logging::error(axis.isApprox(Vec3::UnitX(), Precision(1e-12)),
                   "CutBuilder: only planes normal to global X are supported yet");

    const auto cut_elements = section.cut_element_ids;
    for (ID element_id : cut_elements) {
        ElementPtr& element = model_data.elements[static_cast<std::size_t>(element_id)];
        auto* c3d8 = dynamic_cast<C3D8*>(element.get());
        logging::error(c3d8 != nullptr,
                       "CutBuilder: only C3D8 elements are supported yet");

        split_c3d8_x_plane(model_data, section, *c3d8);
    }

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
