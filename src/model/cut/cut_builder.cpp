/**
 * @file cut_builder.cpp
 * @brief Builds topological cuts for pretension sections.
 */

#include "cut_builder.h"
#include "../element/element.h"

#include <algorithm>
#include <limits>

namespace fem::model {

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

    // The cut is considered to pass through an element only when the element
    // has nodes strictly on both sides of the plane. Nodes lying exactly on
    // the plane alone do not constitute a cut.
    const Precision tolerance = Precision(100) * std::numeric_limits<Precision>::epsilon();

    section.cut_element_ids.clear();

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
        logging::info(
            true,
            "PRETENSIONSECTION '",
            section.name,
            "': detected ",
            section.cut_element_ids.size(),
            " cut element(s)");
    }
}

} // namespace fem::model
