#include "model.h"

namespace fem::model {

Field Model::build_pretension_force_matrix() {
    Field force{"PRETENSION_FORCE", FieldDomain::NODE,
                _data->field_rows(FieldDomain::NODE), 3};
    force.set_zero();
    for (const auto& section : _data->pretension_sections) {
        if (!section || section->state != pretension::State::Loading ||
            section->control != pretension::Control::Force ||
            section->interface_pairs.empty()) continue;
        const Precision norm = section->axis_direction.norm();
        logging::error(norm > Precision(0),
            "Pretension section '", section->name, "' has a zero axis direction");
        const Vec3 axis = section->axis_direction / norm;
        const Precision pair_force = section->prescribed_value /
            static_cast<Precision>(section->interface_pairs.size());
        for (const auto& pair : section->interface_pairs) {
            for (Dim component = 0; component < 3; ++component) {
                force(static_cast<Index>(pair.side_a), component) -= pair_force * axis(component);
                force(static_cast<Index>(pair.side_b), component) += pair_force * axis(component);
            }
        }
    }
    return force;
}

Field Model::build_pretension_gap_matrix(const Field& displacement) {
    logging::error(displacement.domain == FieldDomain::NODE && displacement.components >= 3,
        "Pretension gap output requires a nodal displacement field");
    Field gap{"PRETENSION_GAP", FieldDomain::NODE,
              _data->field_rows(FieldDomain::NODE), 1};
    gap.set_zero();
    for (const auto& section : _data->pretension_sections) {
        if (!section || section->interface_pairs.empty()) continue;
        const Vec3 axis = section->axis_direction.normalized();
        for (const auto& pair : section->interface_pairs) {
            const Vec3 relative =
                displacement.row_vec3(static_cast<Index>(pair.side_b)) -
                displacement.row_vec3(static_cast<Index>(pair.side_a));
            const Precision measured = axis.dot(relative);
            const Precision value = section->control == pretension::Control::Displacement
                ? (section->state == pretension::State::Locked
                    ? section->locked_gap : section->prescribed_value)
                : measured;
            gap(static_cast<Index>(pair.side_a), 0) = value;
        }
    }
    return gap;
}

void Model::capture_pretension_gaps(const Field& displacement) {
    logging::error(displacement.domain == FieldDomain::NODE && displacement.components >= 3,
        "Capturing pretension gaps requires a nodal displacement field");
    for (const auto& section : _data->pretension_sections) {
        if (!section || section->state == pretension::State::Open ||
            section->interface_pairs.empty()) continue;
        const Vec3 axis = section->axis_direction.normalized();
        Precision sum = 0;
        Index count = 0;
        for (const auto& pair : section->interface_pairs) {
            const Vec3 relative =
                displacement.row_vec3(static_cast<Index>(pair.side_b)) -
                displacement.row_vec3(static_cast<Index>(pair.side_a));
            sum += axis.dot(relative);
            ++count;
        }
        logging::error(count > 0,
            "Pretension section '", section->name,
            "' has no valid interface pair for gap capture");
        section->last_solved_gap = sum / static_cast<Precision>(count);
        section->has_solved_gap = true;
        logging::info(true, "PRETENSION '", section->name,
            "': solved mean gap = ", section->last_solved_gap);
    }
}

} // namespace fem::model
