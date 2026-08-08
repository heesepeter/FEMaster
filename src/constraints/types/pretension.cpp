/**
 * @file pretension.cpp
 * @brief Builds axial compatibility equations for pretension sections.
 */

#include "pretension.h"

#include "../../core/logging.h"

#include <cmath>
#include <utility>

namespace fem::constraint {

Equations get_pretension_equations(
    SystemDofIds& system_dof_ids,
    model::ModelData& model_data,
    pretension::PretensionSection& section) {
    (void) model_data;

    Equations equations;
    if (section.state != pretension::State::Locked) {
        return equations;
    }

    logging::error(!section.interface_pairs.empty(),
                   "Pretension section '", section.name,
                   "' has no interface node pairs");

    const Precision axis_norm = section.axis_direction.norm();
    logging::error(axis_norm > Precision(0),
                   "Pretension section '", section.name,
                   "' has a zero axis direction");

    const Vec3 axis = section.axis_direction / axis_norm;
    const Precision weight =
        Precision(1) / static_cast<Precision>(section.interface_pairs.size());
    constexpr Precision coefficient_tolerance = Precision(1e-12);

    std::vector<EquationEntry> entries;
    entries.reserve(section.interface_pairs.size() * 6);

    for (const auto& pair : section.interface_pairs) {
        for (Dim component = 0; component < 3; ++component) {
            const Precision coefficient = weight * axis(component);
            if (std::abs(coefficient) <= coefficient_tolerance) {
                continue;
            }

            const Dim dof = component;
            const bool side_a_active =
                pair.side_a >= 0 &&
                pair.side_a < system_dof_ids.rows() &&
                system_dof_ids(pair.side_a, dof) >= 0;
            const bool side_b_active =
                pair.side_b >= 0 &&
                pair.side_b < system_dof_ids.rows() &&
                system_dof_ids(pair.side_b, dof) >= 0;

            if (side_a_active) {
                entries.push_back({pair.side_a, dof, -coefficient});
            }
            if (side_b_active) {
                entries.push_back({pair.side_b, dof, coefficient});
            }
        }
    }

    if (!entries.empty()) {
        Equation equation(std::move(entries), section.locked_gap);
        equation.source = EquationSourceKind::Manual;
        equations.push_back(std::move(equation));
    }

    return equations;
}

} // namespace fem::constraint
