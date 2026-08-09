/**
 * @file pretension.cpp
 * @brief Builds axial compatibility equations for pretension sections.
 */

#include "pretension.h"

#include "../../core/logging.h"

#include <array>
#include <cmath>
#include <utility>

namespace fem::constraint {

Equations get_pretension_equations(
    SystemDofIds& system_dof_ids,
    model::ModelData& model_data,
    pretension::PretensionSection& section) {
    (void) model_data;

    Equations equations;
    if (section.state == pretension::State::Open) {
        return equations;
    }

    // Force control is represented by equal and opposite nodal loads on the
    // two cut sides. It must not also create a prescribed-gap equation.
    if (section.control == pretension::Control::Force) {
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
    constexpr Precision coefficient_tolerance = Precision(1e-12);

    Vec3 reference = std::abs(axis(0)) < Precision(0.9)
        ? Vec3::UnitX() : Vec3::UnitY();
    Vec3 tangent_1 = axis.cross(reference).normalized();
    Vec3 tangent_2 = axis.cross(tangent_1).normalized();
    const std::array<Vec3, 3> bases = {tangent_1, tangent_2, axis};
    const Precision normal_gap = section.state == pretension::State::Loading
        ? section.prescribed_value
        : section.locked_gap;
    const std::array<Precision, 3> rhs_values = {Precision(0), Precision(0), normal_gap};

    for (Index basis_index = 0; basis_index < bases.size(); ++basis_index) {
        const Vec3& basis = bases[basis_index];

        for (const auto& pair : section.interface_pairs) {
            std::vector<EquationEntry> entries;
            entries.reserve(6);
            for (Dim component = 0; component < 3; ++component) {
                const Precision coefficient = basis(component);
                if (std::abs(coefficient) <= coefficient_tolerance) {
                    continue;
                }

                const bool side_a_active =
                    pair.side_a >= 0 &&
                    pair.side_a < system_dof_ids.rows() &&
                    system_dof_ids(pair.side_a, component) >= 0;
                const bool side_b_active =
                    pair.side_b >= 0 &&
                    pair.side_b < system_dof_ids.rows() &&
                    system_dof_ids(pair.side_b, component) >= 0;

                if (side_a_active) {
                    entries.push_back({pair.side_a, component, -coefficient});
                }
                if (side_b_active) {
                    entries.push_back({pair.side_b, component, coefficient});
                }
            }
            if (!entries.empty()) {
                Equation equation(std::move(entries), rhs_values[basis_index]);
                equation.source = EquationSourceKind::Manual;
                equations.push_back(std::move(equation));
            }
        }
    }

    return equations;
}

} // namespace fem::constraint
