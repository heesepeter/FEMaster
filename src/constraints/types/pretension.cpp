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
    const std::array<Vec3, 2> tangents = {tangent_1, tangent_2};
    const Precision normal_gap = section.state == pretension::State::Loading
        ? section.prescribed_value
        : section.locked_gap;

    const auto append_projection = [&](std::vector<EquationEntry>& entries,
                                       const pretension::InterfacePair& pair,
                                       const Vec3& basis,
                                       Precision scale) {
        for (Dim component = 0; component < 3; ++component) {
            const Precision coefficient = scale * basis(component);
            if (std::abs(coefficient) <= coefficient_tolerance) continue;
            if (pair.side_a >= 0 && pair.side_a < system_dof_ids.rows() &&
                system_dof_ids(pair.side_a, component) >= 0) {
                entries.push_back({pair.side_a, component, -coefficient});
            }
            if (pair.side_b >= 0 && pair.side_b < system_dof_ids.rows() &&
                system_dof_ids(pair.side_b, component) >= 0) {
                entries.push_back({pair.side_b, component, coefficient});
            }
        }
    };
    const auto add_equation = [&](std::vector<EquationEntry> entries,
                                  Precision rhs) {
        if (entries.empty()) return;
        Equation equation(std::move(entries), rhs);
        equation.source = EquationSourceKind::Manual;
        equations.push_back(std::move(equation));
    };

    // The two faces remain tied in both tangential directions for every
    // control mode.
    for (const Vec3& tangent : tangents) {
        for (const auto& pair : section.interface_pairs) {
            std::vector<EquationEntry> entries;
            entries.reserve(6);
            append_projection(entries, pair, tangent, Precision(1));
            add_equation(std::move(entries), Precision(0));
        }
    }

    if (section.control == pretension::Control::Force) {
        // Leave one collective axial gap free and constrain all remaining
        // pair gaps to it. The distributed nodal forces then act on exactly
        // this single generalized pretension displacement.
        const auto& reference_pair = section.interface_pairs.front();
        for (Index i = 1; i < section.interface_pairs.size(); ++i) {
            std::vector<EquationEntry> entries;
            entries.reserve(12);
            append_projection(entries, section.interface_pairs[i], axis, Precision(1));
            append_projection(entries, reference_pair, axis, Precision(-1));
            add_equation(std::move(entries), Precision(0));
        }
    } else {
        // Displacement loading and LOCK prescribe the collective gap.
        for (const auto& pair : section.interface_pairs) {
            std::vector<EquationEntry> entries;
            entries.reserve(6);
            append_projection(entries, pair, axis, Precision(1));
            add_equation(std::move(entries), normal_gap);
        }
    }

    return equations;
}

} // namespace fem::constraint
