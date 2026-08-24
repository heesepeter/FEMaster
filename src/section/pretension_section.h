#pragma once

#include "../core/types_eig.h"
#include "../core/types_num.h"

#include <memory>
#include <string>
#include <vector>

namespace fem::pretension {

enum class Control {
    Force,
    Displacement
};

enum class State {
    Open,
    Loading,
    Locked
};

struct InterfacePair {
    ID side_a = -1;
    ID side_b = -1;
};

struct PretensionSection {
    using Ptr = std::shared_ptr<PretensionSection>;

    std::string name;

    std::string interface_surface_set_a;
    std::string interface_surface_set_b;

    Vec3 axis_origin = Vec3::Zero();
    Vec3 axis_direction = Vec3::UnitX();
    std::vector<ID> side_a_nodes;
    std::vector<ID> side_b_nodes;
    std::vector<InterfacePair> interface_pairs;
    Control control = Control::Displacement;
    State state = State::Open;

    Precision prescribed_value = 0;
    Precision locked_gap = 0;
    Precision last_solved_gap = 0;
    bool has_solved_gap = false;

    bool prepared = false;

};

}
