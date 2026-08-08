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

struct CutFacet {
    std::vector<ID> node_ids;
    Precision area = 0;
};

struct PretensionSection {
    using Ptr = std::shared_ptr<PretensionSection>;

    std::string name;

    std::string cylinder_surface_set;

    Vec3 axis_origin = Vec3::Zero();
    Vec3 axis_direction = Vec3::UnitX();
    Precision cut_coordinate = 0;

    std::vector<ID> side_a_nodes;
    std::vector<ID> side_b_nodes;
    std::vector<ID> cut_element_ids;
    std::vector<CutFacet> facets;

    Control control = Control::Displacement;
    State state = State::Open;

    Precision prescribed_value = 0;
    Precision locked_gap = 0;

    bool prepared = false;
};

}
