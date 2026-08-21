#pragma once

#include "../core/types_eig.h"
#include "../core/types_num.h"

#include <array>
#include <map>
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

struct InterfacePair {
    ID side_a = -1;
    ID side_b = -1;
};

struct PretensionSection {
    using Ptr = std::shared_ptr<PretensionSection>;

    std::string name;

    std::string cylinder_surface_set;
    std::string interface_surface_set_a;
    std::string interface_surface_set_b;

    Vec3 axis_origin = Vec3::Zero();
    Vec3 axis_direction = Vec3::UnitX();
    Precision cut_coordinate = 0;
    Precision snap_ratio = Precision(0.02);

    std::vector<ID> side_a_nodes;
    std::vector<ID> side_b_nodes;
    std::vector<InterfacePair> interface_pairs;
    std::map<std::array<long long, 3>, InterfacePair> interface_node_cache;
    std::vector<ID> cut_element_ids;
    std::map<ID, ID> cut_quality_origins;
    std::map<ID, std::string> cut_quality_types;
    std::vector<CutFacet> facets;

    Control control = Control::Displacement;
    State state = State::Open;

    Precision prescribed_value = 0;
    Precision locked_gap = 0;
    Precision last_solved_gap = 0;
    bool has_solved_gap = false;

    bool prepared = false;

    bool uses_surface_pair() const {
        return !interface_surface_set_a.empty() &&
               !interface_surface_set_b.empty();
    }
};

}
