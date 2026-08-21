#include "model.h"

#include "cut/cut_builder.h"

#include "../section/section.h"
#include "../section/section_beam.h"
#include "../section/section_solid.h"
#include "../section/section_shell_abd.h"
#include "../section/section_shell_integrated.h"
#include "../section/section_truss.h"
#include "../feature/point_mass.h"

#include "../bc/load_collector.h"
#include "../bc/support_collector.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <set>

namespace fem {
namespace model {
void Model::step_begin() {
    try {
        for (auto& elem : _data->elements) {
            if (elem != nullptr) {
                if (auto* structural = elem->as<StructuralElement>()) {
                    structural->step_begin();
                }
            }
        }
    } catch (...) {
        step_end();
        throw;
    }
}

void Model::step_end() {
    for (auto& elem : _data->elements) {
        if (elem != nullptr) {
            if (auto* structural = elem->as<StructuralElement>()) {
                structural->step_end();
            }
        }
    }
}

void Model::add_connector(const std::string& set1,
                          const std::string& set2,
                          const std::string& coordinate_system,
                          constraint::ConnectorType type) {
    logging::error(_data->node_sets.has(set1), "Node set ", set1, " does not exist");
    logging::error(_data->node_sets.has(set2), "Node set ", set2, " does not exist");

    logging::error(_data->coordinate_systems.has(coordinate_system), "Coordinate system ", coordinate_system, " does not exist");

    logging::error(_data->node_sets.get(set1)->size() == 1, "Set 1 must contain exactly one node");
    logging::error(_data->node_sets.get(set2)->size() == 1, "Set 2 must contain exactly one node");

    ID id1 = _data->node_sets.get(set1)->first();
    ID id2 = _data->node_sets.get(set2)->first();

    _data->connectors.emplace_back(id1, id2, _data->coordinate_systems.get(coordinate_system), type);
}

void Model::add_coupling(const std::string &master_set, const std::string &slave_set, Dofs coupled_dofs, constraint::CouplingType type, bool is_surface) {
    logging::error(_data->node_sets.get(master_set)->size() == 1, "Master set must contain exactly one node");

    if (is_surface) {
        logging::error(_data->surface_sets.has(slave_set), "Slave set ", slave_set, " is not a defined surface set");
    } else {
        logging::error(_data->node_sets.has(slave_set), "Slave set ", slave_set, " is not a defined node set");
    }

    ID master_node = _data->node_sets.get(master_set)->first();

    if (is_surface) {
        _data->couplings.emplace_back(master_node, _data->surface_sets.get(slave_set), coupled_dofs, type);
    } else {
        _data->couplings.emplace_back(master_node, _data->node_sets.get(slave_set), coupled_dofs, type);
    }
    // Attach master region pointer for reporting
    if (!_data->couplings.empty()) {
        auto& c = _data->couplings.back();
        c.master_region = _data->node_sets.get(master_set);
    }
}

void Model::add_tie(const std::string& master_set,
                    const std::string& slave_set,
                    Precision distance,
                    bool adjust) {
    // ---------------------------------------------------------------------
    // Resolve slave: prefer node set, otherwise surface set, otherwise error
    // ---------------------------------------------------------------------
    NodeRegion::Ptr    slave_node_ptr    = nullptr;
    SurfaceRegion::Ptr slave_surface_ptr = nullptr;

    if (_data->node_sets.has(slave_set) && _data->node_sets.get(slave_set) &&
        _data->node_sets.get(slave_set)->size() > 0) {
        slave_node_ptr = _data->node_sets.get(slave_set);
    } else if (_data->surface_sets.has(slave_set) && _data->surface_sets.get(slave_set) &&
               _data->surface_sets.get(slave_set)->size() > 0) {
        slave_surface_ptr = _data->surface_sets.get(slave_set);
    } else {
        logging::error(false,
                       "Slave set ", slave_set,
                       " is neither a defined non-empty node set nor a defined non-empty surface set");
    }

    // ---------------------------------------------------------------------
    // Resolve master: surfaces first, then lines (as before)
    // ---------------------------------------------------------------------
    const bool has_surfaces =
        _data->surface_sets.has(master_set) && _data->surface_sets.get(master_set) &&
        _data->surface_sets.get(master_set)->size() > 0;

    const bool has_lines =
        _data->line_sets.has(master_set) && _data->line_sets.get(master_set) &&
        _data->line_sets.get(master_set)->size() > 0;

    if (has_surfaces) {
        SurfaceRegion::Ptr master_ptr = _data->surface_sets.get(master_set);

        if (slave_node_ptr) {
            _data->ties.emplace_back(master_ptr, slave_node_ptr, distance, adjust);
        } else {
            _data->ties.emplace_back(master_ptr, slave_surface_ptr, distance, adjust);
        }
        return;
    }

    if (has_lines) {
        LineRegion::Ptr master_line_ptr = _data->line_sets.get(master_set);

        if (slave_node_ptr) {
            _data->ties.emplace_back(master_line_ptr, slave_node_ptr, distance, adjust);
        } else {
            _data->ties.emplace_back(master_line_ptr, slave_surface_ptr, distance, adjust);
        }
        return;
    }

    logging::error(false, "Master set ", master_set, " contains neither surfaces nor lines");
}

void Model::add_contact(const std::string& master_set,
                        const std::string& slave_set,
                        Precision distance,
                        Precision penalty,
                        Precision clearance,
                        bool flip_normal) {
    NodeRegion::Ptr    slave_node_ptr    = nullptr;
    SurfaceRegion::Ptr slave_surface_ptr = nullptr;

    if (_data->node_sets.has(slave_set) && _data->node_sets.get(slave_set) &&
        _data->node_sets.get(slave_set)->size() > 0) {
        slave_node_ptr = _data->node_sets.get(slave_set);
    } else if (_data->surface_sets.has(slave_set) && _data->surface_sets.get(slave_set) &&
               _data->surface_sets.get(slave_set)->size() > 0) {
        slave_surface_ptr = _data->surface_sets.get(slave_set);
    } else {
        logging::error(false,
                       "CONTACT: slave set ", slave_set,
                       " is neither a defined non-empty node set nor a defined non-empty surface set");
    }

    logging::error(_data->surface_sets.has(master_set) &&
                   _data->surface_sets.get(master_set) &&
                   _data->surface_sets.get(master_set)->size() > 0,
                   "CONTACT: master set ", master_set, " is not a defined non-empty surface set");

    SurfaceRegion::Ptr master_ptr = _data->surface_sets.get(master_set);
    if (slave_node_ptr) {
        _data->contacts.emplace_back(master_ptr,
                                     slave_node_ptr,
                                     distance,
                                     penalty,
                                     clearance,
                                     flip_normal);
    } else {
        _data->contacts.emplace_back(master_ptr,
                                     slave_surface_ptr,
                                     distance,
                                     penalty,
                                     clearance,
                                     flip_normal);
    }
}

void Model::add_pretension_section(
    const std::string& name,
    const std::string& surface_set,
    Vec3 axis,
    const std::string& position,
    Precision snap_ratio) {
    logging::error(!name.empty(),
                   "PRETENSION SECTION: name must not be empty");
    const bool duplicate_name = std::any_of(
        _data->pretension_sections.begin(),
        _data->pretension_sections.end(),
        [&](const pretension::PretensionSection::Ptr& existing) {
            return existing && existing->name == name;
        });
    logging::error(!duplicate_name,
                   "PRETENSION SECTION: duplicate section name ", name);
    logging::error(_data->surface_sets.has(surface_set),
                   "PRETENSION SECTION: surface set ", surface_set,
                   " does not exist");
    logging::error(axis.norm() > Precision(0),
                   "PRETENSION SECTION: axis must not be zero");
    logging::error(position == "MIDDLE",
                   "PRETENSION SECTION: only POSITION=MIDDLE is supported yet");
    logging::error(snap_ratio >= Precision(0),
                   "PRETENSION SECTION: SNAP must not be negative");
    logging::error(snap_ratio <= Precision(0.1),
                   "PRETENSION SECTION: SNAP must not exceed 0.1");

    auto section = std::make_shared<pretension::PretensionSection>();
    section->name = name;
    section->cylinder_surface_set = surface_set;
    section->axis_direction = axis.normalized();
    section->snap_ratio = snap_ratio;
    Vec3 surface_center = Vec3::Zero();
    Precision minimum_axis_coordinate = std::numeric_limits<Precision>::max();
    Precision maximum_axis_coordinate = std::numeric_limits<Precision>::lowest();
    Index surface_node_count = 0;
    std::vector<bool> used_nodes(static_cast<std::size_t>(_data->max_nodes), false);

    const auto surface_region = _data->surface_sets.get(surface_set);
    for (const ID surface_id : *surface_region) {
        if (surface_id < 0 || surface_id >= static_cast<ID>(_data->surfaces.size())) {
            continue;
        }
        const auto& surface = _data->surfaces[static_cast<std::size_t>(surface_id)];
        if (!surface) {
            continue;
        }
        for (Index i = 0; i < surface->n_nodes; ++i) {
            const ID node_id = surface->nodes()[i];
            if (node_id < 0 || node_id >= static_cast<ID>(_data->max_nodes) ||
                used_nodes[static_cast<std::size_t>(node_id)]) {
                continue;
            }
            used_nodes[static_cast<std::size_t>(node_id)] = true;
            const Index node = static_cast<Index>(node_id);
            surface_center(0) += (*_data->positions)(node, 0);
            surface_center(1) += (*_data->positions)(node, 1);
            surface_center(2) += (*_data->positions)(node, 2);
            const Vec3 node_position = _data->positions->row_vec3(node);
            const Precision axis_coordinate =
                section->axis_direction.dot(node_position);
            minimum_axis_coordinate = std::min(
                minimum_axis_coordinate, axis_coordinate);
            maximum_axis_coordinate = std::max(
                maximum_axis_coordinate, axis_coordinate);
            ++surface_node_count;
        }
    }

    logging::error(surface_node_count > 0,
                   "PRETENSION SECTION: surface set ", surface_set,
                   " has no valid surface nodes");
    surface_center /= static_cast<Precision>(surface_node_count);
    const Precision middle_axis_coordinate =
        (minimum_axis_coordinate + maximum_axis_coordinate) * Precision(0.5);
    section->axis_origin = surface_center +
        (middle_axis_coordinate - section->axis_direction.dot(surface_center)) *
            section->axis_direction;
    section->cut_coordinate = Precision(0);

    _data->pretension_sections.push_back(std::move(section));
}

void Model::add_pretension_interface_section(
    const std::string& name,
    const std::string& surface_set_a,
    const std::string& surface_set_b) {
    logging::error(!name.empty(),
                   "PRETENSION SECTION: name must not be empty");
    logging::error(surface_set_a != surface_set_b,
                   "PRETENSION SECTION: SURFACE_A and SURFACE_B must differ");
    logging::error(_data->surface_sets.has(surface_set_a),
                   "PRETENSION SECTION: surface set ", surface_set_a,
                   " does not exist");
    logging::error(_data->surface_sets.has(surface_set_b),
                   "PRETENSION SECTION: surface set ", surface_set_b,
                   " does not exist");
    const bool duplicate_name = std::any_of(
        _data->pretension_sections.begin(),
        _data->pretension_sections.end(),
        [&](const pretension::PretensionSection::Ptr& existing) {
            return existing && existing->name == name;
        });
    logging::error(!duplicate_name,
                   "PRETENSION SECTION: duplicate section name ", name);

    auto section = std::make_shared<pretension::PretensionSection>();
    section->name = name;
    section->interface_surface_set_a = surface_set_a;
    section->interface_surface_set_b = surface_set_b;
    _data->pretension_sections.push_back(std::move(section));
}

namespace {

std::set<ID> surface_set_nodes(
    ModelData& data,
    const std::string& set_name) {
    std::set<ID> nodes;
    const auto region = data.surface_sets.get(set_name);
    for (const ID surface_id : *region) {
        if (surface_id < 0 ||
            surface_id >= static_cast<ID>(data.surfaces.size())) continue;
        const SurfacePtr& surface = data.surfaces[static_cast<std::size_t>(surface_id)];
        if (!surface) continue;
        for (Index i = 0; i < surface->n_nodes; ++i) {
            nodes.insert(surface->nodes()[i]);
        }
    }
    return nodes;
}

std::set<ID> surface_owner_elements(
    ModelData& data,
    const std::string& set_name) {
    std::set<ID> owners;
    const auto region = data.surface_sets.get(set_name);
    for (const ID surface_id : *region) {
        if (surface_id < 0 || surface_id >=
            static_cast<ID>(data.surface_element_ids.size())) continue;
        const ID owner = data.surface_element_ids[static_cast<std::size_t>(surface_id)];
        if (owner >= 0) owners.insert(owner);
    }
    return owners;
}

void prepare_surface_pair(
    ModelData& data,
    pretension::PretensionSection& section) {
    const std::set<ID> nodes_a = surface_set_nodes(
        data, section.interface_surface_set_a);
    const std::set<ID> nodes_b = surface_set_nodes(
        data, section.interface_surface_set_b);
    logging::error(!nodes_a.empty() && !nodes_b.empty(),
                   "PRETENSION SECTION '", section.name,
                   "': both interface surfaces must contain nodes");
    logging::error(nodes_a.size() == nodes_b.size(),
                   "PRETENSION SECTION '", section.name,
                   "': SURFACE_A and SURFACE_B have different node counts");

    Vec3 center_a = Vec3::Zero();
    Vec3 center_b = Vec3::Zero();
    Precision scale = Precision(0);
    for (const ID node : nodes_a) {
        const Vec3 p = data.positions->row_vec3(static_cast<Index>(node));
        center_a += p;
        scale = std::max(scale, p.norm());
    }
    for (const ID node : nodes_b) {
        const Vec3 p = data.positions->row_vec3(static_cast<Index>(node));
        center_b += p;
        scale = std::max(scale, p.norm());
    }
    center_a /= static_cast<Precision>(nodes_a.size());
    center_b /= static_cast<Precision>(nodes_b.size());
    const Precision pair_tolerance = Precision(1e-8) * std::max(scale, Precision(1));

    std::map<ID, ID> pairs;
    std::set<ID> unused_b = nodes_b;
    for (const ID node_a : nodes_a) {
        const Vec3 position_a = data.positions->row_vec3(static_cast<Index>(node_a));
        ID best_node = -1;
        Precision best_distance = std::numeric_limits<Precision>::max();
        for (const ID node_b : unused_b) {
            const Precision distance = (position_a - data.positions->row_vec3(
                static_cast<Index>(node_b))).norm();
            if (distance < best_distance) {
                best_distance = distance;
                best_node = node_b;
            }
        }
        logging::error(best_node >= 0 && best_distance <= pair_tolerance,
                       "PRETENSION SECTION '", section.name,
                       "': circle-face nodes cannot be paired within tolerance ",
                       pair_tolerance);
        pairs.emplace(node_a, best_node);
        unused_b.erase(best_node);
    }
    logging::error(unused_b.empty(),
                   "PRETENSION SECTION '", section.name,
                   "': not all SURFACE_B nodes could be paired");

    std::set<ID> shared_nodes;
    for (const auto& [node_a, node_b] : pairs) {
        if (node_a == node_b) shared_nodes.insert(node_a);
    }

    std::map<ID, ID> duplicated_b;
    if (!shared_nodes.empty()) {
        const std::set<ID> seeds = surface_owner_elements(
            data, section.interface_surface_set_b);
        logging::error(!seeds.empty(),
                       "PRETENSION SECTION '", section.name,
                       "': SURFACE_B has no owning solid elements");
        std::map<ID, std::vector<ID>> node_elements;
        for (const ElementPtr& element : data.elements) {
            if (!element) continue;
            for (const ID node : *element) {
                if (shared_nodes.count(node) == 0) {
                    node_elements[node].push_back(element->elem_id);
                }
            }
        }
        std::set<ID> side_b_elements = seeds;
        std::deque<ID> queue(seeds.begin(), seeds.end());
        while (!queue.empty()) {
            const ID element_id = queue.front();
            queue.pop_front();
            const ElementPtr& element = data.elements[static_cast<std::size_t>(element_id)];
            for (const ID node : *element) {
                if (shared_nodes.count(node) > 0) continue;
                for (const ID neighbor : node_elements[node]) {
                    if (side_b_elements.insert(neighbor).second) queue.push_back(neighbor);
                }
            }
        }
        for (const ID node : shared_nodes) {
            duplicated_b[node] = data.append_node(
                data.positions->row_vec3(static_cast<Index>(node)));
        }
        for (const ID element_id : side_b_elements) {
            ElementPtr& element = data.elements[static_cast<std::size_t>(element_id)];
            for (const auto& [old_node, new_node] : duplicated_b) {
                element->replace_node(old_node, new_node);
            }
        }
        const auto region_b = data.surface_sets.get(section.interface_surface_set_b);
        for (const ID surface_id : *region_b) {
            SurfacePtr& surface = data.surfaces[static_cast<std::size_t>(surface_id)];
            if (!surface) continue;
            for (Index i = 0; i < surface->n_nodes; ++i) {
                const auto found = duplicated_b.find(surface->nodes()[i]);
                if (found != duplicated_b.end()) surface->nodes()[i] = found->second;
            }
        }
    }

    section.interface_pairs.clear();
    section.side_a_nodes.clear();
    section.side_b_nodes.clear();
    for (const auto& [node_a, original_b] : pairs) {
        const auto duplicated = duplicated_b.find(original_b);
        const ID node_b = duplicated == duplicated_b.end()
            ? original_b : duplicated->second;
        section.interface_pairs.push_back({node_a, node_b});
        section.side_a_nodes.push_back(node_a);
        section.side_b_nodes.push_back(node_b);
    }
    Vec3 axis = center_b - center_a;
    if (axis.norm() <= pair_tolerance) {
        const auto region_a = data.surface_sets.get(section.interface_surface_set_a);
        SurfacePtr reference_surface;
        for (const ID surface_id : *region_a) {
            if (surface_id < 0 ||
                surface_id >= static_cast<ID>(data.surfaces.size())) continue;
            reference_surface = data.surfaces[static_cast<std::size_t>(surface_id)];
            if (reference_surface) break;
        }
        logging::error(reference_surface != nullptr,
                       "PRETENSION SECTION '", section.name,
                       "': SURFACE_A contains no valid surface");
        axis = reference_surface->normal(*data.positions, Vec2::Zero());
    }
    logging::error(axis.norm() > Precision(0),
                   "PRETENSION SECTION '", section.name,
                   "': cannot determine interface axis");
    section.axis_direction = axis.normalized();
    section.axis_origin = (center_a + center_b) * Precision(0.5);
    logging::info(true,
                  "PRETENSIONSECTION '", section.name,
                  "': paired ", section.interface_pairs.size(),
                  " circle-face node(s), split ", shared_nodes.size(),
                  " shared node(s)");
}

} // namespace

void Model::prepare_pretension_sections() {
    for (auto& section : _data->pretension_sections) {
        logging::error(section != nullptr,
                       "PRETENSION SECTION: null section encountered");

        if (section->prepared) {
            continue;
        }

        if (section->uses_surface_pair()) {
            prepare_surface_pair(*_data, *section);
        } else {
            CutBuilder::split(*_data, *section);
        }
        section->prepared = true;
    }
}

void Model::set_pretension_load(
    const std::string& name,
    pretension::Control control,
    Precision value) {
    for (auto& section : _data->pretension_sections) {
        if (section && section->name == name) {
            section->control = control;
            section->prescribed_value = value;
            section->state = pretension::State::Loading;
            if (control == pretension::Control::Force) {
                section->has_solved_gap = false;
            }
            return;
        }
    }

    logging::error(false,
                   "PRETENSION: section ", name, " does not exist");
}

void Model::lock_pretension_section(const std::string& name) {
    for (auto& section : _data->pretension_sections) {
        if (section && section->name == name) {
            logging::error(section->state != pretension::State::Open,
                           "PRETENSION LOCK: section ", name,
                           " has not been loaded yet");
            if (section->control == pretension::Control::Force) {
                logging::error(section->has_solved_gap,
                               "PRETENSION LOCK: force-controlled section ", name,
                               " has no solved gap to lock");
                section->locked_gap = section->last_solved_gap;
                section->control = pretension::Control::Displacement;
            } else if (section->state == pretension::State::Loading) {
                section->locked_gap = section->prescribed_value;
            }
            section->state = pretension::State::Locked;
            logging::info(true, "PRETENSION '", name,
                          "': locked gap = ", section->locked_gap);
            return;
        }
    }

    logging::error(false,
                   "PRETENSION LOCK: section ", name, " does not exist");
}

void Model::add_rbm(const std::string& set) {
    logging::error(_data->elem_sets.has(set), "RBM: element set ", set, " not found");
    logging::error(_data->elem_sets.get(set) && _data->elem_sets.get(set)->size() > 0,
                   "RBM: element set ", set, " is empty");
    _data->rbms.emplace_back(_data->elem_sets.get(set));
}

void Model::add_cload(const std::string& nset, Vec6 load, const std::string& orientation, const std::string& amplitude) {
    logging::error(_data->node_sets.has(nset), "Node set ", nset, " does not exist");
    auto region_ptr = _data->node_sets.get(nset);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude), "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_cload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

void Model::add_cload(const ID id, Vec6 load, const std::string& orientation, const std::string& amplitude) {
    auto region_ptr = std::make_shared<NodeRegion>("INTERNAL");
    region_ptr->add(id);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude), "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_cload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

void Model::add_dload(const std::string& sfset, Vec3 load, const std::string& orientation, const std::string& amplitude) {
    logging::error(_data->surface_sets.has(sfset), "Surface set ", sfset, " does not exist");
    auto region_ptr = _data->surface_sets.get(sfset);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude), "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_dload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

void Model::add_dload(ID id, Vec3 load, const std::string& orientation, const std::string& amplitude) {
    auto region_ptr = std::make_shared<SurfaceRegion>("INTERNAL");
    region_ptr->add(id);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude), "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_dload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

void Model::add_pload(const std::string& sfset, Precision load, const std::string& amplitude) {
    logging::error(_data->surface_sets.has(sfset), "Surface set ", sfset, " does not exist");
    auto region_ptr = _data->surface_sets.get(sfset);
    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude), "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }
    _data->load_cols.get()->add_pload(region_ptr, load, amplitude_ptr);
}

void Model::add_pload(ID id, Precision load, const std::string& amplitude) {
    auto region_ptr = std::make_shared<SurfaceRegion>("INTERNAL");
    region_ptr->add(id);
    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude), "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }
    _data->load_cols.get()->add_pload(region_ptr, load, amplitude_ptr);
}

void Model::add_vload(const std::string& elset, Vec3 load, const std::string& orientation, const std::string& amplitude) {
    logging::error(_data->elem_sets.has(elset), "Element set ", elset, " does not exist");
    auto region_ptr = _data->elem_sets.get(elset);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude), "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_vload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

void Model::add_vload(const ID id, Vec3 load, const std::string& orientation, const std::string& amplitude) {
    auto region_ptr = std::make_shared<ElementRegion>("INTERNAL");
    region_ptr->add(id);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude), "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_vload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

void Model::add_inertialload(const std::string& elset,
                             Vec3 center,
                             Vec3 center_acceleration,
                             Vec3 angular_velocity,
                             Vec3 angular_acceleration,
                             bool consider_point_masses) {
    logging::error(_data->elem_sets.has(elset), "Element set ", elset, " does not exist");
    auto region_ptr = _data->elem_sets.get(elset);
    _data->load_cols.get()->add_inertialload(region_ptr,
                                             center,
                                             center_acceleration,
                                             angular_velocity,
                                             angular_acceleration,
                                             consider_point_masses);
}

void Model::define_amplitude(const std::string& name, bc::Interpolation interpolation) {
    auto amplitude = _data->amplitudes.activate(name);
    amplitude->set_interpolation(interpolation);
    amplitude->clear_samples();
}

void Model::add_amplitude_sample(const std::string& name, Precision time, Precision value) {
    logging::error(_data->amplitudes.has(name), "Amplitude ", name, " has not been defined");
    _data->amplitudes.get(name)->add_sample(time, value);
}

void Model::add_tload(std::string& temp_field, Precision ref_temp) {
    auto field = _data->get_field(temp_field);
    logging::error(field != nullptr, "Temperature field ", temp_field, " does not exist");
    logging::error(field->domain == FieldDomain::NODE, "Temperature field ", temp_field, " must be a node field");
    logging::error(field->components == 1, "Temperature field ", temp_field, " must have 1 component");
    _data->load_cols.get()->add_tload(field, ref_temp);

    // // TODO
    // logging::error(_fields_temperature.has(temp_field), "Temperature field ", temp_field, " does not exist");
    //
    // auto temp_ptr = _fields_temperature.get(temp_field);
    // for (ElementPtr& elem : _data->elements) {
    //     if (elem == nullptr) continue;
    //     if (auto sel = elem->as<StructuralElement>())
    //         sel->apply_tload((*_data->load_cols.get()), *temp_ptr, ref_temp);
    // }
}

void Model::add_support(const std::string& nset, const StaticVector<6> constraint, const std::string& orientation) {
    logging::error(_data->supp_cols.has_any(), "No support collectors have been defined");
    logging::error(_data->supp_cols.get() != nullptr, "No support collectors is currently active");
    if (!orientation.empty())
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");

    bc::SupportCollector::Ptr supp_col = _data->supp_cols.get();
    supp_col->add_supp(_data->node_sets.get(nset), constraint, _data->coordinate_systems.get(orientation));
}

void Model::add_support(const ID id, const StaticVector<6> constraint, const std::string& orientation) {
    logging::error(_data->supp_cols.has_any(), "No support collectors have been defined");
    logging::error(_data->supp_cols.get() != nullptr, "No support collectors is currently active");
    if (!orientation.empty())
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");

    // create a new NodeRegion
    NodeRegion::Ptr region = std::make_shared<NodeRegion>("INTERNAL");
    region->add(id);

    bc::SupportCollector::Ptr supp_col = _data->supp_cols.get();
    supp_col->add_supp(region, constraint, _data->coordinate_systems.get(orientation));
}

void Model::solid_section(const std::string& set, const std::string& material, const std::string& orientation) {
    logging::error(_data->elem_sets.has(set), "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material), "Material ", material, " is not a defined material");
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");
    }

    SolidSection::Ptr section = std::make_shared<SolidSection>();
    section->material_    = _data->materials.get(material);
    section->region_      = _data->elem_sets.get(set);
    section->orientation_ = orientation.empty() ? nullptr : _data->coordinate_systems.get(orientation);
    this->_data->sections.push_back(section);
}

void Model::beam_section(const std::string& set, const std::string& material,  const std::string& profile, Vec3 orientation) {
    logging::error(_data->elem_sets.has(set), "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material), "Material ", material, " is not a defined material");
    logging::error(_data->profiles.has(profile), "Profile ", profile, " is not a defined profile");
    BeamSection::Ptr sec = std::make_shared<BeamSection>();
    sec->material_  = _data->materials.get(material);
    sec->region_    = _data->elem_sets.get(set);
    sec->profile_   = _data->profiles.get(profile);
    sec->direction_ = orientation;
    this->_data->sections.push_back(sec);
}

void Model::truss_section(const std::string& set, const std::string& material, Precision area) {
    logging::error(_data->elem_sets.has(set), "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material), "Material ", material, " is not a defined material");

    TrussSection::Ptr sec = std::make_shared<TrussSection>(
        _data->materials.get(material),
        _data->elem_sets.get(set),
        area
    );
    this->_data->sections.push_back(sec);
}

void Model::shell_section(const std::string& set,
                          const std::string& material,
                          Precision          thickness,
                          const std::string& orientation,
                          Index              csys_axis) {
    logging::error(_data->elem_sets.has(set), "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material), "Material ", material, " is not a defined material");
    logging::error(csys_axis >= 0 && csys_axis < 3, "SHELLSECTION: CSYSAXIS must be 1, 2 or 3");
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");
    }

    IntegratedShellSection::Ptr sec = std::make_shared<IntegratedShellSection>(
        _data->materials.get(material),
        _data->elem_sets.get(set),
        thickness,
        orientation.empty() ? nullptr : _data->coordinate_systems.get(orientation),
        csys_axis
    );
    this->_data->sections.push_back(sec);
}

void Model::shell_section_abd(const std::string& set,
                              const std::string& material,
                              Precision          thickness,
                              const Mat6&        abd,
                              const Mat2&        shear,
                              const std::string& orientation,
                              Index              csys_axis) {
    logging::error(_data->elem_sets.has(set), "Element set ", set, " is not a defined element set");
    logging::error(csys_axis >= 0 && csys_axis < 3, "SHELLSECTION: CSYSAXIS must be 1, 2 or 3");
    if (!material.empty()) {
        logging::error(_data->materials.has(material), "Material ", material, " is not a defined material");
    }
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation), "Coordinate system ", orientation, " does not exist");
    }

    ABDShellSection::Ptr sec = std::make_shared<ABDShellSection>(
        material.empty() ? nullptr : _data->materials.get(material),
        _data->elem_sets.get(set),
        thickness,
        abd,
        shear,
        orientation.empty() ? nullptr : _data->coordinate_systems.get(orientation),
        csys_axis
    );
    this->_data->sections.push_back(sec);
}

void Model::add_point_mass_feature(const std::string& nset,
                                   Precision mass,
                                   Vec3 rotary_inertia,
                                   Vec3 spring,
                                   Vec3 rotary_spring) {
    logging::error(_data->node_sets.has(nset), "Node set ", nset, " is not a defined node set");
    auto feat = std::make_shared<feature::PointMass>();
    feat->region_                  = _data->node_sets.get(nset);
    feat->mass_                    = mass;
    feat->rotary_inertia_          = rotary_inertia;
    feat->spring_constants_        = spring;
    feat->rotary_spring_constants_ = rotary_spring;
    this->_data->features.push_back(std::move(feat));
}

std::ostream& operator<<(std::ostream& ostream, const model::Model& model) {
    ostream << "max nodes = " << model._data->max_nodes << '\n';
    ostream << "max elements = " << model._data->max_elems << '\n';
    ostream << "max surfaces = " << model._data->max_surfaces << "\n";

    logging::info(true, "Materials");
    logging::up();
    for (const auto& material : model._data->materials) {
        material.second->info();
    }
    logging::down();

    logging::info(true, "Sections");
    logging::up();
    for (const auto& section : model._data->sections) {
        section->info();
    }
    logging::down();

    logging::info(true, "Profiles");
    logging::up();
    for(const auto &profile : model._data->profiles) {
        profile.second->info();
    }
    logging::down();

    logging::info(true, "Element sets");
    logging::up();
    for (const auto& elem_set : model._data->elem_sets) {
        elem_set.second->info();
    }
    logging::down();

    logging::info(true, "Node sets");
    logging::up();
    for (const auto& node_set : model._data->node_sets) {
        node_set.second->info();
    }
    logging::down();

    return ostream;

    // // print materials
    // ostream << "Materials:\n";
    // for (const auto& material : model._data->materials) {
    //     ostream << *material.second;
    // }
    // return ostream;
}
} // namespace model
} // namespace fem
