#include "model.h"

#include "cut/cut_builder.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <set>

namespace fem::model {

void Model::add_pretension_section(const std::string& name,
                                   const std::string& surface_set,
                                   Vec3 axis,
                                   const std::string& position,
                                   Precision snap_ratio) {
    logging::error(!name.empty(), "PRETENSION SECTION: name must not be empty");
    const bool duplicate = std::any_of(
        _data->pretension_sections.begin(), _data->pretension_sections.end(),
        [&](const auto& section) { return section && section->name == name; });
    logging::error(!duplicate,
        "PRETENSION SECTION: duplicate section name ", name);
    logging::error(_data->surface_sets.has(surface_set),
        "PRETENSION SECTION: surface set ", surface_set, " does not exist");
    logging::error(axis.norm() > Precision(0),
        "PRETENSION SECTION: axis must not be zero");
    logging::error(position == "MIDDLE",
        "PRETENSION SECTION: only POSITION=MIDDLE is supported yet");
    logging::error(snap_ratio >= Precision(0) && snap_ratio <= Precision(0.1),
        "PRETENSION SECTION: SNAP must be between 0 and 0.1");

    auto section = std::make_shared<pretension::PretensionSection>();
    section->name = name;
    section->cylinder_surface_set = surface_set;
    section->axis_direction = axis.normalized();
    section->snap_ratio = snap_ratio;

    Vec3 center = Vec3::Zero();
    Precision minimum = std::numeric_limits<Precision>::max();
    Precision maximum = std::numeric_limits<Precision>::lowest();
    std::set<ID> used_nodes;
    for (const ID surface_id : *_data->surface_sets.get(surface_set)) {
        if (surface_id < 0 ||
            surface_id >= static_cast<ID>(_data->surfaces.size())) continue;
        const auto& surface = _data->surfaces[static_cast<std::size_t>(surface_id)];
        if (!surface) continue;
        for (Index i = 0; i < surface->n_nodes; ++i) {
            const ID node = surface->nodes()[i];
            if (!used_nodes.insert(node).second) continue;
            const Vec3 point = _data->positions->row_vec3(static_cast<Index>(node));
            center += point;
            const Precision coordinate = section->axis_direction.dot(point);
            minimum = std::min(minimum, coordinate);
            maximum = std::max(maximum, coordinate);
        }
    }
    logging::error(!used_nodes.empty(),
        "PRETENSION SECTION: surface set ", surface_set,
        " has no valid surface nodes");
    center /= static_cast<Precision>(used_nodes.size());
    const Precision middle = (minimum + maximum) * Precision(0.5);
    section->axis_origin = center +
        (middle - section->axis_direction.dot(center)) * section->axis_direction;
    section->cut_coordinate = Precision(0);
    _data->pretension_sections.push_back(std::move(section));
}

void Model::add_pretension_interface_section(const std::string& name,
                                             const std::string& surface_set_a,
                                             const std::string& surface_set_b) {
    logging::error(!name.empty(), "PRETENSION SECTION: name must not be empty");
    logging::error(surface_set_a != surface_set_b,
        "PRETENSION SECTION: SURFACE_A and SURFACE_B must differ");
    logging::error(_data->surface_sets.has(surface_set_a),
        "PRETENSION SECTION: surface set ", surface_set_a, " does not exist");
    logging::error(_data->surface_sets.has(surface_set_b),
        "PRETENSION SECTION: surface set ", surface_set_b, " does not exist");
    const bool duplicate = std::any_of(
        _data->pretension_sections.begin(), _data->pretension_sections.end(),
        [&](const auto& section) { return section && section->name == name; });
    logging::error(!duplicate,
        "PRETENSION SECTION: duplicate section name ", name);

    auto section = std::make_shared<pretension::PretensionSection>();
    section->name = name;
    section->interface_surface_set_a = surface_set_a;
    section->interface_surface_set_b = surface_set_b;
    _data->pretension_sections.push_back(std::move(section));
}

namespace {

std::set<ID> surface_nodes(ModelData& data, const std::string& set_name) {
    std::set<ID> nodes;
    for (const ID surface_id : *data.surface_sets.get(set_name)) {
        if (surface_id < 0 ||
            surface_id >= static_cast<ID>(data.surfaces.size())) continue;
        const auto& surface = data.surfaces[static_cast<std::size_t>(surface_id)];
        if (!surface) continue;
        for (Index i = 0; i < surface->n_nodes; ++i) {
            nodes.insert(surface->nodes()[i]);
        }
    }
    return nodes;
}

std::set<ID> surface_owners(ModelData& data, const std::string& set_name) {
    std::set<ID> owners;
    for (const ID surface_id : *data.surface_sets.get(set_name)) {
        if (surface_id < 0 ||
            surface_id >= static_cast<ID>(data.surface_element_ids.size())) continue;
        const ID owner = data.surface_element_ids[static_cast<std::size_t>(surface_id)];
        if (owner >= 0) owners.insert(owner);
    }
    return owners;
}

void prepare_surface_pair(ModelData& data,
                          pretension::PretensionSection& section) {
    const auto nodes_a = surface_nodes(data, section.interface_surface_set_a);
    const auto nodes_b = surface_nodes(data, section.interface_surface_set_b);
    logging::error(!nodes_a.empty() && !nodes_b.empty(),
        "PRETENSION SECTION '", section.name,
        "': both interface surfaces must contain nodes");
    logging::error(nodes_a.size() == nodes_b.size(),
        "PRETENSION SECTION '", section.name,
        "': SURFACE_A and SURFACE_B have different node counts");

    Vec3 center_a = Vec3::Zero();
    Vec3 center_b = Vec3::Zero();
    Precision scale = 0;
    for (const ID node : nodes_a) {
        const Vec3 point = data.positions->row_vec3(static_cast<Index>(node));
        center_a += point;
        scale = std::max(scale, point.norm());
    }
    for (const ID node : nodes_b) {
        const Vec3 point = data.positions->row_vec3(static_cast<Index>(node));
        center_b += point;
        scale = std::max(scale, point.norm());
    }
    center_a /= static_cast<Precision>(nodes_a.size());
    center_b /= static_cast<Precision>(nodes_b.size());
    const Precision tolerance = Precision(1e-8) * std::max(scale, Precision(1));

    std::map<ID, ID> pairs;
    std::set<ID> unused_b = nodes_b;
    for (const ID node_a : nodes_a) {
        const Vec3 point_a = data.positions->row_vec3(static_cast<Index>(node_a));
        ID best = -1;
        Precision distance_best = std::numeric_limits<Precision>::max();
        for (const ID node_b : unused_b) {
            const Precision distance = (point_a - data.positions->row_vec3(
                static_cast<Index>(node_b))).norm();
            if (distance < distance_best) {
                distance_best = distance;
                best = node_b;
            }
        }
        logging::error(best >= 0 && distance_best <= tolerance,
            "PRETENSION SECTION '", section.name,
            "': circle-face nodes cannot be paired within tolerance ", tolerance);
        pairs.emplace(node_a, best);
        unused_b.erase(best);
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
        const auto seeds = surface_owners(data, section.interface_surface_set_b);
        logging::error(!seeds.empty(),
            "PRETENSION SECTION '", section.name,
            "': SURFACE_B has no owning solid elements");

        std::map<ID, std::vector<ID>> node_elements;
        for (const auto& element : data.elements) {
            if (!element) continue;
            for (const ID node : *element) {
                if (!shared_nodes.count(node)) {
                    node_elements[node].push_back(element->elem_id);
                }
            }
        }
        std::set<ID> side_b_elements = seeds;
        std::deque<ID> queue(seeds.begin(), seeds.end());
        while (!queue.empty()) {
            const ID element_id = queue.front();
            queue.pop_front();
            const auto& element = data.elements[static_cast<std::size_t>(element_id)];
            for (const ID node : *element) {
                if (shared_nodes.count(node)) continue;
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
            auto& element = data.elements[static_cast<std::size_t>(element_id)];
            for (const auto& [old_node, new_node] : duplicated_b) {
                element->replace_node(old_node, new_node);
            }
        }
        for (const ID surface_id : *data.surface_sets.get(
                 section.interface_surface_set_b)) {
            auto& surface = data.surfaces[static_cast<std::size_t>(surface_id)];
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
    if (axis.norm() <= tolerance) {
        SurfacePtr reference;
        for (const ID surface_id : *data.surface_sets.get(
                 section.interface_surface_set_a)) {
            if (surface_id >= 0 &&
                surface_id < static_cast<ID>(data.surfaces.size())) {
                reference = data.surfaces[static_cast<std::size_t>(surface_id)];
                if (reference) break;
            }
        }
        logging::error(reference != nullptr,
            "PRETENSION SECTION '", section.name,
            "': SURFACE_A contains no valid surface");
        axis = reference->normal(*data.positions, Vec2::Zero());
    }
    logging::error(axis.norm() > Precision(0),
        "PRETENSION SECTION '", section.name,
        "': cannot determine interface axis");
    section.axis_direction = axis.normalized();
    section.axis_origin = (center_a + center_b) * Precision(0.5);
    logging::info(true, "PRETENSIONSECTION '", section.name,
        "': paired ", section.interface_pairs.size(),
        " circle-face node(s), split ", shared_nodes.size(), " shared node(s)");
}

} // namespace

void Model::prepare_pretension_sections() {
    for (auto& section : _data->pretension_sections) {
        logging::error(section != nullptr,
            "PRETENSION SECTION: null section encountered");
        if (section->prepared) continue;
        if (section->uses_surface_pair()) {
            prepare_surface_pair(*_data, *section);
        } else {
            CutBuilder::split(*_data, *section);
        }
        section->prepared = true;
    }
}

void Model::set_pretension_load(const std::string& name,
                                pretension::Control control,
                                Precision value) {
    for (auto& section : _data->pretension_sections) {
        if (!section || section->name != name) continue;
        section->control = control;
        section->prescribed_value = value;
        section->state = pretension::State::Loading;
        if (control == pretension::Control::Force) section->has_solved_gap = false;
        return;
    }
    logging::error(false, "PRETENSION: section ", name, " does not exist");
}

void Model::lock_pretension_section(const std::string& name) {
    for (auto& section : _data->pretension_sections) {
        if (!section || section->name != name) continue;
        logging::error(section->state != pretension::State::Open,
            "PRETENSION LOCK: section ", name, " has not been loaded yet");
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
    logging::error(false, "PRETENSION LOCK: section ", name, " does not exist");
}

} // namespace fem::model
