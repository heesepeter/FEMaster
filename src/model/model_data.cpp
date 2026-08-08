/**
 * @file model_data.cpp
 * @brief Provides the translation unit companion for `ModelData`.
 *
 * Currently this file only includes the header to satisfy build systems that
 * expect a `.cpp` next to the declaration.
 *
 * @see src/model/model_data.h
 */

#include "model_data.h"
#include "element/element.h"

namespace fem {
namespace model {

ID ModelData::append_node(const Vec3& position) {
    logging::error(positions != nullptr,
                   "ModelData: cannot append a node before POSITION has been initialized");
    logging::error(positions_reference != nullptr,
                   "ModelData: cannot append a node before POSITION_REFERENCE has been initialized");
    logging::error(node_sets.has_all() && node_sets.all() != nullptr,
                   "ModelData: global node set is not initialized");

    const auto all_nodes = node_sets.all();
    const ID node_id = all_nodes->size() == 0 ? ID(0) : all_nodes->last() + ID(1);

    logging::error(node_id >= ID(0) && node_id < max_nodes,
                   "ModelData: node capacity exhausted while appending node ", node_id,
                   " (capacity: ", max_nodes, ")");

    const Index row = static_cast<Index>(node_id);
    for (Index component = 0; component < 3; ++component) {
        (*positions)(row, component) = position(component);
        (*positions_reference)(row, component) = position(component);
    }

    // Appended interface nodes belong to the global node set only. Do not use
    // Sets::add() here because an input NSET may still be active while the
    // topology is being split.
    if (node_sets.has_all() && node_sets.all() != nullptr) {
        node_sets.all()->add(node_id);
    }
    return node_id;
}

ID ModelData::next_free_element_id() const {
    for (ID element_id = 0; element_id < max_elems; ++element_id) {
        if (elements[static_cast<std::size_t>(element_id)] == nullptr) {
            return element_id;
        }
    }

    logging::error(false,
                   "ModelData: element capacity exhausted while inserting a new element");
    return ID(-1);
}

void ModelData::insert_element(ElementPtr element) {
    logging::error(element != nullptr,
                   "ModelData: cannot insert a null element");

    const ID element_id = element->elem_id;
    logging::error(element_id >= 0 && element_id < max_elems,
                   "ModelData: element id ", element_id,
                   " is outside the allocated element capacity ", max_elems);
    logging::error(elements[static_cast<std::size_t>(element_id)] == nullptr,
                   "ModelData: element slot ", element_id,
                   " is already occupied");

    element->_model_data = this;
    elements[static_cast<std::size_t>(element_id)] = std::move(element);

    if (elem_sets.has_all() && elem_sets.all() != nullptr) {
        elem_sets.all()->add(element_id);
    }
}

Index ModelData::field_rows(FieldDomain domain) {
    switch (domain) {
        case FieldDomain::UNKNOWN:
            logging::error(false, "ModelData: cannot allocate UNKNOWN fields");
            return 0;
        case FieldDomain::NODE:
            return static_cast<Index>(max_nodes);
        case FieldDomain::ELEMENT:
            return static_cast<Index>(max_elems);
        case FieldDomain::ELEMENT_NODAL:
            logging::error(element_nodal_offsets != nullptr,
                           "ModelData: element nodal offset field is not initialized");
            logging::error((*element_nodal_offsets)(static_cast<Index>(max_elems)) > 0,
                           "ModelData: no element nodes are available, cannot allocate ELEMENT_NODAL fields");
            return static_cast<Index>((*element_nodal_offsets)(static_cast<Index>(max_elems)));
        case FieldDomain::ELEMENT_IP:
            logging::error(element_ip_offsets != nullptr,
                           "ModelData: element IP offset field is not initialized");
            logging::error((*element_ip_offsets)(static_cast<Index>(max_elems)) > 0,
                           "ModelData: no integration points are available, cannot allocate ELEMENT_IP fields");
            return static_cast<Index>((*element_ip_offsets)(static_cast<Index>(max_elems)));
        case FieldDomain::ELEMENT_MP:
            logging::error(element_mp_offsets != nullptr,
                           "ModelData: element MP offset field is not initialized");
            logging::error((*element_mp_offsets)(static_cast<Index>(max_elems)) > 0,
                           "ModelData: no material points are available, cannot allocate ELEMENT_MP fields");
            return static_cast<Index>((*element_mp_offsets)(static_cast<Index>(max_elems)));
    }
    logging::error(false, "ModelData: unknown field domain");
    return 0;
}

void ModelData::initialize_element_enumeration() {
    logging::error(element_nodal_offsets == nullptr &&
                   element_ip_offsets    == nullptr &&
                   element_mp_offsets    == nullptr,
                   "ModelData: element enumeration has already been initialized");

    element_nodal_offsets = std::make_shared<Field>(
        "ELEMENT_NODAL_OFFSETS", FieldDomain::ELEMENT, static_cast<Index>(max_elems + 1), 1);
    element_ip_offsets = std::make_shared<Field>(
        "ELEMENT_IP_OFFSETS", FieldDomain::ELEMENT, static_cast<Index>(max_elems + 1), 1);
    element_mp_offsets = std::make_shared<Field>(
        "ELEMENT_MP_OFFSETS", FieldDomain::ELEMENT, static_cast<Index>(max_elems + 1), 1);

    Index nodal_offset = 0;
    Index ip_offset    = 0;
    Index mp_offset    = 0;

    for (Index row = 0; row < static_cast<Index>(max_elems); ++row) {
        (*element_nodal_offsets)(row) = static_cast<Precision>(nodal_offset);
        (*element_ip_offsets)(row)    = static_cast<Precision>(ip_offset);
        (*element_mp_offsets)(row)    = static_cast<Precision>(mp_offset);

        const ElementPtr& element = elements[row];
        if (element != nullptr) {
            element->elem_nodal_offset = static_cast<ID>(nodal_offset);
            element->elem_ip_offset    = static_cast<ID>(ip_offset);
            element->elem_mp_offset    = static_cast<ID>(mp_offset);

            nodal_offset += static_cast<Index>(element->n_nodes());
            ip_offset    += static_cast<Index>(element->num_ip());
            mp_offset    += static_cast<Index>(element->num_ip()) * element->num_mp_per_ip();
        }
    }

    const Index sentinel = static_cast<Index>(max_elems);
    (*element_nodal_offsets)(sentinel) = static_cast<Precision>(nodal_offset);
    (*element_ip_offsets)(sentinel)    = static_cast<Precision>(ip_offset);
    (*element_mp_offsets)(sentinel)    = static_cast<Precision>(mp_offset);
    max_integration_points = static_cast<ID>(ip_offset);
    max_material_points    = static_cast<ID>(mp_offset);
}

bool ModelData::has_field(const std::string& name) const {
    return fields.find(name) != fields.end();
}

Field::Ptr ModelData::get_field(const std::string& name) const {
    auto it = fields.find(name);
    if (it == fields.end()) {
        return nullptr;
    }
    return it->second;
}

Field::Ptr ModelData::create_field(const std::string& name, FieldDomain domain, Index components, bool fill_nan, bool reg) {
    logging::error(!name.empty(), "Field name cannot be empty");

    if (!reg) {
        const Index rows = field_rows(domain);
        auto field = std::make_shared<Field>(name, domain, rows, components);
        if (fill_nan) {
            field->fill_nan();
        }
        return field;
    }

    auto it = fields.find(name);
    if (it != fields.end()) {
        Field::Ptr field = it->second;
        logging::error(field->domain == domain, "Field '", name, "': domain mismatch");
        logging::error(field->components == components, "Field '", name, "': components mismatch");
        return field;
    }

    const Index rows = field_rows(domain);
    auto field = std::make_shared<Field>(name, domain, rows, components);
    if (fill_nan) {
        field->fill_nan();
    }
    fields.emplace(name, field);
    return field;
}

Field ModelData::create_field_(const std::string& name, FieldDomain domain, Index components, bool fill_nan) {
    logging::error(!name.empty(), "Field name cannot be empty");

    const Index rows = field_rows(domain);
    auto field = Field(name, domain, rows, components);
    if (fill_nan) {
        field.fill_nan();
    }
    return field;
}

Field ModelData::element_nodal_to_nodal(const Field& element_nodal,
                                        const Field& element_weights,
                                        const std::string& name) const {
    logging::error(element_nodal.domain == FieldDomain::ELEMENT_NODAL,
                   "ModelData: element_nodal_to_nodal requires ELEMENT_NODAL input field '",
                   element_nodal.name, "'");
    logging::error(element_weights.domain == FieldDomain::ELEMENT,
                   "ModelData: element_nodal_to_nodal requires ELEMENT weight field '",
                   element_weights.name, "'");
    logging::error(element_weights.components == 1,
                   "ModelData: element weight field '", element_weights.name,
                   "' must have exactly one component");
    logging::error(element_weights.rows == static_cast<Index>(max_elems),
                   "ModelData: element weight field '", element_weights.name,
                   "' has ", element_weights.rows, " rows, expected ", max_elems);
    logging::error(element_nodal_offsets != nullptr,
                   "ModelData: element nodal offset field is not initialized");
    logging::error(!name.empty(), "ModelData: projected nodal field name cannot be empty");

    const Field& offsets = *element_nodal_offsets;
    const Index expected_rows = static_cast<Index>(offsets(static_cast<Index>(max_elems), 0));
    logging::error(element_nodal.rows == expected_rows,
                   "ModelData: ELEMENT_NODAL field '", element_nodal.name,
                   "' has ", element_nodal.rows, " rows, expected ", expected_rows);

    Field nodal{name, FieldDomain::NODE, static_cast<Index>(max_nodes), element_nodal.components};
    nodal.set_zero();
    std::vector<Precision> weight_sum(static_cast<std::size_t>(max_nodes), Precision(0));

    for (Index elem_idx = 0; elem_idx < static_cast<Index>(max_elems); ++elem_idx) {
        const ElementPtr& element = elements[elem_idx];
        if (!element) {
            continue;
        }

        const ID elem_id = element->elem_id;
        logging::error(elem_id >= 0 && elem_id < max_elems,
                       "ModelData: element id out of range in element_nodal_to_nodal: ", elem_id);

        const Precision weight = element_weights(static_cast<Index>(elem_id), 0);
        if (weight == Precision(0)) {
            continue;
        }

        const Index offset = static_cast<Index>(offsets(static_cast<Index>(elem_id), 0));
        const Index next_offset = static_cast<Index>(offsets(static_cast<Index>(elem_id) + 1, 0));
        logging::error(next_offset >= offset,
                       "ModelData: invalid element nodal offsets for element ", elem_id);
        logging::error(next_offset - offset == static_cast<Index>(element->n_nodes()),
                       "ModelData: element nodal offset span does not match element node count for element ",
                       elem_id);

        for (Index local_node = 0; local_node < static_cast<Index>(element->n_nodes()); ++local_node) {
            const Index element_row = offset + local_node;
            const Index node_id = static_cast<Index>(element->nodes()[local_node]);
            logging::error(node_id < static_cast<Index>(max_nodes),
                           "ModelData: node id out of range in element_nodal_to_nodal: ", node_id);

            for (Index component = 0; component < element_nodal.components; ++component) {
                nodal(node_id, component) += weight * element_nodal(element_row, component);
            }
            weight_sum[static_cast<std::size_t>(node_id)] += weight;
        }
    }

    for (Index node = 0; node < static_cast<Index>(max_nodes); ++node) {
        const Precision weight = weight_sum[static_cast<std::size_t>(node)];
        if (weight == Precision(0)) {
            continue;
        }
        for (Index component = 0; component < element_nodal.components; ++component) {
            nodal(node, component) /= weight;
        }
    }

    return nodal;
}

} // namespace model
} // namespace fem
