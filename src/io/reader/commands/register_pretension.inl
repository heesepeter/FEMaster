// register_pretension.inl - registers *PRETENSIONSECTION

#include <functional>
#include <memory>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader { class Parser; }

namespace fem::io::reader::commands {

using PretensionSectionCountSink = std::function<void()>;

inline void register_pretension_section(
    fem::io::dsl::Registry& registry,
    model::Model& model) {
    registry.command("PRETENSIONSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Define a pretension section using one merged interface face "
            "or two existing face sets.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME")
                    .required()
                    .doc("Pretension section name")
                .key("SURFACE_A")
                    .required()
                    .doc("First existing interface face set")
                .key("SURFACE_B")
                    .optional("")
                    .doc("Second interface face set; omit for a merged interface"));

        auto name = std::make_shared<std::string>();
        auto surface_a = std::make_shared<std::string>();
        auto surface_b = std::make_shared<std::string>();

        command.on_enter([name, surface_a, surface_b](const fem::io::dsl::Keys& keys) {
            *name = keys.raw("NAME");
            *surface_a = keys.raw("SURFACE_A");
            *surface_b = keys.raw("SURFACE_B");
            if (surface_b->empty()) *surface_b = *surface_a;
        });

        command.on_exit([&model, name, surface_a, surface_b](const fem::io::dsl::Keys&) {
            model.add_pretension_interface_section(*name, *surface_a, *surface_b);
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

inline void register_pretension_section_count(
    fem::io::dsl::Registry& registry,
    PretensionSectionCountSink sink) {
    registry.command("PRETENSIONSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").required()
                .key("SURFACE_A").required()
                .key("SURFACE_B").optional(""));

        command.on_exit([sink](const fem::io::dsl::Keys&) {
            sink();
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

inline void register_pretension(
    fem::io::dsl::Registry& registry,
    model::Model& model,
    fem::io::reader::Parser& parser) {
    registry.command("PRETENSION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "LOADCASE"}));
        command.doc("Apply or lock a pretension section.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("SECTION")
                    .required()
                    .doc("Pretension section name")
                .key("ACTION")
                    .optional("LOCK")
                    .allowed({"LOCK", "LOAD"})
                .key("CONTROL")
                    .optional("DISPLACEMENT")
                    .allowed({"FORCE", "DISPLACEMENT"})
                .key("VALUE")
                    .optional("0")
                    .doc("Force or relative displacement value"));

        command.on_enter([&model, &parser](const fem::io::dsl::Keys& keys) {
            const std::string section = keys.raw("SECTION");
            const std::string action = keys.raw("ACTION");

            if (parser.active_loadcase()) {
                parser.queue_pretension_action(
                    section, action, keys.raw("CONTROL"),
                    keys.get<fem::Precision>("VALUE"));
                return;
            }

            if (action == "LOCK") {
                model.lock_pretension_section(section);
                return;
            }

            const std::string control = keys.raw("CONTROL");
            const fem::Precision value = keys.get<fem::Precision>("VALUE");
            const auto mode = control == "FORCE"
                ? fem::pretension::Control::Force
                : fem::pretension::Control::Displacement;

            model.set_pretension_load(section, mode, value);
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
