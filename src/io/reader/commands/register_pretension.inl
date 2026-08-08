// register_pretension.inl - registers *PRETENSIONSECTION

#include <array>
#include <functional>
#include <memory>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

using PretensionSectionCountSink = std::function<void()>;

inline void register_pretension_section(
    fem::io::dsl::Registry& registry,
    model::Model& model) {
    registry.command("PRETENSIONSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Define a pretension section by a cylindrical surface, axis and cut position.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME")
                    .required()
                    .doc("Pretension section name")
                .key("SURFACE")
                    .required()
                    .doc("Selected cylindrical surface set")
                .key("POSITION")
                    .optional("MIDDLE")
                    .allowed({"MIDDLE"})
                    .doc("Cut position along the cylinder axis"));

        auto name = std::make_shared<std::string>();
        auto surface = std::make_shared<std::string>();
        auto position = std::make_shared<std::string>();

        command.on_enter([name, surface, position](const fem::io::dsl::Keys& keys) {
            *name = keys.raw("NAME");
            *surface = keys.raw("SURFACE");
            *position = keys.raw("POSITION");
        });

        command.variant(
            fem::io::dsl::Variant::make()
                .segment(fem::io::dsl::Segment::make()
                    .range(fem::io::dsl::LineRange{}.min(1).max(1))
                    .pattern(fem::io::dsl::Pattern::make()
                        .fixed<fem::Precision, 3>()
                        .name("AXIS")
                        .desc("Cylinder-axis direction"))
                    .bind([&model, name, surface, position](
                              const std::array<fem::Precision, 3>& values) {
                        model.add_pretension_section(
                            *name,
                            *surface,
                            fem::Vec3{values[0], values[1], values[2]},
                            *position);
                    })));
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
                .key("SURFACE").required()
                .key("POSITION").optional("MIDDLE"));

        command.variant(
            fem::io::dsl::Variant::make()
                .segment(fem::io::dsl::Segment::make()
                    .range(fem::io::dsl::LineRange{}.min(1).max(1))
                    .pattern(fem::io::dsl::Pattern::make()
                        .fixed<fem::Precision, 3>()
                        .name("AXIS"))
                    .bind([sink](const std::array<fem::Precision, 3>&) {
                        sink();
                    })));
    });
}

inline void register_pretension(
    fem::io::dsl::Registry& registry,
    model::Model& model) {
    registry.command("PRETENSION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
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

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const std::string section = keys.raw("SECTION");
            const std::string action = keys.raw("ACTION");

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
