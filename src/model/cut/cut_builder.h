#pragma once

#include "../model_data.h"
#include "../../section/pretension_section.h"

namespace fem::model {

class CutBuilder {
public:
    static void split(
        ModelData& model_data,
        pretension::PretensionSection& section);
};

} // namespace fem::model