// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef REPLACEMENT_H_
#define REPLACEMENT_H_

#include "InputOutput.h"

namespace combi_ff {

namespace topology_builder {

void ReplaceMacrosByParameters(
    const std::string& filename_molecules_with_macros,
    const std::string& family_code, const IOFileProperties& io_file_properties);
void CreateMTBWithParameters(
    const std::string& filename_molecules_with_macros,
    const std::string& family_code, const IOFileProperties& io_file_properties,
    const std::vector<std::pair<std::string, std::string>>& replacements);

}  // namespace topology_builder

}  // namespace combi_ff

#endif