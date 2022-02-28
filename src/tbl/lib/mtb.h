// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef MTB_H_
#define MTB_H_

#include "TopologicalProperty.h"
#include "XmlParser.h"

namespace combi_ff {

namespace topology_builder {

class IOFileProperties;

static const std::vector<std::string> parameter_count_comments({"NB", "NBA",
                                                                "NIDA", "NDA"});
static const std::vector<std::string> parameter_comments({"IB", "JB", "KB",
                                                          "LB"});
static const std::vector<PropertyType> property_types_in_mtb(
    {topologicalPropertyMap.find("bond")->second,
     topologicalPropertyMap.find("angle")->second,
     topologicalPropertyMap.find("improper_dihedral")->second,
     topologicalPropertyMap.find("torsional_dihedral")->second});

void CreateMTB(const std::string& filename_molecules_with_macros,
               const std::string& family_code,
               const IOFileProperties& io_file_properties);

}  // namespace topology_builder

}  // namespace combi_ff
#endif