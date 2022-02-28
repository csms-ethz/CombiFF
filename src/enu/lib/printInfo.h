// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef PRINTINFO_H
#define PRINTINFO_H
#include <iomanip>

#include "Range.h"

namespace combi_ff {

namespace enu {

class SubstructureCollection;

void PrintRangedProperty(const Range& r, const std::string whitespaces,
                         const std::string name);
void PrintRestrictions(const std::string whitespaces, const size_t& max_degree,
                       const RangeVector& ranged_properties);
void PrintSubstructures(
    const std::string& indent2, const std::string& indent3,
    const std::vector<SubstructureCollection>& substructures);
void PrintSubstructureRange(const Range& r);

#endif

}  // namespace enu

}  // namespace combi_ff