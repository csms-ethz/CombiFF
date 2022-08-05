// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef RANGE_H
#define RANGE_H

#include <cstdlib>
#include <unordered_map>
#include <utility>
#include <vector>
#include <string>

namespace combi_ff {

// add new restriction here (before num_ranges)
typedef enum {
  range_unsaturations,
  range_bonds,
  range_single_bonds,
  range_double_bonds,
  range_triple_bonds,
  range_quadruple_bonds,
  range_rings,
  num_ranges
} RangeRestriction;

typedef std::unordered_map<std::string, RangeRestriction> RangeMap;

static const RangeMap possible_ranges{
    {"unsaturations", range_unsaturations},
    {"total_bonds", range_bonds},
    {"single_bonds", range_single_bonds},
    {"double_bonds", range_double_bonds},
    {"triple_bonds", range_triple_bonds},
    {"quadruple_bonds", range_quadruple_bonds},
    {"cycles", range_rings}};

typedef std::pair<int, int> Range;
typedef std::vector<Range> RangeVector;

bool IsInRange(size_t i, Range r);
bool AreInRange(const std::vector<size_t>& ranges,
                const RangeVector& ranged_properties);

}  // namespace combi_ff

#endif
