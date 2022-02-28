// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "printInfo.h"

#include <iostream>

#include "Substructure.h"

namespace combi_ff {

namespace enu {

void PrintRestrictions(const std::string whitespaces, const size_t& max_degree,
                       const RangeVector& ranged_properties) {
  std::cout << whitespaces << std::setw(18) << std::left << "max_bond_degree:";

  if (max_degree != 5)
    std::cout << max_degree << '\n';

  else
    std::cout << "unrestricted\n";

  PrintRangedProperty(ranged_properties[range_unsaturations], whitespaces,
                      "unsaturations:");
  PrintRangedProperty(ranged_properties[range_bonds], whitespaces,
                      "total_bonds:");
  PrintRangedProperty(ranged_properties[range_single_bonds], whitespaces,
                      "single_bonds:");
  PrintRangedProperty(ranged_properties[range_double_bonds], whitespaces,
                      "double_bonds:");
  PrintRangedProperty(ranged_properties[range_triple_bonds], whitespaces,
                      "triple_bonds:");
  PrintRangedProperty(ranged_properties[range_quadruple_bonds], whitespaces,
                      "quadruple_bonds");
  PrintRangedProperty(ranged_properties[range_rings], whitespaces, "cycles:");
  // std::cout << "******************************************************\n";
}

void PrintRangedProperty(const Range& r, const std::string whitespaces,
                         const std::string name) {
  std::cout << whitespaces << std::setw(18) << std::left << name;

  if (r != Range({0, -1})) {
    if (r.first != r.second) {
      if (r.second == -1)
        std::cout << ">= " << r.first << '\n';

      else
        std::cout << r.first << "-" << r.second << "\n";

    } else
      std::cout << r.first << '\n';
  }

  else
    std::cout << "unrestricted\n";
}

void PrintSubstructures(
    const std::string& indent2, const std::string& indent3,
    const std::vector<SubstructureCollection>& substructures) {
  for (size_t i = 0; i < substructures.size(); i++) {
    std::cout << indent2 << "-" << substructures[i].GetCode() << " with range ";
    PrintSubstructureRange(substructures[i].GetRange());
    std::cout << indent3 << "atom vecs are: ";

    for (auto&& M : substructures[i].GetSubstructureMatrices())
      std::cout << "(" << M.GetAtomVector() << ") ";

    std::cout << '\n';

    if (substructures[i].GetXOR())
      std::cout << indent3 << "(XOR: all AtomVectors have to be different\n";

    else if (substructures[i].GetAND())
      std::cout << indent3 << "(AND: all AtomVectors have to be the same\n";
  }
}

/**************************************
pint Information on the Fragment ranges
**************************************/
void PrintSubstructureRange(const Range& r) {
  if (r != Range({0, -1})) {
    if (r.first != r.second) {
      if (r.second == -1)
        std::cout << ">=" << std::to_string(r.first) << '\n';

      else
        std::cout << std::to_string(r.first) << "-" << std::to_string(r.second)
                  << '\n';

    } else
      std::cout << std::to_string(r.first) << '\n';
  }

  else
    std::cout << "unrestricted\n" << '\n';
}

}  // namespace enu

}  // namespace combi_ff