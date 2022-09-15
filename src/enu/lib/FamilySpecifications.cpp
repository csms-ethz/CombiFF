// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "FamilySpecifications.h"

#include "InputOutput.h"

namespace combi_ff {

namespace enu {

FamilySpecifications::FamilySpecifications(InputOutput& IO)
    : abstr_substructures(AbstractSubstructureMap(0)),
      pseudoatoms(PseudoatomMap(0)),
      alias_map(AliasMap(0)),
      families(FamilyVector(0)) {
  if (IO.GetEnumSpec().used_families.size()) {
    try {
      CreateSubstructures(abstr_substructures,
                          IO.GetInputFileNamesAt(substructure_file));

    } catch (combi_ff::input_warning& w) {
      std::cout << "?Warning: " << w.what() << std::endl;
    }

    try {
      CreatePseudoatoms(pseudoatoms, IO.GetInputFileNamesAt(pseudoatom_file));

    } catch (combi_ff::input_warning& w) {
      std::cout << "?Warning: " << w.what() << std::endl;
    }

    try {
      CreateAliases(alias_map, IO.GetInputFileNamesAt(alias_file));

    } catch (combi_ff::input_warning& w) {
      std::cout << "?Warning: " << w.what() << std::endl;
    }

    try {
      CreateFamilies(families, IO.GetInputFileNamesAt(family_file),
                     IO.GetEnumSpec().used_families, abstr_substructures,
                     alias_map, pseudoatoms);

    } catch (combi_ff::input_warning& w) {
      std::cout << "?Warning: " << w.what() << std::endl;
    }
  }
}

const FamilyVector& FamilySpecifications::GetFamilies() const {
  return families;
}

const PseudoatomMap& FamilySpecifications::GetPseudoatoms() const {
  return pseudoatoms;
}

}  // namespace enu

}  // namespace combi_ff