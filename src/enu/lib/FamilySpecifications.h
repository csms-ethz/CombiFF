#ifndef FAMILYSPECIFICATIONS_H_
#define FAMILYSPECIFICATIONS_H_

#include "Family.h"

namespace combi_ff {

namespace enu {

class InputOutput;

class FamilySpecifications {

 public:

  FamilySpecifications(InputOutput& IO);
  const FamilyVector& GetFamilies() const;
  const PseudoatomMap& GetPseudoatoms() const ;

 private:

  AbstractSubstructureMap abstr_substructures;
  PseudoatomMap pseudoatoms;
  AliasMap alias_map;
  FamilyVector families;

};

} //namespace enu

} //namespace combi_ff

#endif