// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef TB_ATOM_H_
#define TB_ATOM_H_

#include "AtomTypes.h"

namespace combi_ff {

namespace topology_builder {

typedef enum { core, link } LinkageType;

class TblAtom {
 public:
  TblAtom() = default;
  TblAtom(const std::string& atom_id, const AtomTypeSet& atomtypes,
          const LinkageType linktype);

  void SetAtomID(const std::string& id);
  void SetAtomTypeSet(const AtomTypeSet& atomTypeSet);
  void AddNeighbor(const size_t nbrIdx);
  const AtomTypeSet& GetAtomTypes() const;
  const LinkageType GetLinkageType() const;
  const std::string& GetAtomID() const;
  const NeighborVector& GetNeighbors() const;
  void SetLinkageType(const LinkageType linktype_);
  bool IsCoreAtom() const;
  bool IsLinkAtom() const;

 private:
  std::string atom_id{""};
  AtomTypeSet atom_types{};
  LinkageType link_type{core};
  const TblAtom* core_atom;
  NeighborVector neighbors;
};

}  // namespace topology_builder

}  // namespace combi_ff

#endif