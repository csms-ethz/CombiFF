// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "TblAtom.h"

namespace combi_ff {

namespace topology_builder {

/********************
TblAtom implementation
********************/

/*
constructor
*/
TblAtom::TblAtom(const std::string& atom_id, const AtomTypeSet& atom_types,
                 const LinkageType link_type)
    : atom_id(atom_id), atom_types(atom_types), link_type(link_type) {}

/*
Setter functions
*/
void TblAtom::SetAtomID(const std::string& id) { atom_id = id; }

void TblAtom::SetAtomTypeSet(const AtomTypeSet& atom_type_set) {
  atom_types = atom_type_set;
}

void TblAtom::SetLinkageType(const LinkageType link_type_) {
  link_type = link_type_;
}

void TblAtom::AddNeighbor(const size_t nbr_idx) {
  neighbors.push_back(nbr_idx);
}

/*
Getter functions
*/
const AtomTypeSet& TblAtom::GetAtomTypes() const { return atom_types; }

const LinkageType TblAtom::GetLinkageType() const { return link_type; }

const std::string& TblAtom::GetAtomID() const { return atom_id; }

const NeighborVector& TblAtom::GetNeighbors() const { return neighbors; }

bool TblAtom::IsCoreAtom() const { return link_type == core; }

bool TblAtom::IsLinkAtom() const { return link_type == link; }

}  // namespace topology_builder

}  // namespace combi_ff