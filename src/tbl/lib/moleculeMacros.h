// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef MOLECULEMACROS_H_
#define MOLECULEMACROS_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "MatrixTbl.h"
#include "XmlParser.h"

namespace combi_ff {

namespace topology_builder {

class TblFragment;
class TopologicalPropertyBase;
class AtomTypeSet;
class IOFileProperties;

struct AtomInMolecule {
  size_t idx{};
  std::vector<size_t> neighbors{};
  std::string atom_id{};
  std::string parameter{};
  std::string atomtype{};
  const std::vector<std::string> connected_link_atoms{};

  AtomInMolecule() = default;
  AtomInMolecule(size_t idx, const std::string& par, const std::string& id,
                 const std::vector<std::string>& connected_link_atoms);
};

struct topo {
  std::shared_ptr<TopologicalPropertyBase> property{nullptr};
  TblFragment const* frag{NULL};
  std::string fragment_id{""};

  topo() = default;
};

bool operator<(const std::shared_ptr<TopologicalPropertyBase>& b1,
               const std::shared_ptr<TopologicalPropertyBase>& b2);
bool operator<(const std::pair<std::string, const AtomTypeSet*>& at1,
               const std::pair<std::string, const AtomTypeSet*>& at2);
void CreateMoleculesWithMacros(
    const std::string& filename_molecule_decomposition,
    std::string& filename_molecules_with_macros, const std::string& family_code,
    const std::vector<TblFragment>& tbl_fragments,
    const IOFileProperties& io_file_properties);
void CreateMoleculeWithMacros(XmlElement_ptr molecule_decomposition,
                              XmlElement& molecule_macros,
                              const std::vector<TblFragment>& tbl_fragments,
                              const bool united_atom,
                              const bool third_neighbor_exclusions,
                              const bool unique_torsionals);
const std::string GetParameterNameFromFragment(
    const TblFragment& frag, const std::string& fragment_id,
    const TopologicalPropertyBase& prop,
    const std::unordered_map<std::string, std::pair<std::string, std::string>>&
        core_atom_types,
    const std::vector<std::string>& involved_atoms_codes,
    const std::unordered_map<std::string, std::vector<std::string>>&
        connectedFragments);
const std::string GetParameterNameAutomatically(
    const topology_builder::AdjacencyMatrix& A,
    std::shared_ptr<TopologicalPropertyBase>& cur_property,
    const std::vector<AtomInMolecule>& indexed_atoms,
    const std::unordered_map<std::string, std::pair<std::string, std::string>>&
        core_atom_types);

}  // namespace topology_builder

}  // namespace combi_ff

#endif