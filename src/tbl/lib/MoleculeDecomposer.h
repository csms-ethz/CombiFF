// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef MOLECULEDECOMPOSER_H_
#define MOLECULEDECOMPOSER_H_

#include <fstream>
#include <string>
#include <vector>

#include "InputOutput.h"
#include "Link.h"
#include "MatrixTbl.h"
#include "TblAtom.h"
#include "XmlParser.h"
namespace combi_ff {

namespace topology_builder {

class Link;
class TblFragment;
class IOFileProps;

class FamilyDecomposer {
 public:
  FamilyDecomposer() = default;
  FamilyDecomposer(const std::string& filename_family_isomer_enumeration,
                   const IOFileProperties& io_file_properties,
                   const std::vector<TblFragment>& tbl_fragments);

  void CheckXMLFormat(XmlElement& family_isomer_enumeration_root);
  size_t GetNumIsomers(const XmlElement& family_isomer_enumeration_root);
  void CreateFamilyMoleculeDecompositions();

  const std::string& GetFilename() const;
  const std::string& GetFamilyCode() const;

 private:
  std::string filename_molecule_decomposition{""};
  std::string filename_family_isomer_enumeration{""};
  std::string family_code{""};
  const IOFileProperties& io_file_properties{IOFileProperties()};
  const std::vector<TblFragment>& tbl_fragments{std::vector<TblFragment>(0)};
};

class MoleculeDecomposer : private FamilyDecomposer {
 public:
  MoleculeDecomposer() = default;
  MoleculeDecomposer(const std::string& smiles,
                     const std::vector<TblFragment>& tbl_fragments,
                     XmlElement& molecule_decomposition);

  void CreateMoleculeDecomposition(XmlElement& molecule_decomposition);

  bool CheckValidity(XmlElement& molecule_decomposition);

  size_t GetMatchIndex(const ComparisonMatrix& M, const size_t i);

  bool FindBondLinkingMatch(const TblFragment& frag);

  bool UllmannMatch(ComparisonMatrix& M, std::vector<bool>& matched_columns,
                    std::vector<bool>& matched_rows, int k, const size_t n,
                    const size_t m, const FragmentMatrixTbl& fragment_matrix,
                    const std::vector<TblAtom>& fragment_tbl_atoms);

  bool Refine(combi_ff::ComparisonMatrix& M, int k, const size_t n,
              const size_t m, const FragmentMatrixTbl& fragment_matrix,
              const std::vector<TblAtom>& fragment_tbl_atoms);

 private:
  topology_builder::AdjacencyMatrix A{0};
  size_t N{0};
  const std::vector<TblFragment>& tbl_fragments{std::vector<TblFragment>(0)};
  std::vector<bool> is_assigned_core_atom{std::vector<bool>(0)};
  std::vector<bool> is_assigned_core_atom_all_true{std::vector<bool>(0)};
  std::vector<int> is_assigned_link_atom{std::vector<int>(0)};
  std::vector<const TblAtom*> core_atoms{std::vector<const TblAtom*>(0)};
  std::vector<std::vector<const TblAtom*>> link_atoms{
      std::vector<std::vector<const TblAtom*>>(0)};
  std::vector<const TblFragment*> fragments{std::vector<const TblFragment*>(0)};
  std::vector<HalfLink> half_links{std::vector<HalfLink>(0)};
  std::vector<Link> links{std::vector<Link>(0)};
  bool all_assigned{false};
};

}  // namespace topology_builder

}  // namespace combi_ff

#endif
