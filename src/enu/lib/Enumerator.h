// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef ENUMERATOR_H
#define ENUMERATOR_H

#include <list>
#include <memory>
#include <string>
#include <vector>

#include "AdjacencyMatrixEnu.h"
#include "LambdaVector.h"
#include "MaxFillAlg.h"
#include "Pseudoatom.h"
#include "Range.h"
#include "SmilesGenerator.h"
#include "Substructure.h"

namespace combi_ff {

namespace enu {

class EnumSpecifications;
class Family;

typedef SmilesGenerator<size_t, combi_ff::Atom> SmilesGeneratorEnu;

struct EnumeratorArgs {
  std::shared_ptr<AdjacencyMatrix> A{nullptr};
  std::unique_ptr<const std::string> formula_ptr{nullptr};
  std::unique_ptr<MaxFillAlg> max_fill_alg_ptr{nullptr};
  std::unique_ptr<AdjacencyMatrix> full_matrix_BU_ptr{nullptr};
  std::unique_ptr<AdjacencyMatrix> extended_matrix_BU_ptr{nullptr};
  std::unique_ptr<AdjacencyMatrix> full_matrix_ptr{nullptr};
  std::shared_ptr<AdjacencyMatrix> extended_matrix_ptr{nullptr};
  size_t n_pseudoatoms;
  EnumeratorArgs() = default;
};

class Enumerator {
 public:
  Enumerator(const EnumSpecifications& enum_spec,
             const std::string& output_file_name);

  Enumerator(const Family& family, const PseudoatomMap& pseudoatoms,
             const bool stereo, const std::string& output_file_name);

  size_t GetNumIsomers() const;

  void EnumerateIsomers();

  void EnumerateFormula(const AtomVector<combi_ff::Atom>& atom_types,
                        const LambdaVector& lambda, const size_t max_degree,
                        const bool stereo);

  void GetNImplicitAtoms(const AtomVector<combi_ff::Atom>& non_H_atoms,
                         size_t& num_pseudoatoms_hyd);

  void GetFullAndExtendedMatrixSkeletons(
      const size_t N_full, const AtomVector<combi_ff::Atom>& atoms_hat);

  bool CheckMolecule();

  void WriteMolecule();

  void GetFullMatrix();

  void GetExtendedMatrix(RepresentationSystem*& u0,
                         RepresentationSystem*& u_automorph,
                         RepresentationSystem*& u_id);

  void EnumerateStereoSmiles(const size_t n_rings,
                             const RepresentationSystem& u_automorph,
                             const RepresentationSystem& u_id,
                             const SmilesGeneratorEnu& smiles_gen);

  void PrintState(const size_t num_curr_isomers);

  void PrintConstitutionalIsomer(const std::string& smiles);
  void PrintClosingIsomerTag();

  std::string CreateCanonicalFormula(
      const LambdaVector& l, const AtomVector<combi_ff::Atom>& used_atoms);
  /* void Enumerator::printMolecule(Enu::AdjacencyMatrix& A,const std::string&
   * formula,const RepresentationSystem& u0,const RepresentationSystem&
   * uAutomorph,const RepresentationSystem& id,const size_t nPSA,const
   * Enu::AdjacencyMatrix& fullMatrixBU,const Enu::AdjacencyMatrix&
   * extendedMatrixBU);*/

 private:
  const std::string code;
  const std::vector<std::list<LambdaVector>>& lambda_ranges_vec;
  const std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors;
  const RangeVector& ranged_properties;
  // note: substructures and pseudoatoms are pointers s.t. we can use NULL for
  // direct enumeration. If direct enumeration is ever removed,these can be
  // changed to references
  const std::vector<SubstructureCollection>* substructures = {NULL};
  const PseudoatomMap* pseudoatoms = {NULL};
  const size_t max_degree;
  const bool stereo;
  std::ofstream output_file;
  int deg_unsaturations;
  size_t num_isomers, num_isomers_for_listing, num_isomers_for_formula;
  const bool has_pseudoatom;
  size_t formula_idx;
  EnumeratorArgs enumerator_arguments;
  std::vector<size_t> ranges;
};

}  // namespace enu

}  // namespace combi_ff

#endif
