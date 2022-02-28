// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef FAMILY_H
#define FAMILY_H

#include <iomanip>

#include "Alias.h"
#include "Pseudoatom.h"
#include "Substructure.h"

namespace combi_ff {

namespace enu {

/**************************
DECLARATION OF CLASS Family
**************************/
class Family {
 public:
  Family() = default;
  Family(const std::string& code, const std::string& version,
         const std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
         const std::vector<std::list<LambdaVector>>& lambda_ranges_vec,
         const size_t max_degree, const RangeVector& ranged_properties,
         const std::vector<SubstructureCollection>& substructures,
         const bool has_pseudoatom);

  const std::string& GetCode() const;
  const std::string& GetVersion() const;
  const std::vector<AtomVector<combi_ff::Atom>>& GetUsedAtomVectors() const;
  const std::vector<std::list<LambdaVector>>& GetLambdaRangesVec() const;
  size_t GetMaxDeg() const;
  const RangeVector& GetRangedProperties() const;
  const std::vector<SubstructureCollection>& GetSubstructures() const;
  const PseudoatomMap& GetPseudoatoms() const;
  bool GetHasPseudoatom() const;

  void Print() const;

 private:
  /*MEMBER VARIABLE DECLARATIONS*/
  std::string code{""};
  std::string version{""};
  std::vector<AtomVector<combi_ff::Atom>> used_atom_vectors{
      std::vector<AtomVector<combi_ff::Atom>>(0)};
  std::vector<std::list<LambdaVector>> lambda_ranges_vec{
      std::vector<std::list<LambdaVector>>(0)};
  size_t max_degree{4};
  RangeVector ranged_properties{RangeVector(0)};
  std::vector<SubstructureCollection> substructures{
      std::vector<SubstructureCollection>(0)};
  bool has_pseudoatom{false};
};

typedef std::vector<Family> FamilyVector;

void CreateFamilies(FamilyVector& families,
                    std::list<std::string>& family_definition_file_names,
                    const StringVector& used_families,
                    const AbstractSubstructureMap& abstr_substructures,
                    const AliasMap& alias_map,
                    const PseudoatomMap& pseudoatoms);

void AddUnitedAtom(const std::string& formula, size_t& j,
                   std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
                   std::list<LambdaVector>& lambda_ranges);

void AddPseudoatom(const std::string& formula, size_t& j,
                   std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
                   std::list<LambdaVector>& lambda_ranges,
                   const PseudoatomMap& pseudoatoms,
                   const std::string& version);

void AddAtom(const std::string& formula, size_t& j, const AliasMap& alias_map,
             std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
             std::list<LambdaVector>& lambda_ranges,
             const std::string& version);

void AddAtomType(std::string& name, const AliasMap& alias_map,
                 std::vector<AtomVector<combi_ff::Atom>>& atom_vectors,
                 const std::string& version);

void GetNextFamily(const std::string& version, FamilyVector& families,
                   const XmlElement_ptr family,
                   const StringVector& used_families, const AliasMap& alias_map,
                   const PseudoatomMap& pseudoatoms,
                   const AbstractSubstructureMap& abstr_substructures,
                   const bool all_fam, StringVector& codes,
                   StringVector& acros);

void GetFormula(const std::string& formula, const AliasMap& alias_map,
                const PseudoatomMap& pseudoatoms,
                std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
                std::vector<std::list<LambdaVector>>& lambda_ranges_vec,
                bool& has_pseudoatom, const std::string& version);

void GetRestriction(const std::string& range_value,
                    /*std::string propName, */ Range& r);

void GetSubstructures(const XmlElement_ptr substructures_xml,
                      const AbstractSubstructureMap& abstr_substructures,
                      std::vector<SubstructureCollection>& substructures,
                      const AliasMap& alias_map, const std::string& version);

// void GetSubstruc(std::istringstream& line, const AbstractSubstructureMap&
// Substructures, std::vector<SubstructureCollection>& substructureCols,
//                  const AliasMap& alias_map);

void AddAtomTypes(std::string& frg, size_t& j, std::vector<int>& XORidx,
                  std::vector<int>& ANDidx,
                  std::vector<AtomVector<combi_ff::Atom>>& atom_vectors,
                  const AliasMap& alias_map, const std::string& version);

void PrintFamilies(FamilyVector& families);

void CleanUpAtomVectorsFrag(
    std::vector<AtomVector<combi_ff::Atom>>& atom_vectors,
    std::vector<int>& XORidx, std::vector<int>& ANDidx);

void CleanUpAtomVectorsFam(
    bool hasXOR, std::vector<int>& XORidx,
    std::list<LambdaVector>& lambda_range,
    std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
    std::vector<std::list<LambdaVector>>& lambda_ranges_vec);

}  // namespace enu

}  // namespace combi_ff

#endif
