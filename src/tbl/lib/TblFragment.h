// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef TBFRAGMENT_H_
#define TBFRAGMENT_H_

#include "AtomTypes.h"
#include "MatrixTbl.h"
#include "TopologicalProperty.h"

namespace combi_ff {

namespace topology_builder {

class TblAtom;

typedef std::unordered_map<std::string, double> BondDegreeMap;

static BondDegreeMap bond_degrees{{"single", 1},
                                  {"double", 2},
                                  {"triple", 3},
                                  {"quadruple", 4},
                                  {"aromatic", 1.5}};

class TblFragment {
 public:
  TblFragment() = default;
  TblFragment(const std::string& version, const std::string& code, size_t N);
  TblFragment(const std::string& version, const std::string& code, size_t N,
              const std::vector<TblAtom>& atom_types,
              const std::vector<std::shared_ptr<TopologicalPropertyBase>>&
                  topological_properties,
              const FragmentMatrixTbl& M);
  TblFragment(const TblFragment& frag);

  static const std::string fragment_prefix;

  std::string GetCode() const;
  const std::vector<TblAtom>& GetTblAtoms() const;
  const std::vector<std::shared_ptr<TopologicalPropertyBase>>&
  GetTopologicalProperties() const;
  const std::string& GetVersion() const;
  const FragmentMatrixTbl& GetMatrix() const;

  size_t GetCoreIndex(const size_t idx) const;

  size_t GetNS() const;
  size_t GetNA() const;
  size_t GetND() const;
  size_t GetNT() const;
  size_t GetNQ() const;
  size_t GetNumCoreAtoms() const;
  bool HasCycle() const;
  void DetermineCoreAtoms();
  const TblAtom& GetTblAtom(const size_t idx) const;
  const TblAtom& GetTblAtom(const std::string& id) const;
  const TblAtom* GetCoreAtom(const std::string& id) const;

 private:
  std::string version{""};
  FragmentMatrixTbl F{FragmentMatrixTbl(0)};
  std::string code{""};
  std::vector<TblAtom> atom_types{std::vector<TblAtom>(0)};
  std::vector<const TblAtom*> core_atoms{std::vector<const TblAtom*>(
      0)};  // if an atom in atom_types is a core atom, core_atoms points to it,
            // if it's a link atom it points to the bonded core atom
  std::vector<std::shared_ptr<TopologicalPropertyBase>> topological_properties{
      std::vector<std::shared_ptr<TopologicalPropertyBase>>(0)};
  size_t n_single{0};
  size_t n_aromatic{0};
  size_t n_double{0};
  size_t n_triple{0};
  size_t n_quadruple{0};
  size_t n_core_atoms{0};
  bool cyclic{false};
};

bool GetPriority(const TblFragment& f1, const TblFragment& f2);
void CreateTblFragments(const std::list<std::string>& frag_file_names,
                        std::vector<TblFragment>& tbl_fragments,
                        const std::vector<std::string>& fragments_to_use,
                        const AtomTypeSetMap& atom_type_set);
void AddFragment(std::vector<TblFragment>& tbl_fragments,
                 const std::vector<std::string>& fragments_to_use,
                 XmlElement_ptr fragment, const AtomTypeSetMap& atom_type_set,
                 const std::string& version);
void CheckFragmentXML(const XmlElement_ptr fragment);
void AddAtoms(const XmlElement_ptr atoms, const AtomTypeSetMap& atom_type_set,
              std::vector<TblAtom>& atom_types, const std::string& version);
void AddBonds(const XmlElement_ptr bonds, std::vector<TblAtom>& atom_types,
              FragmentMatrixTbl& M, const std::string& code);
void AddTopologicalProperties(
    const XmlElement_ptr topological_property,
    std::vector<std::shared_ptr<TopologicalPropertyBase>>&
        topological_properties,
    std::vector<TblAtom>& atom_types, const std::string& code);
void SetLinkageType(const XmlElement_ptr link_type, TblAtom& new_atom,
                    bool& is_core_atom);

void CheckAtomTypeXML(const XmlElement_ptr atom_info,
                      const AtomTypeSetMap& atom_type_set,
                      const std::string& version);
void CheckLinkTypeXML(const XmlElement_ptr link_type);
void CheckAtomsXML(const XmlElement_ptr atoms);
void CheckAtomXML(const XmlElement_ptr atom);
void CheckBondXML(const XmlElement_ptr bonds);
void CheckBondsXML(const XmlElement_ptr bond);
void CheckInvolvedAtomsXML(const XmlElement_ptr involved_atoms,
                           const size_t num_involved_atoms);
void CheckBondValidity(TblAtom& atom1, TblAtom& atom2, const double degree);
void CheckTopologicalPropertyXML(const XmlElement_ptr topological_property);
std::vector<size_t> CreateInvolvedAtomIndexVector(
    const XmlElement_ptr involved_atoms, std::vector<TblAtom>& atom_types,
    const std::string& code);

}  // namespace topology_builder

}  // namespace combi_ff

#endif
