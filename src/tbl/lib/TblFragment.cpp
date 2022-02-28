// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "TblFragment.h"

#include "TblAtom.h"
#include "TopologicalProperty.h"
#include "version.h"

namespace combi_ff {

namespace topology_builder {

/***********************************
IMPLEMENTATION OF TblFragment METHODS
***********************************/

/***********
CONSTRUCTORS
***********/

TblFragment::TblFragment(const std::string& version, const std::string& code,
                         size_t N)
    : version(version), F(FragmentMatrixTbl(N)), code(code), atom_types(N) {}

TblFragment::TblFragment(
    const std::string& version, const std::string& code, size_t N,
    const std::vector<TblAtom>& atom_types,
    const std::vector<std::shared_ptr<TopologicalPropertyBase>>&
        topological_properties,
    const FragmentMatrixTbl& M)
    : version(version),
      F(M),
      code(code),
      atom_types(atom_types),
      topological_properties(topological_properties),
      n_single(0),
      n_aromatic(0),
      n_double(0),
      n_triple(0),
      n_quadruple(0),
      n_core_atoms(0),
      cyclic(false) {
  // calculate the number of single, double, triple and quadruple bonds in the
  // Fragment
  for (size_t i = 0; i < N; i++) {
    for (size_t j = i + 1; j < N; j++) {
      if (F.GetElement(i, j) == 1)
        n_single++;

      else if (F.GetElement(i, j) == 2)
        n_double++;

      else if (F.GetElement(i, j) == 3)
        n_triple++;

      else if (F.GetElement(i, j) == 4)
        n_quadruple++;

      else if (F.GetElement(i, j) == 1.5)
        n_aromatic++;

      else if (F.GetElement(i, j) != 0)
        throw input_error("fragment has unrecognized bond degree " +
                          std::to_string(F.GetElement(i, j)) + '\n');
    }
  }

  for (auto&& atom : atom_types) {
    if (atom.GetLinkageType() == core) n_core_atoms++;
  }

  // Check if the fragment has a cyclicle
  cyclic = F.HasCycle();
}

TblFragment::TblFragment(const TblFragment& frag)
    : version(frag.version),
      F(frag.F),
      code(frag.code),
      atom_types(frag.atom_types),
      topological_properties(frag.topological_properties),
      n_single(frag.n_single),
      n_aromatic(frag.n_aromatic),
      n_double(frag.n_double),
      n_triple(frag.n_triple),
      n_quadruple(frag.n_quadruple),
      n_core_atoms(frag.n_core_atoms),
      cyclic(frag.cyclic) {}

/*
static fragment prefix
*/
const std::string TblFragment::fragment_prefix = "~";

/*************
GetTER METHODS
*************/
std::string TblFragment::GetCode() const { return code; }

const std::vector<TblAtom>& TblFragment::GetTblAtoms() const {
  return atom_types;
}

const FragmentMatrixTbl& TblFragment::GetMatrix() const { return F; }

size_t TblFragment::GetNS() const { return n_single; }

size_t TblFragment::GetNA() const { return n_aromatic; }

size_t TblFragment::GetND() const { return n_double; }

size_t TblFragment::GetNT() const { return n_triple; }

size_t TblFragment::GetNQ() const { return n_quadruple; }

size_t TblFragment::GetNumCoreAtoms() const { return n_core_atoms; }

bool TblFragment::HasCycle() const { return cyclic; }

const std::vector<std::shared_ptr<TopologicalPropertyBase>>&
TblFragment::GetTopologicalProperties() const {
  return topological_properties;
}
const std::string& TblFragment::GetVersion() const { return version; }

// Get the core index that belongs to the atom with the index idx. If the atom
// itself is a core atom, idx is returned, and otherwise, the connections of
// atom with index idx are searched until a core atom is found
size_t TblFragment::GetCoreIndex(const size_t idx) const {
  if (atom_types[idx].GetLinkageType() == core)
    return idx;

  else {
    for (size_t i = 0; i < F.GetN(); i++) {
      if (F.GetElement(i, idx) && atom_types[i].IsCoreAtom()) return i;
    }
  }

  throw std::runtime_error("atom " + std::to_string(idx) + " in fragment " +
                           code + " is not connected to a core atom\n");
}

void TblFragment::DetermineCoreAtoms() {
  core_atoms = std::vector<const TblAtom*>(0);

  for (auto&& atom : atom_types) {
    if (atom.GetLinkageType() == core)
      core_atoms.push_back(&atom);

    else {
      bool found(false);

      for (auto&& nbr : atom.GetNeighbours()) {
        if (atom_types[nbr].GetLinkageType() == core) {
          core_atoms.push_back(&atom_types[nbr]);
          found = true;
          break;
        }
      }

      if (!found)
        throw input_error("link atom " + atom.GetAtomID() +
                          " is not connected to a core atom in fragment " +
                          code);
    }
  }
}

const TblAtom& TblFragment::GetTblAtom(const size_t idx) const {
  return atom_types[idx];
}

const TblAtom& TblFragment::GetTblAtom(const std::string& id) const {
  for (auto&& atom : atom_types) {
    if (atom.GetAtomID() == id) return atom;
  }

  throw std::runtime_error("no atom with id " + id + " in fragment " + code +
                           "\n");
}

const TblAtom* TblFragment::GetCoreAtom(const std::string& id) const {
  for (size_t i = 0; i < atom_types.size(); i++) {
    if (atom_types[i].GetAtomID() == id) return core_atoms[i];
  }

  throw std::runtime_error("no atom with id " + id + " in fragment " + code +
                           "\n");
}

/********************************************************************************
CREATE THE TblFragment GIVEN IN THE INPUT FILE fragFile AND SAVE THEM TO
tbl_fragments
********************************************************************************/
void CreateTblFragments(const std::list<std::string>& frag_file_names,
                        std::vector<TblFragment>& tbl_fragments,
                        const std::vector<std::string>& fragments_to_use,
                        const AtomTypeSetMap& atom_type_set) {
  std::string blockName, s;
  std::cout << "creating TblFragments\n";
  tbl_fragments = std::vector<TblFragment>(fragments_to_use.size());

  if (!frag_file_names.size())
    throw std::runtime_error("didn't specify any fragment files\n");

  for (auto&& fragFileName : frag_file_names) {
    XmlParserIn parser(fragFileName, XmlParserIn::read_all);
    // std::cout << parser << std::endl;
    const XmlTree& tree = parser.GetTree_const();
    const XmlElement& root = tree.GetRoot_const();
    tree.CheckRootTagName("fragments");
    root.CheckAttribute("version");
    const std::string& version = root.attributes.find("version")->second;

    if (version != combi_ff::current_version)
      std::cout << "?Warning: currently running combi_ff version "
                << combi_ff::current_version << " but atom type file "
                << fragFileName << " is version " << version << "\n";

    tbl_fragments.reserve(root.children.size());

    for (auto&& child : root.children) {
      child->CheckTagName("fragment");
      AddFragment(tbl_fragments, fragments_to_use, child, atom_type_set,
                  version);
    }
  }

  // Check that every fragment given in fragments_to_use was actually found and
  // added to tbl_fragments
  if (fragments_to_use.size()) {
    for (size_t i = 0; i < tbl_fragments.size(); i++) {
      std::cout << tbl_fragments[i].GetCode() << std::endl;

      if (tbl_fragments.at(i).GetCode() == "") {
        throw input_error("fragment " + fragments_to_use.at(i) +
                          " which is listed in fragments_to_use was not found "
                          "in any of the specified fragment files\n");
      }
    }
  }

  std::sort(tbl_fragments.begin(), tbl_fragments.end(), GetPriority);

  for (auto&& frag : tbl_fragments) frag.DetermineCoreAtoms();

  std::cout << "fragments: \n";

  for (auto&& frag : tbl_fragments) {
    std::cout << frag.GetCode() << '\n';

    for (auto&& at : frag.GetTblAtoms())
      std::cout << at.GetAtomTypes().GetName() << " ";

    std::cout << std::endl;
    frag.GetMatrix().Print();
    std::cout << std::endl;
  }
}

// determine the priority between two fragments
bool GetPriority(const TblFragment& f1, const TblFragment& f2) {
  if (f1.HasCycle() && !(f2.HasCycle()))
    return true;

  else if (!(f1.HasCycle()) && f2.HasCycle())
    return false;

  else if (f1.GetNT() > f2.GetNT())
    return true;

  else if (f1.GetNT() < f2.GetNT())
    return false;

  else if (f1.GetNA() > f2.GetNA())
    return true;

  else if (f1.GetNA() < f2.GetNA())
    return false;

  else if (f1.GetND() > f2.GetND())
    return true;

  else if (f1.GetND() < f2.GetND())
    return false;

  else if (f1.GetNS() > f2.GetNS())
    return true;

  else if (f1.GetNS() < f2.GetNS())
    return false;

  else if (f1.GetMatrix().GetN() > f2.GetMatrix().GetN())
    return true;

  else if (f1.GetMatrix().GetN() < f2.GetMatrix().GetN())
    return false;

  else if (f1.GetNumCoreAtoms() > f2.GetNumCoreAtoms())
    return true;

  else if (f1.GetNumCoreAtoms() < f2.GetNumCoreAtoms())
    return false;

  return false;
}

void AddFragment(std::vector<TblFragment>& tbl_fragments,
                 const std::vector<std::string>& fragments_to_use,
                 XmlElement_ptr fragment, const AtomTypeSetMap& atom_type_set,
                 const std::string& version) {
  CheckFragmentXML(fragment);
  std::string code = fragment->attributes["code"];

  if (code.substr(0, TblFragment::fragment_prefix.size()) !=
      TblFragment::fragment_prefix)
    throw input_error("Please use prefix " + TblFragment::fragment_prefix +
                      " in front of fragment codes. Missing prefix for " +
                      code + '\n');

  std::vector<TblAtom> atom_types;
  size_t N(0);
  FragmentMatrixTbl M;
  std::vector<std::shared_ptr<TopologicalPropertyBase>> topological_properties(
      0);

  for (auto&& topological_propertyList : fragment->children) {
    if (topological_propertyList->tag == "atoms") {
      AddAtoms(topological_propertyList, atom_type_set, atom_types, version);
      N = atom_types.size();
      M = FragmentMatrixTbl(N);

    } else if (topological_propertyList->tag == "bonds") {
      if (!N)
        throw input_error("defining bonds, but no atoms found in fragment " +
                          code + "\n");

      AddBonds(topological_propertyList, atom_types, M, code);
      AddTopologicalProperties(topological_propertyList, topological_properties,
                               atom_types, code);

    } else
      AddTopologicalProperties(topological_propertyList, topological_properties,
                               atom_types, code);
  }

  if (fragments_to_use.size()) {
    auto&& it =
        std::find(fragments_to_use.begin(), fragments_to_use.end(), code);

    if (it != fragments_to_use.end())
      tbl_fragments.at(std::distance(fragments_to_use.begin(), it)) =
          (TblFragment(version, code, N, atom_types, topological_properties,
                       M));

  } else
    tbl_fragments.push_back(
        TblFragment(version, code, N, atom_types, topological_properties, M));
}

void AddAtoms(const XmlElement_ptr atoms, const AtomTypeSetMap& atom_type_set,
              std::vector<TblAtom>& atom_types, const std::string& version) {
  CheckAtomsXML(atoms);
  bool hasCoreAtom(false);

  for (auto&& atom : atoms->children) {
    CheckAtomXML(atom);
    TblAtom new_atom;
    new_atom.SetAtomID(atom->attributes.find("id")->second);

    for (auto&& atom : atom_types) {
      if (atom.GetAtomID() == new_atom.GetAtomID())
        throw input_error("listed more than one atom with id " +
                          new_atom.GetAtomID() + " in a fragment\n");
    }

    for (auto&& atom_info : atom->children) {
      if (atom_info->tag == "atomtype") {
        CheckAtomTypeXML(atom_info, atom_type_set, version);
        new_atom.SetAtomTypeSet(atom_type_set.find(atom_info->value)->second);

      } else if (atom_info->tag == "linktype") {
        CheckLinkTypeXML(atom_info);
        SetLinkageType(atom_info, new_atom, hasCoreAtom);

      } else
        throw unexpected_tag_error(atom_info->tag, "atomtype or linktype");
    }

    atom_types.push_back(new_atom);
  }

  if (!hasCoreAtom)
    throw input_error("fragment must have at least one core atom!\n");
}

void SetLinkageType(const XmlElement_ptr link_type, TblAtom& new_atom,
                    bool& is_core_atom) {
  std::string linkaGetype = link_type->attributes.find("type")->second;

  if (linkaGetype == "core") {
    new_atom.SetLinkageType(core);
    is_core_atom = true;

    if (new_atom.GetAtomTypes().GetAtomSet().size() != 1)
      throw input_error(
          "core atom " + new_atom.GetAtomID() +
          " can only have an att as \"atomtype\", but has the atomtypeSet " +
          new_atom.GetAtomTypes().GetName() + "\n");

  } else if (linkaGetype == "link")
    new_atom.SetLinkageType(link);

  else
    throw input_error("unknown link_type " + linkaGetype);
}

void CheckLinkTypeXML(const XmlElement_ptr link_type) {
  link_type->CheckAttributeSize(1);
  link_type->CheckAttribute("type");
  link_type->CheckNumberOfChildren_equal(0);
}

void CheckAtomTypeXML(const XmlElement_ptr atomtype,
                      const AtomTypeSetMap& atom_type_set,
                      const std::string& version) {
  atomtype->CheckAttributeSize(0);
  atomtype->CheckValue();
  auto&& it = atom_type_set.find(atomtype->value);

  if (it == atom_type_set.end())
    throw input_error("couldn't find atom type (Set) " + atomtype->value +
                      " in atom type Set\n");

  if (it->second.GetVersion() != version)
    std::cout << "?Warning: current fragment file is version " << version
              << ", but atom type (Set) " << atomtype->value << " is version "
              << it->second.GetVersion() << "\n";
}

void CheckAtomsXML(const XmlElement_ptr atoms) {
  atoms->CheckAttributeSize(0);
  atoms->CheckNumberOfChildren_atLeast(1);
}

void CheckAtomXML(const XmlElement_ptr atom) {
  atom->CheckTagName("atom");
  atom->CheckAttributeSize(1);
  atom->CheckAttribute("id");
  atom->CheckNumberOfChildren_equal(2);
  const std::string& id = atom->attributes.find("id")->second;

  if (!(id.front() >= 'a' && id.front() <= 'z'))
    throw input_error(
        "Please start atom ids with lower case letter in <atom> tag. Not the "
        "case for atom " +
        id + '\n');
}

void AddBonds(const XmlElement_ptr bonds, std::vector<TblAtom>& atom_types,
              FragmentMatrixTbl& M, const std::string& code) {
  CheckBondsXML(bonds);

  for (auto&& bond : bonds->children) {
    CheckBondXML(bond);
    std::string degreeString = bond->attributes["degree"];
    double degree = bond_degrees[degreeString];
    const XmlElement_ptr involved_atoms = bond->GetFirstChild();
    CheckInvolvedAtomsXML(involved_atoms, 2);
    std::vector<size_t> involved_atomsIdx =
        CreateInvolvedAtomIndexVector(involved_atoms, atom_types, code);
    auto&& atom1 = atom_types[involved_atomsIdx[0]];
    auto&& atom2 = atom_types[involved_atomsIdx[1]];
    CheckBondValidity(atom1, atom2, degree);
    atom1.AddNeighbour(involved_atomsIdx[1]);
    atom2.AddNeighbour(involved_atomsIdx[0]);
    M.SetElement(involved_atomsIdx[0], involved_atomsIdx[1], degree);
  }
}

std::vector<size_t> CreateInvolvedAtomIndexVector(
    const XmlElement_ptr involved_atoms, std::vector<TblAtom>& atom_types,
    const std::string& code) {
  std::vector<size_t> involved_atomsIdx(0);
  involved_atomsIdx.reserve(involved_atoms->children.size());

  for (auto&& involved_atom : involved_atoms->children) {
    std::string atomID = involved_atom->value;
    bool found = false;

    for (size_t i = 0; i < atom_types.size(); i++) {
      if (atomID == atom_types[i].GetAtomID()) {
        found = true;
        involved_atomsIdx.push_back(i);
        break;
      }
    }

    if (!found)
      throw input_error("atom with id " + atomID +
                        " not found in atoms list in fragment " + code + "\n");
  }

  return involved_atomsIdx;
}

void CheckBondValidity(TblAtom& atom1, TblAtom& atom2, const double degree) {
  if (atom1.GetLinkageType() == link && atom2.GetLinkageType() == link)
    throw input_error("not allowed to link 2 link atoms toGether: " +
                      atom1.GetAtomID() + " and " + atom2.GetAtomID() + '\n');

  // should I Check if bond degree is not larger than deg of att for core atoms?
  /*if(atom1.GetLinkageType() == core){

  }

  if(atom2.GetLinkageType() == core){

  }*/
}

void CheckInvolvedAtomsXML(const XmlElement_ptr involved_atoms,
                           const size_t num_involved_atoms) {
  involved_atoms->CheckTagName("involved_atoms");
  involved_atoms->CheckAttributeSize(0);
  involved_atoms->CheckNumberOfChildren_equal(num_involved_atoms);

  for (auto&& involved_atom : involved_atoms->children) {
    involved_atom->CheckTagName("involved_atom");
    involved_atom->CheckAttributeSize(0);
    involved_atom->CheckValue();
    involved_atom->CheckNumberOfChildren_equal(0);
  }
}

void CheckBondsXML(const XmlElement_ptr bonds) {
  bonds->CheckAttributeSize(0);
  bonds->CheckNoValue();
}

void CheckBondXML(const XmlElement_ptr bond) {
  bond->CheckTagName("bond");
  bond->CheckAttributeSize(1);
  bond->CheckAttribute("degree");
  bond->CheckNumberOfChildren_equal(1);
  std::string degree = bond->attributes.find("degree")->second;

  if (bond_degrees.find(degree) == bond_degrees.end())
    throw input_error(
        "unknown bond degree " + degree +
        " (options are single, double, triple, quadruple, aromatic)\n");
}

void AddTopologicalProperties(
    const XmlElement_ptr topological_propertyList,
    std::vector<std::shared_ptr<TopologicalPropertyBase>>&
        topological_properties,
    std::vector<TblAtom>& atom_types, const std::string& code) {
  for (auto&& topological_property : topological_propertyList->children) {
    CheckTopologicalPropertyXML(topological_property);
    topological_properties.push_back(CreateNewTopologicalProperty(
        topologicalPropertyMap.find(topological_property->tag)->second));
    std::shared_ptr<TopologicalPropertyBase>& curType =
        topological_properties.back();
    const XmlElement_ptr involved_atoms = topological_property->GetFirstChild();
    CheckInvolvedAtomsXML(involved_atoms, curType->GetNumInvolvedAtoms());

    for (auto&& att : topological_property->attributes)
      curType->AddAttribute(att.first, att.second);

    curType->SetInvolvedAtoms(
        CreateInvolvedAtomIndexVector(involved_atoms, atom_types, code));
  }
}

void CheckTopologicalPropertyXML(const XmlElement_ptr topological_property) {
  // topological_property can have 0 or 1 attributes
  if (topological_property->attributes.size() &&
      topological_property->tag != "bond") {
    topological_property->CheckAttributeSize(1);
    topological_property->CheckAttribute("parameter");
  }

  if (topologicalPropertyMap.find(topological_property->tag) ==
      topologicalPropertyMap.end())
    throw input_error(topological_property->tag +
                      " is not a known topological property\n");

  topological_property->CheckNumberOfChildren_equal(1);
}

void CheckFragmentXML(const XmlElement_ptr fragment) {
  fragment->CheckAttributeSize(1);
  fragment->CheckAttribute("code");
}

}  // namespace topology_builder

}  // namespace combi_ff