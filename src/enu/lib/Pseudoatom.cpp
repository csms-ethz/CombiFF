#include "Pseudoatom.h"
#include "printInfo.h"
#include "exceptions.h"
#include "ContainerOperators.h"
#include "version.h"

#include <sstream>

namespace combi_ff {

namespace enu {

/***********************************************
IMPLEMENTATION OF AbstractSubstructure METHODS
***********************************************/

/***********
CONSTRUCTORS
***********/

Pseudoatom::Pseudoatom(const FragmentMatrix& M,
                       const std::string& code,
                       const std::string& version,
                       size_t degree)
  :	M(M), code(code), version(version), degree(degree) {}


/***************
GetTER FUNCTIONS
***************/
const FragmentMatrix& Pseudoatom::GetMatrix() const {
  return M;
}

const std::string& Pseudoatom::GetCode() const {
  return code;
}

const std::string& Pseudoatom::GetVersion() const {
  return version;
}

size_t Pseudoatom::GetDegree() const {
  return degree;
}

const AtomVector<combi_ff::Atom>& Pseudoatom::GetAtoms() const {
  return M.GetAtomVector();
}


/**********************
CREATE FAMILY FRAGMENTS
**********************/
void CreatePseudoatoms(PseudoatomMap& pseudoatoms,
                       const std::list<std::string>& pseudoatom_file_names) {
  if (!pseudoatom_file_names.size()) {
    throw combi_ff::input_warning("no pseudoatom files");
    return;
  }

  std::cout << "creating pseudoatoms\n";

  for (auto && pseudoatom_file_name : pseudoatom_file_names) {
    XmlParserIn parser(pseudoatom_file_name, XmlParserIn::read_all);
    const XmlElement& pseudoatoms_root = parser.GetTree().GetRoot();
    pseudoatoms_root.CheckTagName("pseudoatoms");
    pseudoatoms_root.CheckAttribute("version");
    pseudoatoms_root.CheckNumberOfChildren_atLeast(1);
    const std::string& version = pseudoatoms_root.attributes.find("version")->second;

    if (version != combi_ff::current_version)
      std::cout << "?Warning: currently running combi_ff version " << combi_ff::current_version
                << " but pseudoatom file " << pseudoatom_file_name
                << " is version " << version << "\n";

    for (auto && pseudoatom : pseudoatoms_root.children)
      GetNextPseudoatom(pseudoatoms, pseudoatom, version);
  }

  std::cout << "******************************************************\n";
}

/****************
create pseudoatom
****************/
void GetNextPseudoatom(PseudoatomMap& pseudoatoms,
                       const XmlElement_ptr pseudoatom,
                       const std::string& version) {
  FragmentMatrix M;
  size_t degree;
  pseudoatom->CheckTagName("pseudoatom");
  pseudoatom->CheckAttributeSize(1);
  pseudoatom->CheckNumberOfChildren_equal(2);
  pseudoatom->CheckNoValue();
  pseudoatom->CheckAttribute("code");
  const std::string& code = "'" + pseudoatom->attributes.find("code")->second + "'";

  if (pseudoatoms.find(code) != pseudoatoms.end())
    throw input_error("pseudoatom with code " + code + " occurs more than once\n");

  auto&& pseudoatom_property_it = pseudoatom->children.begin();
  const XmlElement_ptr atoms_xml = *pseudoatom_property_it;
  atoms_xml->CheckTagName("atoms");
  atoms_xml->CheckNumberOfChildren_atLeast(1);
  atoms_xml->CheckNoValue();
  atoms_xml->CheckAttributeSize(0);
  AtomVector<combi_ff::Atom> atoms(0);

  for (auto && atom_xml : atoms_xml->children) {
    atom_xml->CheckTagName("atom");
    atom_xml->CheckValue();

    if (atom_xml->value == "H")
      throw combi_ff::input_error("In " + code +
                                  ": Pseudoatoms cannot contain explicit hydrogens. Please use implicit hydrogen notation if you wish to include hydrogen atoms in the pseudoatom");

    atoms.push_back(Atom(atom_xml->value));
  }

  const XmlElement_ptr adjacency_stack = *(++pseudoatom_property_it);
  adjacency_stack->CheckTagName("adjacency_stack");
  //adjacency_stack->checkValue();
  M = FragmentMatrix(atoms.size());
  std::istringstream adjacency_stack_stream(adjacency_stack->value);
  AdjacencyVector v(0);
  int value;

  while (adjacency_stack_stream >> value)
    v.push_back(value);

  size_t idx = 0;

  for (size_t i = 0; i < M.GetN() - 1; i++) {
    for (size_t j = i + 1; j < M.GetN(); j++) {
      if (idx >= v.size())
        throw combi_ff::input_error("stack too small for number of atoms given for pseudoatom " + code);

      M.SetElement(i, j, v[idx++]);
    }
  }

  M.SetAtomVector(atoms);

  if (idx != v.size())
    throw combi_ff::input_error("stack too big for number of atoms given for pseudoaotm " + code);

  if (M.GetAtomVector().front().GetDegree() <= M.AccumulateRow(0))
    throw combi_ff::input_error("first atom of the pseudoatom is already fully linked for pseudoatom" +
                                code);

  degree = M.GetAtomVector().front().GetDegree() - (int)M.AccumulateRow(0);

  for (size_t row = 1; row < M.GetN(); row++) {
    if (M.GetAtomVector()[row].GetDegree() != M.AccumulateRow(row))
      throw combi_ff::input_error("atom " + std::to_string(row) +
                                  " is not fully linked (either too many or not enough bonds) for pseudoatom " + code);
  }

  pseudoatoms.emplace(code, Pseudoatom(M, code, version, degree));
  std::cout << code << " (degree: " << degree << ")" << std::endl;
  std::cout << M << std::endl;
}

} //namespace enu

} //namespace combi_ff