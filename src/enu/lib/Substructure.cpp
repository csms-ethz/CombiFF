#include "Substructure.h"
#include "exceptions.h"
#include "version.h"
#include <sstream>

namespace combi_ff {

namespace enu {

/***********
CONSTRUCTORS
***********/
AbstractSubstructure::AbstractSubstructure(const size_t N, const AdjacencyVector& v,
                                           const std::string& code,
                                           const std::string& version)
  : M(FragmentMatrix(N, AtomVector<combi_ff::Atom>(N, Atom("*")), v)),
    code(code),
    version(version) {}

AbstractSubstructure::AbstractSubstructure(const SymmetricalMatrix<size_t>& M,
                                           const std::string& code,
                                           const std::string& version)
  : M(M),
    code(code),
    version(version) {}

AbstractSubstructure::AbstractSubstructure(const AbstractSubstructure& FF)
  : M(FF.GetMatrix().GetN(), (FF.GetMatrix().GetElements())),
    code(FF.GetCode()),
    version(FF.GetVersion()) {}

/***************
GetTER FUNCTIONS
***************/
const SymmetricalMatrix<size_t>& AbstractSubstructure::GetMatrix() const {
  return M;
}

const std::string& AbstractSubstructure::GetCode() const {
  return code;
}

const std::string& AbstractSubstructure::GetVersion() const {
  return version;
}

/***********
CONSTRUCTORS
***********/
SubstructureCollection::SubstructureCollection(const AbstractSubstructure& F,
                                               const Range& r,
                                               const AtomVector<combi_ff::Atom>& atoms,
                                               const bool AND,
                                               const  bool XOR)
  : AND(bool(AND)),
    XOR(bool(XOR)),
    code(F.GetCode()),
    r(r) {
  FragmentMatrix M(F.GetMatrix().GetN(), atoms, (F.GetMatrix().GetElements()));
  substructure_matrices = (std::vector<FragmentMatrix>(1, M));
}

SubstructureCollection::SubstructureCollection(const AbstractSubstructure& F,
                                               const Range& r,
                                               const  FragmentMatrix& M,
                                               const bool AND,
                                               const  bool XOR)
  : AND(AND),
    XOR(XOR),
    code(F.GetCode()),
    r(r) {
  substructure_matrices = std::vector<FragmentMatrix>(1, M);
}

SubstructureCollection::SubstructureCollection(const AbstractSubstructure& F,
                                               const Range& r,
                                               std::vector<AtomVector<combi_ff::Atom>>& AtomVectors,
                                               const  bool AND,
                                               const  bool XOR)
  : substructure_matrices(0),
    AND(AND),
    XOR(XOR),
    code(F.GetCode()),
    r(r) {
  substructure_matrices.reserve(AtomVectors.size());

  for (auto && av : AtomVectors) {
    size_t n = av.size();
    size_t num_hyd(n);

    for (auto && a : av)
      num_hyd += a.GetNumFixedHydrogens();

    //if there are united atoms in the current atom vector, fully exend the matrix to explicitly include them
    if (num_hyd > n) {
      const SymmetricalMatrix<size_t> matrix_ = F.GetMatrix();
      FragmentMatrix substructure_matrix(num_hyd);
      AtomVector<combi_ff::Atom> atoms_tmp = av;
      size_t next_H_index(n);

      for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++)
          substructure_matrix.SetElement(i, j, matrix_.GetElement(i, j));

        for (size_t k = 0; k < av[i].GetNumFixedHydrogens(); k++) {
          substructure_matrix.SetElement(i, next_H_index++, 1);
          atoms_tmp.push_back(Atom("H"));
        }

        //atoms_tmp[i].SetNH(0);
        atoms_tmp[i].SetNumFixedHydrogens(0);
      }

      substructure_matrix.SetAtomVector(atoms_tmp);
      av = atoms_tmp;
      n = num_hyd;
      substructure_matrices.push_back(substructure_matrix);

    } else
      substructure_matrices.push_back(FragmentMatrix(F.GetMatrix().GetN(),
                                                     av,
                                                     F.GetMatrix().GetElements()));

    substructure_matrices.back().ResetAtomNeighbours();
  }
}

/*************
GetTER METHODS
*************/
const std::string& SubstructureCollection::GetCode() const {
  return code;
}

const Range& SubstructureCollection::GetRange() const {
  return r;
}

const std::vector<FragmentMatrix>& SubstructureCollection::GetSubstructureMatrices() const {
  return substructure_matrices;
}

bool SubstructureCollection::GetXOR() const {
  return XOR;
}

bool SubstructureCollection::GetAND() const {
  return AND;
}


/**********************
CREATE FAMILY FRAGMENTS
**********************/
void CreateSubstructures(AbstractSubstructureMap& abstract_substructures,
                         const std::list<std::string>& substructure_file_names) {
  if (!substructure_file_names.size()) {
    throw combi_ff::input_warning("no substructure input file");
    return;
  }

  std::cout << "creating prototype Substructures\n";

  for (auto && substructure_file_name : substructure_file_names) {
    XmlParserIn parser(substructure_file_name, XmlParserIn::read_all);
    const XmlElement& substructures_root = parser.GetTree().GetRoot();
    substructures_root.CheckTagName("substructures");
    substructures_root.CheckNumberOfChildren_atLeast(1);
    substructures_root.CheckAttribute("version");
    const std::string version = substructures_root.attributes.find("version")->second;

    if (version != combi_ff::current_version)
      std::cout << "?Warning: currently running combi_ff version " << combi_ff::current_version
                << " but substructure file " << substructure_file_name
                << " is version " << version << "\n";

    //read substrFile and add all the abstract fragments
    for (auto && substructure : substructures_root.children)
      GetNextAbstractSubstructure(abstract_substructures, substructure, version);
  }

  for (auto && s : abstract_substructures) {
    std::cout << s.second.GetCode() << std::endl;
    std::cout << s.second.GetMatrix() << std::endl;
  }

  std::cout << "******************************************************\n";
}

/***************************************
CREATE AND ADD AN AbstractSubstructure
***************************************/
void GetNextAbstractSubstructure(AbstractSubstructureMap& abstract_substructures,
                                 const XmlElement_ptr substructure,
                                 const std::string& version) {
  substructure->CheckTagName("substructure");
  substructure->CheckNumberOfChildren_equal(2);
  substructure->CheckAttributeSize(1);
  substructure->CheckAttribute("code");
  substructure->CheckNoValue();
  std::string code = substructure->attributes.find("code")->second;

  if (abstract_substructures.find(code) != abstract_substructures.end())
    throw combi_ff::input_error("Substructure code " + code + " occurs more than once\n");

  auto&& substructurePropertyIt = substructure->children.begin();
  const XmlElement_ptr num_atoms = *substructurePropertyIt;
  num_atoms->CheckTagName("num_atoms");
  num_atoms->CheckValue();
  num_atoms->CheckNumberOfChildren_equal(0);
  num_atoms->CheckAttributeSize(0);
  size_t n_atoms = std::stoi(num_atoms->value);

  if (!n_atoms)
    throw input_error("substructure must consist of at least one atom. Not the case for " + code +
                      "\n");

  FragmentMatrix M(n_atoms, std::vector<Atom>(n_atoms, Atom("*")));
  const XmlElement_ptr adjacency_stack = *(++substructurePropertyIt);
  adjacency_stack->CheckTagName("adjacency_stack");
  adjacency_stack->CheckValue();
  adjacency_stack->CheckNumberOfChildren_equal(0);
  adjacency_stack->CheckAttributeSize(0);
  AdjacencyVector values(0);
  std::istringstream adjacency_stackStream(adjacency_stack->value);
  int value;

  while (adjacency_stackStream >> value)
    values.push_back(value);

  size_t idx = 0;

  for (size_t i = 0; i < n_atoms - 1; i++) {
    for (size_t j = i + 1; j < n_atoms; j++) {
      if (idx >= values.size())
        throw combi_ff::input_error("stack too small for the given number of atoms");

      M.SetElement(i, j, values[idx++]);
    }
  }

  if (idx != values.size())
    throw combi_ff::input_error("stack too big for the given number of atoms");

  abstract_substructures.emplace(code, AbstractSubstructure(M, code, version));
  return;
}

} //namespace enu

} //namespace combi_ff