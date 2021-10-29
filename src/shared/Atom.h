#ifndef ATOM_H
#define ATOM_H

#include <cassert>
#include "Permutation.h"
#include <vector>
#include <list>
#include "LambdaVector.h"
#include "ElementTypes.h"
#include "StringVector.h"


namespace combi_ff {

typedef std::vector<combi_ff::Connectivity> ConnectivityVec;
typedef ConnectivityVec NeighborVector;

class Atom {
 public:
  /* constructors */
  Atom() = default;
  Atom(combi_ff::ElementSymbol symbol, combi_ff::Connectivity degree);
  Atom(combi_ff::ElementSymbol symbol);
  Atom(combi_ff::ElementSymbol name_, combi_ff::ElementSymbol symbol);
  Atom(const Atom& a);

  /* operators */
  bool operator==(const Atom& b) const;
  bool operator!=(const Atom& b) const;
  bool operator<(const Atom& b) const;
  bool operator>(const Atom& b) const;

  /* getter functions */
  const combi_ff::ElementSymbol& GetUnitedAtomSymbol() const;
  const combi_ff::ElementSymbol& GetElementSymbol() const;
  combi_ff::Connectivity GetDegree() const;
  combi_ff::ElementNumber GetElementNumber() const;
  const combi_ff::NeighborVector& GetNeighbours() const ;
  size_t GetNumNeighbors() const ;
  combi_ff::ElementPriority GetElementPriority() const;
  combi_ff::Connectivity GetNumConnections() const ;
  combi_ff::Connectivity GetNumFixedHydrogens() const ;
  combi_ff::Connectivity GetNumTotalHydrogens() const ;
  combi_ff::Connectivity GetNumHydrogens() const ;
  bool GetHasFixedHydrogens() const;
  bool IsPseudoatom() const;
  double GetMass() const ;
  bool GetHydrogenInSmiles() const;
  const std::string GetFormalCharge() const;

  /* setter functions */
  void SetNumFixedHydrogens(size_t i);
  void SetUnitedAtomSymbol(std::string s);
  void SetDegree(size_t i) ;
  void SetNumHydrogens(size_t i);
  void SetHasFixedHydrogens(bool b);
  void SetPseudoatom(bool b);
  void SetHydrogensInSmiles(bool b);

  /* other member functions */
  void AddNeighbour(size_t i);
  void RemoveNeighbour(size_t i) ;
  void SetNeighbours(combi_ff::NeighborVector n);
  void EraseNeighbours();
  void RemoveHydrogens() ;

 private:
  combi_ff::element_types::ElementProperties properties;
  combi_ff::ElementSymbol symbol_united_atom {""};
  combi_ff::Connectivity num_hydrogens {0};
  combi_ff::Connectivity num_fixed_hydrogens {0};
  bool has_fixed_hydrogens {false};
  bool is_pseudoatom {false};
  combi_ff::NeighborVector neighbors {NeighborVector(0)};
};

class CnvAtom final : public Atom {

 public:
  CnvAtom() = default;
  CnvAtom(combi_ff::ElementSymbol symbol);
  CnvAtom(combi_ff::ElementSymbol name_, combi_ff::ElementSymbol symbol);
  CnvAtom(const CnvAtom& a);
  void SetFormalCharge(const std::string& f);
  const std::string& GetFormalCharge() const;


 private:
  std::string formal_charge {""};

};

std::ostream& operator<<(std::ostream& stream, const Atom& a) ;

template<typename AtomClass>
using AtomVector = std::vector<AtomClass>;
typedef std::vector<combi_ff::Connectivity> ConnectivityVector;

void SortAtoms(AtomVector<Atom>& atoms, std::list<combi_ff::LambdaVector>& lambda) ;
void SortAtoms(AtomVector<Atom>& atoms, std::list<combi_ff::LambdaVector>& lambda,
               std::vector<int>& XORidx) ;
const combi_ff::LambdaVector Sort(AtomVector<Atom>& atoms) ;

template<typename AtomClass>
std::string CreateCanonicalFormulaFromAtomVector(const AtomVector<AtomClass>& atoms) {
  AtomVector<AtomClass> atoms_full(atoms);

  for (size_t i = 0; i < atoms_full.size(); i++) {
    if (atoms_full[i].GetNumHydrogens()) {
      for (uint j = 0; j < atoms_full[i].GetNumHydrogens(); j++)
        atoms_full.push_back(AtomClass("H"));

      atoms_full[i].SetNumHydrogens(0);
    }

    if (atoms_full[i].GetNumFixedHydrogens()) {
      for (uint j = 0; j < atoms_full[i].GetNumFixedHydrogens(); j++)
        atoms_full.push_back(AtomClass("H"));

      atoms_full[i].SetNumFixedHydrogens(0);
    }
  }

  for (size_t i = 0; i < atoms_full.size(); i++) {
    for (size_t j = i + 1; j < atoms_full.size(); j++) {
      if (atoms_full[i].GetElementPriority() > atoms_full[j].GetElementPriority())
        std::swap(atoms_full[i], atoms_full[j]);
    }
  }

  size_t nH(0);
  combi_ff::LambdaVector lambda_types(atoms_full.size(), 1);
  combi_ff::StringVector atom_types(atoms_full.size());

  for (size_t i = 0; i < atoms_full.size(); i++) {
    atom_types[i] = (atoms_full[i].GetElementSymbol());
    lambda_types[i] = 1;
    nH += atoms_full[i].GetNumFixedHydrogens();

    for (size_t j = i + 1; j < atoms_full.size(); j++) {
      if (atoms_full[i].GetElementSymbol() == atoms_full[j].GetElementSymbol()) {
        lambda_types[i]++;
        nH += atoms_full[j].GetNumFixedHydrogens();
        atoms_full.erase(atoms_full.begin() + j);
        j--;
      }
    }
  }

  lambda_types.resize(atoms_full.size());
  atom_types.resize(atoms_full.size());

  if (nH) {
    //case that formula only consists of hydrogen atoms
    if (atom_types[0] == "H")
      lambda_types[0] += nH;

    //case that hydrogen is the atom type with highest priority
    else if (atom_types.size() > 1 && atom_types[1] == "H")
      lambda_types[1] += nH;

    //case that hydrogen is atom type with 2nd highest priority
    //(i.e. there is at least one carbon atom in the formula as carbon is the only atom with a higher priority than hydrogen)
    else {
      if (atom_types[0] == "C") {
        atom_types.insert(atom_types.begin() + 1, "H");
        lambda_types.insert(lambda_types.begin() + 1, nH);

      } else {
        atom_types.insert(atom_types.begin(), "H");
        lambda_types.insert(lambda_types.begin(), nH);
      }
    }
  }

  std::string formula("");

  for (size_t i = 0; i < atom_types.size(); i++)
    formula += atom_types[i] + std::to_string(lambda_types[i]);

  return formula;
}

void PermuteVector(AtomVector<Atom>& atoms,
                   const combi_ff::Permutations& permutations);

} //namespace combi_ff


#endif
