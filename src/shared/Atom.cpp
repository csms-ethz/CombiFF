// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "Atom.h"

#include <algorithm>

#include "ContainerOperators.h"
#include "exceptions.h"

namespace combi_ff {

Atom::Atom(combi_ff::ElementSymbol symbol, combi_ff::Connectivity degree)
    : symbol_united_atom(symbol),
      num_hydrogens(0),
      num_fixed_hydrogens(0),
      has_fixed_hydrogens(0),
      is_pseudoatom(0),
      neighbors(NeighborVector(0)) {
  properties.symbol = symbol;
  properties.degree = degree;
  properties.element_nr = 100;
  properties.priority = 100;
  properties.hydrogen_in_smiles = false;
  properties.mass = 0;
  neighbors.reserve(degree);
}

Atom::Atom(combi_ff::ElementSymbol name_, combi_ff::ElementSymbol symbol)
    : Atom(symbol) {
  symbol_united_atom = name_;
}

Atom::Atom(combi_ff::ElementIdentifier identifier_)
    : num_hydrogens(0),
      num_fixed_hydrogens(0),
      has_fixed_hydrogens(0),
      is_pseudoatom(0),
      neighbors(NeighborVector(0)) {
  int num_hydrogens_(-1);

  if (std::isdigit(identifier_.back())) {
    num_hydrogens_ = (std::stoi(std::string(1, identifier_.back())));

    if (*(identifier_.end() - 2) != 'H')
      throw input_error("united atom only possible with hydrogen, but got " +
                        identifier_ + "\n");

    identifier_ = identifier_.substr(0, identifier_.size() - 2);
  }

  auto element = element_types::element_property_map.find(identifier_);

  if (element == element_types::element_property_map.end())
    throw input_error("unknown element identifier (atom type) \"" +
                      identifier_ + "\".");

  properties = element->second;
  symbol_united_atom = properties.symbol;

  if (num_hydrogens_ != -1) SetNumFixedHydrogens((size_t)num_hydrogens_);

  neighbors.reserve(properties.degree);
}

CnvAtom::CnvAtom(combi_ff::ElementIdentifier identifier_)
    : Atom(identifier_), formal_charge("") {}

CnvAtom::CnvAtom(combi_ff::ElementSymbol name_, combi_ff::ElementSymbol symbol)
    : Atom(name_, symbol), formal_charge("") {}

Atom::Atom(const Atom& a)
    : properties(a.properties),
      symbol_united_atom(a.GetUnitedAtomSymbol()),
      num_hydrogens(a.GetNumHydrogens()),
      num_fixed_hydrogens(a.GetNumFixedHydrogens()),
      has_fixed_hydrogens(a.GetHasFixedHydrogens()),
      is_pseudoatom(a.is_pseudoatom),
      neighbors(a.neighbors) {}

CnvAtom::CnvAtom(const CnvAtom& a) : Atom(a), formal_charge(a.formal_charge) {}

bool Atom::operator==(const Atom& b) const {
  if (GetElementSymbol() == b.GetElementSymbol() &&
      GetDegree() == b.GetDegree() &&
      GetNumConnections() == b.GetNumConnections() &&
      GetNumTotalHydrogens() == b.GetNumTotalHydrogens())
    return true;

  else
    return false;
}

bool Atom::operator!=(const Atom& b) const { return !(*this == b); }

bool Atom::operator<(const Atom& b) const {
  if (GetDegree() < b.GetDegree())
    return false;

  else if (GetDegree() > b.GetDegree())
    return true;

  else if (GetNumConnections() > b.GetNumConnections())
    return false;

  else if (GetNumConnections() < b.GetNumConnections())
    return true;

  if (GetElementSymbol() == "H") return false;

  if (b.GetElementSymbol() == "H") return true;

  // return(a.GetTypeName() < b.GetTypeName());
  return (GetElementNumber() < b.GetElementNumber());
}

bool Atom::operator>(const Atom& b) const {
  if (GetDegree() < b.GetDegree())
    return true;

  else if (GetDegree() > b.GetDegree())
    return false;

  else if (GetNumConnections() >
           b.GetNumConnections())  // higher degree and lower num connections
                                   // means that the bond degree of the
                                   // connections is higher. could still be e.g.
                                   // 1 + 3 or 2 + 2, but this is easy to
                                   // evaluate, i.e. good balance between
                                   // speedup for smaller partitions and cost to
                                   // evaluate criterion
    return true;

  else if (GetNumConnections() < b.GetNumConnections())
    return false;

  if (GetElementSymbol() == "H") return true;

  if (b.GetElementSymbol() == "H") return false;

  // return(a.GetTypeName() > b.GetTypeName());
  return (GetElementNumber() > b.GetElementNumber());
}

const combi_ff::ElementSymbol& Atom::GetElementSymbol() const {
  return properties.symbol;
}

const combi_ff::ElementSymbol& Atom::GetUnitedAtomSymbol() const {
  return symbol_united_atom;
}

combi_ff::Connectivity Atom::GetDegree() const { return properties.degree; }

combi_ff::ElementNumber Atom::GetElementNumber() const {
  return properties.element_nr;
}

void Atom::AddNeighbour(size_t i) {
  /*auto&& range = std::equal_range(neighbors.begin(), neighbors.end(), i);

  if(range.first == range.second)
    neighbors.insert(range.first, i);*/

  auto el = std::upper_bound(neighbors.begin(), neighbors.end(), i);
  if (el == neighbors.begin() || *(el - 1) != i) neighbors.insert(el, i);
}

void Atom::RemoveNeighbour(size_t i) {
  auto el = std::find(neighbors.begin(), neighbors.end(), i);
  if (el != neighbors.end()) neighbors.erase(el);
}

void Atom::SetNeighbours(NeighborVector n) { neighbors = n; }

void Atom::EraseNeighbours() { neighbors.clear(); }

double Atom::GetMass() const {
  return properties.mass + (double)GetNumTotalHydrogens() *
                               combi_ff::element_types::hydrogen_mass;
}

void Atom::SetNumFixedHydrogens(size_t i) {
  // assert(nH == 0);
  properties.degree += num_fixed_hydrogens;
  properties.degree -= i;
  num_fixed_hydrogens = i;
  has_fixed_hydrogens = true;

  if (symbol_united_atom.size() > 2 &&
      symbol_united_atom[symbol_united_atom.size() - 2] == 'H' &&
      symbol_united_atom.back() != '>') {
    symbol_united_atom.erase(symbol_united_atom.size() - 1);
    symbol_united_atom += std::to_string(i);

  } else
    symbol_united_atom += "H" + std::to_string(i);
}

void Atom::SetHasFixedHydrogens(bool b) { has_fixed_hydrogens = b; }

bool Atom::GetHasFixedHydrogens() const { return has_fixed_hydrogens; }

const std::string Atom::GetFormalCharge() const { return ""; }

const std::string& CnvAtom::GetFormalCharge() const { return formal_charge; }

const combi_ff::NeighborVector& Atom::GetNeighbours() const {
  return neighbors;
}

size_t Atom::GetNumNeighbors() const {
  return neighbors.size() + GetNumTotalHydrogens();
}

size_t Atom::GetNumConnections() const { return neighbors.size(); }

size_t Atom::GetNumFixedHydrogens() const { return num_fixed_hydrogens; }

size_t Atom::GetNumTotalHydrogens() const {
  return num_fixed_hydrogens + num_hydrogens;
}

Connectivity Atom::GetNumHydrogens() const { return num_hydrogens; }

bool Atom::GetHydrogenInSmiles() const { return properties.hydrogen_in_smiles; }

void Atom::SetHydrogensInSmiles(bool b) { properties.hydrogen_in_smiles = b; }

void CnvAtom::SetFormalCharge(const std::string& f) {
  if (formal_charge.size()) {
    if (formal_charge.front() == '+') {
      if (std::isdigit(formal_charge.back()))
        SetDegree((size_t)((int)GetDegree() -
                           std::stoi(std::string(1, formal_charge.back()))));

      else
        SetDegree(GetDegree() - 1);

    } else {  // formal charge negative
      if (std::isdigit(formal_charge.back()))
        SetDegree((size_t)((int)GetDegree() +
                           std::stoi(std::string(1, formal_charge.back()))));

      else
        SetDegree(GetDegree() + 1);
    }
  }

  formal_charge = f;

  if (formal_charge.size()) {
    if (formal_charge.front() == '+') {
      if (std::isdigit(formal_charge.back()))
        SetDegree((size_t)((int)GetDegree() +
                           std::stoi(std::string(1, formal_charge.back()))));

      else
        SetDegree(GetDegree() + 1);

    } else {  // formal charge negative
      if (std::isdigit(formal_charge.back()))
        SetDegree((size_t)((int)GetDegree() -
                           std::stoi(std::string(1, formal_charge.back()))));

      else
        SetDegree(GetDegree() - 1);
    }
  }
}

combi_ff::ElementPriority Atom::GetElementPriority() const {
  return properties.priority;
}

void Atom::SetUnitedAtomSymbol(std::string s) { symbol_united_atom = s; }

void Atom::SetDegree(size_t i) { properties.degree = i; }

bool Atom::IsPseudoatom() const { return is_pseudoatom; }

void Atom::SetPseudoatom(bool b) { is_pseudoatom = b; }

void Atom::SetNumHydrogens(size_t i) {
  properties.degree += num_hydrogens;
  properties.degree -= i;
  num_hydrogens = i;
  // assert(fixedNH == 0 && !hasFixedNH);

  if (symbol_united_atom.size() > 2 &&
      symbol_united_atom[symbol_united_atom.size() - 2] == 'H' &&
      symbol_united_atom.back() != '>') {
    symbol_united_atom.erase(symbol_united_atom.size() - 1);
    symbol_united_atom += std::to_string(i);

  } else
    symbol_united_atom += "H" + std::to_string(i);
}

void Atom::RemoveHydrogens() {
  num_hydrogens = 0;
  num_fixed_hydrogens = 0;
  symbol_united_atom = properties.symbol;
}

std::ostream& operator<<(std::ostream& stream, const Atom& a) {
  // stream << "Atom " << a.GetUnitedAtomSymbol() << " (connectivity " <<
  // a.GetDegree() << ")\n";
  stream << a.GetUnitedAtomSymbol();
  return stream;
}

/*std::ostream& operator<<(std::ostream& stream, const AtomVector& atoms) {
  for(auto && a : atoms)
    stream << a.GetUnitedAtomSymbol() << " ";

  return stream;
}*/

void SortAtoms(AtomVector<Atom>& atoms,
               std::list<combi_ff::LambdaVector>& lambda) {
  for (size_t i = 0; i < atoms.size(); i++) {
    for (size_t j = i + 1; j < atoms.size(); j++) {
      if (atoms[j].GetDegree() > atoms[i].GetDegree()) {
        std::swap(atoms[i], atoms[j]);

        for (auto&& l : lambda) std::swap(l[i], l[j]);

      } else if (atoms[j].GetDegree() == atoms[i].GetDegree()) {
        if (atoms[j].GetElementSymbol() != "H" &&
            atoms[j].GetElementSymbol() < atoms[i].GetElementSymbol()) {
          std::swap(atoms[i], atoms[j]);

          for (auto&& l : lambda) std::swap(l[i], l[j]);

        } else if (atoms[i].GetElementSymbol() == "H" &&
                   atoms[j].GetElementSymbol() != "H") {
          std::swap(atoms[i], atoms[j]);

          for (auto&& l : lambda) std::swap(l[i], l[j]);
        }
      }
    }
  }
}

void SortAtoms(AtomVector<Atom>& atoms,
               std::list<combi_ff::LambdaVector>& lambda,
               std::vector<int>& XORidx) {
  for (size_t i = 0; i < atoms.size(); i++) {
    for (size_t j = i + 1; j < atoms.size(); j++) {
      if (atoms[j].GetDegree() > atoms[i].GetDegree()) {
        std::swap(atoms[i], atoms[j]);
        std::swap(XORidx[i], XORidx[j]);

        for (auto&& l : lambda) std::swap(l[i], l[j]);

      } else if (atoms[j].GetDegree() == atoms[i].GetDegree()) {
        if (atoms[j].GetElementSymbol() != "H" &&
            atoms[j].GetElementSymbol() < atoms[i].GetElementSymbol()) {
          std::swap(atoms[i], atoms[j]);
          std::swap(XORidx[i], XORidx[j]);

          for (auto&& l : lambda) std::swap(l[i], l[j]);

        } else if (atoms[i].GetElementSymbol() == "H" &&
                   atoms[j].GetElementSymbol() != "H") {
          std::swap(atoms[i], atoms[j]);
          std::swap(XORidx[i], XORidx[j]);

          for (auto&& l : lambda) std::swap(l[i], l[j]);
        }
      }
    }
  }
}

const LambdaVector Sort(AtomVector<Atom>& atoms) {
  // sort atoms by bond capacity (i.e. less hydrogens) s.t. atoms with the same
  // number of NH and fixedNH are next to each other, as these should be
  // indistinguishable in the enumeration
  for (size_t i = 0; i < atoms.size(); i++) {
    for (size_t j = i + 1; j < atoms.size(); j++) {
      if (atoms[j].GetDegree() > atoms[i].GetDegree())
        std::swap(atoms[i], atoms[j]);

      else if (atoms[j].GetDegree() == atoms[i].GetDegree() &&
               atoms[j] > atoms[i])
        std::swap(atoms[i], atoms[j]);
    }
  }

  combi_ff::LambdaVector lambda(atoms.size(), 1);
  std::vector<AtomVector<Atom>> atom_blocks(1, AtomVector<Atom>(1, atoms[0]));
  atom_blocks.reserve(atoms.size());
  size_t idx(0);

  for (size_t i = 1; i < atoms.size(); i++) {
    if (atoms[i - 1] != atoms[i]) {
      atom_blocks.push_back(AtomVector<Atom>(1, atoms[i]));
      atom_blocks.back().reserve(atoms.size() - i);
      idx++;

    } else {
      atom_blocks.back().push_back(atoms[i]);
      lambda[idx]++;
    }
  }

  lambda.resize(idx + 1);

  for (size_t i = 0; i < lambda.size(); i++) {
    for (size_t j = i + 1; j < lambda.size(); j++) {
      if (atom_blocks[i].front().GetDegree() <
          atom_blocks[j].front().GetDegree()) {
        std::swap(atom_blocks[i], atom_blocks[j]);
        std::swap(lambda[i], lambda[j]);

      } else if (atom_blocks[i].front().GetDegree() ==
                     atom_blocks[j].front().GetDegree() &&
                 lambda[i] > lambda[j]) {
        std::swap(atom_blocks[i], atom_blocks[j]);
        std::swap(lambda[i], lambda[j]);

      } else if (atom_blocks[i].front().GetDegree() ==
                     atom_blocks[j].front().GetDegree() &&
                 lambda[i] == lambda[j] &&
                 atom_blocks[i].front() > atom_blocks[j].front()) {
        std::swap(atom_blocks[i], atom_blocks[j]);
        std::swap(lambda[i], lambda[j]);
      }
    }
  }

  atoms = AtomVector<Atom>(0);

  for (auto&& ab : atom_blocks) atoms.insert(atoms.end(), ab.begin(), ab.end());

  return lambda;
}

void PermuteVector(AtomVector<Atom>& atoms, const Permutations& permutations) {
  for (auto&& p : permutations) {
    // idx.swap(p.first, p.second);
    std::swap(atoms[p.first], atoms[p.second]);
  }
}

}  // namespace combi_ff
