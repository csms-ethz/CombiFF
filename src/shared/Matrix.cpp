#include "Matrix.h"

namespace combi_ff {

/*
 ComparisonMatrix class members
*/

ComparisonMatrix::ComparisonMatrix(size_t N, size_t M) 
: Matrix<bool>(N, M, false) {}

void ComparisonMatrix::Print() const {
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < M; j++)
      std::cout << GetElement(i, j) << " ";

    std::cout << std::endl;
  }
}

size_t ComparisonMatrix::AccumulateColumn(const size_t col) const {
  size_t sum(0);
  size_t lim = N * M;

  for (size_t i = col; i < lim; i += M)
    sum += elements[i];

  return sum;
}

size_t ComparisonMatrix::AccumulateRow(const size_t row) const {
  size_t sum(0);
  size_t lim = row * M + M;

  for (size_t j = row * M; j < lim; j++)
    sum += elements[j];

  return sum;
}

/*
 FragmentMatrix class members
*/
FragmentMatrix::FragmentMatrix() : SymmetricalMatrix(0), atoms(0) {}

FragmentMatrix::FragmentMatrix(size_t N)
  : SymmetricalMatrix(N), atoms(combi_ff::AtomVector<combi_ff::Atom>(N)) {}

FragmentMatrix::FragmentMatrix(size_t N,
                               combi_ff::AtomVector<combi_ff::Atom> atoms)
  : SymmetricalMatrix(N), atoms(atoms) {}

FragmentMatrix::FragmentMatrix(size_t N, AdjacencyVector v)
  : SymmetricalMatrix<size_t>(N), atoms(N, combi_ff::Atom("*")) {
  for (size_t i = 0; i < N_minus_one; i++) {
    for (size_t j = i + 1; j < N; j++) {
      if (v[i * N + j])
        SetElement(i, j, v[i * N + j]); //to Set neighbours correctly
    }
  }
}

FragmentMatrix::FragmentMatrix(size_t N,
                               combi_ff::AtomVector<combi_ff::Atom> atoms,
                               AdjacencyVector v) : SymmetricalMatrix<size_t>(N), atoms(atoms) {
  for (size_t i = 0; i < N_minus_one; i++) {
    for (size_t j = i + 1; j < N; j++) {
      if (v[i * N + j])
        SetElement(i, j, v[i * N + j]); //to Set neighbours correctly
    }
  }
}

FragmentMatrix::FragmentMatrix(SymmetricalMatrix<size_t>& M,
                               combi_ff::AtomVector<combi_ff::Atom>& atoms) : SymmetricalMatrix<size_t>
  (M.GetN()),
  atoms(atoms) {
  for (size_t i = 0; i < N_minus_one; i++) {
    for (size_t j = i + 1; j < N; j++) {
      if (M.GetElements()[i * N + j ])
        SetElement(i, j, M.GetElements()[i * N + j]); //to Set neighbours correctly
    }
  }
}
void FragmentMatrix::SetElement(const size_t i, const size_t j,
                                const size_t v) {
  if (v != 0) {
    atoms[i].AddNeighbour(j);
    atoms[j].AddNeighbour(i);

  } else {
    if (GetElement(i, j)) {
      atoms[i].RemoveNeighbour(j);
      atoms[j].RemoveNeighbour(i);
    }
  }

  SymmetricalMatrix<size_t>::SetElement(i, j, v);
}

void FragmentMatrix::SetAtomVector(const combi_ff::AtomVector<Atom>& a) {
  atoms.assign(a.begin(), a.end());

  for (size_t i = 0; i < N_minus_one; i++) {
    for (size_t j = i + 1; j < N; j++) {
      if (GetElement(i, j)) {
        atoms[i].AddNeighbour(j);
        atoms[j].AddNeighbour(i);
      }
    }
  }
}

const combi_ff::AtomVector<combi_ff::Atom>& FragmentMatrix::GetAtomVector()
const {
  return atoms;
}

void FragmentMatrix::print() const {
  std::cout << "    ";

  for (size_t i = 0; i < atoms.size(); i++)
    std::cout << atoms[i].GetUnitedAtomSymbol() << " ";

  std::cout << std::endl;

  for (size_t i = 0; i < N; i++) {
    std::cout << std::setw(3) << atoms[i].GetUnitedAtomSymbol() << " ";

    for (size_t j = 0; j < N; j++)
      std::cout << std::setw((int)atoms[j].GetUnitedAtomSymbol().size()) <<
                GetElement(i, j) << " ";

    std::cout << std::endl;
  }
}

void FragmentMatrix::PrintIndented() const {
  std::cout << '\t';

  for (size_t i = 0; i < N; i++)
    std::cout << atoms[i].GetUnitedAtomSymbol() << " ";

  std::cout << '\n';

  for (size_t i = 0; i < N; i++) {
    std::cout << '\t';

    for (size_t j = 0; j < N; j++)
      std::cout << Matrix<size_t>::GetElement(i, j) << " ";

    std::cout << std::endl;
  }
}

void FragmentMatrix::ResetAtomNeighbours() {
  for (auto && a : atoms)
    a.EraseNeighbours();

  for (size_t i = 0; i < N_minus_one; i++) {
    for (size_t j = i + 1; j < N; j++) {
      if (GetElement(i, j)) {
        atoms[i].AddNeighbour(j);
        atoms[j].AddNeighbour(i);
      }
    }
  }
}

std::string FragmentMatrix::GetPrintIndented() const {
  std::string s("\t\t\t");

  for (size_t i = 0; i < N; i++)
    s += atoms[i].GetUnitedAtomSymbol() + " ";

  s +=  '\n';

  for (size_t i = 0; i < N; i++) {
    s += "\t\t\t";

    for (size_t j = 0; j < N; j++)
      s += std::to_string(Matrix<size_t>::GetElement(i, j)) + " ";

    s += "\n";
  }

  return s;
}

void FragmentMatrix::GetNumMultipleBonds(size_t& n_single, size_t& n_double,
                                         size_t& n_triple,
                                         size_t& n_quadruple) const {
  for (size_t i = 0; i < N_minus_one; i++) {
    n_single += atoms[i].GetNumHydrogens() + atoms[i].GetNumFixedHydrogens();

    for (size_t j = i + 1; j < N; j++) {
      size_t bond = GetElement(i, j);

      if (bond == 1)
        n_single++;

      else if (bond == 2)
        n_double++;

      else if (bond == 3)
        n_triple++;

      else if (bond == 4)
        n_quadruple++;
    }
  }
}

bool FragmentMatrix::HasCycle() const {
  for (size_t i = 0; i < N; i++) {
    const combi_ff::NeighborVector& neighbors = atoms[i].GetNeighbours();

    for (auto neighbor : neighbors) {
      std::vector<bool> visited_(N, false);
      SymmetricalMatrix<bool> visited(N);
      visited_[neighbor] = true;
      visited.SetElement(i, neighbor, true);
      std::stack<size_t> stack;
      stack.push(neighbor);

      while (!stack.empty()) {
        size_t curr_element = stack.top();
        stack.pop();

        for (size_t ii = 0; ii < N; ii++) {
          if (!visited.GetElement(curr_element, ii) &&
              GetElement(curr_element, ii) != 0) {
            if (ii == i)
              return true;

            visited.SetElement(curr_element, ii, true);
            visited_[ii] = true;
            stack.push(ii);
          }
        }
      }
    }
  }

  return false;
}

} // namespace combi_ff