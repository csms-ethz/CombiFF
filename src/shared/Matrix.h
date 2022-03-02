// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef Matrix_H
#define Matrix_H

#include <algorithm>
#include <iomanip>
#include <numeric>
#include <stack>

#include "Atom.h"
#include "ContainerOperators.h"
#include "Permutation.h"
#include "Range.h"
#include "exceptions.h"

namespace combi_ff {

typedef std::vector<size_t> DegreeVector;
typedef std::vector<size_t> TypeVector;
typedef std::vector<size_t> AdjacencyVector;

/**** MATRIX CLASS ****/
template <typename T>
class Matrix {
 public:
  Matrix(size_t N, size_t M, T init);
  Matrix(size_t N, size_t M, std::vector<T> v);

  T GetElement(const size_t i, const size_t j) const;
  const std::vector<T>& GetElements() const;
  size_t GetN() const;
  size_t GetM() const;

  void SetElement(const size_t i, const size_t j, const T value);
  void SetElements(const std::vector<T> v);
  void Print() const;

  T AccumulateColumn(const size_t col) const;
  T AccumulateRow(const size_t row) const;
  T AccumulateColumn(const size_t imin, const size_t imax,
                     const size_t j) const;
  T AccumulateRow(const size_t i, const size_t jmin, const size_t jmax) const;

 protected:
  size_t N, M;  // matrix dimensions. N = number of columns, M = number of rows
  std::vector<T> elements;  // store matrix elements in 1D vector
};

/**** COMPARISON MATRIX CLASS ****/
class ComparisonMatrix final : public Matrix<bool> {
 public:
  ComparisonMatrix(size_t N, size_t M);

  void Print() const;
  size_t AccumulateColumn(const size_t col) const;
  size_t AccumulateRow(const size_t row) const;
};

/**** SYMMETRICAL MATRIX CLASS ****/
template <typename T>
class SymmetricalMatrix : public Matrix<T> {
 public:
  SymmetricalMatrix();
  SymmetricalMatrix(size_t N);
  SymmetricalMatrix(size_t N, std::vector<T> v);

  void SetElement(const size_t i, const size_t j, const T v);
  void FillColumn(const std::vector<T> d, const size_t j);
  void FillRow(const std::vector<T> d, const size_t i);

  std::string GetStack() const;

  void Permute(const combi_ff::Permutations& permutations);
  void Permute(const combi_ff::Permutation permutation);

  void PrintIndented() const;

  bool HasCycle() const;
  bool IsConnected() const;

  bool operator==(const SymmetricalMatrix& B) const;
  bool operator>(const SymmetricalMatrix& B) const;
  bool operator<(const SymmetricalMatrix& B) const;

 protected:
  size_t N_minus_one;
};

/**** FRAGMENT MATRIX CLASS ****/
class FragmentMatrix final : public SymmetricalMatrix<size_t> {
 public:
  FragmentMatrix();
  FragmentMatrix(size_t N);
  FragmentMatrix(size_t N, combi_ff::AtomVector<combi_ff::Atom> atoms);
  FragmentMatrix(size_t N, AdjacencyVector v);
  FragmentMatrix(size_t N, combi_ff::AtomVector<combi_ff::Atom> atoms,
                 AdjacencyVector v);
  FragmentMatrix(SymmetricalMatrix<size_t>& M,
                 combi_ff::AtomVector<combi_ff::Atom>& atoms);

  void SetElement(const size_t i, const size_t j, const size_t v);
  void SetAtomVector(const combi_ff::AtomVector<Atom>& a);

  const combi_ff::AtomVector<combi_ff::Atom>& GetAtomVector() const;

  void print() const;
  void PrintIndented() const;
  std::string GetPrintIndented() const;

  void ResetAtomNeighbors();

  void GetNumMultipleBonds(size_t& n_single, size_t& n_double, size_t& n_triple,
                           size_t& n_quadruple) const;
  bool HasCycle() const;

 private:
  combi_ff::AtomVector<combi_ff::Atom> atoms;
};

/**** ADJACENCY MATRIX CLASS ****/
template <typename T, typename AtomClass>
class AdjacencyMatrix final : public SymmetricalMatrix<T> {
 public:
  AdjacencyMatrix();
  AdjacencyMatrix(size_t N);
  AdjacencyMatrix(size_t N, combi_ff::LambdaVector& lambda_,
                  combi_ff::AtomVector<AtomClass>& atoms_);
  AdjacencyMatrix(const AdjacencyMatrix& B);

  void SetElements(std::vector<T> v);
  void SetElement(const size_t i, const size_t j, const T v);
  void SetLambda(combi_ff::LambdaVector& l);
  void SetAtomVector(const combi_ff::AtomVector<AtomClass>& a);
  void SetIsAromaticCarbon(const std::vector<bool>& a);
  AdjacencyMatrix<T, AtomClass> ExtendToFullMatrix() const;

  double GetMass() const;
  const combi_ff::AtomVector<AtomClass>& GetAtomVector() const;
  combi_ff::AtomVector<AtomClass>& GetAtomVectorNonConst();
  AtomClass& GetAtom(size_t i);
  const combi_ff::LambdaVector& GetLambda() const;
  const TypeVector& GetTypes() const;
  const std::vector<size_t>& GetIndices() const;
  size_t GetTypeNr(const size_t i) const;
  size_t GetMaxRow(size_t i) const;
  size_t GetMinRow(size_t i) const;
  size_t GetNumDoubleBonds() const;
  void GetNumMultipleBonds(size_t& n_single, size_t& n_double, size_t& n_triple,
                           size_t& n_quadruple, size_t& n_aromatic) const;
  void GetNumMultipleBonds(std::vector<size_t>& ranges) const;
  std::vector<bool>& GetIsAromaticCarbon();
  const std::vector<bool>& GetIsAromaticCarbonConst() const;
  bool IsLastRowOfALambdaBlock(size_t i);

  void Print() const;
  void PrintToFile(std::ostream& outputFile) const;

  void ResetAtomNeighbors();

  bool hasCycle() const;
  bool IsConnected() const;
  std::vector<bool> IsInCycle() const;
  bool AreInSameCycle(size_t atom1, size_t atom2) const;
  size_t GetNumPiAtoms() const;

  bool SmallerThanPermuted(const size_t curr_row,
                           const std::vector<size_t>& permuted_indices,
                           size_t& i_pos, size_t& j_pos, size_t& i_pos_perm,
                           size_t& j_pos_perm, size_t& i0, size_t& j0,
                           bool& diff, const size_t smallest_diff_index,
                           const size_t highest_permuted_index);
  bool EqualTo(const std::vector<size_t>& permuted_indices, size_t& j_pos);

  void MakeCanonical(const combi_ff::RepresentationSystem& u);
  void MakeCanonical(const combi_ff::RepresentationSystem& u,
                     std::vector<size_t>& index, const size_t limit,
                     Permutations& matrix_permutations);

  void SortAtomVector();
  std::vector<size_t> SortAtomVector(std::vector<bool>& is_aromatic);
  Permutations SortAtomVector_(std::vector<bool>& is_aromatic);

 private:
  combi_ff::LambdaVector lambda;
  TypeVector
      types;  // e.g. for atoms = [C, C, C, H, H] => types = [0, 0, 0, 1, 1]
  combi_ff::AtomVector<AtomClass> atoms;
  std::vector<size_t> min_row, max_row;
  std::vector<size_t> indices;
  std::vector<bool> is_aromatic_carbon;
  std::vector<bool> lambda_block_end;
};

/**** upper triangular matrix class ****/
template <typename T>
class TriangularMatrix final : public Matrix<T> {
 public:
  TriangularMatrix() : Matrix<T>(0, 0, 0) {}
  TriangularMatrix(size_t N) : Matrix<T>(N, N, 0) {}
  void Print() const {
    for (size_t i = 0; i < this->N; i++) {
      for (size_t j = 0; j < i; j++) std::cout << "  ";

      for (size_t j = i; j < this->N; j++)
        std::cout << this->GetElement(i, j) << " ";

      std::cout << std::endl;
    }
  }
};

/*
Matrix class member functions
*/
template <typename T>
Matrix<T>::Matrix(size_t N, size_t M, T init)
    : N(N), M(M), elements(std::vector<T>(N * M, init)) {}

template <typename T>
Matrix<T>::Matrix(size_t N, size_t M, std::vector<T> v)
    : N(N), M(M), elements(v) {}

template <typename T>
T Matrix<T>::GetElement(const size_t i, const size_t j) const {
  return elements[i * M + j];
}
template <typename T>
const std::vector<T>& Matrix<T>::GetElements() const {
  return elements;
}
template <typename T>
size_t Matrix<T>::GetN() const {
  return N;
}
template <typename T>
size_t Matrix<T>::GetM() const {
  return M;
}
template <typename T>
void Matrix<T>::SetElement(const size_t i, const size_t j, const T value) {
  elements[i * M + j] = value;
}
template <typename T>
void Matrix<T>::SetElements(const std::vector<T> v) {
  elements = v;
}
template <typename T>
void Matrix<T>::Print() const {
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < M; j++) std::cout << GetElement(i, j) << " ";

    std::cout << std::endl;
  }
}
template <typename T>
T Matrix<T>::AccumulateColumn(const size_t col) const {
  T sum(0);
  size_t lim = N * M;

  for (size_t i = col; i < lim; i += M) sum += elements[i];

  return sum;
}
template <typename T>
T Matrix<T>::AccumulateRow(const size_t row) const {
  T sum(0);
  size_t lim = row * M + M;

  for (size_t j = row * M; j < lim; j++) sum += elements[j];

  return sum;
}
template <typename T>
T Matrix<T>::AccumulateColumn(const size_t imin, const size_t imax,
                              const size_t j) const {
  T sum(0);
  size_t lim = imax * M + j;

  for (size_t i = imin * M + j; i < lim; i += M) sum += elements[i];

  return sum;
}
template <typename T>
T Matrix<T>::AccumulateRow(const size_t i, const size_t jmin,
                           const size_t jmax) const {
  T sum(0);
  size_t lim = i * M + jmax;

  for (size_t j = i * M + jmin; j < lim; j++) sum += elements[j];

  return sum;
  // return std::accumulate(A.GetElements().begin() + i * A.GetN() + jmin,
  // A.GetElements().begin() + i * A.GetN() + jmax ,0);
}

/*
SymmetricalMatrix class member functions
*/
template <typename T>
SymmetricalMatrix<T>::SymmetricalMatrix()
    : Matrix<T>(0, 0, T(0)), N_minus_one(0) {}

template <typename T>
SymmetricalMatrix<T>::SymmetricalMatrix(size_t N)
    : Matrix<T>(N, N, T(0)), N_minus_one(N - 1) {}

template <typename T>
SymmetricalMatrix<T>::SymmetricalMatrix(size_t N, std::vector<T> v)
    : Matrix<T>(N, N, v), N_minus_one(N - 1) {
  for (size_t i = 0; i < N_minus_one; i++) {
    if (this->GetElement(i, i) != 0) {
      PrintIndented();
      throw input_error("symmetry matrix does not have zero diagonal\n");
    }

    for (size_t j = 0; j < N; j++) {
      if (this->GetElement(i, j) != this->GetElement(j, i)) {
        PrintIndented();
        throw input_error("symmetry matrix is not symmetrical\n");
      }
    }
  }
}

template <typename T>
void SymmetricalMatrix<T>::SetElement(const size_t i, const size_t j,
                                      const T v) {
  Matrix<T>::SetElement(i, j, v);
  Matrix<T>::SetElement(j, i, v);
}
template <typename T>
void SymmetricalMatrix<T>::FillColumn(const std::vector<T> d, const size_t j) {
  for (size_t i = 0; i < this->N; i++) Matrix<T>::SetElement(i, j, d[i]);
}
template <typename T>
void SymmetricalMatrix<T>::FillRow(const std::vector<T> d, const size_t i) {
  for (size_t j = 0; j < this->N; j++) Matrix<T>::SetElement(i, j, d[j]);
}
template <typename T>
void SymmetricalMatrix<T>::Permute(const combi_ff::Permutations& permutations) {
  for (size_t i = 0; i < permutations.size(); i++) {
    if (permutations[i].first != permutations[i].second)
      Permute(permutations[i]);
  }
}
template <typename T>
void SymmetricalMatrix<T>::Permute(const combi_ff::Permutation permutation) {
  size_t p1 = permutation.first;
  size_t p2 = permutation.second;
  std::vector<T> row1(0);
  std::vector<T> row2(0);

  for (size_t i = 0; i < this->N; i++) {
    row1.push_back(Matrix<T>::GetElement(p1, i));
    row2.push_back(Matrix<T>::GetElement(p2, i));
  }

  FillRow(row1, p2);
  FillRow(row2, p1);
  std::vector<T> col1(0);
  std::vector<T> col2(0);

  for (size_t j = 0; j < this->N; j++) {
    col1.push_back(Matrix<T>::GetElement(j, p1));
    col2.push_back(Matrix<T>::GetElement(j, p2));
  }

  FillColumn(col1, p2);
  FillColumn(col2, p1);
}
template <typename T>
void SymmetricalMatrix<T>::PrintIndented() const {
  for (size_t i = 0; i < this->N; i++) {
    std::cout << '\t';

    for (size_t j = 0; j < this->N; j++)
      std::cout << Matrix<T>::GetElement(i, j) << " ";

    std::cout << std::endl;
  }

  std::cout << std::endl;
}
template <typename T>
bool SymmetricalMatrix<T>::IsConnected() const {
  std::vector<bool> visited(this->N, false);
  std::stack<size_t> stack;
  stack.push(0);
  visited[0] = true;

  while (!stack.empty()) {
    size_t curr_element = stack.top();
    stack.pop();

    for (size_t i = 0; i < this->N; i++) {
      if (!visited[i] && Matrix<T>::GetElement(curr_element, i) != 0) {
        visited[i] = true;
        stack.push(i);
      }
    }
  }

  return (std::find(visited.begin(), visited.end(), false) == visited.end());
}
template <typename T>
bool SymmetricalMatrix<T>::operator==(const SymmetricalMatrix& B) const {
  for (size_t i = 0; i < N_minus_one; i++) {
    for (size_t j = i + 1; j < this->N; j++) {
      if (B.GetElement(i, j) != Matrix<T>::GetElement(i, j)) return false;
    }
  }

  return true;
}

template <typename T>
bool SymmetricalMatrix<T>::operator>(const SymmetricalMatrix& B) const {
  for (size_t i = 0; i < N_minus_one; i++) {
    for (size_t j = i + 1; j < this->N; j++) {
      if (B.GetElement(i, j) != Matrix<T>::GetElement(i, j)) {
        if (B.GetElement(i, j) > Matrix<T>::GetElement(i, j))
          return false;

        else
          return true;
      }
    }
  }

  return false;
}
template <typename T>
bool SymmetricalMatrix<T>::operator<(const SymmetricalMatrix& B) const {
  for (size_t i = 0; i < N_minus_one; i++) {
    for (size_t j = i + 1; j < this->N; j++) {
      if (B.GetElement(i, j) != Matrix<T>::GetElement(i, j)) {
        if (B.GetElement(i, j) < Matrix<T>::GetElement(i, j))
          return false;

        else
          return true;
      }
    }
  }

  return false;
}
template <typename T>
std::string SymmetricalMatrix<T>::GetStack() const {
  std::string stack("");

  for (size_t i = 0; i < N_minus_one; i++) {
    for (size_t j = i + 1; j < this->N; j++)
      stack += std::to_string(Matrix<T>::GetElement(i, j));
  }

  return stack;
}
template <typename T>
bool SymmetricalMatrix<T>::HasCycle() const {
  for (size_t i = 0; i < this->N; i++) {
    for (size_t j = 0; j < this->N; j++) {
      if (Matrix<T>::GetElement(i, j)) {
        std::vector<bool> visited_(this->N, false);
        SymmetricalMatrix<bool> visited(this->N);
        visited_[j] = true;
        visited.SetElement(i, j, true);
        std::stack<size_t> stack;
        stack.push(j);

        while (!stack.empty()) {
          size_t currEl = stack.top();
          stack.pop();

          for (size_t ii = 0; ii < this->N; ii++) {
            if (!visited.GetElement(currEl, ii) &&
                Matrix<T>::GetElement(currEl, ii) != 0) {
              if (ii == i) return true;

              visited.SetElement(currEl, ii, true);
              visited_[ii] = true;
              stack.push(ii);
            }
          }
        }
      }
    }
  }

  return false;
}

/*
AdjacencyMatrix class member functions
*/
template <typename T, typename AtomClass>
AdjacencyMatrix<T, AtomClass>::AdjacencyMatrix() : SymmetricalMatrix<T>(0) {}

template <typename T, typename AtomClass>
AdjacencyMatrix<T, AtomClass>::AdjacencyMatrix(size_t N)
    : SymmetricalMatrix<T>(N),
      atoms(AtomVector<AtomClass>(N)),
      is_aromatic_carbon(N, false),
      lambda_block_end(N, false) {
  min_row.reserve(N);
  max_row.reserve(N);
  indices = std::vector<size_t>(N);
  std::iota(
      indices.begin(), indices.end(),
      0);  // fills indices with increasing indices [0, 1, 2, 3, ..., N - 1]
  lambda.reserve(N);
  types.reserve(N);
  atoms.reserve(N);
}
template <typename T, typename AtomClass>
AdjacencyMatrix<T, AtomClass>::AdjacencyMatrix(
    size_t N, combi_ff::LambdaVector& lambda_,
    combi_ff::AtomVector<AtomClass>& atoms_)
    : AdjacencyMatrix<T, AtomClass>(N) {
  lambda.assign(lambda_.begin(), lambda_.end());
  atoms.assign(atoms_.begin(), atoms_.end());
  size_t sum_lambda = 0;

  for (size_t i = 0; i < lambda.size(); i++) {
    sum_lambda += lambda[i];
    lambda_block_end[sum_lambda - 1] = true;

    for (size_t j = 0; j < lambda[i]; j++) types.push_back(i);
  }

  if (N > 1) lambda_block_end[N - 2] = true;

  for (size_t i = 0; i < this->N_minus_one; i++) {
    for (size_t j = i + 1; j < N; j++) {
      if (this->GetElement(i, j)) {
        atoms[i].AddNeighbor(j);
        atoms[j].AddNeighbor(i);
      }
    }
  }

  for (size_t ii = 0; ii < N; ii++) {
    size_t r = GetTypeNr(ii);
    min_row.push_back(std::accumulate(lambda.begin(), lambda.begin() + r, 0));
    max_row.push_back(
        std::accumulate(lambda.begin(), lambda.begin() + r + 1, 0));
  }

  max_row.back()--;
}
template <typename T, typename AtomClass>
AdjacencyMatrix<T, AtomClass>::AdjacencyMatrix(const AdjacencyMatrix& B)
    : SymmetricalMatrix<T>(B.N),
      lambda(B.lambda),
      types(B.types),
      atoms(B.atoms),
      min_row(B.min_row),
      max_row(B.max_row),
      indices(B.indices),
      is_aromatic_carbon(B.is_aromatic_carbon),
      lambda_block_end(B.lambda_block_end) {
  const std::vector<T> v = B.GetElements();
  SetElements(v);
}

template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::SetElements(std::vector<T> v) {
  SymmetricalMatrix<T>::SetElements(v);

  for (size_t i = 0; i < this->N_minus_one; i++) {
    for (size_t j = i + 1; j < this->N; j++) {
      if (this->GetElement(i, j)) {
        atoms[i].AddNeighbor(j);
        atoms[j].AddNeighbor(i);
      }
    }
  }
}
template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::Print() const {
  std::cout << "    " << std::setfill(' ');

  for (size_t i = 0; i < atoms.size(); i++)
    std::cout << std::setw(3) << atoms[i].GetUnitedAtomSymbol() << " ";

  std::cout << std::endl;

  for (size_t i = 0; i < this->N; i++) {
    std::cout << std::setw(3) << atoms[i].GetUnitedAtomSymbol() << " ";

    for (size_t j = 0; j < this->N; j++)
      std::cout << std::setw(3) << this->GetElement(i, j) << " ";

    std::cout << std::endl;
  }
}

template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::SetElement(const size_t i, const size_t j,
                                               const T v) {
  if (v != 0) {
    atoms[i].AddNeighbor(j);
    atoms[j].AddNeighbor(i);

  } else {
    if (this->GetElement(i, j)) {
      atoms[i].RemoveNeighbor(j);
      atoms[j].RemoveNeighbor(i);
    }
  }

  SymmetricalMatrix<T>::SetElement(i, j, v);
}

template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::SetAtomVector(
    const combi_ff::AtomVector<AtomClass>& a) {
  atoms.assign(a.begin(), a.end());

  for (size_t i = 0; i < this->N_minus_one; i++) {
    for (size_t j = i + 1; j < this->N; j++) {
      if (this->GetElement(i, j)) {
        atoms[i].AddNeighbor(j);
        atoms[j].AddNeighbor(i);
      }
    }
  }
}
template <typename T, typename AtomClass>
double AdjacencyMatrix<T, AtomClass>::GetMass() const {
  double mass = 0;

  for (auto&& a : atoms) mass += a.GetMass();

  return mass;
}
template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::SetLambda(combi_ff::LambdaVector& l) {
  lambda = l;

  for (size_t i = 0; i < lambda.size(); i++) {
    for (size_t j = 0; j < lambda[i]; j++) types.push_back(i);
  }

  for (size_t i = 0; i < this->N_minus_one; i++) {
    for (size_t j = i + 1; j < this->N; j++) {
      if (this->GetElement(i, j)) {
        atoms[i].AddNeighbor(j);
        atoms[j].AddNeighbor(i);
      }
    }
  }

  min_row.reserve(this->N);
  max_row.reserve(this->N);

  for (size_t ii = 0; ii < this->N; ii++) {
    size_t r = GetTypeNr(ii);
    min_row.push_back(std::accumulate(lambda.begin(), lambda.begin() + r, 0));
    max_row.push_back(
        std::accumulate(lambda.begin(), lambda.begin() + r + 1, 0));
  }

  max_row.back()--;
}
template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::ResetAtomNeighbors() {
  for (auto&& a : atoms) a.EraseNeighbors();

  for (size_t i = 0; i < this->N_minus_one; i++) {
    for (size_t j = i + 1; j < this->N; j++) {
      if (this->GetElement(i, j)) {
        atoms[i].AddNeighbor(j);
        atoms[j].AddNeighbor(i);
      }
    }
  }
}
template <typename T, typename AtomClass>
const combi_ff::AtomVector<AtomClass>&
AdjacencyMatrix<T, AtomClass>::GetAtomVector() const {
  return atoms;
}
template <typename T, typename AtomClass>
combi_ff::AtomVector<AtomClass>&
AdjacencyMatrix<T, AtomClass>::GetAtomVectorNonConst() {
  return atoms;
}
template <typename T, typename AtomClass>
AtomClass& AdjacencyMatrix<T, AtomClass>::GetAtom(size_t i) {
  return atoms[i];
}
template <typename T, typename AtomClass>
const combi_ff::LambdaVector& AdjacencyMatrix<T, AtomClass>::GetLambda() const {
  return lambda;
}
template <typename T, typename AtomClass>
const TypeVector& AdjacencyMatrix<T, AtomClass>::GetTypes() const {
  return types;
}
template <typename T, typename AtomClass>
const std::vector<size_t>& AdjacencyMatrix<T, AtomClass>::GetIndices() const {
  return indices;
}
template <typename T, typename AtomClass>
size_t AdjacencyMatrix<T, AtomClass>::GetTypeNr(const size_t i) const {
  return types[i];
}
template <typename T, typename AtomClass>
size_t AdjacencyMatrix<T, AtomClass>::GetMaxRow(size_t i) const {
  return max_row[i];
}
template <typename T, typename AtomClass>
size_t AdjacencyMatrix<T, AtomClass>::GetMinRow(size_t i) const {
  return min_row[i];
}
template <typename T, typename AtomClass>
bool AdjacencyMatrix<T, AtomClass>::IsLastRowOfALambdaBlock(size_t i) {
  return lambda_block_end[i];
}
template <typename T, typename AtomClass>
bool AdjacencyMatrix<T, AtomClass>::hasCycle() const {
  for (size_t i = 0; i < this->N; i++) {
    const combi_ff::NeighborVector& neighbors = atoms[i].GetNeighbors();

    for (auto neighbor : neighbors) {
      std::vector<bool> visited_(this->N, false);
      SymmetricalMatrix<bool> visited(this->N);
      visited_[neighbor] = true;
      visited.SetElement(i, neighbor, true);
      std::stack<size_t> stack;
      stack.push(neighbor);

      while (!stack.empty()) {
        size_t curr_element = stack.top();
        stack.pop();

        for (size_t ii = 0; ii < this->N; ii++) {
          if (!visited.GetElement(curr_element, ii) &&
              this->GetElement(curr_element, ii) != 0) {
            if (ii == i) return true;

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
template <typename T, typename AtomClass>
bool AdjacencyMatrix<T, AtomClass>::IsConnected() const {
  std::vector<bool> visited(this->N, false);
  std::stack<size_t> stack;
  visited[0] = true;

  // add neighbors of 1st atom to stack manually to save time, since we know
  // that visited[neighbors[0]] is false for all the neighbors
  for (auto&& nbr : atoms[0].GetNeighbors()) {
    visited[nbr] = true;
    stack.push(nbr);
  }

  while (!stack.empty()) {
    size_t currEl = stack.top();
    stack.pop();

    for (auto&& nbr : atoms[currEl].GetNeighbors()) {
      if (!visited[nbr]) {
        visited[nbr] = true;

        if (atoms[nbr].GetDegree() >
            1)  // if the degree is 1, the only neighbor is the current one
                // that we already visited
          stack.push(nbr);
      }
    }
  }

  return (std::find(visited.begin(), visited.end(), false) == visited.end());
}
template <typename T, typename AtomClass>
bool AdjacencyMatrix<T, AtomClass>::SmallerThanPermuted(
    const size_t curr_row, const std::vector<size_t>& permuted_indices,
    size_t& i_pos, size_t& j_pos, size_t& i_pos_perm, size_t& j_pos_perm,
    size_t& i0, size_t& j0, bool& diff, const size_t smallest_diff_index,
    const size_t highest_permuted_index) {
  // permutedIndices.assign(indices.begin(), indices.end());
  // PermuteVector(permutedIndices, permutations);
  size_t lim = std::min(max_row[curr_row], highest_permuted_index + 1);

  for (size_t i = min_row[curr_row]; i < lim; i++) {
    size_t i2 = permuted_indices[i];

    for (size_t j = std::max(i + 1, smallest_diff_index); j < this->N; j++) {
      size_t j2 = permuted_indices[j];

      if (this->GetElement(i, j) != this->GetElement(i2, j2)) {
        diff = true;

        if (this->GetElement(i, j) >
            this->GetElement(
                i2,
                j2)) {  // matrix is larger than permuted one in current block
          j_pos = j;
          return false;

        } else {  // matrix is smaller than permuted one in current block -> not
                  // canonical
          i_pos = i;
          j_pos = j;
          i_pos_perm = i2;
          j_pos_perm = j2;
          return true;
        }
      }

      if (this->GetElement(i2, j2)) {  // permuted value is larger than 0
        if (i2 < j2) {
          if (i2 > i0 || (i2 == i0 && j2 > j0)) {
            i0 = i2;
            j0 = j2;
          }

        } else {
          if (j2 > i0 || (j2 == i0 && i2 > j0)) {
            i0 = j2;
            j0 = i2;
          }
        }
      }
    }
  }

  // matrix is equal to permuted one in current block
  return false;
}

template <typename T, typename AtomClass>
bool AdjacencyMatrix<T, AtomClass>::EqualTo(
    const std::vector<size_t>& permuted_indices, size_t& j_pos) {
  // permutedIndices.assign(indices.begin(), indices.end());
  // PermuteVector(permutedIndices, permutations);
  for (size_t i = 0; i < this->N; i++) {
    size_t i2 = permuted_indices[i];

    for (size_t j = i + 1; j < this->N_minus_one; j++) {
      size_t j2 = permuted_indices[j];

      if (this->GetElement(i, j) != this->GetElement(i2, j2)) {
        j_pos = j;
        return false;
      }
    }
  }

  return true;
}
template <typename T, typename AtomClass>
size_t AdjacencyMatrix<T, AtomClass>::GetNumDoubleBonds() const {
  size_t n_double(0);

  for (size_t i = 0; i < this->N_minus_one; i++) {
    for (size_t j = i + 1; j < this->N; j++) {
      if (this->GetElement(i, j) == 2) n_double++;
    }
  }

  return n_double;
}
template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::GetNumMultipleBonds(
    size_t& n_single, size_t& n_double, size_t& n_triple, size_t& n_quadruple,
    size_t& n_aromatic) const {
  for (size_t i = 0; i < this->N_minus_one; i++) {
    for (size_t j = i + 1; j < this->N; j++) {
      T deg = this->GetElement(i, j);

      if (deg == 1)
        n_single++;

      else if (deg == 2)
        n_double++;

      else if (deg == 3)
        n_triple++;

      else if (deg == 4)
        n_quadruple++;

      else if (deg == 1.5)
        n_aromatic++;
    }

    n_single += atoms[i].GetNumHydrogens() + atoms[i].GetNumFixedHydrogens();
  }
}
template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::GetNumMultipleBonds(
    std::vector<size_t>& ranges) const {
  for (size_t i = 0; i < this->N_minus_one; i++) {
    ranges[range_single_bonds] +=
        atoms[i].GetNumHydrogens() + atoms[i].GetNumFixedHydrogens();

    for (size_t j = i + 1; j < this->N; j++) {
      size_t deg = this->GetElement(i, j);

      if (deg == 1)
        ranges[range_single_bonds]++;

      else if (deg == 2)
        ranges[range_double_bonds]++;

      else if (deg == 3)
        ranges[range_triple_bonds]++;

      else if (deg == 4)
        ranges[range_quadruple_bonds]++;
    }
  }

  ranges[range_bonds] =
      std::accumulate(ranges.begin() + range_single_bonds,
                      ranges.begin() + range_quadruple_bonds, 0);
}

template <typename T, typename AtomClass>
std::vector<bool> AdjacencyMatrix<T, AtomClass>::IsInCycle() const {
  std::vector<bool> in_cylce(this->N);

  for (size_t j = 0; j < this->N; j++) {
    bool cyc = false;
    const combi_ff::NeighborVector& neighbors = atoms[j].GetNeighbors();

    for (auto neighbor : neighbors) {
      if (cyc) break;

      std::vector<bool> visited_(this->N, false);
      SymmetricalMatrix<bool> visited(this->N);
      visited_[neighbor] = true;
      visited.SetElement(j, neighbor, true);
      std::stack<size_t> stack;
      stack.push(neighbor);

      while (!cyc && !stack.empty()) {
        size_t curr_element = stack.top();
        stack.pop();

        for (size_t i = 0; i < this->N; i++) {
          if (!visited.GetElement(curr_element, i) &&
              this->GetElement(curr_element, i) != 0) {
            visited.SetElement(curr_element, i, true);
            visited_[i] = true;
            stack.push(i);

            if (i == j) cyc = true;
          }
        }
      }
    }

    if (cyc) in_cylce[j] = true;
  }

  return in_cylce;
}
template <typename T, typename AtomClass>
bool AdjacencyMatrix<T, AtomClass>::AreInSameCycle(size_t atom1,
                                                   size_t atom2) const {
  const combi_ff::NeighborVector& neighbors = atoms[atom1].GetNeighbors();
  SymmetricalMatrix<bool> visited(this->N);
  std::stack<size_t> stack;
  bool visit1 = false;
  bool visit2 = false;

  for (auto neighbor : neighbors) {
    if (neighbor != atom2) {
      visited.SetElement(atom1, neighbor, true);
      stack.push(neighbor);

      while (!stack.empty()) {
        size_t currEl = stack.top();
        stack.pop();

        for (auto&& nbr : atoms[currEl].GetNeighbors()) {
          if (!visited.GetElement(currEl, nbr)) {
            visited.SetElement(currEl, nbr, true);

            if (nbr != atom1) stack.push(nbr);

            if (nbr == atom2)
              visit2 = true;

            else if (nbr == atom1)
              visit1 = true;
          }
        }
      }
    }
  }

  return (visit1 && visit2);
}
template <typename T, typename AtomClass>
size_t AdjacencyMatrix<T, AtomClass>::GetNumPiAtoms() const {
  const std::vector<bool> inCycle = IsInCycle();
  size_t n_pi(0);

  for (size_t i = 0; i < this->N; i++) {
    if (inCycle[i]) {
      const combi_ff::NeighborVector& neighbors = atoms[i].GetNeighbors();
      int n_double(0);

      for (auto neighbor : neighbors) {
        if (inCycle[neighbor]) {
          if (this->GetElement(i, neighbor) == 2) {
            // nPi++;
            n_double++;
          }

          //!!! For now, assume that a molecule that has a triple bond in *any*
          //! cycle
          // is not aromatic, but this is not actually true !!!
          else if (this->GetElement(i, neighbor) == 3)
            return 0;
        }

        // !!! For now, assume that a molecucle that has an atom that's
        // connected with a double bond to an atom outside the ring is not
        // aromatic, but this is not actually true !!!
        else if (this->GetElement(i, neighbor) > 1)
          return 0;
      }

      // !!! For now, assume that a molecucle that has an atom that's connected
      // to more than one atom outside the ring is not aromatic, but this is not
      // actually true !!!
      if ((neighbors.size() + atoms[i].GetNumHydrogens() +
           atoms[i].GetNumFixedHydrogens()) > 3)
        return 0;

      if (n_double == 1) n_pi++;

      // !!! for now, assume that a molecule that has an atom with two double
      // bonds in a cycle is not aromatic, but this is not actually true !!!
      else if (n_double > 1)
        return 0;

      else if (n_double == 0 &&
               (atoms[i].GetDegree() + atoms[i].GetNumHydrogens() +
                atoms[i].GetNumFixedHydrogens()) > 1 &&
               (atoms[i].GetDegree() + atoms[i].GetNumHydrogens() +
                atoms[i].GetNumFixedHydrogens()) < 4) {
        if (neighbors.size() <= 3) n_pi += 2;
      }
    }
  }

  return n_pi;
}
template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::MakeCanonical(
    const combi_ff::RepresentationSystem& u) {
  AdjacencyMatrix B(*this);
  combi_ff::PermutationIterator it_perm(u);
  combi_ff::Permutations maximizing_automorphic_permutation;
  const std::vector<size_t>* permuted_indices;

  while (it_perm.GetNextPermutation()) {
    bool diff(false);
    permuted_indices = it_perm.GetPermutedIndices();

    for (size_t i = 0; i < this->N; i++) {
      for (size_t j = std::max(i + 1, it_perm.GetSmallestDiffIndex());
           j < this->N; j++) {
        if (B.GetElement(i, j) <
            this->GetElement((*permuted_indices)[i], (*permuted_indices)[j])) {
          B = AdjacencyMatrix(*this);
          maximizing_automorphic_permutation =
              *(it_perm.GetCombinedPermutation());
          B.Permute(maximizing_automorphic_permutation);
          diff = true;
          break;

        } else if (B.GetElement(i, j) >
                   this->GetElement((*permuted_indices)[i],
                                    (*permuted_indices)[j])) {
          if (it_perm.GetCurrentIndex() > j) it_perm.SetCurrentIndex(j);

          diff = true;
          break;
        }
      }

      if (diff) break;
    }
  }

  *this = AdjacencyMatrix(B);
  combi_ff::PermuteVector(is_aromatic_carbon,
                          maximizing_automorphic_permutation);
  combi_ff::PermuteVector(atoms, maximizing_automorphic_permutation);
  ResetAtomNeighbors();
}

template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::MakeCanonical(
    const combi_ff::RepresentationSystem& u, std::vector<size_t>& index,
    const size_t limit, Permutations& matrix_permutations) {
  AdjacencyMatrix B(*this);
  combi_ff::PermutationIterator it_perm(u);
  combi_ff::Permutations maximizing_automorphic_permutation;
  const std::vector<size_t>* permuted_indices;

  size_t iteration = 0;

  while (it_perm.GetNextPermutation() && ++iteration < limit) {
    bool diff(false);
    permuted_indices = it_perm.GetPermutedIndices();

    for (size_t i = 0; i < this->N; i++) {
      for (size_t j = std::max(i + 1, it_perm.GetSmallestDiffIndex());
           j < this->N; j++) {
        if (B.GetElement(i, j) <
            this->GetElement((*permuted_indices)[i], (*permuted_indices)[j])) {
          B = AdjacencyMatrix(*this);
          maximizing_automorphic_permutation =
              *(it_perm.GetCombinedPermutation());
          B.Permute(maximizing_automorphic_permutation);
          diff = true;
          break;

        } else if (B.GetElement(i, j) >
                   this->GetElement((*permuted_indices)[i],
                                    (*permuted_indices)[j])) {
          if (it_perm.GetCurrentIndex() > j) it_perm.SetCurrentIndex(j);

          diff = true;
          break;
        }
      }

      if (diff) break;
    }
  }

  if (iteration == limit)
    std::cerr
        << "?Warning: didn't finish matrix canonicalization due to iteration "
           "limit. You can change the limit in src/cnv/Handler.h\n";

  *this = AdjacencyMatrix(B);
  matrix_permutations.insert(matrix_permutations.end(),
                             maximizing_automorphic_permutation.begin(),
                             maximizing_automorphic_permutation.end());
  combi_ff::PermuteVector(index, maximizing_automorphic_permutation);
  combi_ff::PermuteVector(is_aromatic_carbon,
                          maximizing_automorphic_permutation);
  combi_ff::PermuteVector(atoms, maximizing_automorphic_permutation);
  ResetAtomNeighbors();
}

template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::SortAtomVector() {
  //   std::cout << "sort1\n";
  combi_ff::LambdaVector lambda_new(0);
  //   A.print();

  for (size_t i = 0; i < this->N; i++) {
    for (size_t j = 0; j < this->N; j++) {
      if (atoms[i] < atoms[j]) {
        std::swap(atoms[i], atoms[j]);
        this->Permute(combi_ff::Permutation({i, j}));
      }
    }
  }

  for (auto&& a : atoms) a.SetNeighbors(combi_ff::NeighborVector(0));

  SetAtomVector(atoms);
  lambda_new.push_back(1);
  std::vector<combi_ff::AtomVector<AtomClass>> atom_blocks(0);
  atom_blocks.push_back(combi_ff::AtomVector<AtomClass>(1, atoms[0]));

  for (size_t i = 1; i < this->N; i++) {
    if (atoms[i - 1] != atoms[i]) {
      lambda_new.push_back(1);
      atom_blocks.push_back(combi_ff::AtomVector<AtomClass>(1, atoms[i]));

    } else {
      atom_blocks.back().push_back(atoms[i]);
      lambda_new.back()++;
    }
  }

  for (size_t i = 0; i < lambda_new.size(); i++) {
    for (size_t j = i + 1; j < lambda_new.size(); j++) {
      if (atom_blocks[i].front().GetDegree() ==
              atom_blocks[j].front().GetDegree() &&
          lambda_new[i] > lambda_new[j]) {
        //  std::cout << "b have to swap " << i << " and " << j << std::endl;
        //  A.print();
        size_t ii =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + i, 0);
        size_t jj =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + j + 1, -1);
        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";
        //  std::cout << std::endl;
        std::swap(atom_blocks[i], atom_blocks[j]);

        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";

        //  std::cout << std::endl;

        for (size_t k = 0; k < std::max(lambda_new[i], lambda_new[j]); k++) {
          // std::cout << "  b " << A.GetN() << " " << ii << " " << jj <<
          // std::endl;
          this->Permute(combi_ff::Permutation({ii, jj}));
          // std::cout << &isAromatic << std::endl;
          std::swap(atoms[ii], atoms[jj]);
          // A.print();
          ii++;
          jj--;
        }

        std::swap(lambda_new[i], lambda_new[j]);

      } else if (atom_blocks[i].front().GetDegree() ==
                     atom_blocks[j].front().GetDegree() &&
                 lambda_new[i] == lambda_new[j] &&
                 atom_blocks[i].front() > atom_blocks[j].front()) {
        //  std::cout << "b have to swap " << i << " and " << j << std::endl;
        //  A.print();
        size_t ii =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + i, 0);
        size_t jj =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + j + 1, -1);
        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";
        //  std::cout << std::endl;
        std::swap(atom_blocks[i], atom_blocks[j]);

        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";

        //  std::cout << std::endl;

        for (size_t k = 0; k < std::max(lambda_new[i], lambda_new[j]); k++) {
          // std::cout << "  b " << A.GetN() << " " << ii << " " << jj <<
          // std::endl;
          this->Permute(combi_ff::Permutation({ii, jj}));
          // std::cout << &isAromatic << std::endl;
          std::swap(atoms[ii], atoms[jj]);
          // A.print();
          ii++;
          jj--;
        }

        std::swap(lambda_new[i], lambda_new[j]);
      }
    }
  }

  ResetAtomNeighbors();
  SetLambda(lambda_new);
}

template <typename T, typename AtomClass>
std::vector<size_t> AdjacencyMatrix<T, AtomClass>::SortAtomVector(
    std::vector<bool>& is_aromatic) {
  //  A.print();
  std::vector<size_t> idx(GetIndices());
  combi_ff::LambdaVector lambda_new(0);
  Permutations matrix_permutations(0);

  for (uint i = 0; i < this->N; i++) {
    for (uint j = 0; j < this->N; j++) {
      if (atoms[i] < atoms[j]) {
        std::swap(atoms[i], atoms[j]);
        this->Permute(combi_ff::Permutation({i, j}));
        matrix_permutations.push_back(combi_ff::Permutation({i, j}));
        std::swap(is_aromatic[i], is_aromatic[j]);
        std::swap(idx[i], idx[j]);
      }
    }
  }

  for (auto&& a : atoms) a.SetNeighbors(combi_ff::NeighborVector(0));

  SetAtomVector(atoms);
  lambda_new.push_back(1);
  std::vector<combi_ff::AtomVector<AtomClass>> atom_blocks(0);
  atom_blocks.push_back(combi_ff::AtomVector<AtomClass>(1, atoms[0]));

  for (uint i = 1; i < this->N; i++) {
    if (atoms[i - 1] != atoms[i]) {
      lambda_new.push_back(1);
      atom_blocks.push_back(combi_ff::AtomVector<AtomClass>(1, atoms[i]));

    } else {
      atom_blocks.back().push_back(atoms[i]);
      lambda_new.back()++;
    }
  }

  for (uint i = 0; i < lambda_new.size(); i++) {
    for (uint j = i + 1; j < lambda_new.size(); j++) {
      // assert(atomBlocks[i].front().GetDegree() >=
      // atomBlocks[j].front().GetDegree());
      if (atom_blocks[i].front().GetDegree() ==
              atom_blocks[j].front().GetDegree() &&
          lambda_new[i] > lambda_new[j]) {
        //  std::cout << "b have to swap " << i << " and " << j << std::endl;
        //  A.print();
        uint ii =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + i, 0);
        uint jj =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + j + 1, -1);
        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";
        //  std::cout << std::endl;
        std::swap(atom_blocks[i], atom_blocks[j]);

        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";

        //  std::cout << std::endl;

        for (uint k = 0; k < std::max(lambda_new[i], lambda_new[j]); k++) {
          // std::cout << "  b " << A.GetN() << " " << ii << " " << jj <<
          // std::endl;
          this->Permute(combi_ff::Permutation({ii, jj}));
          matrix_permutations.push_back(combi_ff::Permutation({i, j}));
          std::swap(is_aromatic[ii], is_aromatic[jj]);
          std::swap(idx[ii], idx[jj]);
          // std::cout << &isAromatic << std::endl;
          std::swap(atoms[ii], atoms[jj]);
          // A.print();
          ii++;
          jj--;
        }

        std::swap(lambda_new[i], lambda_new[j]);

      } else if (atom_blocks[i].front().GetDegree() ==
                     atom_blocks[j].front().GetDegree() &&
                 lambda_new[i] == lambda_new[j] &&
                 atom_blocks[i].front() > atom_blocks[j].front()) {
        //  std::cout << "b have to swap " << i << " and " << j << std::endl;
        //  A.print();
        uint ii =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + i, 0);
        uint jj =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + j + 1, -1);
        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";
        //  std::cout << std::endl;
        std::swap(atom_blocks[i], atom_blocks[j]);

        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";

        //  std::cout << std::endl;

        for (uint k = 0; k < std::max(lambda_new[i], lambda_new[j]); k++) {
          // std::cout << "  b " << A.GetN() << " " << ii << " " << jj <<
          // std::endl;
          this->Permute(combi_ff::Permutation({ii, jj}));
          matrix_permutations.push_back(combi_ff::Permutation({i, j}));
          std::swap(is_aromatic[ii], is_aromatic[jj]);
          std::swap(idx[ii], idx[jj]);
          // std::cout << &isAromatic << std::endl;
          std::swap(atoms[ii], atoms[jj]);
          // A.print();
          ii++;
          jj--;
        }

        std::swap(lambda_new[i], lambda_new[j]);
      }
    }
  }

  ResetAtomNeighbors();
  SetLambda(lambda_new);
  return idx;
}

template <typename T, typename AtomClass>
Permutations AdjacencyMatrix<T, AtomClass>::SortAtomVector_(
    std::vector<bool>& is_aromatic) {
  //  A.print();
  std::vector<size_t> idx(GetIndices());
  combi_ff::LambdaVector lambda_new(0);
  Permutations matrix_permutations(0);

  for (uint i = 0; i < this->N; i++) {
    for (uint j = 0; j < this->N; j++) {
      if (atoms[i] < atoms[j]) {
        std::swap(atoms[i], atoms[j]);
        this->Permute(combi_ff::Permutation({i, j}));
        matrix_permutations.push_back(combi_ff::Permutation({i, j}));
        std::swap(is_aromatic[i], is_aromatic[j]);
        std::swap(idx[i], idx[j]);
      }
    }
  }

  for (auto&& a : atoms) {
    a.SetNeighbors(combi_ff::NeighborVector(0));
  }

  SetAtomVector(atoms);

  lambda_new.push_back(1);
  std::vector<combi_ff::AtomVector<AtomClass>> atom_blocks(0);
  atom_blocks.push_back(combi_ff::AtomVector<AtomClass>(1, atoms[0]));

  for (uint i = 1; i < this->N; i++) {
    if (atoms[i - 1] != atoms[i]) {
      lambda_new.push_back(1);
      atom_blocks.push_back(combi_ff::AtomVector<AtomClass>(1, atoms[i]));

    } else {
      atom_blocks.back().push_back(atoms[i]);
      lambda_new.back()++;
    }
  }

  for (uint i = 0; i < lambda_new.size(); i++) {
    for (uint j = i + 1; j < lambda_new.size(); j++) {
      // assert(atomBlocks[i].front().GetDegree() >=
      // atomBlocks[j].front().GetDegree());
      if (atom_blocks[i].front().GetDegree() ==
              atom_blocks[j].front().GetDegree() &&
          lambda_new[i] > lambda_new[j]) {
        //  std::cout << "b have to swap " << i << " and " << j << std::endl;
        //  A.print();
        uint ii =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + i, 0);
        uint jj =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + j + 1, -1);
        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";
        //  std::cout << std::endl;
        std::swap(atom_blocks[i], atom_blocks[j]);

        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";

        //  std::cout << std::endl;

        for (uint k = 0; k < std::max(lambda_new[i], lambda_new[j]); k++) {
          // std::cout << "  b " << A.GetN() << " " << ii << " " << jj <<
          // std::endl;
          this->Permute(combi_ff::Permutation({ii, jj}));
          matrix_permutations.push_back(combi_ff::Permutation({ii, jj}));
          std::swap(is_aromatic[ii], is_aromatic[jj]);
          std::swap(idx[ii], idx[jj]);
          // std::cout << &isAromatic << std::endl;
          std::swap(atoms[ii], atoms[jj]);

          // A.print();
          ii++;
          jj--;
        }

        std::swap(lambda_new[i], lambda_new[j]);

      } else if (atom_blocks[i].front().GetDegree() ==
                     atom_blocks[j].front().GetDegree() &&
                 lambda_new[i] == lambda_new[j] &&
                 atom_blocks[i].front() > atom_blocks[j].front()) {
        //  std::cout << "b have to swap " << i << " and " << j << std::endl;
        //  A.print();
        uint ii =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + i, 0);
        uint jj =
            std::accumulate(lambda_new.begin(), lambda_new.begin() + j + 1, -1);
        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";
        //  std::cout << std::endl;
        std::swap(atom_blocks[i], atom_blocks[j]);

        //  for(auto && a : atomBlocks)
        //      std::cout << &a << " ";

        //  std::cout << std::endl;

        for (uint k = 0; k < std::max(lambda_new[i], lambda_new[j]); k++) {
          // std::cout << "  b " << A.GetN() << " " << ii << " " << jj <<
          // std::endl;
          this->Permute(combi_ff::Permutation({ii, jj}));
          matrix_permutations.push_back(combi_ff::Permutation({ii, jj}));
          std::swap(is_aromatic[ii], is_aromatic[jj]);
          std::swap(idx[ii], idx[jj]);
          std::swap(atoms[ii], atoms[jj]);
          // std::cout << &isAromatic << std::endl;
          // A.print();
          ii++;
          jj--;
        }

        std::swap(lambda_new[i], lambda_new[j]);
      }
    }
  }

  ResetAtomNeighbors();
  SetLambda(lambda_new);
  return matrix_permutations;
}

template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::PrintToFile(
    std::ostream& outputFile) const {
  outputFile << "    ";

  for (size_t i = 0; i < atoms.size(); i++)
    outputFile << atoms[i].GetUnitedAtomSymbol() << atoms[i].GetFormalCharge()
               << " ";

  outputFile << std::endl;

  for (size_t i = 0; i < this->N; i++) {
    std::string symbol =
        atoms[i].GetUnitedAtomSymbol() + atoms[i].GetFormalCharge();
    outputFile << std::setw(3) << symbol << " ";

    for (size_t j = 0; j < this->N; j++)
      outputFile << std::setw((int)atoms[j].GetUnitedAtomSymbol().size() +
                              (int)atoms[j].GetFormalCharge().size())
                 << this->GetElement(i, j) << " ";

    outputFile << std::endl;
  }
}
template <typename T, typename AtomClass>
std::vector<bool>& AdjacencyMatrix<T, AtomClass>::GetIsAromaticCarbon() {
  return is_aromatic_carbon;
}
template <typename T, typename AtomClass>
const std::vector<bool>&
AdjacencyMatrix<T, AtomClass>::GetIsAromaticCarbonConst() const {
  return is_aromatic_carbon;
}
template <typename T, typename AtomClass>
void AdjacencyMatrix<T, AtomClass>::SetIsAromaticCarbon(
    const std::vector<bool>& a) {
  is_aromatic_carbon = a;
}
template <typename T, typename AtomClass>
AdjacencyMatrix<T, AtomClass>
AdjacencyMatrix<T, AtomClass>::ExtendToFullMatrix() const {
  size_t n_implicit_H(0);

  for (auto&& a : GetAtomVector()) n_implicit_H += a.GetNumTotalHydrogens();

  if (!n_implicit_H) return *this;

  AdjacencyMatrix<T, AtomClass> A_full(this->N + n_implicit_H);
  AtomVector<AtomClass> atoms_full(GetAtomVector());

  for (auto&& a : atoms_full) a.RemoveHydrogens();

  AtomVector<AtomClass> atmp(n_implicit_H, AtomClass("H"));
  atoms_full.insert(atoms_full.end(), atmp.begin(), atmp.end());
  A_full.SetAtomVector(atoms_full);
  size_t Hindex(this->N);

  for (size_t i = 0; i < this->N; i++) {
    for (size_t j = 0; j < this->N; j++)
      A_full.SetElement(i, j, this->GetElement(i, j));

    for (size_t j = 0; j < GetAtomVector()[i].GetNumTotalHydrogens(); j++) {
      // assert(Hindex < A.GetAtomVector().size());
      A_full.SetElement(i, Hindex++, 1);
    }
  }

  return A_full;
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& m) {
  size_t N = m.GetN();
  size_t M = m.GetM();

  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < M; j++) stream << m.GetElement(i, j) << " ";

    stream << std::endl;
  }

  return stream;
}

}  // namespace combi_ff

#endif
