// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "MaxFillAlg.h"

#include <chrono>

namespace combi_ff {

namespace enu {

MaxFillAlg::MaxFillAlg(enu::AdjacencyMatrix& A, const bool stereo,
                       const size_t max_degree)
    : N(A.GetN()),
      stereo(stereo),
      i(0),
      j(0),
      A(A),
      degree_vec(std::vector<int>(N)),
      found(false),
      accumulated_row_matrix(TriangularMatrix<size_t>(N + 1)),
      accumulated_column_matrix(TriangularMatrix<size_t>(N + 1)),
      perm_it(N),
      next_step_ptr(
          &MaxFillAlg::MaxFillAlg::ForwardStep) {  // Set fw to true initially
  const LambdaVector& lambda = A.GetLambda();
  // num_perms[i] indicates, with how many other atom indices the index of atom
  // i can be permuted note: index permutations are only allowed within the same
  // lambda partition, and with higher indices e.g. for lambdaHat = [1, 3, 1, 4,
  // 6], we have num_perms = [1, 3, 2, 1, 1, 4, 3, 2, 1, 6, 5, 4, 3, 2, 1]
  std::vector<size_t> num_perms(
      std::accumulate(lambda.begin(), lambda.end(), 0));
  size_t idx(0);

  for (size_t ii = 0; ii < lambda.size(); ii++) {
    for (size_t jj = 0; jj < lambda[ii]; jj++)
      num_perms[idx++] = (lambda[ii] - jj);
  }

  // a RepresentationSystem contains all the possible index permutations for the
  // different atom indices
  //  e.g. for lambdaHat = [2, 1, 4], the RepresentationSystem is
  //  (0,0), (0,1)                 (for idx 0)
  //  (1,1)                        (for idx 1)
  //  (2,2)                        (for idx 2)
  //  (3,3), (3,4), (3,5), (3,6)   (for idx 3)
  //  (4,4), (4,5), (4,6)          (for idx 4)
  //  (5,5), (5,6)                 (for idx 5)
  //  (6,6)                        (for idx 6)
  // representation system for the identical permutations
  id = RepresentationSystem(N - 1);

  for (size_t ii = 0; ii < id.size(); ii++) {
    id[ii] = StabilizerVector(1, Permutations(1, {ii, ii}));
  }

  // representation system for the Automorphism Group U0
  RepresentationSystem u0(N - 1);

  for (size_t ii = 0; ii < u0.size(); ii++) {
    u0[ii].resize(num_perms[ii]);
    size_t idx = 0;

    for (size_t jj = 0; jj < num_perms[ii]; jj++) {
      u0[ii][idx++] = (Permutations(1, {ii, ii + jj}));
    }
  }

  row_stabilizer_rep = std::vector<RepresentationSystem>(1);

  row_stabilizer_rep[0] = u0;

  for (size_t ii = 1; ii < row_stabilizer_rep.size(); ii++) {
    row_stabilizer_rep[ii].resize(u0.size());
    row_stabilizer_rep[ii].assign(id.begin(), id.end());

    for (size_t jj = 0; jj < row_stabilizer_rep[ii].size(); jj++)
      row_stabilizer_rep[ii][jj].reserve(u0[jj].size());
  }

  // lambda_prime is used for the semicanonicity test, storing the refined
  // partitions for each row

  for (size_t ii = 0; ii < N; ii++)
    degree_vec[ii] = (int)A.GetAtomVector()[ii].GetDegree();

  M = SymmetricalMatrix<int>(N);
  L = TriangularMatrix<int>(N);
  C = TriangularMatrix<int>(N);
  // fill the supplementary matrices M, L, C
  GetMLC(max_degree);
}

/*******************************************************
FIND THE LEXICOGRAPHICALLY NEXT SMALLER ADJACENCY MATRIX
*******************************************************/
bool MaxFillAlg::GetNextCanonicalMatrix() {
  found = false;

  // loop is necessary to avoid call stack overflow due to recursion!
  do {
    (this->*next_step_ptr)();

    if (found) return true;
  } while (j > 0);

  return false;
}

/************************************
DETERMINE THE VALUES l_hat, c_hat and m
************************************/
std::pair<int, int> MaxFillAlg::GetHat() {
  return {degree_vec[i] - (int)accumulated_row_matrix.GetElement(i, j),
          degree_vec[j] - (int)accumulated_column_matrix.GetElement(i, j)};
}

/*****************************
Set (i,j) TO THE NEXT POSITION
*****************************/
void MaxFillAlg::SetNextIndex() {
  // if we're at the end of a row, increase the row number and Set the column
  // number s.t. we're one value off the diagonal
  if (j == N - 1) {
    i++;
    j = i + 1;
  }

  // otherwise just increase j
  else
    j++;
}

/*********************************
Set (i,j) TO THE PREVIOUS POSITION
*********************************/
void MaxFillAlg::SetPreviousIndex() {
  // if we're at a position that's one element off the diagonal, decrease the
  // row number and Set j to the last column number
  if (j == i + 1) {
    i--;
    j = N - 1;
  }

  // otherwise just decrease j
  else
    j--;
}

/********************************************
FORWARD STEP OF THE MAXIMUM FILLING ALGORITHM
********************************************/
void MaxFillAlg::ForwardStep() {
  SetNextIndex();
  int l_hat, c_hat;
  std::tie(l_hat, c_hat) = GetHat();
  int m = M.GetElement(i, j);
  // calculate the largest possible value that the current A can take for
  // element (i,j) and the smallest possible value it has to take
  size_t x =
      std::min({l_hat, c_hat,
                m});  // note: l_hat, c_hat, m are all positive by definition!
  size_t min = std::max(
      {0, l_hat - L.GetElement(i, j),
       c_hat -
           C.GetElement(i, j)});  // note: min is at least 0, i.e. non-negative

  // if the maximum possible entry is smaller than the minimum necessary entry
  // -> backstep
  if (x < min) next_step_ptr = &MaxFillAlg::BackwardStep;

  // otherwise, accept the maximum possible entry as the new element at position
  // (i,j)
  else {
    A.SetElement(i, j, x);
    accumulated_row_matrix.SetElement(
        i, j + 1, accumulated_row_matrix.GetElement(i, j) + x);
    accumulated_column_matrix.SetElement(
        i + 1, j, accumulated_column_matrix.GetElement(i, j) + x);

    // if we're at the end of a row, calculate the next lambdas
    if (i == N - 2) {
      // off-diagonal element
      accumulated_row_matrix.SetElement(
          i + 1, i + 2, accumulated_column_matrix.GetElement(i + 1, i + 1));

      // if all the rows are filled OR all the rows of the current atom type are
      // filled -> test for canonicity
      CanonicityTest();

    } else if (j == N - 1) {
      accumulated_row_matrix.SetElement(
          i + 1, i + 2, accumulated_column_matrix.GetElement(i + 1, i + 1));
    }
  }
}

/*********************************************
BACKWARD STEP OF THE MAXIMUM FILLING ALGORITHM
*********************************************/
void MaxFillAlg::BackwardStep() {
  // if we're back at the first entry of the matrix, the algorithm is finished
  if (j == 1) {
    j--;  // set j to 0 to stop while loop in GetNextCanonicalMatrix
    return;

  } else {
    // decrease the index pair (i,j)
    SetPreviousIndex();
    // calculate l_hat, c_hat, and m
    int l_hat, c_hat;
    std::tie(l_hat, c_hat) = GetHat();
    int x_minus_one = (int)A.GetElement(i, j) - 1;

    // check if entry x at (i,j) can be decreased
    if ((x_minus_one >= 0) && (l_hat - x_minus_one <= L.GetElement(i, j)) &&
        (c_hat - x_minus_one <= C.GetElement(i, j))) {
      A.SetElement(i, j, x_minus_one);
      accumulated_row_matrix.SetElement(
          i, j + 1, accumulated_row_matrix.GetElement(i, j + 1) - 1);
      accumulated_column_matrix.SetElement(
          i + 1, j, accumulated_column_matrix.GetElement(i + 1, j) - 1);
      next_step_ptr = &MaxFillAlg::ForwardStep;

    } else
      next_step_ptr = &MaxFillAlg::BackwardStep;
  }
}

void MaxFillAlg::CanonicityTest() {
  if (A.IsConnected() && IsCanonical()) found = true;

  SetNextIndex();
  next_step_ptr = &MaxFillAlg::BackwardStep;
}

bool MaxFillAlg::IsCanonical() {
  perm_it.Reset(&row_stabilizer_rep[0], N - 2);
  const std::vector<size_t>* permuted_indices;

  while (perm_it.GetNextPermutation()) {
    permuted_indices = perm_it.GetPermutedIndices();
    if (A.SmallerThanPermuted(*permuted_indices)) {
      return false;  // not canonical
    }
  }

  return true;  // canonical
}

const RepresentationSystem& MaxFillAlg::GetUAutomorph() const {
  return row_stabilizer_rep.back();
}
const RepresentationSystem& MaxFillAlg::GetUId() const { return id; }

const RepresentationSystem& MaxFillAlg::GetU0() const {
  return row_stabilizer_rep.front();
}

RepresentationSystem& MaxFillAlg::GetUAutomorph_non_const() {
  return row_stabilizer_rep.back();
}
RepresentationSystem& MaxFillAlg::GetU0_non_const() {
  return row_stabilizer_rep.front();
}
RepresentationSystem& MaxFillAlg::GetUId_non_const() { return id; }

const std::vector<RepresentationSystem>& MaxFillAlg::GetStabilizers() const {
  return row_stabilizer_rep;
}

/*********************************************************************************
CALCULATE THE SUPPLEMENTARY MATRICES M, L, AND C FOR THE MAXIMUM FILLING
ALGORITHM
- the element (i,j) of M contains the maximum possible bond degree between atoms
i and j
- the element (i,j) of L corresponds to the maximum possible row capacity after
(i,j)
- the element (i,j) of C corresponds to the maximum possible column capacity
after (i,j)
*********************************************************************************/
void MaxFillAlg::GetMLC(const size_t max_degree) {
  size_t N = M.GetN();

  if (N == 2)  //(see e.g. case {CH3}2)
    M.SetElement(0, 1, degree_vec[1]);

  // fill M by filling the upper triangle for indices (i,j)
  //  - diagonal entries are 0, as there are no bonds between an atom and itself
  //  - if the degree of atoms i and j is the same, Set the value to the minimum
  //  of max_degree and deg -1
  //  - otherwise, Set the value to the minimum of max_degree and the minimal
  //  value between deg(atom i) and deg(atom j)
  else {
    for (size_t ii = 0; ii < N; ii++) {
      for (size_t jj = ii + 1; jj < N; jj++) {
        size_t di(degree_vec[ii]), dj(degree_vec[jj]);

        if (di == dj)
          M.SetElement(ii, jj, (int)std::min(max_degree, di - 1));

        else
          M.SetElement(ii, jj, (int)std::min(max_degree, std::min(di, dj)));
      }
    }
  }

  // fill L by going through indices (i,j)
  for (size_t ii = 0; ii < N; ii++) {
    for (size_t jj = ii + 1; jj < N; jj++)
      L.SetElement(ii, jj,
                   std::min(degree_vec[ii], M.AccumulateRow(ii, jj + 1, N)));
  }

  // fill C by going through indices (i,j)
  for (size_t ii = 0; ii < N; ii++) {
    for (size_t jj = ii + 1; jj < N; jj++)
      C.SetElement(ii, jj,
                   std::min(degree_vec[jj], M.AccumulateColumn(ii + 1, N, jj)));
  }
}

}  // namespace enu

}  // namespace combi_ff