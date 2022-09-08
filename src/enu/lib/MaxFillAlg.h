// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef MAXFILLALG_H
#define MAXFILLALG_H

#include "AdjacencyMatrixEnu.h"

namespace combi_ff {

namespace enu {

class MaxFillAlg {
 public:
  MaxFillAlg(enu::AdjacencyMatrix& A, const bool stereo,
             const size_t max_degree);

  bool GetNextCanonicalMatrix();
  std::pair<int, int> GetHat();
  void SetNextIndex();
  void SetPreviousIndex();
  void ForwardStep();
  void BackwardStep();
  void CanonicityTest();

  bool IsCanonical();
  bool IsCanonical(size_t row);
  void DetermineNewIJ(const size_t i0, const size_t j0, const size_t i_pos,
                      const size_t j_pos, const size_t i_pos_perm,
                      const size_t j_pos_perm);
  bool FindStabilizer(const std::vector<size_t>& permuted_indices,
                      RepresentationSystem& u_next,
                      const size_t highest_permuted_idx);

  const RepresentationSystem& GetUAutomorph() const;
  const RepresentationSystem& GetU0() const;
  const RepresentationSystem& GetUId() const;
  RepresentationSystem& GetUAutomorph_non_const();
  RepresentationSystem& GetU0_non_const();
  RepresentationSystem& GetUId_non_const();
  const std::vector<RepresentationSystem>& GetStabilizers() const;
  void GetMLC(const size_t max_degree);

 private:
  const size_t N;
  const bool stereo;
  size_t i, j;
  enu::AdjacencyMatrix& A;
  std::vector<int> degree_vec;
  SymmetricalMatrix<int> M;
  TriangularMatrix<int> L;
  TriangularMatrix<int> C;
  std::vector<RepresentationSystem> row_stabilizer_rep;
  RepresentationSystem id;
  bool found;
  TriangularMatrix<size_t> accumulated_row_matrix;
  TriangularMatrix<size_t> accumulated_column_matrix;
  PermutationIterator perm_it;
  void (combi_ff::enu::MaxFillAlg::*next_step_ptr)();
};

}  // namespace enu

}  // namespace combi_ff

#endif
