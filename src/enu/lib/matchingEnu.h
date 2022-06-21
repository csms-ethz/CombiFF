// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef MATCHINGENU_H_
#define MATCHINGENU_H_

#include "AdjacencyMatrixEnu.h"
#include "Permutation.h"

namespace combi_ff {

namespace enu {

class SubstructureCollection;

bool FindBenzMatch(enu::AdjacencyMatrix& A, bool& canonical,
                   const RepresentationSystem& u0);

bool UllmannMatchBenz(ComparisonMatrix& M,
                      std::vector<bool>& matched_molecule_atoms, int k,
                      const size_t n, const size_t m,
                      const enu::AdjacencyMatrix& A,
                      const FragmentMatrix& fragment_matrix, bool& found,
                      size_t& num_matches,
                      std::vector<bool>& is_aromatic_carbon, bool& canonical,
                      const RepresentationSystem& u0);

bool FindFragMatches(const enu::AdjacencyMatrix& A,
                     const std::vector<SubstructureCollection>& substructures);

int FindFragMatch(const enu::AdjacencyMatrix& A, const FragmentMatrix& frag,
                  std::vector<std::vector<bool>>& involved_atoms);

bool UllmannMatch(ComparisonMatrix& M,
                  std::vector<bool>& matched_molecule_atoms, int k,
                  const size_t n, const size_t m, const enu::AdjacencyMatrix& A,
                  const FragmentMatrix& fragment_matrix, int& num_matches,
                  std::vector<std::vector<bool>>& involved_atoms);

bool Refine(combi_ff::ComparisonMatrix& M, int k, const size_t n,
            const size_t m, const enu::AdjacencyMatrix& A,
            const combi_ff::FragmentMatrix& fragment_matrix);

}  // namespace enu

}  // namespace combi_ff

#endif