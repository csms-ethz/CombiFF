// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef PERMUTATION_H
#define PERMUTATION_H

#include <cassert>
#include <fstream>
#include <utility>
#include <vector>

namespace combi_ff {

typedef std::pair<size_t, size_t> Permutation;
typedef std::vector<Permutation> Permutations;
typedef std::vector<Permutations> StabilizerVector;
typedef std::vector<StabilizerVector> RepresentationSystem;
typedef std::vector<size_t> RepresentationIterator;

class PermutationIterator {
 public:
  /* constructor */
  explicit PermutationIterator(const combi_ff::RepresentationSystem& u);
  explicit PermutationIterator(const size_t N);
  PermutationIterator(const PermutationIterator& it_p);

  /* getter functions */
  const std::vector<size_t>* GetPermutedIndices();
  const std::vector<size_t>* GetPermutedIndices(size_t& highest_permuted_index);
  bool GetNextPermutation();
  const Permutations* GetCombinedPermutation();
  size_t GetSmallestDiffIndex() const;

  /* setter functions */
  void SetCurrentIndex(const size_t current_new);
  size_t GetCurrentIndex() const;
  void SetCurrentIndexToSmallestDiffIndex();
  void Reset();
  void Reset(const combi_ff::RepresentationSystem* u_,
             const size_t current_new);

 private:
  size_t current_index;
  size_t smallest_diff_index;
  size_t size_u_minus_one;
  size_t size_p;
  std::vector<size_t> permuted_indices;
  Permutations permutations;
  const combi_ff::RepresentationSystem* u;
  RepresentationIterator p;
};

bool operator==(const Permutation& a, const Permutation& b);
bool operator!=(const Permutation& a, const Permutation& b);
bool operator>=(const Permutation& a, const Permutation& b);
bool operator<=(const Permutation& a, const Permutation& b);
std::ostream& operator<<(std::ostream& stream, const Permutations& a);
std::ostream& operator<<(std::ostream& stream, const Permutation& a);

size_t Min(const Permutation& a);
size_t Max(const Permutation& a);
Permutation& Max(Permutation& a, Permutation& b);
Permutation& Max(Permutation& a, Permutation& b, Permutation& c);

template <typename T>
void PermuteVector(std::vector<T>& index, const Permutations& permutations) {
  for (auto&& p : permutations) {
    // idx.swap(p.first, p.second);
    std::swap(index[p.first], index[p.second]);
  }
}

}  // namespace combi_ff

#endif
