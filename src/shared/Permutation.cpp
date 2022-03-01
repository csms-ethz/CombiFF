// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "Permutation.h"

#include <iostream>
#include <numeric>

#include "ContainerOperators.h"

namespace combi_ff {

PermutationIterator::PermutationIterator(const size_t N)
    : current_index(N - 2),
      smallest_diff_index(current_index),
      size_u_minus_one(current_index),
      permuted_indices(N),
      permutations(Permutations(0)),
      u(NULL),
      p(RepresentationIterator(N - 1, 0)) {}

PermutationIterator::PermutationIterator(
    const combi_ff::RepresentationSystem& u)
    : current_index(u.size() - 1),
      smallest_diff_index(current_index),
      size_u_minus_one(current_index),
      permuted_indices(u.size() + 1),
      permutations(Permutations(0)),
      u(&u),
      p(RepresentationIterator(u.size(), 0)) {}

PermutationIterator::PermutationIterator(const PermutationIterator& it_p)
    : current_index(it_p.current_index),
      smallest_diff_index(it_p.smallest_diff_index),
      size_u_minus_one(it_p.size_u_minus_one),
      permuted_indices(it_p.permuted_indices),
      permutations(it_p.permutations),
      u(it_p.u),
      p(it_p.p) {}

const std::vector<size_t>* PermutationIterator::GetPermutedIndices() {
  std::iota(permuted_indices.begin(), permuted_indices.end(), 0);

  for (size_t ii = smallest_diff_index; ii < p.size(); ii++) {
    if (p[ii]) {
      for (auto&& perm : (*u)[ii][p[ii]])
        std::swap(permuted_indices[perm.first], permuted_indices[perm.second]);
    }
  }

  return &permuted_indices;
}

const std::vector<size_t>* PermutationIterator::GetPermutedIndices(
    size_t& highest_permuted_index) {
  std::iota(permuted_indices.begin(), permuted_indices.end(), 0);
  size_t hi = 0;

  for (size_t ii = smallest_diff_index; ii < p.size(); ii++) {
    if (p[ii]) {
      for (auto&& perm : (*u)[ii][p[ii]]) {
        std::swap(permuted_indices[perm.first], permuted_indices[perm.second]);
        hi = std::max(hi, perm.second);
      }
    }
  }

  highest_permuted_index = hi;
  return &permuted_indices;
}

bool PermutationIterator::GetNextPermutation() {
  while (current_index <
         p.size()) {  // note: current_index is unsigned -> when current_index =
                      // -1, modular arithmetics leads to current_index > p.size
    if (p[current_index] + 1 < (*u)[current_index].size()) {
      p[current_index]++;
      smallest_diff_index = std::min(smallest_diff_index, current_index);
      current_index = size_u_minus_one;
      return true;

    } else {
      p[current_index] = 0;
      current_index--;
    }
  }

  return false;
}

const Permutations* PermutationIterator::GetCombinedPermutation() {
  permutations.clear();

  for (size_t ii = smallest_diff_index; ii < p.size(); ii++) {
    if (p[ii]) {
      permutations.insert(permutations.end(), (*u)[ii][p[ii]].begin(),
                          (*u)[ii][p[ii]].end());
    }
  }

  return &permutations;
}

size_t PermutationIterator::GetSmallestDiffIndex() const {
  return smallest_diff_index;
}

size_t PermutationIterator::GetCurrentIndex() const { return current_index; }

void PermutationIterator::SetCurrentIndex(const size_t current_new) {
  current_index = current_new;
  std::fill(p.begin() + current_index + 1, p.end(), 0);
}

void PermutationIterator::SetCurrentIndexToSmallestDiffIndex() {
  current_index = smallest_diff_index;
  std::fill(p.begin() + current_index + 1, p.end(), 0);
}

void PermutationIterator::Reset() {
  current_index = u->size() - 1;
  smallest_diff_index = current_index;
  std::fill(p.begin(), p.end(), 0);
}

void PermutationIterator::Reset(const combi_ff::RepresentationSystem* u_,
                                const size_t current_new) {
  current_index = current_new;
  smallest_diff_index = current_index;
  std::fill(p.begin(), p.end(), 0);
  u = u_;
  size_u_minus_one = u->size() - 1;
}

bool operator==(const Permutation& a, const Permutation& b) {
  if (a.first == b.first && a.second == b.second) return true;

  return false;
}

bool operator!=(const Permutation& a, const Permutation& b) {
  return !(a == b);
}

bool operator>=(const Permutation& a, const Permutation& b) {
  // if(min(a) > min(b) || (min(a) == min(b) && max(a) >= max(b)))
  if (a.first > b.first || (a.first == b.first && a.second >= b.second))
    return true;

  else
    return false;
}

bool operator<=(const Permutation& a, const Permutation& b) {
  // if(min(a) < min(b) || (min(a) == min(b) && max(a) <= max(b)))
  if (a.first < b.first || (a.first == b.first && a.second <= b.second))
    return true;

  else
    return false;
}

std::ostream& operator<<(std::ostream& stream, const Permutations& a) {
  for (size_t i = 0; i < a.size(); i++)
    stream << "(" << a[i].first << "," << a[i].second << ")";

  return stream;
}

std::ostream& operator<<(std::ostream& stream, const Permutation& a) {
  stream << "(" << a.first << "," << a.second << ")";
  return stream;
}

size_t Min(const Permutation& a) { return std::min(a.first, a.second); }
size_t Max(const Permutation& a) { return std::max(a.first, a.second); }

Permutation& Max(Permutation& a, Permutation& b) {
  if (a >= b)
    return a;

  else
    return b;
}

Permutation& Max(Permutation& a, Permutation& b, Permutation& c) {
  return Max(a, Max(b, c));
}

size_t NumDiff(const std::vector<size_t>& original,
               const std::vector<size_t>& permutated) {
  size_t diff(0);

  for (size_t i = 0; i < original.size(); i++) {
    if (original[i] != permutated[i]) diff++;
  }

  return diff;
}

size_t NumPerm(const std::vector<size_t>& original,
               const std::vector<size_t>& permutated) {
  const size_t diff = NumDiff(original, permutated);

  switch (diff) {
    // e.g. 1,0,2,3 and 1,0,2,3
    case (0):
      // std::cout << 0 << std::endl;
      return 0;

    // e.g. 1,0,2,3 and 0,1,2,3
    case (2):
      // std::cout << 1 << std::endl;
      return 1;

    // e.g. 1,0,2,3 and 0,2,1,3
    case (3):
      // std::cout << 2 << std::endl;
      return 2;

    // not clear yet how many swaps -> e.g. 1,0,2,3 to 3,2,0,1 are only 2, but
    // to 0,2,3,1 are 3
    case (4): {
      for (size_t i = 0; i < 4; i++) {
        if (permutated[i] == original[0]) {
          // if zeroth value of original is swapped with the same number in both
          // vectors, even permutations (2 swaps), e.g. 1,0,2,3 and 3,2,0,1
          // note: doesn't necessarily have to be element zero, can also be any
          // of the other numbers
          if (permutated[0] == original[i]) {
            // std::cout << 2 << std::endl;
            return 2;
          }

          // if zeroth value of permutated is swapped with different number, odd
          // permutations (3 swaps), e.g. 1,0,2,3 and 0,2,3,1 note: doesn't
          // necessarily have to be element zero, can also be any of the other
          // numbers

          else {
            // std::cout << 3 << std::endl;
            return 3;
          }
        }
      }

      throw std::runtime_error("error in StereoGenerator::NumPerm for case 4");
    }

    default:
      throw std::runtime_error("diff should be 0, 2, 3, or 4, but it is " +
                               std::to_string(diff));
  }
}

void NeighborOrder(const size_t idx, const size_t idx_permuted,
                   const std::vector<size_t>& permuted_indices,
                   std::vector<size_t>& nbrs_original_order,
                   std::vector<size_t>& nbrs_permuted_order,
                   const std::vector<size_t>& coming_from,
                   const std::vector<std::vector<size_t>>& going_to,
                   const std::vector<std::vector<size_t>>& ring_connections) {
  nbrs_original_order.resize(1);
  nbrs_original_order[0] = coming_from[idx];
  nbrs_original_order.insert(nbrs_original_order.end(),
                             ring_connections[idx].begin(),
                             ring_connections[idx].end());
  nbrs_original_order.insert(nbrs_original_order.end(), going_to[idx].begin(),
                             going_to[idx].end());
  nbrs_permuted_order.resize(1);
  nbrs_permuted_order[0] =
      std::distance(permuted_indices.begin(),
                    std::find(permuted_indices.begin(), permuted_indices.end(),
                              coming_from[idx_permuted]));

  for (const auto& rc : ring_connections[idx_permuted])
    nbrs_permuted_order.push_back(std::distance(
        permuted_indices.begin(),
        std::find(permuted_indices.begin(), permuted_indices.end(), rc)));

  for (const auto& gt : going_to[idx_permuted])
    nbrs_permuted_order.push_back(std::distance(
        permuted_indices.begin(),
        std::find(permuted_indices.begin(), permuted_indices.end(), gt)));
}

}  // namespace combi_ff
