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
      size_p(N - 1),
      permuted_indices(N),
      permutations(Permutations(0)),
      u(NULL),
      p(RepresentationIterator(N - 1, 0)) {}

PermutationIterator::PermutationIterator(
    const combi_ff::RepresentationSystem& u)
    : current_index(u.size() - 1),
      smallest_diff_index(current_index),
      size_u_minus_one(current_index),
      size_p(u.size()),
      permuted_indices(u.size() + 1),
      permutations(Permutations(0)),
      u(&u),
      p(RepresentationIterator(u.size(), 0)) {}

PermutationIterator::PermutationIterator(const PermutationIterator& it_p)
    : current_index(it_p.current_index),
      smallest_diff_index(it_p.smallest_diff_index),
      size_u_minus_one(it_p.size_u_minus_one),
      size_p(it_p.p.size()),
      permuted_indices(it_p.permuted_indices),
      permutations(it_p.permutations),
      u(it_p.u),
      p(it_p.p) {}

const std::vector<size_t>* PermutationIterator::GetPermutedIndices() {
  std::iota(permuted_indices.begin(), permuted_indices.end(), 0);

  for (size_t ii = smallest_diff_index; ii < size_p; ii++) {
    if (p[ii]) {
      for (const Permutation& perm : (*u)[ii][p[ii]])
        std::swap(permuted_indices[perm.first], permuted_indices[perm.second]);
    }
  }

  return &permuted_indices;
}

const std::vector<size_t>* PermutationIterator::GetPermutedIndices(
    size_t& highest_permuted_index) {
  std::iota(permuted_indices.begin(), permuted_indices.end(), 0);
  highest_permuted_index = 0;

  for (size_t ii = smallest_diff_index; ii < size_p; ii++) {
    if (p[ii]) {
      for (const Permutation& perm : (*u)[ii][p[ii]]) {
        std::swap(permuted_indices[perm.first], permuted_indices[perm.second]);
        if (perm.second > highest_permuted_index)
          highest_permuted_index = perm.second;
      }
    }
  }

  return &permuted_indices;
}

bool PermutationIterator::GetNextPermutation() {
  // note: current_index is unsigned -> when current_index =
  // -1, modular arithmetics leads to current_index > p.size
  while (current_index < size_p) {
    if (p[current_index] + 1 < (*u)[current_index].size()) {
      p[current_index]++;
      if (current_index < smallest_diff_index)
        smallest_diff_index = current_index;
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

  for (size_t ii = smallest_diff_index; ii < size_p; ii++) {
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

  u = u_;
  size_u_minus_one = u->size() - 1;
  p.resize(u->size());
  std::fill(p.begin(), p.end(), 0);
  size_p = p.size();
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

}  // namespace combi_ff
