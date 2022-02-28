// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef CONTAINEROPERATORS_H
#define CONTAINEROPERATORS_H

#include <fstream>
#include <list>
#include <stdexcept>
#include <vector>

namespace combi_ff {

template <class T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& v) {
  if (v.size()) {
    for (size_t i = 0; i < v.size() - 1; i++) stream << v[i] << " ";

    stream << v.back();
  }

  return stream;
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const std::list<T>& list) {
  for (auto&& list_element : list) stream << list_element << " ";

  return stream;
}

template <class T, class R>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<R>& b) {
  if (a.size() != b.size())
    throw std::logic_error("can't subtract vectors of different size");

  std::vector<T> c(a.size());

  for (size_t i = 0; i < a.size(); i++) c[i] = a[i] - b[i];

  return c;
}

template <class T>
bool operator<(const std::vector<T>& a, const std::vector<T>& b) {
  if (a.size() != b.size())
    throw std::logic_error("can't compare vectors of different size");

  for (size_t i = 0; i < a.size(); i++) {
    if (a[i] < b[i])
      return true;

    else if (a[i] > b[i])
      return false;
  }

  return false;
}

}  // namespace combi_ff

#endif
