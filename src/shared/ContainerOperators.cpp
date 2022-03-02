// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "ContainerOperators.h"

namespace combi_ff {

std::vector<int>& operator++(std::vector<int>& b) {
  auto&& it = std::find(b.rbegin(), b.rend(), 0);

  if (it != b.rend()) *it = 1;

  std::fill(b.rbegin(), it, 0);
  /*for (int i = (int)b.size() - 1; i >= 0; i--) {
    if (b[i] == 0) {
      b[i] = 1;
      return b;

    } else if (b[i] == 1)
      b[i] = 0;
  }*/
  return b;
}

}  // namespace combi_ff