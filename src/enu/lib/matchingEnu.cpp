#include "Range.h"
#include "matchingEnu.h"
#include "Substructure.h"
#include <iomanip>

namespace combi_ff {

namespace enu {


bool FindFragMatches(const enu::AdjacencyMatrix& A,
                     const std::vector<SubstructureCollection>& substructures) {
  std::vector<int> num_matches(0);

  if (!substructures.size())
    return true;

  std::vector<std::vector<bool>> involved_atoms;

  for (auto && F : substructures) {
    num_matches.push_back(0);

    if (F.GetAND()) {
      bool found(false);

      for (auto && M : F.GetSubstructureMatrices()) {
        num_matches.back() = FindFragMatch(A, M, involved_atoms);

        if (IsInRange(num_matches.back(), F.GetRange())) {
          found = true;
          break;
        }
      }

      if (!found)
        return false;

    } else if (F.GetXOR()) {
      size_t num_found(0);

      for (auto && M : F.GetSubstructureMatrices()) {
        if (FindFragMatch(A, M, involved_atoms))
          num_found++;
      }

      if (!IsInRange(num_found, F.GetRange()))
        return false;

    } else {
      for (auto && M : F.GetSubstructureMatrices())
        num_matches.back() += FindFragMatch(A, M, involved_atoms);

      if (!IsInRange(num_matches.back(), F.GetRange()))
        return false;
    }
  }

  //for(size_t i = 0; i < substructures.size(); i++)
  //std::cout << "number of matches for " << FFragments.at(i).GetCode() << " is "  << num_matches[i] << " (required range is " << *(FFragments.at(i).GetRange()) << " ) " <<'\n';
  return true;
}


int FindFragMatch(const enu::AdjacencyMatrix& A, const FragmentMatrix& fragment_matrix,
                  std::vector<std::vector<bool>>& involved_atoms) {
  int num_matches(0);
  const AtomVector<combi_ff::Atom>& atoms = A.GetAtomVector();
  const AtomVector<combi_ff::Atom>& fragAtoms = fragment_matrix.GetAtomVector();
  const size_t m = A.GetN();
  const size_t n = fragment_matrix.GetN();
  ComparisonMatrix M(n, m);

  if (n > m)
    return false;

  for (size_t i = 0; i < n; i ++) {
    bool potential_match(false);

    for (size_t j = 0; j < m; j++) {
      if ((fragAtoms[i].GetElementSymbol() == atoms[j].GetElementSymbol() &&
           fragAtoms[i].GetNumConnections() <=
           atoms[j].GetNumConnections()/*fragAtoms[i].GetNumHydrogens() + fragAtoms[i].GetNumFixedHydrogens() == atoms[j].GetNumHydrogens() + atoms[j].GetNumFixedHydrogens()*/)
          || fragAtoms[i].GetUnitedAtomSymbol() == "*") {
        potential_match = true;
        M.SetElement(i, j, true);
      }
    }

    if (!potential_match)
      return false;
  }

  std::vector<bool> matched_cols(m, false);
  std::vector<bool> matched_rows(n, false);
  int k = -1;
  ComparisonMatrix M_save = M;
  std::vector<bool> matched_rows_save(n);
  std::vector<bool> matched_cols_save(m);
  UllmannMatch(M, matched_cols, matched_rows, k, n, m, A, fragment_matrix, num_matches,
               involved_atoms);
  return num_matches;
}


bool UllmannMatch(ComparisonMatrix& M, std::vector<bool>& matched_cols,
                  std::vector<bool>& matched_rows,
                  int k, const size_t n, const size_t m, const enu::AdjacencyMatrix& A,
                  const FragmentMatrix& fragment_matrix,
                  int& num_matches, std::vector<std::vector<bool>>& involved_atoms) {
  if (k == (int)n - 1) {
    for (size_t i = 0; i < n; i++) {
      if (M.AccumulateRow(i) != 1)
        return false;
    }

    for (size_t j = 0; j < m; j++) {
      if (M.AccumulateColumn(j) > 1)
        return false;
    }

    std::vector<bool> involved_atoms_local(m, false);

    for (size_t j = 0; j < m; j++) {
      for (size_t i = 0; i < n; i++) {
        if (M.GetElement(i, j))
          involved_atoms_local[j] = true;
      }
    }

    for (auto && atomList : involved_atoms) {
      int nO(0); //number of overlaps

      for (size_t j = 0; j < m; j++) {
        if (atomList[j] && involved_atoms_local[j])
          nO++;

        //one overlap is allowed, i.e. two same involved atoms
        if (nO > 2) {
          //if(nO >= fragMat.GetN()) {
          //std::cout << "found again:\n" << M << '\n';
          return false;
        }
      }
    }

    //return true;
    //std::cout << "found M:\n" << M << '\n';
    involved_atoms.push_back(involved_atoms_local);
    num_matches ++;
    return false;
  }

  ComparisonMatrix M_save(n, m);
  std::vector<bool> matched_rows_save(n);
  std::vector<bool> matched_cols_save(m);

  for (size_t l = 0; l < m; l++) {
    if (M.GetElement(k + 1, l) == true && matched_cols[l] == false && matched_rows[k + 1] == false) {
      M_save = M;
      matched_rows_save = matched_rows;
      matched_cols_save = matched_cols;

      for (size_t j = 0; j < m; j++)
        M.SetElement(k + 1, j, false);

      for (size_t i = 0; i < n; i++)
        M.SetElement(i, l, false);

      M.SetElement(k + 1, l, true);
      matched_cols[l] = true;
      matched_rows[k + 1] = true;

      if (Refine(M, k + 1, n, m, A, fragment_matrix)) {
        if (UllmannMatch(M, matched_cols, matched_rows, k + 1, n, m, A, fragment_matrix, num_matches,
                         involved_atoms))
          return true;
      }

      M = M_save;
      matched_rows = matched_rows_save;
      matched_cols = matched_cols_save;
    }
  }

  return false;
}




bool FindBenzMatch(enu::AdjacencyMatrix& A, bool& canonical, const RepresentationSystem& u0) {
  //AtomVector a = Atoms({"*", "*", "*", "*", "*", "*"});
  FragmentMatrix benzene(6,
                         AdjacencyVector({0, 2, 0, 0, 0, 1,
                                          2, 0, 1, 0, 0, 0,
                                          0, 1, 0, 2, 0, 0,
                                          0, 0, 2, 0, 1, 0,
                                          0, 0, 0, 1, 0, 2,
                                          1, 0, 0, 0, 2, 0
                                         }));
  const AtomVector<combi_ff::Atom>& atoms = A.GetAtomVector();
  const FragmentMatrix& benzene_matrix = benzene;
  const AtomVector<combi_ff::Atom>& fragAtoms = benzene_matrix.GetAtomVector();
  const size_t m = A.GetN();
  const size_t n = benzene_matrix.GetN();
  ComparisonMatrix M(n, m);
  std::vector<bool>& is_aromatic = A.GetIsAromaticCarbon();
  std::fill(is_aromatic.begin(), is_aromatic.end(), 0);

  if (n > m)
    return false;

  bool foundPotentialMatchForFragmentAtom(false);

  for (size_t i = 0; i < n; i ++) {
    for (size_t j = 0; j < m; j++) {
      if ((!isdigit(fragAtoms[i].GetUnitedAtomSymbol().back()) &&
           atoms[j].GetElementSymbol() == fragAtoms[i].GetElementSymbol())
          || atoms[j].GetUnitedAtomSymbol() == fragAtoms[i].GetUnitedAtomSymbol()
          || (fragAtoms[i].GetUnitedAtomSymbol() == "*")) {
        foundPotentialMatchForFragmentAtom = true;
        M.SetElement(i, j, true);
      }
    }

    if (!foundPotentialMatchForFragmentAtom)
      return false;
  }

  std::vector<bool> matched_cols(m, false);
  std::vector<bool> matched_rows(n, false);
  int k = -1;
  ComparisonMatrix M_save = M;
  std::vector<bool> matched_rows_save(n);
  std::vector<bool> matched_cols_save(m);
  bool found = false;
  size_t num_matches = 0;
  UllmannMatchBenz(M, matched_cols, matched_rows, k, n, m, A, benzene_matrix, found, num_matches,
                   is_aromatic, canonical, u0);
  return found;
}


bool UllmannMatchBenz(ComparisonMatrix& M,
                      std::vector<bool>& matched_cols,
                      std::vector<bool>& matched_rows,
                      int k,
                      const size_t n,
                      const size_t m,
                      const enu::AdjacencyMatrix& A,
                      const FragmentMatrix& fragment_matrix,
                      bool& found,
                      size_t& num_matches,
                      std::vector<bool>& is_aromatic,
                      bool& canonical,
                      const RepresentationSystem& u0) {
  if (k == (int)n - 1) {
    for (size_t i = 0; i < n; i++) {
      if (M.AccumulateRow(i) != 1)
        return false;
    }

    for (size_t i = 0; i < m; i++) {
      if (M.AccumulateColumn(i) > 1)
        return false;
    }

    //return true;
    std::vector<bool> is_aromatic_save = is_aromatic;

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < m ; j++) {
        if (M.GetElement(i, j))
          (is_aromatic_save)[j] = true;
      }
    }

    size_t num_with_2_plus_aromatic_neighbors = 0;

    for (size_t i = 0; i < m; i++) {
      if ((is_aromatic_save)[i] == true) {
        size_t num_aro_nbr = 0;

        for (auto && nbr : A.GetAtomVector()[i].GetNeighbours()) {
          if (is_aromatic_save[nbr])
            num_aro_nbr++;
        }

        if (num_aro_nbr > 2)
          num_with_2_plus_aromatic_neighbors++;
      }
    }

    if (num_with_2_plus_aromatic_neighbors == 0)
      is_aromatic = is_aromatic_save;

    else
      return false;

    num_matches ++;
    found = true;
    std::vector<bool> in_benzene_ring(A.GetN(), false);

    for (size_t j = 0; j < m; j++) {
      for (size_t i = 0; i < n; i++) {
        if (M.GetElement(i, j))
          in_benzene_ring[j] = true;
      }
    }

    enu::AdjacencyMatrix tmp(A);

    for (size_t i = 0; i < A.GetN(); i++) {
      for (size_t j = 0; j < A.GetN(); j++) {
        if (in_benzene_ring[i] && in_benzene_ring[j]) {
          if (A.GetElement(i, j) == 1)
            tmp.SetElement(i, j, 2);

          else if (A.GetElement(i, j) == 2)
            tmp.SetElement(i, j, 1);
        }
      }
    }

    tmp.MakeCanonical(u0);

    if (tmp > A)
      canonical = false;

    return false;
  }

  ComparisonMatrix M_save(n, m);
  std::vector<bool> matched_rows_save(n);
  std::vector<bool> matched_cols_save(m);

  for (size_t l = 0; l < m; l++) {
    if (M.GetElement(k + 1, l) == true && matched_cols[l] == false && matched_rows[k + 1] == false) {
      M_save = M;
      matched_rows_save = matched_rows;
      matched_cols_save = matched_cols;

      for (size_t j = 0; j < m; j++)
        M.SetElement(k + 1, j, false);

      for (size_t i = 0; i < n; i++)
        M.SetElement(i, l, false);

      M.SetElement(k + 1, l, true);
      matched_cols[l] = true;
      matched_rows[k + 1] = true;

      if (Refine(M, k + 1, n, m, A, fragment_matrix)) {
        if (canonical &&
            UllmannMatchBenz(M, matched_cols, matched_rows, k + 1, n, m, A, fragment_matrix, found, num_matches,
                             is_aromatic, canonical, u0))
          return true;
      }

      M = M_save;
      matched_rows = matched_rows_save;
      matched_cols = matched_cols_save;
    }
  }

  return false;
}

bool Refine(combi_ff::ComparisonMatrix& M, int k, const size_t n,
            const size_t m, const enu::AdjacencyMatrix& A,
            const combi_ff::FragmentMatrix& fragment_matrix) {
  bool changed(true), valid(false), found(false);

  while (changed) {
    changed = false;

    for (size_t i = k + 1; i < n; i++) {
      for (size_t j = 0; j < m; j++) {
        if (M.GetElement(i, j)) {
          valid = true;

          for (auto && i2 : fragment_matrix.GetAtomVector()[i].GetNeighbours()) {
            found = false;

            for (auto && j2 :  A.GetAtomVector()[j].GetNeighbours()) {
              if (M.GetElement(i2, j2)) {
                if ((fragment_matrix.GetElement(i, i2) == A.GetElement(j, j2))) {
                  found = true;
                  break;
                }
              }
            }

            if (!found) {
              valid = false;
              break;
            }
          }

          if (!valid) {
            M.SetElement(i, j, false);
            changed = true;
            bool nonzero = false;

            for (size_t h = 0; h < m; h++) {
              if (M.GetElement(i, h))
                nonzero = true;
            }

            if (!nonzero)
              return false;
          }
        }
      }
    }
  }

  return true;
}

} //namespace enu

} //namespace combi_ff