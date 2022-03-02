// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "StereoGenerator.h"

namespace combi_ff {

namespace enu {

StereoGenerator::StereoGenerator(const enu::AdjacencyMatrix& A,
                                 const size_t n_rings,
                                 const RepresentationSystem& u_automorph,
                                 const RepresentationSystem& u_id,
                                 const SmilesGeneratorEnu& smiles_gen)
    : A(A),
      N(A.GetN()),
      u_automorph(u_automorph),
      u_id(u_id),
      num_stereo_centers(0),
      num_ct_bonds(0),
      visited_indices(smiles_gen.GetVisitedIndices()),
      smiles_blocks(smiles_gen.GetSmilesBlocks()),
      coming_from(smiles_gen.GetComingFrom()),
      going_to(smiles_gen.GetGoingTo()),
      ring_connections(smiles_gen.GetRingConnections()),
      is_in_cycle(A.GetN(), false),
      potential_stereo_idx(std::vector<size_t>(0)),
      true_stereo(std::vector<bool>(N, false)),
      potential_cis_trans_stereo(std::vector<bool>(N, false)),
      potential_cis_trans_stereo_idx(std::vector<size_t>(0)),
      true_cis_trans_stereo_idx(std::vector<size_t>(0)),
      double_bond_partner(std::vector<int>(N, -1)),
      single_bond_partner(std::vector<int>(N, -1)),
      true_stereo_idx(std::vector<size_t>(0)),
      potential_para_stereo_idx(std::vector<size_t>(0)),
      potential_para_cis_trans_stereo_idx(std::vector<size_t>(0)),
      true_cis_trans_pairs(std::vector<CisTransPair>(0)),
      para_cis_trans_pairs(std::vector<CisTransPair>(0)),
      position_in_true_stereo_idx(std::vector<int>(N, -1)),
      true_cis_trans_stereo(std::vector<bool>(N, false)),
      visited_idx(std::vector<size_t>(N)),
      stereo_smiles(
          std::vector<std::tuple<std::string, int, std::pair<int, int>>>(0)) {
  for (size_t i = 0; i < visited_indices.size(); i++)
    visited_idx[visited_indices[i]] =
        i;  // visited_idx is inverse of visited_indices

  if (n_rings) is_in_cycle = A.IsInCycle();
}

std::vector<std::tuple<std::string, int, std::pair<int, int>>>&
StereoGenerator::GetStereoSmiles() {
  return stereo_smiles;
}

void StereoGenerator::GenerateStereoSmiles() {
  FindPotentialStereoCenters();

  if (!(potential_stereo_idx.size() || true_stereo_idx.size() ||
        potential_cis_trans_stereo_idx.size()))
    return;

  FindTrueStereoCenters();

  if (!(true_stereo_idx.size() || true_cis_trans_stereo_idx.size())) return;

  FindPotentialParaStereoCenters();
  std::vector<Configuration> valid_configurations;
  valid_configurations.reserve(true_stereo_idx.size() *
                               true_cis_trans_stereo_idx.size());
  FindValidConfigurations(valid_configurations);
  WriteStereoSmiles(valid_configurations);
}

void StereoGenerator::FindPotentialStereoCenters() {
  const Atom* a;

  for (size_t i = 0; i < N; i++) {
    a = &A.GetAtomVector()[i];
    size_t n_neighbors = a->GetNeighbors().size();

    // test if an atom is a potential tetrahedral stereo center
    if (n_neighbors == 4 ||
        (n_neighbors == 3 &&
         (a->GetNumFixedHydrogens() + a->GetNumHydrogens() == 1)))
      TestIfPotentialTetrahedralStereoCenter(*a, i);

    // test if an atom is a potential cis trans stereo center
    else if ((n_neighbors == 3 ||
              (n_neighbors == 2 &&
               (a->GetNumFixedHydrogens() + a->GetNumHydrogens() == 1))) &&
             !is_in_cycle[i])
      TestIfPotentialCisTransStereoCenter(*a, i);

    // atoms are ordered by bond degree, so as soon as an atom in the atom vec
    // has a bond degree of less than 2, all the following atoms will also have
    // a bond degree of less than 2 and we can stop
    else if (n_neighbors < 2)
      break;
  }
}

void StereoGenerator::TestIfPotentialTetrahedralStereoCenter(const Atom& a,
                                                             const size_t idx) {
  bool same_first_nbrs(false);
  bool same_first_nbrs_single_bond(false);
  const Atom *nbr1, *nbr2;

  for (size_t i = 0; i < a.GetNeighbors().size(); i++) {
    const size_t n1 = a.GetNeighbors()[i];
    nbr1 = &A.GetAtomVector()[n1];

    for (size_t j = i + 1; j < a.GetNeighbors().size(); j++) {
      const size_t n2 = a.GetNeighbors()[j];
      nbr2 = &A.GetAtomVector()[n2];

      if (nbr1->GetUnitedAtomSymbol() == nbr2->GetUnitedAtomSymbol() &&
          A.GetElement(idx, n1) == A.GetElement(idx, n2) &&
          nbr1->GetNumNeighbors() == nbr2->GetNumNeighbors() &&
          (nbr1->GetNumFixedHydrogens() + nbr1->GetNumHydrogens() ==
           nbr2->GetNumFixedHydrogens() + nbr2->GetNumHydrogens())) {
        same_first_nbrs = true;

        if (nbr1->GetDegree() == 1) {
          same_first_nbrs_single_bond = true;
          break;
        }
      }
    }

    if (same_first_nbrs_single_bond) break;
  }

  if (!same_first_nbrs_single_bond) {
    if (same_first_nbrs)
      potential_stereo_idx.push_back(idx);

    else {
      true_stereo[idx] = true;
      true_stereo_idx.push_back(idx);
    }
  }
}

void StereoGenerator::TestIfPotentialCisTransStereoCenter(const Atom& a,
                                                          const size_t idx) {
  bool has_double_bond(false);
  size_t double_bond_nbr(0);
  bool same_first_nbrs(false);
  bool same_first_nbrs_single_bond(false);

  for (auto&& nbr : a.GetNeighbors()) {
    if (A.GetElement(idx, nbr) == 2) {
      has_double_bond = true;
      double_bond_nbr = nbr;
      double_bond_partner[idx] = (int)nbr;

    } else if (A.GetElement(idx, nbr) == 1 &&
               (single_bond_partner[idx] == -1 ||
                visited_idx[nbr] < visited_idx[single_bond_partner[idx]]))
      single_bond_partner[idx] = (int)nbr;
  }

  if (has_double_bond) {
    const Atom *nbr1, *nbr2;

    for (size_t i = 0; i < a.GetNeighbors().size(); i++) {
      size_t n1 = a.GetNeighbors()[i];

      if (n1 != double_bond_nbr) {
        nbr1 = &A.GetAtomVector()[n1];

        for (size_t j = i + 1; j < a.GetNeighbors().size(); j++) {
          size_t n2 = a.GetNeighbors()[j];

          if (n2 != double_bond_nbr) {
            nbr2 = &A.GetAtomVector()[n2];

            if (nbr1->GetUnitedAtomSymbol() == nbr2->GetUnitedAtomSymbol() &&
                A.GetElement(idx, n1) == A.GetElement(idx, n2) &&
                nbr1->GetNumNeighbors() == nbr2->GetNumNeighbors() &&
                (nbr1->GetNumFixedHydrogens() + nbr1->GetNumHydrogens() ==
                 nbr2->GetNumFixedHydrogens() + nbr2->GetNumHydrogens())) {
              same_first_nbrs = true;

              if (nbr1->GetDegree() == 1) {
                same_first_nbrs_single_bond = true;
                break;
              }
            }
          }
        }
      }

      // not stereo
      if (same_first_nbrs_single_bond) break;
    }

    if (same_first_nbrs) {
      if (!same_first_nbrs_single_bond) {
        potential_cis_trans_stereo[idx] = 1;
        potential_cis_trans_stereo_idx.push_back(idx);
      }

    } else {
      potential_cis_trans_stereo[idx] = 1;
      potential_cis_trans_stereo_idx.push_back(idx);
    }
  }
}

void StereoGenerator::FindTrueStereoCenters() {
  PermutationIterator perm_it(u_automorph);
  const Atom* a;
  const std::vector<size_t>* permuted_indices;

  while (
      perm_it.GetNextPermutation() &&
      (potential_stereo_idx.size() || potential_cis_trans_stereo_idx.size())) {
    permuted_indices = perm_it.GetPermutedIndices();

    for (size_t i = 0; i < potential_stereo_idx.size(); i++) {
      size_t idx = potential_stereo_idx[i];
      a = &A.GetAtomVector()[idx];

      // only look at the nbr permutations if the permutation leaves the current
      // atom identical
      if ((*permuted_indices)[idx] == idx) {
        if (NeighborsArePermuted(*a, *permuted_indices)) {
          true_stereo[idx] = false;
          potential_stereo_idx.erase(potential_stereo_idx.begin() + i);
          potential_para_stereo_idx.push_back(idx);
          i--;
        }
      }
    }

    for (size_t i = 0; i < potential_cis_trans_stereo_idx.size(); i++) {
      size_t idx = potential_cis_trans_stereo_idx[i];
      a = &A.GetAtomVector()[idx];

      if ((*permuted_indices)[idx] == idx) {
        if (NeighborsArePermuted(*a, *permuted_indices)) {
          potential_cis_trans_stereo[idx] = false;
          potential_cis_trans_stereo_idx.erase(
              potential_cis_trans_stereo_idx.begin() + i);
          potential_para_cis_trans_stereo_idx.push_back(idx);
          i--;
        }
      }
    }
  }

  for (auto&& idx : potential_stereo_idx) {
    true_stereo[idx] = true;
    true_stereo_idx.push_back(idx);
  }

  /*
  remove from potentical cis trans stereo and true stereo, if the neighbor is
  not a true or potential  cis trans stereo
  */
  bool removed_ct_stereo(false);

  do {
    removed_ct_stereo = false;

    for (size_t i = 0; i < potential_cis_trans_stereo_idx.size(); i++) {
      size_t idx = potential_cis_trans_stereo_idx[i];
      size_t nbr = double_bond_partner[idx];

      if (!potential_cis_trans_stereo[nbr]) {
        potential_cis_trans_stereo[idx] = false;
        potential_cis_trans_stereo_idx.erase(
            potential_cis_trans_stereo_idx.begin() + i);
        potential_para_cis_trans_stereo_idx.push_back(idx);
        i--;
        removed_ct_stereo = true;
      }
    }
  } while (removed_ct_stereo && potential_cis_trans_stereo_idx.size());

  for (auto&& idx : potential_cis_trans_stereo_idx) {
    true_cis_trans_stereo[idx] = true;
    true_cis_trans_stereo_idx.push_back(idx);
  }

  std::sort(true_stereo_idx.begin(), true_stereo_idx.end());

  for (size_t i = 0; i < true_stereo_idx.size(); i++)
    position_in_true_stereo_idx[true_stereo_idx[i]] = (int)i;

  std::sort(true_cis_trans_stereo_idx.begin(), true_cis_trans_stereo_idx.end());

  for (size_t i = 0; i < true_cis_trans_stereo_idx.size(); i++) {
    if (double_bond_partner[true_cis_trans_stereo_idx[i]] >
        (int)true_cis_trans_stereo_idx[i]) {
      CisTransPair ct(true_cis_trans_stereo_idx[i],
                      double_bond_partner[true_cis_trans_stereo_idx[i]],
                      single_bond_partner[true_cis_trans_stereo_idx[i]],
                      single_bond_partner
                          [double_bond_partner[true_cis_trans_stereo_idx[i]]]);

      if (visited_idx[true_cis_trans_stereo_idx[i]] >
          visited_idx[double_bond_partner[true_cis_trans_stereo_idx[i]]]) {
        std::swap(ct.idx1, ct.idx2);
        std::swap(ct.nbr1, ct.nbr2);
      }

      true_cis_trans_pairs.push_back(ct);
    }
  }

  for (size_t i = 0; i < true_cis_trans_pairs.size(); i++) {
    for (size_t j = i + 1; j < true_cis_trans_pairs.size(); j++) {
      if (visited_idx[true_cis_trans_pairs[i].idx1] >
          visited_idx[true_cis_trans_pairs[j].idx1])
        std::swap(true_cis_trans_pairs[i], true_cis_trans_pairs[j]);

      else if (visited_idx[true_cis_trans_pairs[i].idx1] ==
                   visited_idx[true_cis_trans_pairs[j].idx1] &&
               visited_idx[true_cis_trans_pairs[i].idx2] >
                   visited_idx[true_cis_trans_pairs[j].idx2])
        std::swap(true_cis_trans_pairs[i], true_cis_trans_pairs[j]);
    }
  }
}

void StereoGenerator::FindPotentialParaStereoCenters() {
  if (potential_para_stereo_idx.size() ||
      potential_para_cis_trans_stereo_idx.size()) {
    RepresentationSystem local_stabilizers = u_id;  // = stabilizers;
    PermutationIterator perm_it(u_automorph);
    const Permutations* stab;

    while (perm_it.GetNextPermutation()) {
      bool add(true);
      stab = perm_it.GetCombinedPermutation();

      // if the permutation swaps two true stereocenters, it might not be an
      // automorphism of the adj matrix when the two stereocenters have a
      // different configuration
      //  (note if one of them , e.g. the first, is a true s.c., then both are,
      //  i.e. if true_stereo[pp.first] is true it must hold that
      //  trueStereo[pp.second] is also true, otherwise they couldn't be
      //  swapped!
      for (auto&& pp : *stab) {
        if (true_stereo[pp.first] || true_cis_trans_stereo[pp.first]) {
          add = false;
          break;
        }
      }

      if (add) {
        local_stabilizers[perm_it.GetSmallestDiffIndex()].push_back(*stab);
        perm_it.SetCurrentIndexToSmallestDiffIndex();
      }
    }

    PermutationIterator perm_it_local(local_stabilizers);
    const std::vector<size_t>* permuted_indices;

    while (perm_it_local.GetNextPermutation() &&
           (potential_para_stereo_idx.size() ||
            potential_para_cis_trans_stereo_idx.size())) {
      permuted_indices = perm_it_local.GetPermutedIndices();

      for (size_t i = 0; i < potential_para_stereo_idx.size(); i++) {
        size_t idx = potential_para_stereo_idx[i];

        if ((*permuted_indices)[idx] == idx) {
          if (NeighborsArePermuted(A.GetAtomVector()[idx], *permuted_indices)) {
            potential_para_stereo_idx.erase(potential_para_stereo_idx.begin() +
                                            i);
            i--;
          }
        }
      }

      for (size_t i = 0; i < potential_para_cis_trans_stereo_idx.size(); i++) {
        size_t idx = potential_para_cis_trans_stereo_idx[i];

        if ((*permuted_indices)[idx] == idx) {
          if (NeighborsArePermuted(A.GetAtomVector()[idx], *permuted_indices)) {
            potential_para_cis_trans_stereo_idx.erase(
                potential_para_cis_trans_stereo_idx.begin() + i);
            i--;
          }
        }
      }
    }

    if (potential_para_stereo_idx.size() ||
        potential_para_cis_trans_stereo_idx.size()) {
      std::vector<bool> para_tetra_stereo(N, false);
      std::vector<bool> para_cis_trans_stereo(N, false);

      for (auto&& idx : potential_para_stereo_idx)
        para_tetra_stereo[idx] = true;

      for (auto&& idx : potential_para_cis_trans_stereo_idx)
        para_cis_trans_stereo[idx] = true;

      // remove CT para stereo centers that are not connected to another para or
      // true CT stereocenter
      bool removed_para_ct_stereo(false);

      do {
        removed_para_ct_stereo = false;

        for (size_t i = 0; i < potential_para_cis_trans_stereo_idx.size();
             i++) {
          size_t idx = potential_para_cis_trans_stereo_idx[i];
          size_t nbr = double_bond_partner[idx];

          if (!true_cis_trans_stereo[nbr] && !para_cis_trans_stereo[nbr]) {
            para_cis_trans_stereo[idx] = false;
            potential_para_cis_trans_stereo_idx.erase(
                potential_para_cis_trans_stereo_idx.begin() + i);
            i--;
            removed_para_ct_stereo = true;
          }
        }
      } while (removed_para_ct_stereo &&
               potential_para_cis_trans_stereo_idx.size());

      std::vector<Permutation> identical_stereo_center_pairs(0);
      std::vector<Permutation> identical_ct_stereo_center_pairs(0);
      PermutationIterator perm_it(u_automorph);

      while (perm_it.GetNextPermutation()) {
        permuted_indices = perm_it.GetPermutedIndices();

        for (size_t i = 0; i < true_stereo_idx.size(); i++) {
          size_t idx = true_stereo_idx[i];

          if ((*permuted_indices)[idx] > idx &&
              std::find(identical_stereo_center_pairs.begin(),
                        identical_stereo_center_pairs.end(),
                        Permutation({idx, (*permuted_indices)[idx]})) ==
                  identical_stereo_center_pairs.end()) {
            identical_stereo_center_pairs.push_back(
                Permutation({idx, (*permuted_indices)[idx]}));
          }
        }

        for (size_t i = 0; i < true_cis_trans_stereo_idx.size(); i++) {
          size_t idx = true_cis_trans_stereo_idx[i];

          if ((*permuted_indices)[idx] > idx &&
              std::find(identical_ct_stereo_center_pairs.begin(),
                        identical_ct_stereo_center_pairs.end(),
                        Permutation({idx, (*permuted_indices)[idx]})) ==
                  identical_ct_stereo_center_pairs.end()) {
            identical_ct_stereo_center_pairs.push_back(
                Permutation({idx, (*permuted_indices)[idx]}));
          }
        }
      }

      std::vector<std::vector<size_t>> identical_stereo_centers(0);
      std::vector<int> position_in_identical_stereo_centers(N, -1);
      std::vector<std::vector<size_t>> identical_ct_stereo_centers(0);
      std::vector<int> position_in_identical_ct_stereo_centers(N, -1);
      identical_stereo_centers.reserve(true_stereo_idx.size() / 2);
      identical_stereo_centers.reserve(true_cis_trans_stereo_idx.size() / 2);

      for (auto&& is : identical_stereo_center_pairs) {
        if (position_in_identical_stereo_centers[is.first] != -1) {
          identical_stereo_centers
              [position_in_identical_stereo_centers[is.first]]
                  .push_back(is.first);
          identical_stereo_centers
              [position_in_identical_stereo_centers[is.first]]
                  .push_back(is.second);
          position_in_identical_stereo_centers[is.second] =
              position_in_identical_stereo_centers[is.first];

        } else if (position_in_identical_stereo_centers[is.second] != -1) {
          identical_stereo_centers
              [position_in_identical_stereo_centers[is.second]]
                  .push_back(is.first);
          identical_stereo_centers
              [position_in_identical_stereo_centers[is.second]]
                  .push_back(is.second);
          position_in_identical_stereo_centers[is.first] =
              position_in_identical_stereo_centers[is.second];

        } else {
          identical_stereo_centers.emplace_back(
              std::vector<size_t>{is.first, is.second});
          position_in_identical_stereo_centers[is.first] =
              (int)identical_stereo_centers.size() - 1;
          position_in_identical_stereo_centers[is.second] =
              (int)identical_stereo_centers.size() - 1;
        }
      }

      for (auto&& is : identical_ct_stereo_center_pairs) {
        if (position_in_identical_ct_stereo_centers[is.first] != -1) {
          identical_ct_stereo_centers
              [position_in_identical_ct_stereo_centers[is.first]]
                  .push_back(is.first);
          identical_ct_stereo_centers
              [position_in_identical_ct_stereo_centers[is.first]]
                  .push_back(is.second);
          position_in_identical_ct_stereo_centers[is.second] =
              position_in_identical_ct_stereo_centers[is.first];

        } else if (position_in_identical_ct_stereo_centers[is.second] != -1) {
          identical_ct_stereo_centers
              [position_in_identical_ct_stereo_centers[is.second]]
                  .push_back(is.first);
          identical_ct_stereo_centers
              [position_in_identical_ct_stereo_centers[is.second]]
                  .push_back(is.second);
          position_in_identical_ct_stereo_centers[is.first] =
              position_in_identical_ct_stereo_centers[is.second];

        } else {
          identical_ct_stereo_centers.emplace_back(
              std::vector<size_t>{is.first, is.second});
          position_in_identical_ct_stereo_centers[is.first] =
              (int)identical_ct_stereo_centers.size() - 1;
          position_in_identical_ct_stereo_centers[is.second] =
              (int)identical_ct_stereo_centers.size() - 1;
        }
      }

      std::vector<size_t> symmetry_center_atoms(0);
      symmetry_center_atoms.reserve(N);
      std::vector<size_t> minimal_path;

      for (auto&& id : identical_stereo_centers) {
        std::sort(id.begin(), id.end());
        auto&& it = std::unique(id.begin(), id.end());
        id.resize(std::distance(id.begin(), it));

        for (size_t i = 0; i < id.size(); i++) {
          for (size_t j = i + 1; j < id.size(); j++) {
            FindShortestPath(minimal_path, id[i], id[j]);

            if (minimal_path.size() % 2 != 0)
              symmetry_center_atoms.push_back(
                  minimal_path[minimal_path.size() / 2]);
          }
        }
      }

      for (auto&& id : identical_ct_stereo_centers) {
        std::sort(id.begin(), id.end());
        auto&& it = std::unique(id.begin(), id.end());
        id.resize(std::distance(id.begin(), it));

        for (size_t i = 0; i < id.size(); i++) {
          for (size_t j = i + 1; j < id.size(); j++) {
            FindShortestPath(minimal_path, id[i], id[j]);

            if (minimal_path.size() % 2 != 0)
              symmetry_center_atoms.push_back(
                  minimal_path[minimal_path.size() / 2]);
          }
        }
      }

      std::sort(symmetry_center_atoms.begin(), symmetry_center_atoms.end());
      auto&& it1 = std::unique(symmetry_center_atoms.begin(),
                               symmetry_center_atoms.end());
      symmetry_center_atoms.resize(
          std::distance(symmetry_center_atoms.begin(), it1));
      std::vector<bool> is_symmetry_atom(N, false);

      for (auto&& idx : symmetry_center_atoms) is_symmetry_atom[idx] = true;

      potential_para_stereo_idx.resize(0);
      potential_para_cis_trans_stereo_idx.resize(0);
      potential_para_stereo_idx.reserve(N);
      potential_para_cis_trans_stereo_idx.reserve(N);

      for (size_t i = 0; i < N; i++) {
        if (!is_in_cycle[i]) {
          if (is_symmetry_atom[i]) {
            if (para_tetra_stereo[i])
              potential_para_stereo_idx.push_back(i);

            else if (para_cis_trans_stereo[i])
              potential_para_cis_trans_stereo_idx.push_back(i);

          } else if (para_cis_trans_stereo[i] &&
                     is_symmetry_atom[double_bond_partner[i]] &&
                     para_cis_trans_stereo[double_bond_partner[i]])
            potential_para_cis_trans_stereo_idx.push_back(i);

        } else {
          if (para_tetra_stereo[i])
            potential_para_stereo_idx.push_back(i);

          else if (para_cis_trans_stereo[i])
            potential_para_cis_trans_stereo_idx.push_back(i);
        }
      }

      for (size_t i = 0; i < potential_para_cis_trans_stereo_idx.size(); i++) {
        if (double_bond_partner[potential_para_cis_trans_stereo_idx[i]] >
            (int)potential_para_cis_trans_stereo_idx[i]) {
          if (visited_idx[potential_para_cis_trans_stereo_idx[i]] <
              visited_idx[double_bond_partner
                              [potential_para_cis_trans_stereo_idx[i]]]) {
            para_cis_trans_pairs.emplace_back(
                potential_para_cis_trans_stereo_idx[i],
                double_bond_partner[potential_para_cis_trans_stereo_idx[i]],
                single_bond_partner[potential_para_cis_trans_stereo_idx[i]],
                single_bond_partner
                    [double_bond_partner
                         [potential_para_cis_trans_stereo_idx[i]]]);

          } else {
            para_cis_trans_pairs.emplace_back(
                double_bond_partner[potential_para_cis_trans_stereo_idx[i]],
                potential_para_cis_trans_stereo_idx[i],
                single_bond_partner
                    [double_bond_partner
                         [potential_para_cis_trans_stereo_idx[i]]],
                single_bond_partner[potential_para_cis_trans_stereo_idx[i]]);
          }
        }
      }

      for (size_t i = 0; i < para_cis_trans_pairs.size(); i++) {
        for (size_t j = i + 1; j < para_cis_trans_pairs.size(); j++) {
          if (visited_idx[para_cis_trans_pairs[i].idx1] >
              visited_idx[para_cis_trans_pairs[j].idx1])
            std::swap(para_cis_trans_pairs[i], para_cis_trans_pairs[j]);

          else if (visited_idx[para_cis_trans_pairs[i].idx1] ==
                       visited_idx[para_cis_trans_pairs[j].idx1] &&
                   visited_idx[para_cis_trans_pairs[i].idx2] >
                       visited_idx[para_cis_trans_pairs[j].idx2])
            std::swap(para_cis_trans_pairs[i], para_cis_trans_pairs[j]);
        }
      }
    }
  }
}

void StereoGenerator::FindValidConfigurations(
    std::vector<Configuration>& valid_configurations) {
  std::vector<Config> valid_configurations_tetra(0);
  std::vector<Config> valid_configurations_ct(0);

  if (true_stereo_idx.size())
    FindValidTrueTetrahedralConfigurations(valid_configurations_tetra);

  if (true_cis_trans_pairs.size())
    FindValidTrueCTConfigurations(valid_configurations_ct);

  AddValidTrueConfigurations(valid_configurations, valid_configurations_tetra,
                             valid_configurations_ct);

  if ((valid_configurations_tetra.size() || valid_configurations_ct.size()) &&
      (potential_para_stereo_idx.size() ||
       potential_para_cis_trans_stereo_idx.size()))
    FindValidParaConfigurations(valid_configurations);
}

void StereoGenerator::FindValidTrueTetrahedralConfigurations(
    std::vector<Config>& valid_configurations_tetra) {
  const size_t n = true_stereo_idx.size();
  Config configuration_all(N);
  Config conf(n, true);
  Config all_true(n, true);
  std::vector<size_t> nbrs_original_order;
  std::vector<size_t> nbrs_permuted_order;
  Config configuration_permuted(conf.size());
  std::fill(configuration_all.begin(), configuration_all.end(), -1);
  const std::vector<size_t>* permuted_indices;

  do {
    ++conf;
    bool ignore(false);

    for (size_t i = 0; i < true_stereo_idx.size(); i++)
      configuration_all[true_stereo_idx[i]] = conf[i];

    PermutationIterator perm_it(u_automorph);

    while (perm_it.GetNextPermutation()) {
      permuted_indices = perm_it.GetPermutedIndices();
      // std::cout << *(perm_it.GetCombinedPermutation()) << std::endl;

      for (size_t j = 0; j < true_stereo_idx.size(); j++) {
        size_t idx = true_stereo_idx[j];
        size_t idx_permuted = (*permuted_indices)[idx];

        // if the stereocenter was swapped, determine new stereo configuration
        if (idx != idx_permuted) {
          NeighborOrder(idx, idx_permuted, *permuted_indices,
                        nbrs_original_order, nbrs_permuted_order, coming_from,
                        going_to, ring_connections);

          if (NumPerm(nbrs_original_order, nbrs_permuted_order) % 2 == 0)
            configuration_permuted[j] = configuration_all[idx_permuted];

          else
            configuration_permuted[j] = !configuration_all[idx_permuted];

        } else
          configuration_permuted[j] = configuration_all[idx];
      }

      if (configuration_permuted < conf) {
        ignore = true;
        break;
      }
    }

    if (!ignore) valid_configurations_tetra.push_back(conf);
  } while (conf != all_true);
}

void StereoGenerator::FindValidTrueCTConfigurations(
    std::vector<Config>& valid_configurations_ct) {
  const size_t num_ct = true_cis_trans_pairs.size();
  Config configuration_ct(num_ct, true);
  Config all_true_ct(num_ct, true);
  Config configuration_permuted(num_ct);
  const std::vector<size_t>* permuted_indices;

  do {
    ++configuration_ct;
    bool ignore(false);
    PermutationIterator perm_it(u_automorph);

    while (perm_it.GetNextPermutation()) {
      permuted_indices = perm_it.GetPermutedIndices();

      for (size_t j = 0; j < true_cis_trans_pairs.size(); j++) {
        size_t idx1 = true_cis_trans_pairs[j].idx1;
        size_t idx2 = true_cis_trans_pairs[j].idx2;
        size_t nbr1 = true_cis_trans_pairs[j].nbr1;
        size_t nbr2 = true_cis_trans_pairs[j].nbr2;
        size_t idx_permuted1 = (*permuted_indices)[idx1];
        size_t idx_permuted2 = (*permuted_indices)[idx2];
        size_t nbr_permuted1 = (*permuted_indices)[nbr1];
        size_t nbr_permuted2 = (*permuted_indices)[nbr2];

        if (idx1 != idx_permuted1 || idx2 != idx_permuted2) {
          // find the permuted cistrans pair of (idx1,idx2), i.e.
          // (idx_permuted1, idx_permuted2)
          for (size_t ii = 0; ii < true_cis_trans_pairs.size(); ii++) {
            // case that idx1 is swapped with idx1 of the ii-th
            // true_cis_trans_pairs
            //      and idx2 is swapped with idx2
            if (true_cis_trans_pairs[ii].idx1 == idx_permuted1 &&
                true_cis_trans_pairs[ii].idx2 == idx_permuted2) {
              if (true_cis_trans_pairs[ii].nbr1 == nbr_permuted1) {
                // even permutation (i.e. 0)
                if (true_cis_trans_pairs[ii].nbr2 == nbr_permuted2)
                  configuration_permuted[j] = configuration_ct[ii];

                // odd permutation (i.e. 1)
                else
                  configuration_permuted[j] = !configuration_ct[ii];
              }

              // odd permutation (i.e. 1)
              else if (true_cis_trans_pairs[ii].nbr2 == nbr_permuted2)
                configuration_permuted[j] = !configuration_ct[ii];

              // even permutation (i.e. 2)
              else
                configuration_permuted[j] = configuration_ct[ii];

              break;
            }

            // case that idx1 is swapped with idx2 of the ii-th
            // true_cis_trans_pairs
            //      and idx2 is swapped with idx1
            else if (true_cis_trans_pairs[ii].idx1 == idx_permuted2 &&
                     true_cis_trans_pairs[ii].idx2 == idx_permuted1) {
              if (true_cis_trans_pairs[ii].nbr1 == nbr_permuted2) {
                // even permutation (i.e. 0)
                if (true_cis_trans_pairs[ii].nbr2 == nbr_permuted1)
                  configuration_permuted[j] = configuration_ct[ii];

                // odd permutation (i.e. 1)
                else
                  configuration_permuted[j] = !configuration_ct[ii];
              }

              // odd permutation (i.e. 1)
              else if (true_cis_trans_pairs[ii].nbr2 == nbr_permuted1)
                configuration_permuted[j] = !configuration_ct[ii];

              // even permutation (i.e. 2)
              else
                configuration_permuted[j] = configuration_ct[ii];

              break;
            }
          }

        } else
          configuration_permuted[j] = configuration_ct[j];
      }

      if (configuration_permuted < configuration_ct) {
        ignore = true;
        break;
      }
    }

    if (!ignore) valid_configurations_ct.push_back(configuration_ct);
  } while (configuration_ct != all_true_ct);
}

void StereoGenerator::AddValidTrueConfigurations(
    std::vector<Configuration>& valid_configurations,
    std::vector<Config>& valid_configurations_tetra,
    std::vector<Config>& valid_configurations_ct) {
  if (valid_configurations_tetra.size()) {
    if (valid_configurations_ct.size()) {
      valid_configurations.resize(valid_configurations_tetra.size() *
                                  valid_configurations_ct.size());
      size_t idx(0);

      for (size_t i = 0; i < valid_configurations_tetra.size(); i++) {
        for (size_t j = 0; j < valid_configurations_ct.size(); j++)
          valid_configurations[idx++] = (Configuration(
              valid_configurations_tetra[i], valid_configurations_ct[j]));
      }

    } else {
      valid_configurations.resize(valid_configurations_tetra.size());

      for (size_t i = 0; i < valid_configurations_tetra.size(); i++)
        valid_configurations[i] =
            (Configuration(valid_configurations_tetra[i], Config(0)));
    }

  } else if (valid_configurations_ct.size()) {
    valid_configurations.resize(valid_configurations_ct.size());

    for (size_t j = 0; j < valid_configurations_ct.size(); j++)
      valid_configurations[j] =
          (Configuration(Config(0), valid_configurations_ct[j]));
  }
}

void StereoGenerator::FindValidParaConfigurations(
    std::vector<Configuration>& valid_configurations) {
  std::vector<int> position_in_potential_para_stereo_idx(N, -1);

  for (size_t i = 0; i < potential_para_stereo_idx.size(); i++)
    position_in_potential_para_stereo_idx[potential_para_stereo_idx[i]] =
        (int)i;

  size_t vcs = valid_configurations.size();
  std::vector<size_t> nbrs_original_order, nbrs_permuted_order;
  nbrs_original_order.reserve(4);
  nbrs_permuted_order.reserve(4);
  const Atom* a;
  const std::vector<size_t>* permuted_indices;

  // go through valid configurations
  for (size_t k = 0; k < vcs; k++) {
    auto vc_tet = valid_configurations[k].valid_config_tetra;
    auto vc_ct = valid_configurations[k].valid_config_ct;
    RepresentationSystem local_stabilizers = u_id;  // = stabilizers;
    PermutationIterator perm_it(u_automorph);

    // go through permutations to find stabilizers, taking into account the
    // configurations of the true stereocenters
    while (perm_it.GetNextPermutation()) {
      bool add(true);
      permuted_indices = perm_it.GetPermutedIndices();

      for (size_t i = 0; i < true_stereo_idx.size() && add; i++) {
        size_t idx = true_stereo_idx[i];
        size_t idx_permuted = (*permuted_indices)[idx];

        if (idx != idx_permuted) {
          NeighborOrder(idx, idx_permuted, (*permuted_indices),
                        nbrs_original_order, nbrs_permuted_order, coming_from,
                        going_to, ring_connections);

          if (NumPerm(nbrs_original_order, nbrs_permuted_order) % 2 == 0) {
            if (vc_tet[i] !=
                vc_tet[position_in_true_stereo_idx[(*permuted_indices)[idx]]]) {
              add = false;
              break;
            }

          } else {
            if (vc_tet[i] ==
                vc_tet[position_in_true_stereo_idx[(*permuted_indices)[idx]]]) {
              add = false;
              break;
            }
          }
        }
      }

      for (size_t i = 0; i < true_cis_trans_pairs.size() && add; i++) {
        size_t idx1 = true_cis_trans_pairs[i].idx1;
        size_t idx2 = true_cis_trans_pairs[i].idx2;
        size_t nbr1 = true_cis_trans_pairs[i].nbr1;
        size_t nbr2 = true_cis_trans_pairs[i].nbr2;
        size_t idx_permuted1 = (*permuted_indices)[idx1];
        size_t idx_permuted2 = (*permuted_indices)[idx2];
        size_t nbr_permuted1 = (*permuted_indices)[nbr1];
        size_t nbr_permuted2 = (*permuted_indices)[nbr2];

        if (idx1 != idx_permuted1 || idx2 != idx_permuted2) {
          for (size_t ii = 0; ii < true_cis_trans_pairs.size(); ii++) {
            if (true_cis_trans_pairs[ii].idx1 == idx_permuted1 &&
                true_cis_trans_pairs[ii].idx2 == idx_permuted2) {
              if (true_cis_trans_pairs[ii].nbr1 == nbr_permuted1) {
                if (true_cis_trans_pairs[ii].nbr2 == nbr_permuted2) {
                  if (vc_ct[i] != vc_ct[ii]) {
                    add = false;
                    break;
                  }

                } else {
                  if (vc_ct[i] == vc_ct[ii]) {
                    add = false;
                    break;
                  }
                }

              } else if (true_cis_trans_pairs[ii].nbr2 == nbr_permuted2) {
                if (vc_ct[i] == vc_ct[ii]) {
                  add = false;
                  break;
                }

              } else {
                if (vc_ct[i] != vc_ct[ii]) {
                  add = false;
                  break;
                }
              }

              break;

            } else if (true_cis_trans_pairs[ii].idx1 == idx_permuted2 &&
                       true_cis_trans_pairs[ii].idx2 == idx_permuted1) {
              if (true_cis_trans_pairs[ii].nbr1 == nbr_permuted2) {
                if (true_cis_trans_pairs[ii].nbr2 == nbr_permuted1) {
                  if (vc_ct[i] != vc_ct[ii]) {
                    add = false;
                    break;
                  }

                } else {
                  if (vc_ct[i] == vc_ct[ii]) {
                    add = false;
                    break;
                  }
                }

              } else if (true_cis_trans_pairs[ii].nbr2 == nbr_permuted1) {
                if (vc_ct[i] == vc_ct[ii]) {
                  add = false;
                  break;
                }

              } else {
                if (vc_ct[i] != vc_ct[ii]) {
                  add = false;
                  break;
                }
              }

              break;
            }
          }
        }
      }

      if (add) {
        local_stabilizers[perm_it.GetSmallestDiffIndex()].push_back(
            *(perm_it.GetCombinedPermutation()));
        perm_it.SetCurrentIndexToSmallestDiffIndex();
      }
    }

    Config para_configuration(potential_para_stereo_idx.size(), 1);
    Config para_configuration_ct(para_cis_trans_pairs.size(), 1);
    const Config all_true(para_configuration.size(), 1);
    const Config all_not_para(para_configuration.size(), -1);
    const Config all_true_ct(para_configuration_ct.size(), 1);
    const Config all_not_para_ct(para_configuration_ct.size(), -1);
    bool found_new_configuration(false);

    do {
      ++para_configuration;
      Config para_configuration_orig(para_configuration);

      do {
        ++para_configuration_ct;
        Config para_configuration_ct_orig(para_configuration_ct);
        bool changed_configuration;
        PermutationIterator perm_it_local(local_stabilizers);
        const std::vector<size_t>* permuted_indices;

        do {
          changed_configuration = false;
          perm_it_local.Reset();

          while (perm_it_local.GetNextPermutation()) {
            permuted_indices = perm_it_local.GetPermutedIndices();
            bool skip(false);

            for (size_t i = 0; i < potential_para_stereo_idx.size(); i++) {
              size_t idx = potential_para_stereo_idx[i];
              size_t idx_permuted = (*permuted_indices)[idx];

              if (idx != idx_permuted) {
                NeighborOrder(idx, idx_permuted, (*permuted_indices),
                              nbrs_original_order, nbrs_permuted_order,
                              coming_from, going_to, ring_connections);

                if (NumPerm(nbrs_original_order, nbrs_permuted_order) % 2 ==
                    0) {
                  if (para_configuration[i] !=
                      para_configuration[position_in_potential_para_stereo_idx[(
                          *permuted_indices)[idx]]]) {
                    skip = true;
                    break;
                  }

                } else {
                  if (para_configuration[i] ==
                      para_configuration[position_in_potential_para_stereo_idx[(
                          *permuted_indices)[idx]]]) {
                    skip = true;
                    break;
                  }
                }
              }
            }

            for (size_t i = 0; i < para_cis_trans_pairs.size() && !skip; i++) {
              size_t idx1 = para_cis_trans_pairs[i].idx1;
              size_t idx2 = para_cis_trans_pairs[i].idx2;
              size_t nbr1 = para_cis_trans_pairs[i].nbr1;
              size_t nbr2 = para_cis_trans_pairs[i].nbr2;
              size_t idx_permuted1 = (*permuted_indices)[idx1];
              size_t idx_permuted2 = (*permuted_indices)[idx2];
              size_t nbr_permuted1 = (*permuted_indices)[nbr1];
              size_t nbr_permuted2 = (*permuted_indices)[nbr2];

              if (idx1 != idx_permuted1 || idx2 != idx_permuted2) {
                for (size_t ii = 0; ii < para_cis_trans_pairs.size(); ii++) {
                  if (para_cis_trans_pairs[ii].idx1 == idx_permuted1 &&
                      para_cis_trans_pairs[ii].idx2 == idx_permuted2) {
                    if (para_cis_trans_pairs[ii].nbr1 == nbr_permuted1) {
                      if (para_cis_trans_pairs[ii].nbr2 == nbr_permuted2) {
                        if (para_configuration_ct[i] !=
                            para_configuration_ct[ii]) {
                          skip = true;
                          break;
                        }

                      } else {
                        if (para_configuration_ct[i] ==
                            para_configuration_ct[ii]) {
                          skip = true;
                          break;
                        }
                      }

                    } else if (para_cis_trans_pairs[ii].nbr2 == nbr_permuted2) {
                      if (para_configuration_ct[i] ==
                          para_configuration_ct[ii]) {
                        skip = true;
                        break;
                      }

                    } else {
                      if (para_configuration_ct[i] !=
                          para_configuration_ct[ii]) {
                        skip = true;
                        break;
                      }
                    }

                    break;

                  } else if (para_cis_trans_pairs[ii].idx1 == idx_permuted2 &&
                             para_cis_trans_pairs[ii].idx2 == idx_permuted1) {
                    if (para_cis_trans_pairs[ii].nbr1 == nbr_permuted2) {
                      if (para_cis_trans_pairs[ii].nbr2 == nbr_permuted1) {
                        if (para_configuration_ct[i] !=
                            para_configuration_ct[ii]) {
                          skip = true;
                          break;
                        }

                      } else {
                        if (para_configuration_ct[i] ==
                            para_configuration_ct[ii]) {
                          skip = true;
                          break;
                        }
                      }

                    } else if (para_cis_trans_pairs[ii].nbr2 == nbr_permuted1) {
                      if (para_configuration_ct[i] ==
                          para_configuration_ct[ii]) {
                        skip = true;
                        break;
                      }

                    } else {
                      if (para_configuration_ct[i] !=
                          para_configuration_ct[ii]) {
                        skip = true;
                        break;
                      }
                    }

                    break;
                  }
                }
              }
            }

            if (!skip) {
              for (size_t i = 0; i < potential_para_stereo_idx.size(); i++) {
                size_t idx = potential_para_stereo_idx[i];
                a = &A.GetAtomVector()[idx];

                if ((*permuted_indices)[idx] == idx) {
                  if (NeighborsArePermuted(*a, *permuted_indices) &&
                      para_configuration[i] != -1) {
                    para_configuration[i] = -1;
                    changed_configuration = true;
                  }
                }
              }

              for (size_t i = 0; i < para_cis_trans_pairs.size(); i++) {
                size_t idx1 = para_cis_trans_pairs[i].idx1;
                size_t idx2 = para_cis_trans_pairs[i].idx2;
                size_t idx_permuted1 = (*permuted_indices)[idx1];
                size_t idx_permuted2 = (*permuted_indices)[idx2];
                size_t nbr1 = para_cis_trans_pairs[i].nbr1;
                size_t nbr2 = para_cis_trans_pairs[i].nbr2;
                size_t nbr_permuted1 = (*permuted_indices)[nbr1];
                size_t nbr_permuted2 = (*permuted_indices)[nbr2];

                if (((idx1 == idx_permuted1 && nbr1 != nbr_permuted1) ||
                     (idx2 == idx_permuted2 && nbr2 != nbr_permuted2)) &&
                    para_configuration_ct[i] != -1) {
                  para_configuration_ct[i] = -1;
                  changed_configuration = true;
                }
              }
            }
          }
        } while (changed_configuration);

        bool ignore(false);

        if (para_configuration != all_not_para ||
            para_configuration_ct != all_not_para_ct) {
          Config configuration_all(N, -1);

          for (size_t i = 0; i < true_stereo_idx.size(); i++)
            configuration_all[true_stereo_idx[i]] = vc_tet[i];

          for (size_t i = 0; i < potential_para_stereo_idx.size(); i++)
            configuration_all[potential_para_stereo_idx[i]] =
                para_configuration[i];

          perm_it_local.Reset();

          while (perm_it_local.GetNextPermutation()) {
            permuted_indices = perm_it_local.GetPermutedIndices();
            Config configuration_permuted(para_configuration.size());
            Config configuration_permuted_ct(para_configuration_ct.size(), -1);

            for (size_t j = 0; j < potential_para_stereo_idx.size(); j++) {
              if (para_configuration[j] == -1)
                configuration_permuted[j] = 0;

              else {
                size_t idx = potential_para_stereo_idx[j];
                size_t idx_permuted = (*permuted_indices)[idx];
                NeighborOrder(idx, idx_permuted, *permuted_indices,
                              nbrs_original_order, nbrs_permuted_order,
                              coming_from, going_to, ring_connections);

                if (NumPerm(nbrs_original_order, nbrs_permuted_order) % 2 == 0)
                  configuration_permuted[j] = configuration_all[(
                      *permuted_indices)[potential_para_stereo_idx[j]]];

                else
                  configuration_permuted[j] = !configuration_all[(
                      *permuted_indices)[potential_para_stereo_idx[j]]];
              }
            }

            if (IsSmaller(configuration_permuted, para_configuration_orig)) {
              ignore = true;
              break;
            }

            for (size_t j = 0; j < para_cis_trans_pairs.size(); j++) {
              if (para_configuration_ct[j] == -1)
                configuration_permuted_ct[j] = 0;

              else {
                size_t idx1 = para_cis_trans_pairs[j].idx1;
                size_t idx2 = para_cis_trans_pairs[j].idx2;
                size_t nbr1 = para_cis_trans_pairs[j].nbr1;
                size_t nbr2 = para_cis_trans_pairs[j].nbr2;
                size_t idx_permuted1 = (*permuted_indices)[idx1];
                size_t idx_permuted2 = (*permuted_indices)[idx2];
                size_t nbr_permuted1 = (*permuted_indices)[nbr1];
                size_t nbr_permuted2 = (*permuted_indices)[nbr2];

                if (idx1 != idx_permuted1 || idx2 != idx_permuted2) {
                  for (size_t ii = 0; ii < para_cis_trans_pairs.size(); ii++) {
                    if (para_cis_trans_pairs[ii].idx1 == idx_permuted1 &&
                        para_cis_trans_pairs[ii].idx2 == idx_permuted2) {
                      if (para_cis_trans_pairs[ii].nbr1 == nbr_permuted1) {
                        if (para_cis_trans_pairs[ii].nbr2 == nbr_permuted2)
                          configuration_permuted_ct[j] =
                              para_configuration_ct[ii];

                        else
                          configuration_permuted_ct[j] =
                              !para_configuration_ct[ii];

                      } else if (para_cis_trans_pairs[ii].nbr2 == nbr_permuted2)
                        configuration_permuted_ct[j] =
                            !para_configuration_ct[ii];

                      else
                        configuration_permuted_ct[j] =
                            para_configuration_ct[ii];

                      break;

                    } else if (para_cis_trans_pairs[ii].idx1 == idx_permuted2 &&
                               para_cis_trans_pairs[ii].idx2 == idx_permuted1) {
                      if (para_cis_trans_pairs[ii].nbr1 == nbr_permuted2) {
                        if (para_cis_trans_pairs[ii].nbr2 == nbr_permuted1)
                          configuration_permuted_ct[j] =
                              para_configuration_ct[ii];

                        else
                          configuration_permuted_ct[j] =
                              !para_configuration_ct[ii];

                      } else if (para_cis_trans_pairs[ii].nbr2 == nbr_permuted1)
                        configuration_permuted_ct[j] =
                            !para_configuration_ct[ii];

                      else
                        configuration_permuted_ct[j] =
                            para_configuration_ct[ii];

                      break;
                    }
                  }

                } else
                  configuration_permuted_ct[j] = para_configuration_ct[j];
              }
            }

            if (IsSmaller(configuration_permuted_ct,
                          para_configuration_ct_orig)) {
              ignore = true;
              break;
            }
          }

          if (!ignore) {
            found_new_configuration = true;
            Config new_configuration(0);
            new_configuration.reserve(vc_tet.size() +
                                      para_configuration.size());
            new_configuration.insert(new_configuration.end(), vc_tet.begin(),
                                     vc_tet.end());
            new_configuration.insert(new_configuration.end(),
                                     para_configuration.begin(),
                                     para_configuration.end());
            Config new_configuration_ct(0);
            new_configuration_ct.reserve(vc_ct.size() +
                                         para_configuration_ct.size());
            new_configuration_ct.insert(new_configuration_ct.end(),
                                        vc_ct.begin(), vc_ct.end());
            new_configuration_ct.insert(new_configuration_ct.end(),
                                        para_configuration_ct.begin(),
                                        para_configuration_ct.end());
            valid_configurations.push_back(
                Configuration(new_configuration, new_configuration_ct));
          }
        }

        para_configuration = para_configuration_orig;
        para_configuration_ct = para_configuration_ct_orig;
      } while (para_configuration_ct < all_true_ct);
    } while (para_configuration < all_true);

    if (found_new_configuration) {
      valid_configurations.erase(valid_configurations.begin() + k);
      k--;
      vcs--;

    } else {
      valid_configurations[k].valid_config_tetra.insert(
          valid_configurations[k].valid_config_tetra.end(),
          all_not_para.begin(), all_not_para.end());
      valid_configurations[k].valid_config_ct.insert(
          valid_configurations[k].valid_config_ct.end(),
          all_not_para_ct.begin(), all_not_para_ct.end());
    }
  }
}

void StereoGenerator::WriteStereoSmiles(
    std::vector<Configuration>& valid_configurations) {
  stereo_smiles.reserve(valid_configurations.size());
  size_t num_true_stereo_centers = true_stereo_idx.size();
  num_stereo_centers = num_true_stereo_centers;
  num_ct_bonds = true_cis_trans_pairs.size();
  std::vector<CisTransPair> all_cis_trans_pairs;
  all_cis_trans_pairs.reserve(true_cis_trans_pairs.size() +
                              para_cis_trans_pairs.size());
  all_cis_trans_pairs.insert(all_cis_trans_pairs.end(),
                             true_cis_trans_pairs.begin(),
                             true_cis_trans_pairs.end());
  all_cis_trans_pairs.insert(all_cis_trans_pairs.end(),
                             para_cis_trans_pairs.begin(),
                             para_cis_trans_pairs.end());

  for (size_t i = 0; i < all_cis_trans_pairs.size(); i++) {
    for (size_t j = i + 1; j < all_cis_trans_pairs.size(); j++) {
      if (visited_idx[all_cis_trans_pairs[i].idx1] >
          visited_idx[all_cis_trans_pairs[j].idx1]) {
        for (auto&& valid_conf : valid_configurations)
          std::swap(valid_conf.valid_config_ct[i],
                    valid_conf.valid_config_ct[j]);

        std::swap(all_cis_trans_pairs[i], all_cis_trans_pairs[j]);

      } else if (visited_idx[all_cis_trans_pairs[i].idx1] ==
                     visited_idx[all_cis_trans_pairs[j].idx1] &&
                 visited_idx[all_cis_trans_pairs[i].idx2] >
                     visited_idx[all_cis_trans_pairs[j].idx2]) {
        for (auto&& valid_conf : valid_configurations)
          std::swap(valid_conf.valid_config_ct[i],
                    valid_conf.valid_config_ct[j]);

        std::swap(all_cis_trans_pairs[i], all_cis_trans_pairs[j]);
      }
    }
  }

  Config vc;
  Config vc_flipped(valid_configurations[0].valid_config_tetra.size());
  auto smiles_blocks_local = smiles_blocks;

  for (size_t pos = 0; pos < valid_configurations.size(); pos++) {
    auto&& valid_conf = valid_configurations[pos];
    vc.assign(valid_conf.valid_config_tetra.begin(),
              valid_conf.valid_config_tetra.end());
    std::fill(vc_flipped.begin(), vc_flipped.end(), -1);

    for (size_t i = 0; i < smiles_blocks_local.size(); i++) {
      // only Reset the element_name, as it is the only thing that's potentially
      // changed!
      smiles_blocks_local[i].element_name = smiles_blocks[i].element_name;
    }

    for (size_t i = 0; i < true_stereo_idx.size(); i++) {
      smiles_blocks_local[true_stereo_idx[i]].element_name.insert(0, "[");

      if (vc[i]) {
        vc_flipped[i] = 0;

        if (A.GetAtomVector()[true_stereo_idx[i]].GetNumHydrogens() +
            A.GetAtomVector()[true_stereo_idx[i]].GetNumFixedHydrogens())
          smiles_blocks_local[true_stereo_idx[i]].element_name += "@H]";

        else
          smiles_blocks_local[true_stereo_idx[i]].element_name += "@]";

      } else {
        vc_flipped[i] = 1;

        if (A.GetAtomVector()[true_stereo_idx[i]].GetNumHydrogens() +
            A.GetAtomVector()[true_stereo_idx[i]].GetNumFixedHydrogens())
          smiles_blocks_local[true_stereo_idx[i]].element_name += "@@H]";

        else
          smiles_blocks_local[true_stereo_idx[i]].element_name += "@@]";
      }
    }

    size_t num_stereo_centers_local(num_true_stereo_centers);

    for (size_t i = 0; i < potential_para_stereo_idx.size(); i++) {
      if (vc[i + true_stereo_idx.size()] != -1) {
        num_stereo_centers_local++;
        smiles_blocks_local[potential_para_stereo_idx[i]].element_name.insert(
            0, "[");

        if (vc[i + true_stereo_idx.size()]) {
          vc_flipped[i + true_stereo_idx.size()] = 0;

          if (A.GetAtomVector()[potential_para_stereo_idx[i]]
                  .GetNumHydrogens() +
              A.GetAtomVector()[potential_para_stereo_idx[i]]
                  .GetNumFixedHydrogens())
            smiles_blocks_local[potential_para_stereo_idx[i]].element_name +=
                "@H]";

          else
            smiles_blocks_local[potential_para_stereo_idx[i]].element_name +=
                "@]";

        } else {
          vc_flipped[i + true_stereo_idx.size()] = 1;

          if (A.GetAtomVector()[potential_para_stereo_idx[i]]
                  .GetNumHydrogens() +
              A.GetAtomVector()[potential_para_stereo_idx[i]]
                  .GetNumFixedHydrogens())
            smiles_blocks_local[potential_para_stereo_idx[i]].element_name +=
                "@@H]";

          else
            smiles_blocks_local[potential_para_stereo_idx[i]].element_name +=
                "@@]";
        }
      }
    }

    vc = valid_conf.valid_config_ct;
    size_t num_ct_bonds_local(0);

    for (size_t i = 0; i < all_cis_trans_pairs.size(); i++) {
      if (vc[i] != -1) {
        num_ct_bonds_local++;
        bool first_bond_down(false);

        if (smiles_blocks_local[all_cis_trans_pairs[i].idx1].element_name[0] ==
            '\\')
          first_bond_down = true;

        else if (smiles_blocks_local[all_cis_trans_pairs[i].idx1]
                     .element_name[0] != '/')
          smiles_blocks_local[all_cis_trans_pairs[i].idx1].element_name.insert(
              0, "/");

        // case trans
        if (!vc[i]) {
          if (first_bond_down)
            smiles_blocks_local[all_cis_trans_pairs[i].nbr2]
                .element_name.insert(0, "\\");

          else
            smiles_blocks_local[all_cis_trans_pairs[i].nbr2]
                .element_name.insert(0, "/");
        }

        // case cis
        else {
          if (first_bond_down)
            smiles_blocks_local[all_cis_trans_pairs[i].nbr2]
                .element_name.insert(0, "/");

          else
            smiles_blocks_local[all_cis_trans_pairs[i].nbr2]
                .element_name.insert(0, "\\");
        }
      }
    }

    int enantiomer(-1);

    for (size_t i = 0; i < valid_configurations.size(); i++) {
      if (valid_configurations[i].valid_config_tetra == vc_flipped &&
          pos != i &&
          valid_configurations[i].valid_config_ct ==
              valid_conf.valid_config_ct) {
        enantiomer = (int)i;
        break;
      }
    }

    num_stereo_centers = std::max(num_stereo_centers, num_stereo_centers_local);
    num_ct_bonds = std::max(num_ct_bonds, num_ct_bonds_local);
    std::string smiles("");
    size_t length(0);

    for (auto&& smi : smiles_blocks_local) length += smi.size();

    smiles.reserve(length);

    for (auto&& v : visited_indices) {
      auto&& n = smiles_blocks_local[v];
      smiles += n.opening_braces + n.bond_type + n.element_name;

      for (auto&& ri : n.ring_indices) {
        smiles += ri.first;

        if (ri.second <= 9)
          smiles += std::to_string(ri.second);

        else
          smiles += "%" + std::to_string(ri.second);
      }

      smiles += n.closing_braces;
    }

    stereo_smiles.emplace_back(
        std::tuple<std::string, int, std::pair<int, int>>{
            smiles, enantiomer,
            std::pair<int, int>(num_stereo_centers_local, num_ct_bonds_local)});
  }
}

void StereoGenerator::FindShortestPath(std::vector<size_t>& minimal_path,
                                       size_t i, size_t j) const {
  minimal_path.resize(0);
  std::vector<size_t> dist(N, 10000);
  std::vector<size_t> prev(N);
  dist[i] = 0;
  prev[i] = i;
  std::vector<size_t> q;

  for (size_t ii = 0; ii < N; ii++) q.push_back(ii);

  while (q.size() != 0) {
    size_t vertex_min_dist(q[0]);
    size_t min(dist[q[0]]);
    size_t pos(0);

    for (size_t ii = 0; ii < q.size(); ii++) {
      auto&& node = q[ii];

      if (dist[node] < min) {
        min = dist[node];
        vertex_min_dist = node;
        pos = ii;
      }
    }

    q.erase(q.begin() + pos);
    size_t cur = vertex_min_dist;

    if (cur == j) break;

    for (auto&& nbr : A.GetAtomVector()[cur].GetNeighbors()) {
      if (dist[cur] + 1 < dist[nbr]) {
        dist[nbr] = dist[cur] + 1;
        prev[nbr] = cur;
      }
    }
  }

  size_t cur = j;

  while (prev[cur] != i) {
    minimal_path.push_back(prev[cur]);
    cur = prev[cur];
  }
}

bool StereoGenerator::IsSmaller(Config& perm, Config& orig) const {
  for (size_t i = 0; i < orig.size(); i++) {
    if (perm[i] > orig[i])
      return false;

    else if (perm[i] < orig[i])
      return true;
  }

  return false;
}

Config& operator++(Config& b) {
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
  
bool StereoGenerator::NeighborsArePermuted(
    const Atom& a, const std::vector<size_t>& permuted_indices) {
  for (auto neighbor : a.GetNeighbors()) {
    if (neighbor != permuted_indices[neighbor]) return true;
  }

  return false;
}

size_t StereoGenerator::GetNumStereoCenters() const {
  return num_stereo_centers;
}
size_t StereoGenerator::GetNumCTBonds() const { return num_ct_bonds; }

}  // namespace enu

}  // namespace combi_ff