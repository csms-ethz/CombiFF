// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "SmilesHandler.h"

#include <cmath>
#include <sstream>

#include "smilesToMatrix.h"

namespace combi_ff {

namespace cnv {

void SmilesHandler::Run() {
  cnv::Handler::Run();
  SmilesHandler::PrintFirstLine();

  // determine the different smiles strings from the input_list
  for (auto smiles_it = input_list.begin(); smiles_it != input_list.end();
       ++smiles_it) {
    for (size_t i = 0; i < smiles_it->second.size(); i++) {
      if ((smiles_it->second)[i] == ' ' || (smiles_it->second)[i] == ',') {
        auto next = ++smiles_it;
        smiles_it--;
        input_list.insert(
            next, {smiles_it->first, smiles_it->second.substr(
                                         i + 1, smiles_it->second.size() - i)});
        smiles_it->second = smiles_it->second.substr(0, i);
      }
    }
  }

  for (auto&& smiles_orig : input_list) {
    if (smiles_orig.second.size()) {
      std::cerr << smiles_orig.second << std::endl;

      if (smiles_orig.second.find('.') != std::string::npos) {
        size_t i = 0;
        std::string smiles_canon("");
        std::string fmi("");
        cnv::AdjacencyMatrix A;

        while (i < smiles_orig.second.size()) {
          std::string smiles("");

          for (; i < smiles_orig.second.size(); i++) {
            if (smiles_orig.second[i] == '.') {
              i++;
              break;
            }

            smiles += smiles_orig.second[i];
          }

          max_size = std::min((unsigned long)column_width * 3,
                              std::max(max_size, smiles_orig.second.size()));
          std::string fmi_tmp("");
          std::string smiles_canon_tmp("");
          ConvertSmiles(smiles, smiles_canon_tmp, A);

          if (print_options[cnv::print_family_enumeration])
            FindFamilyIdentifier(smiles_canon, fmi_tmp);

          smiles_canon += smiles_canon_tmp;

          if (i != smiles_orig.second.size()) smiles_canon += '.';

          fmi += fmi_tmp;
        }

        PrintOutput(smiles_orig.second, smiles_canon, fmi, A);

      } else {
        cnv::AdjacencyMatrix A;
        std::string smiles_canon("");
        max_size = std::min((unsigned long)column_width * 3,
                            std::max(max_size, smiles_orig.second.size()));
        ConvertSmiles(smiles_orig.second, smiles_canon, A);
        std::string fmi("");

        for (auto&& a : A.GetAtomVector()) {
          std::cout << a << " " << a.GetTetraStereo() << std::endl;
        }

        if (print_options[cnv::print_family_enumeration])
          FindFamilyIdentifier(smiles_canon, fmi);

        PrintOutput(smiles_orig.second, smiles_canon, fmi, A);
      }
    }
  }
}

void SmilesHandler::FindFamilyIdentifier(const std::string& smiles_canon,
                                         std::string& fmi) {
  bool found(false);
  std::string line, code, formula, smiles_cur;

  for (auto&& family_file_name : fie_file_names) {
    std::ifstream familyFile(family_file_name);

    if (!familyFile.is_open()) {
      throw combi_ff::input_error("family isomer enumeration file " +
                                  family_file_name + " could not be opened");
    }

    while (std::getline(familyFile, line)) {
      if (line[0] != '#') {
        std::istringstream s(line);
        s >> code >> formula >> smiles_cur;

        if (smiles_cur == smiles_canon) {
          found = true;
          break;
        }
      }
    }

    familyFile.close();

    if (found) break;
  }

  if (!found) {
    for (auto&& ffn : fie_file_names) std::cerr << ffn << " ";

    std::cerr << '\n';
    throw combi_ff::input_error(
        "smiles " + smiles_canon +
        " not found in any of the above family isomer enumeration files ");

  } else
    fmi = code;
}

void SmilesHandler::PrintFirstLine() {
  cnv::Handler::PrintFirstLine("# originalSmiles ");
}

void SmilesHandler::ConvertSmiles(const std::string& smiles_orig,
                                  std::string& smiles_canon,
                                  cnv::AdjacencyMatrix& A) {
  if (smiles_orig.size() && smiles_orig != "%") {
    try {
      A = (combi_ff::ConvertToMatrix(smiles_orig));
    }

    // catch unhandled user input
    catch (combi_ff::input_warning& e) {
      std::cerr << "?Warning: " << e.what() << std::endl;
      return;
    }

    if (!A.GetN()) return;

    // convert aromatic bonds 1.5 to alternating double and single bonds
    bool set_last_arom_to_2(false);

    for (int i = 0; i < (int)A.GetN(); i++) {
      for (int j = 0; j < (int)A.GetN(); j++) {
        if (A.GetElement(i, j) == 1.5) {
          if (!set_last_arom_to_2) {
            A.SetElement(i, j, 2);
            set_last_arom_to_2 = true;

          } else {
            if (std::find(A.GetElements().begin() + i * A.GetN(),
                          A.GetElements().begin() + i * A.GetN() + A.GetN(),
                          2) !=
                A.GetElements().begin() + i * A.GetN() + A.GetN())
              A.SetElement(i, j, 1);

            else
              A.SetElement(i, j, 2);
          }
        }
      }
    }

    smiles_canon = "";

    if (print_options[cnv::print_canon_smiles] ||
        print_options[cnv::print_stack] || print_options[cnv::print_matrix] ||
        print_options[cnv::print_family_enumeration] ||
        print_options[cnv::print_canon_atom_vector]) {
      cnv::SmilesGeneratorCnv smiles_gen_orig(A);
      smiles_gen_orig.GenerateSmilesNonCanon();
      std::string smiles_orig_ = smiles_gen_orig.GetSmiles();
      const auto& visited_indices(smiles_gen_orig.GetVisitedIndices());
      std::cout << "original smiles (reconstructed): " << smiles_orig_
                << std::endl;
      std::cout << "visited indices: " << visited_indices << std::endl;

      A.Print();
      for (auto&& a : A.GetAtomVector())
        std::cout << a << " " << a.GetTetraStereo() << " " << a.GetStereo()
                  << std::endl;

      Permutations matrix_permutations =
          A.SortAtomVector_(A.GetIsAromaticCarbon());
      combi_ff::LambdaVector lambda = A.GetLambda();
      size_t N = A.GetN();
      // num_perms[i] indicates, with how many other atom indices the index of
      // atom i can be permuted note: index permutations are only allowed within
      // the same lambda partition, and with higher indices e.g. for lambda =
      // [1, 3, 1, 4, 6], we have num_perms = [1, 3, 2, 1, 1, 4, 3, 2, 1, 6, 5,
      // 4, 3, 2, 1]
      std::vector<size_t> num_perms;
      int ind(0);

      for (size_t i = 0; i < lambda.size(); i++) {
        for (size_t j = 0; j < lambda[i]; j++)
          num_perms.push_back(
              std::accumulate(lambda.begin(), lambda.begin() + i + 1, 0) -
              ind++);
      }

      // a RepresentationSystem contains all the possible index permutations for
      // the different atom indices
      //  e.g. for lambda = [2, 1, 4], the RepresentationSystem is
      //  (0,0), (0,1)         (for idx 0)
      //  (1,1)            (for idx 1)
      //  (2,2)            (for idx 2)
      //  (3,3), (3,4), (3,5), (3,6) (for idx 3)
      //  (4,4), (4,5), (4,6)      (for idx 4)
      //  (5,5), (5,6)         (for idx 5)
      //  (6,6)            (for idx 6)
      combi_ff::RepresentationSystem u(N);

      for (size_t i = 0; i < N; i++) {
        u[i].reserve(num_perms[i]);

        for (size_t j = 0; j < num_perms[i]; j++) {
          u[i].push_back(Permutations(1, {i, i + j}));
        }
      }

      std::vector<size_t> idx = A.GetIndices();
      A.MakeCanonical(u, idx, canon_iteration_limit, matrix_permutations);
      std::cout << matrix_permutations << std::endl;
      std::iota(idx.begin(), idx.end(), 0);
      for (auto p : matrix_permutations) std::swap(idx[p.first], idx[p.second]);
      /*if(output_file.is_open())
        output_file << "canonicalizing permutation is " << idx << '\n';
      else
        std::cout << "canonicalizing permutation is " << idx << '\n';*/

      A.Print();
      cnv::SmilesGeneratorCnv smiles_gen(A);
      smiles_gen.GenerateSmiles();
      smiles_canon = smiles_gen.GetSmiles();
      for (auto&& a : A.GetAtomVector())
        std::cout << a << " " << a.GetTetraStereo() << " " << a.GetStereo()
                  << std::endl;
      const std::vector<size_t>& coming_from_orig =
          smiles_gen_orig.GetComingFrom();
      const std::vector<std::vector<size_t>>& going_to_orig =
          smiles_gen_orig.GetGoingTo();
      const std::vector<std::vector<size_t>>& ring_connections_orig =
          smiles_gen_orig.GetRingConnections();
      const std::vector<size_t>& coming_from_new = smiles_gen.GetComingFrom();
      const std::vector<std::vector<size_t>>& going_to_new =
          smiles_gen.GetGoingTo();
      const std::vector<std::vector<size_t>>& ring_connections_new =
          smiles_gen.GetRingConnections();
      std::cout << idx << std::endl;
      bool stereo = false;
      for (size_t i = 0; i < A.GetN(); i++) {
        size_t idx_original = idx[i];
        const CnvAtom& a_orig = A.GetAtomVector()[idx_original];
        CnvAtom& a_new = A.GetAtomVectorNonConst()[i];
        std::vector<size_t> nbrs_original_order(1);
        nbrs_original_order.reserve(4);
        std::vector<size_t> nbrs_new_order(1);
        nbrs_new_order.reserve(4);

        if (a_new.GetStereo()) {
          stereo = true;
          std::cout << "at " << i << std::endl;
          std::cout << "orig positoin: " << idx_original << std::endl;

          nbrs_original_order[0] = coming_from_orig[idx_original];
          nbrs_original_order.insert(
              nbrs_original_order.end(),
              ring_connections_orig[idx_original].begin(),
              ring_connections_orig[idx_original].end());
          nbrs_original_order.insert(nbrs_original_order.end(),
                                     going_to_orig[idx_original].begin(),
                                     going_to_orig[idx_original].end());
          std::cout << a_orig << " " << a_orig.GetTetraStereo() << std::endl;
          std::cout << a_new << " " << a_new.GetTetraStereo() << std::endl;

          std::cout << nbrs_original_order << std::endl;

          nbrs_new_order[0] = coming_from_new[i];
          nbrs_new_order.insert(nbrs_new_order.end(),
                                ring_connections_new[i].begin(),
                                ring_connections_new[i].end());
          nbrs_new_order.insert(nbrs_new_order.end(), going_to_new[i].begin(),
                                going_to_new[i].end());

          for (auto&& n : nbrs_new_order) n = idx[n];
          std::cout << nbrs_new_order << std::endl;

          std::cout << NumPerm(nbrs_original_order, nbrs_new_order)
                    << std::endl;
          size_t n = NumPerm(nbrs_original_order, nbrs_new_order);
          if ((n % 2)) {
            a_new.SwapTetrahedralStereo();
            // a_new.SwapCTStereo();
          }  // odd
        }
      }
      if (stereo) {
        smiles_canon = smiles_gen.UpdateStereo();
      }

      const auto& visited_indices_(smiles_gen.GetVisitedIndices());
      std::cout << "new order of smiles: ";
      for (auto ii : visited_indices_) std::cout << idx[ii] << " ";
      std::cout << std::endl;
      std::cout << smiles_gen.UpdateStereo() << std::endl;
    }
  }
}

void SmilesHandler::PrintOutput(const std::string& smiles_orig,
                                const std::string& smiles_canon,
                                const std::string& fmi,
                                const cnv::AdjacencyMatrix& A) {
  size_t nUnsat(0), nMB(0), nDB(0), nAB(0), nTB(0), nSB(0), nQB(0), nCyc(0);

  if (print_options[cnv::print_n_double_bonds] ||
      print_options[cnv::print_n_aromatic_bonds] ||
      print_options[cnv::print_n_triple_bonds] ||
      print_options[cnv::print_n_unsaturations] ||
      print_options[cnv::print_n_cycles] ||
      print_options[cnv::print_n_multiple_bonds] ||
      print_options[cnv::print_n_quadruple_bonds]) {
    A.GetNumMultipleBonds(nSB, nDB, nTB, nQB, nAB);
    double DoU_(0);
    assert(!(nAB % 2));
    nMB = nDB + nTB + nAB / 2;

    for (auto&& atom : A.GetAtomVector()) DoU_ += (double)atom.GetDegree() - 2.;

    nUnsat = (size_t)floor((DoU_ / 2.) + 1);
    nCyc = nUnsat - nMB;
  }

  std::ostream* out;

  if (output_file.is_open())
    out = &output_file;

  else
    out = &std::cout;

  *out << /*"#" <<*/ std::setw((int)max_size) << std::left << smiles_orig
       << " ";

  if (print_options[cnv::print_canon_smiles])
    *out << std::setw(column_width) << std::left << smiles_canon << " ";

  if (print_options[cnv::print_canon_formula])
    *out << std::setw(column_width) << std::left
         << CreateCanonicalFormulaFromAtomVector(A.GetAtomVector()) << " ";

  if (print_options[cnv::print_mass])
    *out << std::setw(column_width) << std::left << A.GetMass() << " ";

  if (print_options[cnv::print_family_enumeration])
    *out << std::setw(column_width) << std::left << fmi << " ";

  if (print_options[cnv::print_n_unsaturations])
    *out << std::setw(column_width) << std::left << nUnsat << " ";

  if (print_options[cnv::print_n_multiple_bonds])
    *out << std::setw(column_width) << std::left << nMB << " ";

  if (print_options[cnv::print_n_double_bonds])
    *out << std::setw(column_width) << std::left << nDB << " ";

  if (print_options[cnv::print_n_aromatic_bonds])
    *out << std::setw(column_width) << std::left << nAB << " ";

  if (print_options[cnv::print_n_triple_bonds])
    *out << std::setw(column_width) << std::left << nTB << " ";

  if (print_options[cnv::print_n_quadruple_bonds])
    *out << std::setw(column_width) << std::left << nQB << " ";

  if (print_options[cnv::print_n_cycles])
    *out << std::setw(column_width) << std::left << nCyc << " ";

  if (print_options[cnv::print_canon_atom_vector])
    *out << A.GetAtomVector() << " ";

  if (print_options[cnv::print_stack]) *out << A.GetStack();

  if (print_options[cnv::print_matrix]) {
    *out << '\n';
    A.PrintToFile(*out);
  }

  *out << '\n';
}

}  // namespace cnv

}  // namespace combi_ff