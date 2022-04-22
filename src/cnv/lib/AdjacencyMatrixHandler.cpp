// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "AdjacencyMatrixHandler.h"

#include <cmath>
#include <sstream>

namespace combi_ff {

namespace cnv {

void AdjacencyMatrixHandler::Run() {
  cnv::Handler::Run();
  std::string line;
  std::vector<cnv::AdjacencyMatrix> matrices(0);
  cnv::AdjacencyMatrix A;
  std::vector<std::string> input_file_names(0);

  for (auto it = input_list.begin(); it != input_list.end(); ++it) {
    std::istringstream s(it->second);
    AtomVector<combi_ff::CnvAtom> atoms(0);
    std::vector<double> v(0);

    if (it->first == "command_line") {
      bool isalpha = true;

      while (isalpha) {
        if (it->second.size()) {
          std::string atom_name = (it++)->second;
          std::string atom_type("");

          for (size_t i = 0; i < atom_name.size(); i++) {
            if (std::isalpha(atom_name[i]))
              atom_type += atom_name[i];

            else
              break;
          }

          if (atom_type.size())
            atoms.push_back(combi_ff::CnvAtom(atom_name, atom_type));

          else
            isalpha = false;
        }
      }

      --it;
      double degree;

      for (size_t i = 0; i < atoms.size() * atoms.size(); i++) {
        if (it == input_list.end())
          throw combi_ff::input_error(
              "not enough matrix elements " + std::to_string(v.size()) +
              " for number of atoms " + std::to_string(atoms.size()));

        degree = std::stod((it++)->second);

        if (degree == 1.5)
          throw combi_ff::input_error(
              "Please specify aromatic bonds in adjacency matrices as "
              "alternating single and double bonds for now.");

        v.push_back(degree);
      }

      it--;

    } else {
      if (it->second.size()) {
        std::string atom_name;

        while (s >> atom_name) {
          std::string atom_type("");

          for (size_t i = 0; i < atom_name.size(); i++) {
            if (isalpha(atom_name[i]))
              atom_type += atom_name[i];

            else
              break;
          }

          atoms.push_back(combi_ff::CnvAtom(atom_name, atom_type));
        }

        for (size_t j = 0; j < atoms.size(); j++) {
          std::istringstream ss((++it)->second);
          double degree;

          while (ss >> degree) {
            if (degree == 1.5)
              throw combi_ff::input_error(
                  "Please specify aromatic bonds in adjacency matrices as "
                  "alternating single and double bonds for now.");

            v.push_back(degree);
          }
        }
      }
    }

    if (v.size() != atoms.size() * atoms.size())
      throw combi_ff::input_error("number of matrix elements " +
                                  std::to_string(v.size()) +
                                  " not compatible with number of atoms " +
                                  std::to_string(atoms.size()));

    for (size_t i = 0; i < atoms.size(); i++) {
      if (v[i * atoms.size() + i] != 0)
        throw combi_ff::input_error(
            "adjacency matrix diagonal has non-zero elements");

      for (size_t j = i + 1; j < atoms.size(); j++) {
        if (v[i * atoms.size() + j] != v[j * atoms.size() + i])
          throw combi_ff::input_error(
              "adjacency matrix elements are not symmetrical");
      }
    }

    input_file_names.push_back(it->first);
    matrices.push_back(cnv::AdjacencyMatrix(atoms.size()));
    matrices.back().SetElements(v);
    matrices.back().SetAtomVector(atoms);
  }

  for (size_t i = 0; i < matrices.size(); i++) {
    const std::string& fileName = input_file_names[i];
    cnv::AdjacencyMatrix& AdjacencyMatrixOrig = matrices[i];
    std::string canon_smiles("");
    // std::cout << "\n\n";
    // std::cout << "original matrix:\n";
    // AdjacencyMatrixOrig.print();
    cnv::AdjacencyMatrix A = AdjacencyMatrixHandler::ConvertMatrix(
        AdjacencyMatrixOrig, canon_smiles);
    std::string fmi("");

    if (print_options[print_family_enumeration])
      AdjacencyMatrixHandler::findFamilyIdentifier(canon_smiles, fmi);

    std::cout << "original matrix:\n";
    AdjacencyMatrixOrig.Print();
    std::cout << '\n';
    std::cout << "canonical matrix:\n";
    A.Print();
    std::cout << "input_file:            " << fileName << '\n';
    AdjacencyMatrixHandler::PrintOutput(AdjacencyMatrixOrig, canon_smiles, fmi,
                                        A);
  }
}

void AdjacencyMatrixHandler::findFamilyIdentifier(
    const std::string& canon_smiles, std::string& fmi) {
  bool found(false);
  std::string line, code, formula, smiles_cur;

  for (auto&& familyFileName : fie_file_names) {
    std::ifstream familyFile(familyFileName);

    if (!familyFile.is_open()) {
      std::cout << "!Error: family isomer enumeration file " << familyFileName
                << " could not be opened\n";
      exit(-1);
    }

    while (std::getline(familyFile, line)) {
      if (line[0] != '#') {
        std::istringstream s(line);
        s >> code >> formula >> smiles_cur;

        if (smiles_cur == canon_smiles) {
          found = true;
          break;
        }
      }
    }

    familyFile.close();

    if (found) break;
  }

  if (!found) {
    std::cout << "!Error: smiles " << canon_smiles
              << " not found in any of the following family isomer enumeration "
                 "files ";

    for (auto&& ffn : fie_file_names) std::cout << ffn << " ";

    std::cout << '\n';
    exit(-1);

  } else
    fmi = code;
}

void AdjacencyMatrixHandler::PrintFirstLine() {
  cnv::Handler::PrintFirstLine("# originalMatrix");
}

AdjacencyMatrix AdjacencyMatrixHandler::ConvertMatrix(
    const AdjacencyMatrix& AdjacencyMatrixOrig, std::string& canon_smiles) {
  cnv::AdjacencyMatrix A(AdjacencyMatrixOrig);
  canon_smiles = "";

  if (print_options[print_canon_smiles] || print_options[print_stack] ||
      print_options[print_matrix] || print_options[print_family_enumeration] ||
      print_options[print_canon_atom_vector]) {
    std::vector<size_t> idx = A.SortAtomVector(A.GetIsAromaticCarbon());
    combi_ff::LambdaVector lambda = A.GetLambda();
    size_t N = A.GetN();
    // num_perms[i] indicates, with how many other atom indices the index of
    // atom i can be permuted note: index permutations are only allowed within
    // the same lambda partition, and with higher indices e.g. for lambda = [1,
    // 3, 1, 4, 6], we have num_perms = [1, 3, 2, 1, 1, 4, 3, 2, 1, 6, 5, 4, 3,
    // 2, 1]
    std::vector<size_t> num_perms;
    int ind(0);

    for (size_t i = 0; i < lambda.size(); i++) {
      for (size_t j = 0; j < lambda[i]; j++)
        num_perms.push_back(
            std::accumulate(lambda.begin(), lambda.begin() + i + 1, 0) - ind++);
    }

    // a RepresentationSystem contains all the possible index permutations for
    // the different atom indices
    //  e.g. for lambda = [2, 1, 4], the RepresentationSystem is
    //  (0,0), (0,1)        (for idx 0)
    //  (1,1)           (for idx 1)
    //  (2,2)           (for idx 2)
    //  (3,3), (3,4), (3,5), (3,6)  (for idx 3)
    //  (4,4), (4,5), (4,6)     (for idx 4)
    //  (5,5), (5,6)          (for idx 5)
    //  (6,6)             (for idx 6)
    combi_ff::RepresentationSystem u(N);

    for (size_t i = 0; i < N; i++) {
      u[i].reserve(num_perms[i]);

      for (size_t j = 0; j < num_perms[i]; j++) {
        u[i].push_back(Permutations(1, {i, i + j}));
      }
    }

    A.MakeCanonical(u, idx, canon_iteration_limit);
    combi_ff::cnv::SmilesGeneratorCnv smiles_gen(A);
    smiles_gen.GenerateSmiles();
    canon_smiles = smiles_gen.GetSmiles();
  }

  return A;
}

void AdjacencyMatrixHandler::PrintOutput(
    const cnv::AdjacencyMatrix& AdjacencyMatrixOrig,
    const std::string& canon_smiles, const std::string& fmi,
    const cnv::AdjacencyMatrix& A) {
  int nUnsat(0);
  size_t nB(0), nMB(0), nDB(0), nAB(0), nTB(0), nSB(0), nQB(0), nCyc(0);

  if (print_options[cnv::print_n_bonds] ||
      print_options[cnv::print_n_single_bonds] ||
      print_options[cnv::print_n_double_bonds] ||
      print_options[cnv::print_n_triple_bonds] ||
      print_options[cnv::print_n_unsaturations] ||
      print_options[cnv::print_n_cycles] ||
      print_options[cnv::print_n_multiple_bonds]) {
    A.GetNumMultipleBonds(nSB, nDB, nTB, nQB, nAB);
    double DoU_(0);
    // assert(!(nAB % 2));
    nMB = nDB + nTB + nAB / 2 + nQB;

    for (auto&& atom : A.GetAtomVector()) DoU_ += (double)atom.GetDegree() - 2.;

    nUnsat = (int)floor((DoU_ / 2.) + 1);
    nCyc = nUnsat - nMB;
    nB = nSB + nDB + nTB + nAB / 2 + nQB;
  }

  std::ostream* out;

  if (output_file.is_open())
    out = &output_file;

  else
    out = &std::cout;

  // AdjacencyMatrixOrig.printToFile(*out);

  if (print_options[cnv::print_canon_smiles])
    *out << std::setw(column_width) << std::left << canon_smiles << " ";

  if (print_options[cnv::print_canon_formula])
    *out << std::setw(column_width) << std::left
         << CreateCanonicalFormulaFromAtomVector<combi_ff::CnvAtom>(
                A.GetAtomVector())
         << " ";

  if (print_options[cnv::print_mass]) {
    double mass = 0;

    for (auto&& a : A.GetAtomVector()) mass += a.GetMass();

    *out << std::setw(column_width) << std::left << mass << " ";
  }

  if (print_options[cnv::print_family_enumeration])
    *out << std::setw(column_width) << std::left << fmi << " ";

  if (print_options[cnv::print_n_unsaturations])
    *out << std::setw(column_width) << std::left << nUnsat << " ";

  if (print_options[cnv::print_n_bonds])
    *out << std::setw(column_width) << std::left << nB << " ";

  if (print_options[cnv::print_n_single_bonds])
    *out << std::setw(column_width) << std::left << nSB << " ";

  if (print_options[cnv::print_n_multiple_bonds])
    *out << std::setw(column_width) << std::left << nMB << " ";

  if (print_options[cnv::print_n_double_bonds])
    *out << std::setw(column_width) << std::left << nDB << " ";

  if (print_options[cnv::print_n_triple_bonds])
    *out << std::setw(column_width) << std::left << nTB << " ";

  if (print_options[cnv::print_n_cycles])
    *out << std::setw(column_width) << std::left << nCyc << " ";

  if (print_options[cnv::print_canon_atom_vector])
    *out << "canonical_atom_vector: " << A.GetAtomVector() << " ";

  if (print_options[cnv::print_stack]) *out << A.GetStack();

  if (print_options[cnv::print_matrix]) {
    *out << '\n';
    A.PrintToFile(*out);
  }

  *out << '\n';
}

}  // namespace cnv

}  // namespace combi_ff