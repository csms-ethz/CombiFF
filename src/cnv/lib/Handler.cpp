// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "Handler.h"

#include <iomanip>
#include <iostream>
#include <numeric>

namespace combi_ff {

namespace cnv {

void Handler::Run(const cnv::PrintOption defaultOption) {
  if (!std::accumulate(print_options.begin(), print_options.end(), 0)) {
    std::cout << "# no output arguments -O found. choosing to print ";

    if (defaultOption == print_canon_name)
      std::cout << "the canonical name";

    else if (defaultOption == print_canon_formula)
      std::cout << "the canonical formula";

    else  // default
      std::cout << "the canonical smiles string";

    std::cout << " as a default\n";
    print_options[defaultOption] = true;
  }
}

void Handler::PrintFirstLine(const std::string& first) {
  std::ostream* out;

  if (output_file.is_open())
    out = &output_file;

  else
    out = &std::cout;

  *out << std::setw(column_width) << std::left << first << " ";

  if (print_options[cnv::print_canon_smiles])
    *out << std::setw(column_width + 1) << std::left << "canonicalSmiles ";

  if (print_options[cnv::print_canon_formula])
    *out << std::setw(column_width + 1) << std::left << "canonicalFormula ";

  if (print_options[cnv::print_canon_name])
    *out << std::setw(column_width + 1) << std::left << "canonicalName ";

  if (print_options[cnv::print_mass])
    *out << std::setw(column_width + 1) << std::left << "molecularMass ";

  if (print_options[cnv::print_family_enumeration])
    *out << std::setw(column_width + 1) << std::left << "familyIdentifier ";

  if (print_options[cnv::print_n_unsaturations])
    *out << std::setw(column_width + 1) << std::left << "numUnsaturations ";

  if (print_options[cnv::print_n_multiple_bonds])
    *out << std::setw(column_width + 1) << std::left << "numMultipleBonds ";

  if (print_options[cnv::print_n_double_bonds])
    *out << std::setw(column_width + 1) << std::left << "numDoubleBonds ";

  if (print_options[cnv::print_n_aromatic_bonds])
    *out << std::setw(column_width + 1) << std::left << "numAromaticBonds ";

  if (print_options[cnv::print_n_triple_bonds])
    *out << std::setw(column_width + 1) << std::left << "numTripleBonds ";

  if (print_options[cnv::print_n_quadruple_bonds])
    *out << std::setw(column_width + 1) << std::left << "numQuadrupleBonds ";

  if (print_options[cnv::print_n_cycles])
    *out << std::setw(column_width + 1) << std::left << "numCycles ";

  if (print_options[cnv::print_canon_atom_vector])
    *out << std::setw(column_width + 1) << std::left << "canonicalAtomVector ";

  if (print_options[cnv::print_stack])
    *out << std::setw(column_width + 1) << std::left << "canonicalStack ";

  if (print_options[cnv::print_matrix])
    *out << " (canonical adjacency matrix is printed in next lines)";

  *out << '\n';
}

}  // namespace cnv

}  // namespace combi_ff