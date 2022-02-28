// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "FormulaHandler.h"

#include <cmath>
#include <iomanip>

namespace combi_ff {

namespace cnv {

void FormulaHandler::Run() {
  cnv::Handler::Run(print_canon_formula);

  // determine the different formulas from the input_list
  for (auto it = input_list.begin(); it != input_list.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); i++) {
      if ((it->second)[i] == ' ' || (it->second)[i] == ',') {
        auto next = ++it;
        it--;
        input_list.insert(
            next, {it->first, it->second.substr(i + 1, it->second.size() - i)});
        it->second = it->second.substr(0, i);
      }
    }
  }

  for (auto&& formula : input_list) {
    if (formula.second.size()) FormulaHandler::ConvertFormula(formula.second);
  }
}

void FormulaHandler::ConvertFormula(const std::string& formula) {
  combi_ff::AtomVector<combi_ff::CnvAtom> atoms(0);
  size_t i = 0;

  while (i < formula.size()) {
    std::string atom_name("");
    std::string num("");
    // assert(std::isalpha(formula[i]));

    // while(i < formula.size() && std::isalpha(formula[i]))
    //	atom_name += formula[i++];

    if (formula[i] == 'C' || formula[i] == 'c') {
      if (i + 1 < formula.size() &&
          (formula[i + 1] == 'l' || formula[i + 1] == 'L')) {
        i++;
        atom_name = "Cl";

      } else
        atom_name = "C";

    } else if (formula[i] == 'H' || formula[i] == 'h')
      atom_name = "H";

    else if (formula[i] == 'O' || formula[i] == 'o')
      atom_name = "O";

    else if (formula[i] == 'N' || formula[i] == 'n')
      atom_name = "N";

    else if (formula[i] == 'S' || formula[i] == 's')
      atom_name = "S";

    else if (formula[i] == 'P' || formula[i] == 'p')
      atom_name = "P";

    else if (formula[i] == 'B' || formula[i] == 'r') {
      if (formula[++i] != 'r' && formula[i] != 'R')
        throw combi_ff::input_error("unrecognized atom type in " + formula);

      atom_name = "Br";

    } else if (formula[i] == 'F' || formula[i] == 'f')
      atom_name = "F";

    else if (formula[i] == 'I' || formula[i] == 'i')
      atom_name = "I";

    else
      throw combi_ff::input_error("unrecognized atom type in " + formula);

    i++;

    if (!std::isdigit(formula[i]))
      num = "1";

    else {
      while (i < formula.size() && std::isdigit(formula[i]))
        num += formula[i++];
    }

    AtomVector<combi_ff::CnvAtom> atmp(std::stoi(num), CnvAtom(atom_name));
    atoms.insert(atoms.end(), atmp.begin(), atmp.end());
  }

  FormulaHandler::PrintOutput(formula, atoms);
}

void FormulaHandler::PrintOutput(const std::string& formula_orig,
                                 const AtomVector<combi_ff::CnvAtom>& atoms) {
  size_t n_unsaturations(0);

  if (print_options[cnv::print_n_unsaturations]) {
    double unsat(0);

    for (auto&& atom : atoms) unsat += (double)atom.GetDegree() - 2.;

    n_unsaturations = (size_t)floor((unsat / 2.) + 1);
  }

  std::ostream* out;

  if (output_file.is_open())
    out = &output_file;

  else
    out = &std::cout;

  *out << std::setw(column_width) << std::left << formula_orig << " ";

  if (print_options[cnv::print_canon_formula])
    *out << std::setw(column_width) << std::left
         << CreateCanonicalFormulaFromAtomVector<combi_ff::CnvAtom>(atoms)
         << " ";

  if (print_options[cnv::print_mass]) {
    double mass = 0;

    for (auto&& a : atoms) mass += a.GetMass();

    *out << std::setw(column_width) << std::left << mass << " ";
  }

  if (print_options[cnv::print_n_unsaturations])
    *out << std::setw(column_width) << std::left << n_unsaturations << " ";

  *out << '\n';
}

void FormulaHandler::PrintFirstLine() {
  cnv::Handler::PrintFirstLine("# originalFormula");
}

}  // namespace cnv

}  // namespace combi_ff