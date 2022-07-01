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
  std::string num("");
  std::string atom_name("");

  while (i < formula.size()) {
    if (std::isalpha(formula[i])) {
      if (std::isupper(formula[i])) {
        atom_name = formula[i];

        if (i + 1 < formula.size() && std::islower(formula[i + 1]))
          atom_name += formula[++i];

      } else if (std::islower(formula[i]))
        atom_name = std::string(1, (char)std::toupper(formula[i]));

      i++;

      if (i == formula.size() || !std::isdigit(formula[i]))
        num = "1";

      else {
        num = "";
        while (i < formula.size() && std::isdigit(formula[i]))
          num += formula[i++];
      }

      AtomVector<combi_ff::CnvAtom> atmp(std::stoi(num), CnvAtom(atom_name));
      atoms.insert(atoms.end(), atmp.begin(), atmp.end());

    } else if (formula[i] == '(') {
      combi_ff::AtomVector<combi_ff::CnvAtom> atoms_in_parentheses(0);
      i++;

      while (i < formula.size() && formula[i] != ')') {
        if (std::isalpha(formula[i])) {
          if (std::isupper(formula[i])) {
            atom_name = formula[i];

            if (i + 1 < formula.size() && std::islower(formula[i + 1]))
              atom_name += formula[++i];

          } else if (std::islower(formula[i]))
            atom_name = std::string(1, (char)std::toupper(formula[i]));

        } else
          throw combi_ff::input_error("unexpected character" +
                                      std::string(1, formula[i]) + " in " +
                                      formula);

        i++;

        if (i == formula.size() || !std::isdigit(formula[i]))
          num = "1";

        else {
          num = "";
          while (i < formula.size() && std::isdigit(formula[i]))
            num += formula[i++];
        }

        AtomVector<combi_ff::CnvAtom> atmp(std::stoi(num), CnvAtom(atom_name));
        atoms_in_parentheses.insert(atoms_in_parentheses.end(), atmp.begin(),
                                    atmp.end());
      }

      if (i == formula.size() || formula[i] != ')')
        throw combi_ff::input_error("unclosed opening parenthesis ( in " +
                                    formula);

      i++;

      if (i == formula.size() || !std::isdigit(formula[i]))
        num = "1";

      else {
        num = "";
        while (i < formula.size() && std::isdigit(formula[i]))
          num += formula[i++];
      }

      for (int j = 0; j < std::stoi(num); j++)
        atoms.insert(atoms.end(), atoms_in_parentheses.begin(),
                     atoms_in_parentheses.end());

    } else
      throw combi_ff::input_error(
          "unexpected character" + std::string(1, formula[i]) + " in " +
          formula + ". Please add it in src/cnv/lib/FormulaHandler.cpp.");
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