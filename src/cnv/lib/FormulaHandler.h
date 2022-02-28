// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef FORMULAHANDLER_H_
#define FORMULAHANDLER_H_

#include "Atom.h"
#include "Handler.h"

namespace combi_ff {

namespace cnv {

class FormulaHandler : public cnv::Handler {
 public:
  FormulaHandler(std::ofstream& outputFile, StringVector& fieFileNames,
                 const int columnWidth, std::vector<bool>& printOptions,
                 std::list<std::pair<std::string, std::string>>& inputList)
      : cnv::Handler(outputFile, fieFileNames, columnWidth, printOptions,
                     inputList) {}

  void Run();
  void ConvertFormula(const std::string& formula);
  void PrintOutput(const std::string& formula,
                   const combi_ff::AtomVector<combi_ff::CnvAtom>& atoms);
  void PrintFirstLine();
};

}  // namespace cnv

}  // namespace combi_ff

#endif