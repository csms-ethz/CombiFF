// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef FAMILYIDENTIFIERHANDLER_H_
#define FAMILYIDENTIFIERHANDLER_H_

#include "AdjacencyMatrixCnv.h"
#include "Handler.h"

namespace combi_ff {

namespace cnv {

class FamilyIdentifierHandler : public cnv::Handler {
 public:
  FamilyIdentifierHandler(
      std::ofstream& outputFile, StringVector& fieFileNames,
      const int columnWidth, std::vector<bool>& printOptions,
      std::list<std::pair<std::string, std::string>>& inputList)
      : cnv::Handler(outputFile, fieFileNames, columnWidth, printOptions,
                     inputList) {}

  void Run();
  void PrintFirstLine();
  void ConvertFamilyIdentifier(const std::string& fmi);
  void PrintOutput(const std::string& fmi, const std::string& formula,
                   const std::string& smiles, cnv::AdjacencyMatrix& A);
  void SearchFileForIsomer(const std::string& fmi,
                           const std::string& familyFileName);
};

}  // namespace cnv

}  // namespace combi_ff

#endif