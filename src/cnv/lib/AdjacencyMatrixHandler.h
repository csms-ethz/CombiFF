// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef ADJACENCYMATRIXHANDLER_H_
#define ADJACENCYMATRIXHANDLER_H_

#include "AdjacencyMatrixCnv.h"
#include "Handler.h"

namespace combi_ff {

namespace cnv {

class AdjacencyMatrixHandler : public cnv::Handler {
 public:
  AdjacencyMatrixHandler(
      std::ofstream& outputFile, combi_ff::StringVector& fieFileNames,
      const int columnWidth, std::vector<bool>& printOptions,
      std::list<std::pair<std::string, std::string>>& inputList)
      : cnv::Handler(outputFile, fieFileNames, columnWidth, printOptions,
                     inputList) {}

  void Run();
  void PrintFirstLine();
  cnv::AdjacencyMatrix ConvertMatrix(
      const cnv::AdjacencyMatrix& AdjacencyMatrixOrig, std::string& smilesCan);
  void PrintOutput(const cnv::AdjacencyMatrix& AdjacencyMatrixOrig,
                   const std::string& smilesCan, const std::string& fmi,
                   const AdjacencyMatrix& A);
  void CreateNameVector(std::vector<combi_ff::SmilesBlock>& names,
                        const std::string& smiles);
  void CreateMatrixFromNameVector(
      cnv::AdjacencyMatrix& A, const std::vector<combi_ff::SmilesBlock>& names,
      const size_t nR);
  void findFamilyIdentifier(const std::string& smilesCan, std::string& fmi);
};

}  // namespace cnv

}  // namespace combi_ff

#endif