// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef NAMEHANDLER_H_
#define NAMEHANDLER_H_

#include "Handler.h"

namespace combi_ff {

namespace cnv {

class NameHandler : public cnv::Handler {
 public:
  NameHandler(std::ofstream& outputFile, StringVector& fieFileNames,
              const int columnWidth, std::vector<bool>& printOptions,
              std::list<std::pair<std::string, std::string>>& inputList)
      : cnv::Handler(outputFile, fieFileNames, columnWidth, printOptions,
                     inputList) {}

  void Run();
  void ConvertName(const std::string& name);
  void PrintOutput(const std::string& nameOrig, const std::string& nameCan);
  void PrintFirstLine();
};

}  // namespace cnv

}  // namespace combi_ff

#endif