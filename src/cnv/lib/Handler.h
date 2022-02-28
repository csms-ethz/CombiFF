// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef HANDLER_H
#define HANDLER_H

#include <fstream>
#include <list>
#include <vector>

#include "InputOutput.h"
#include "SmilesGenerator.h"

namespace combi_ff {

namespace cnv {

typedef combi_ff::SmilesGenerator<double, combi_ff::CnvAtom> SmilesGeneratorCnv;

class Handler {
 protected:
  std::ofstream& output_file;
  StringVector& fie_file_names;

  const int column_width;
  const size_t canon_iteration_limit = 100000;
  std::vector<bool>& print_options;
  std::list<std::pair<std::string, std::string>>
      input_list;  // pair: (input source, input)

 public:
  Handler(std::ofstream& output_file, StringVector& fie_file_names,
          const int column_width, std::vector<bool>& print_options,
          std::list<std::pair<std::string, std::string>>& input_list)
      : output_file(output_file),
        fie_file_names(fie_file_names),
        column_width(column_width),
        print_options(print_options),
        input_list(input_list) {}
  void Run(const cnv::PrintOption defaultOption = print_canon_smiles);
  void PrintFirstLine(const std::string& first);
};

}  // namespace cnv

}  // namespace combi_ff

#endif