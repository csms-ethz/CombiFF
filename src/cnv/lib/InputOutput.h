// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef INPUTOUTPUTCVN_H
#define INPUTOUTPUTCVN_H

#include <fstream>
#include <list>
#include <vector>

#include "StringVector.h"

namespace combi_ff {

namespace cnv {

typedef enum {
  smiles,
  formula,
  name,
  family_enumeration,
  matrix,
  not_set
} InputOption;
typedef enum {
  print_canon_smiles,
  print_canon_name,
  print_canon_formula,
  print_mass,
  print_family_enumeration,
  print_n_unsaturations,
  print_n_multiple_bonds,
  print_n_double_bonds,
  print_n_aromatic_bonds,
  print_n_triple_bonds,
  print_n_quadruple_bonds,
  print_n_cycles,
  print_canon_atom_vector,
  print_stack,
  print_matrix,
  num_print_options
} PrintOption;

class InputOutput {
 public:
  InputOutput(std::list<std::pair<std::string, std::string>>& input_list,
              std::vector<bool>& print_options, InputOption& input,
              std::string& output_file_name,
              combi_ff::StringVector& family_enumeration_file_names, int& argc,
              char* argv[]);

  void ReadArguments();
  void PrintInputOptions();
  void GetInputFromFile(const size_t pos);
  void ReadInputOption(const size_t pos,
                       const combi_ff::StringVector& arguments);
  void ReadOutputOption(const size_t pos,
                        const combi_ff::StringVector& arguments);
  void ReadFieFileNamesFromSetupFile(const std::string& setup_file_name);
  const combi_ff::StringVector& GetInputFileNames() const;
  const std::string IncompatibilityMessage(const std::string& input,
                                           const std::string& output) const;

 private:
  std::list<std::pair<std::string, std::string>>& input_list;
  combi_ff::StringVector arguments;
  std::vector<bool>& print_options;
  InputOption& input;
  std::string& output_file_name;
  StringVector& family_enumeration_file_names;
  StringVector input_file_names;
};

}  // namespace cnv

}  // namespace combi_ff

#endif
