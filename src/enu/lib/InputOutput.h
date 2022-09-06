// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

#include <fstream>

#include "EnumSpecification.h"
#include "InputFiles.h"  //define new input file types here

namespace combi_ff {

namespace enu {

typedef std::vector<std::list<std::string>> FileNameVector;

struct IOFileProperties {
  static const std::string output_file_name_default;
  const size_t rnd;  // random number to distinguish IO files if program is run
                     // for different inputs simultaneously in same directory
  FileNameVector input_file_names;
  std::string output_dir{""};
  std::string output_file_name{output_file_name_default};
  const std::string file_name_tmp{".enu_temp_" /*+ std::to_string(rnd)*/};
  StringVector input_file_categories;
  StringVector input_file_name_building_blocks;

  IOFileProperties();
};

class InputOutput {
 public:
  explicit InputOutput(const StringVector& arguments);
  ~InputOutput();

  void ReadInputArguments();
  void GetNextInputOption(size_t& i);

  void AddInputFile(size_t& i);
  void AddAtoms(size_t& i);
  void AddFamilies(size_t& i);
  void AddOutputFile(size_t& i);
  void AddOutputDir(size_t& i);
  void AddMaxDegree(size_t& i);
  void AddRestriction(size_t& i, const size_t position_in_range_vec);
  void AddStereo();
  void AddCountOnly();
  void ReadFileNames(size_t& i, const InputFileType input_type);
  void AddAtom(size_t& j, std::string& formula,
               AtomVector<combi_ff::Atom>& used_atoms);
  void AddUnitedAtom(size_t& j, std::string& formula,
                     AtomVector<combi_ff::Atom>& used_atoms);
  std::list<std::string>& GetInputFileNamesAt(size_t idx);
  const IOFileProperties& GetIOFilProps() const;
  const std::string& GetOutputFileName() const;
  const EnumSpecifications& GetEnumSpec() const;

  void PrintInputOptions();

 private:
  IOFileProperties io_file_properties;
  EnumSpecifications enum_spec;
  StringVector arguments;
};

}  // namespace enu

}  // namespace combi_ff

#endif
