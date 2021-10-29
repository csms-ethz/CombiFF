#ifndef INPUTOUTPUTTBL_H
#define INPUTPUTPUTTBL_H

#include <fstream>
#include <vector>
#include <unordered_map>
#include "InputFiles.h" //define new input file types here
#include "StringVector.h"
#include <list>
#include <iostream>


namespace combi_ff {

namespace topology_builder {

typedef std::vector<std::list<std::string>> FileNameVector;

struct IOFileProperties {
  StringVector fie_families {StringVector(0)};
  FileNameVector input_file_names {FileNameVector(0)};
  std::string output_dir {""};
  std::string output_dir_molecule_decompositions {""};
  std::string output_dir_molecules_with_macros {""};
  std::string output_dir_mtb {""};
  StringVector input_file_categories {StringVector(0)};
  StringVector input_file_name_building_blocks {StringVector(0)};
  bool united_atom {false};
  bool third_neighbor_exclusions {false};
  bool unique_torsionals {false};


  IOFileProperties();
};

class InputOutput {
 private:
  IOFileProperties io_file_properties;
  std::vector<std::string> fragments_to_use;
  std::vector<std::string> arguments;
  size_t i;

 public:
  InputOutput(const std::vector<std::string>& arguments);

  void ReadInput();
  void ReadArguments() ;
  void ReadFileNames(size_t& i, const InputFileType input_type);
  void PrintInputOptions();
  void GetInputFile();
  const std::vector<std::string> GetFragmentsToUse() const;
  const IOFileProperties& GetIOFilProps() const ;
};

} //namespace topology_builder

} //namespace combi_ff

#endif
