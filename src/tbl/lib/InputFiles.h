// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef INPUTFILESTBL_H_
#define INPUTFILESTBL_H_

namespace combi_ff {

namespace topology_builder {

// add new input files here (before num_input_files!)
typedef enum { fragment_file, atomtypes_file, num_input_files } InputFileType;

// compatibility with g++ 4.8
class InputFileTypeHash {
 public:
  size_t operator()(const InputFileType& type) const {
    return std::hash<size_t>()(type);
  }
};

struct inputFileProps {
  std::string category, name_block;

  inputFileProps(const std::string& category, const std::string& name_block)
      : category(category), name_block(name_block) {}
};

typedef std::unordered_map<std::string, InputFileType> InputFileMap;
typedef std::unordered_map<InputFileType, inputFileProps, InputFileTypeHash>
    InputFilePropsMap;

// add new input files here
static const InputFilePropsMap possible_input_file_props{
    {fragment_file, {"fragment file (frg)", "fragmentFile"}},
    {atomtypes_file, {"atomtypes file (att/s)", "atomTypesFile"}}};

static const InputFileMap possible_input_files{{"-fragments", fragment_file},
                                               {"-atomtypes", atomtypes_file}};
//{"@fieFamilies", fiefamilyfile}};
}  // namespace topology_builder

}  // namespace combi_ff

#endif
