#ifndef INPUTFILES_H_
#define INPUTFILES_H_

namespace combi_ff {

namespace enu {

//add new input files here (before num_input_files!)
typedef enum {family_file, substructure_file, alias_file, pseudoatom_file, num_input_files} InputFileType;

//compatibility with g++ 4.8
class InputFileTypeHash {
 public:
  size_t operator()(const InputFileType& type) const{
    return std::hash<size_t>()(type);
  }
};

struct InputFileProperties {
	std::string category, name_block;

	InputFileProperties(const std::string& category, const std::string& name_block) :
		category(category), name_block(name_block) {}
};

typedef std::unordered_map<std::string, InputFileType> InputFileMap;
typedef std::unordered_map<InputFileType, InputFileProperties, InputFileTypeHash> InputFilePropertiesMap;

//add new input files here
static const InputFilePropertiesMap possible_input_file_properties{
	{family_file, {"family Set file (fst)", "familySetFile"}},
	{alias_file, {"alias Set file (ast)", "aliasSetFile"}},
	{substructure_file, {"substructure Set file (sst)", "substructureSetFile"}},
	{pseudoatom_file, {"pseudo atom Set file (pst)", "pseudoAtomSetFile"}}};

static const InputFileMap possible_input_files{
	{"-family_files", family_file},
	{"-element_alias_files", alias_file},
	{"-substructure_files", substructure_file},
	{"-pseudoatom_files", pseudoatom_file}};

} //namespace enu

} //namespace combi_ff



#endif
