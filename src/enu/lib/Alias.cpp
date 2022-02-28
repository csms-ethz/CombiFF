// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "Alias.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

#include "ContainerOperators.h"
#include "exceptions.h"
#include "version.h"

namespace combi_ff {

namespace enu {

AtomTypeAlias::AtomTypeAlias(const std::string& alias_name,
                             const std::string& version)
    : alias_name(alias_name),
      contained_atom_types(std::vector<combi_ff::ElementSymbol>(0)),
      version(version) {}

const std::string& AtomTypeAlias::GetVersion() const { return version; }
const std::string& AtomTypeAlias::GetAliasName() const { return alias_name; }
std::vector<combi_ff::ElementSymbol>& AtomTypeAlias::GetContainedAtomTypes() {
  return contained_atom_types;
}

const std::vector<combi_ff::ElementSymbol>&
AtomTypeAlias::GetContainedAtomTypesConst() const {
  return contained_atom_types;
}
void AtomTypeAlias::SetVersion(const std::string& v) { version = v; }
void AtomTypeAlias::SetAliasName(const std::string& a) { alias_name = a; }
void AtomTypeAlias::AddAtomType(const combi_ff::ElementSymbol& e) {
  contained_atom_types.push_back(e);
}

void CreateAliases(AliasMap& alias_map,
                   const std::list<std::string>& alias_file_names) {
  if (!alias_file_names.size()) {
    throw combi_ff::input_warning("no alias file");
  }

  std::cout << "creating aliases\n";

  for (auto&& alias_file_name : alias_file_names) {
    XmlParserIn parser(alias_file_name, XmlParserIn::read_all);
    const XmlElement& aliases_root = parser.GetTree().GetRoot();
    alias_map.reserve(aliases_root.GetNumberOfChildren());
    aliases_root.CheckNumberOfChildren_atLeast(1);
    aliases_root.CheckAttribute("version");
    const std::string& version =
        aliases_root.attributes.find("version")->second;

    if (version != combi_ff::current_version)
      std::cout << "?Warning: currently running combi_ff version "
                << combi_ff::current_version << " but aliases file "
                << alias_file_name << " is version " << version << "\n";

    for (auto&& alias : aliases_root.children)
      GetNextAlias(alias_map, alias, version);
  }

  std::cout << "******************************************************\n";
}

/**************************************************************
READ THE ALIASES FROM THE BEGINNING OF THE FAMILY LIBRARY FILE
**************************************************************/
void GetNextAlias(AliasMap& alias_map, const XmlElement_ptr alias,
                  const std::string& version) {
  alias->CheckAttribute("name");
  alias->CheckAttributeSize(1);
  alias->CheckNumberOfChildren_atLeast(1);
  std::string map_name(alias->attributes.find("name")->second);

  if (alias_map.find(map_name) != alias_map.end())
    throw input_error("alias " + map_name + " occurs more than once\n");

  alias_map[map_name] = AtomTypeAlias(map_name, version);
  alias_map[map_name].GetContainedAtomTypes().reserve(
      alias->GetNumberOfChildren());

  for (auto&& member : alias->children) {
    member->CheckNumberOfChildren_equal(0);
    member->CheckAttributeSize(0);
    member->CheckValue();

    if (alias_map.find(member->value) != alias_map.end())
      alias_map[map_name].GetContainedAtomTypes().insert(
          alias_map[map_name].GetContainedAtomTypes().end(),
          alias_map[member->value].GetContainedAtomTypesConst().begin(),
          alias_map[member->value].GetContainedAtomTypesConst().end());

    else
      alias_map[map_name].AddAtomType(member->value);
  }

  // sortElements(alias_map[map_name]);
  std::cout << map_name << " -> ";

  for (auto&& a : alias_map[map_name].GetContainedAtomTypesConst())
    std::cout << a << " ";

  std::cout << std::endl;
  /* std::istringstream line(combi_ff::read(aliasesFile));
   std::string s(""), map("");

   //addd whole line to map, ignoring whitespaces, tabs, etc
   while(line >> s)
       map += s;

   size_t j = 0;

   while(j < map.size()) {
       std::string map_name("");

       //everything before '[' belongs to the map name
       while(map[j] != '[')
           map_name += map[j++];

       combi_ff::StringVector elements(0);
       std::string ele("");

       //go through map and read the elements, until ']' is encountered
       do {
           //case that map[j] belongs to an element name
           if(isalpha(map[j]))
               ele += map[j];

           //case that the next element is a united atom type
           else if(map[j] == '{') {
               if(ele != "")
                   throw combi_ff::input_error("didn't expect an alphabetical
   character before \'{\' in " + map);

               while(map[j] != '}')
                   ele += map[j++];

               ele += '}';
               elements.push_back(ele);
               ele = "";
           }

           //case that a comma or ']' is encountered between two element types
           else if(ele != "") {
               //Atom a(ele); //used to test if valid element type
               elements.push_back(ele);
               ele = "";
           }
       }
       while(map[j++] != ']');

       if(!elements.size())
           throw combi_ff::input_error("no alias found for " + map_name);

       //sort the element names alphabetically, and add a new alias
       Enu::sortElements(elements);
       alias_map.emplace(map_name, elements);
       std::cout << "  -" << std::setw(5) << std::left << map_name << " [" <<
   elements << "]\n";
   }*/
}

/************************************************************************************
SORT THE ELEMENTS IN THE VECTOR OF A ALIASES FIRST BY DEGREE, AND THEN
ALPHABETICALLY
************************************************************************************/
void SortElements(combi_ff::StringVector& elements) {
  std::sort(elements.begin(), elements.end(), enu::CompareElems);
  auto&& it = std::unique(elements.begin(), elements.end());
  elements.resize(std::distance(elements.begin(), it));
}

bool CompareElems(std::string e1, std::string e2) {
  if (e1 == "H") return false;

  if (e2 == "H") return false;

  if (e1[0] == '{') e1 = e1.substr(1, e1.size() - 1);

  if (e2[0] == '{') e2 = e2.substr(1, e2.size() - 1);

  if (e1 < e2)
    return true;

  else
    return false;
}

}  // namespace enu

}  // namespace combi_ff
