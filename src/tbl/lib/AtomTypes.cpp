// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "AtomTypes.h"

#include <sstream>

#include "version.h"

namespace combi_ff {

namespace topology_builder {

AtomType::AtomType(const std::string& version,
                   const std::string& atom_type_name,
                   const std::string& element_type_name)
    : version(version),
      atom_type_name(atom_type_name),
      element_type(element_type_name) {}

const Atom& AtomType::GetElementType() const { return element_type; }

const std::string& AtomType::GetAtomTypeName() const { return atom_type_name; }

const std::string& AtomType::GetVersion() const { return version; }

size_t AtomHash::operator()(const Atom& atom) const {
  return std::hash<std::string>()(atom.GetUnitedAtomSymbol());
}

AtomTypeSet::AtomTypeSet() : atom_type_set_name(""), contained_atoms(0) {}

AtomTypeSet::AtomTypeSet(const std::string& version,
                         const std::string& set_name, size_t index)
    : version(version),
      atom_type_set_name(set_name),
      contained_atoms(0),
      index(index) {}

bool AtomTypeSet::operator==(const AtomTypeSet& at) const {
  return atom_type_set_name == at.GetName();
}

bool AtomTypeSet::operator<(const AtomTypeSet& at) const {
  return index < at.GetIndex();
}

const AtomSet& AtomTypeSet::GetAtomSet() const { return contained_atoms; }

const std::string& AtomTypeSet::GetName() const { return atom_type_set_name; }

const std::string& AtomTypeSet::GetVersion() const { return version; }

size_t AtomTypeSet::GetIndex() const { return index; }

void AtomTypeSet::SetName(const std::string& name) {
  atom_type_set_name = name;
}

void AtomTypeSet::SetVersion(const std::string& v) { version = v; }

void AtomTypeSet::AddAtom(const Atom& atom) { contained_atoms.emplace(atom); }

void AtomTypeSet::AddAtomType(const std::string& name) {
  atom_type_names.emplace(name);
}

void CreateAtomTypeSets(const std::list<std::string>& atom_types_file_names,
                        AtomTypeSetMap& atom_type_set) {
  if (!atom_types_file_names.size()) {
    std::cout << "?Warning: not using any atomTypes files\n";
    return;
  }

  atom_type_set["any"].SetName("any");
  atom_type_set["any"].SetVersion(combi_ff::current_version);

  for (auto&& file_name : atom_types_file_names) {
    XmlParserIn parser(file_name, XmlParserIn::read_all);
    // std::cout << parser << std::endl;
    const XmlTree& tree = parser.GetTree_const();
    const XmlElement& root = tree.GetRoot_const();
    tree.CheckRootTagName("atomtypes");
    root.CheckAttribute("version");
    const std::string& version = root.attributes.find("version")->second;

    if (version != combi_ff::current_version)
      std::cout << "?Warning: currently running combi_ff version "
                << combi_ff::current_version << " but atom type file "
                << file_name << " is version " << version << "\n";

    size_t idx = 0;

    for (auto&& child : root.children) {
      if (child->tag == "simpleatomtypes") {
        std::cout << "IDX " << idx << std::endl;
        AddSimpleAtomTypes(atom_type_set, child, version, idx);

      } else if (child->tag == "atomtypesets")
        AddAtomTypeSets(atom_type_set, child, version, idx);

      else
        throw unexpected_tag_error(child->tag,
                                   "simpleatomtypes or atomtypesets");
    }
  }
}

void AddSimpleAtomTypes(AtomTypeSetMap& atom_type_set,
                        const XmlElement_ptr simpleatomtypes,
                        const std::string& version, size_t& idx) {
  for (auto&& atomtype : simpleatomtypes->children) {
    CheckAtomTypeXML(atomtype);
    std::string id = atomtype->attributes["id"];
    std::string element_type = atomtype->children.front()->value;
    std::string atom_name = element_type;
    std::string num_H("");

    if (id == "any")
      throw input_error(
          "atom type Set 'any' is reserved keyword, cannot be explicitly "
          "listed in input\n");

    if (std::isdigit(atom_name.back()) && atom_name.size() > 2 &&
        atom_name[atom_name.size() - 2] == 'H') {
      num_H = atom_name.back();
      atom_name = atom_name.substr(0, atom_name.size() - 2);
    }

    Atom atom(atom_name);

    if (num_H.size()) atom.SetNumFixedHydrogens(std::stoi(num_H));

    if (atom_type_set.find(id) != atom_type_set.end())
      throw input_error("listed atom type " + id + " more than once\n");

    // simpleAtomTypes.emplace(id, atom);
    atom_type_set[id] = AtomTypeSet(version, id, idx++);
    // atom_type_set[att.first].SetName(att.first);
    atom_type_set[id].AddAtom(atom);
    atom_type_set[id].AddAtomType(id);
    atom_type_set["any"].AddAtom(atom);
    atom_type_set["any"].AddAtomType(id);
  }

  // add the simple atom types to the atom_type_set
  // for(auto && att : simpleAtomTypes) {
  // }
  /*std::cout << "simpleatomtypes:\n";

  for(auto && att : simpleAtomTypes)
        std::cout << att.first << " " << att.second << "\n";
  */
}

void AddAtomTypeSets(AtomTypeSetMap& atom_type_set,
                     const XmlElement_ptr atomtype_set_ptrs,
                     const std::string& version, size_t& idx) {
  for (auto&& atomtype_set_ptr : atomtype_set_ptrs->children) {
    CheckAtomTypeSetXML(atomtype_set_ptr);
    std::string set_name = atomtype_set_ptr->attributes["name"];

    if (atom_type_set.find(set_name) != atom_type_set.end())
      throw input_error("listed atom type name or atom type Set name " +
                        set_name + " more than once\n");

    atom_type_set[set_name] = AtomTypeSet(version, set_name, idx++);

    for (auto&& member : atomtype_set_ptr->children) {
      std::string member_id = member->attributes["id"];
      auto it = atom_type_set.find(member_id);

      if (it != atom_type_set.end()) {
        if (it->second.GetVersion() != version)
          std::cout << "?Warning: atom type (Set) " << member_id
                    << " has version " << it->second.GetVersion()
                    << " and is added to atom type Set " << set_name
                    << " which has version " << version << std::endl;

        for (auto&& at : it->second.GetAtomSet()) {
          atom_type_set[set_name].AddAtom(at);
          atom_type_set["any"].AddAtom(at);
        }

        for (auto&& att : it->second.GetAtomTypeNames()) {
          atom_type_set[set_name].AddAtomType(att);
          atom_type_set["any"].AddAtomType(att);
        }

      } else
        throw input_error(member_id +
                          " is neither an att nor an ats (or not defined "
                          "before the current ats " +
                          set_name + "\n");
    }
  }

  /*for(auto && at : atom_type_set) {
        std::cout << at.first << " " ;

        for(auto && att : at.second.GetAtomSet())
                std::cout  << " " << att << " ";

        std::cout << " types: ";

        for(auto&& att : at.second.GetAtomTypeNames())
                std::cout << " " << att << " ";

        std::cout << std::endl;
  }*/
}

void CheckAtomTypeXML(const XmlElement_ptr atomtype) {
  atomtype->CheckTagName("atomtype");
  atomtype->CheckAttributeSize(1);
  atomtype->CheckAttribute("id");
  atomtype->CheckNumberOfChildren_equal(1);
  atomtype->children.front()->CheckTagName("element_type");
  atomtype->children.front()->CheckValue();
  const std::string& id = atomtype->attributes.find("id")->second;

  if (!(id.front() >= 'A' && id.front() <= 'Z'))
    throw input_error(
        "Please start atom type ids with upper case letter in <atomtype> tag. "
        "Not the case for atomtype " +
        id + '\n');
}

void CheckAtomTypeSetXML(const XmlElement_ptr atomtype_set_ptr) {
  atomtype_set_ptr->CheckTagName("atomtypeset");
  atomtype_set_ptr->CheckAttributeSize(1);
  atomtype_set_ptr->CheckAttribute("name");
  atomtype_set_ptr->CheckNumberOfChildren_atLeast(1);
  const std::string& name = atomtype_set_ptr->attributes.find("name")->second;

  if (!(name.front() >= 'A' && name.front() <= 'Z'))
    throw input_error(
        "Please start atom type Set names with upper case letter in "
        "<atomtypeset> tag. Not the case for atomtype Set " +
        name + '\n');

  for (auto&& member : atomtype_set_ptr->children) {
    member->CheckTagName("member");
    member->CheckAttributeSize(1);
    member->CheckAttribute("id");
    member->CheckNoValue();
  }
}

std::ostream& operator<<(std::ostream& stream, const AtomTypeSet& type_set) {
  stream << type_set.GetName() << " contains: ";

  for (auto&& atom : type_set.GetAtomSet())
    stream << atom.GetUnitedAtomSymbol() << " ";

  stream << " atom types: ";

  for (auto&& tt : type_set.GetAtomTypeNames()) stream << tt << " ";

  stream << "idx: " << type_set.GetIndex();
  return stream;
}

}  // namespace topology_builder

}  // namespace combi_ff