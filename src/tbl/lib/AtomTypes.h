// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef ATOMTYPES_H_
#define ATOMTYPES_H_

#include <unordered_set>

#include "Atom.h"
#include "XmlParser.h"

namespace combi_ff {

namespace topology_builder {

class AtomHash;
class AtomTypeSet;

typedef std::unordered_map<std::string, Atom> AtomTypeMap;
typedef std::unordered_set<Atom, AtomHash> AtomSet;
typedef std::unordered_map<std::string, AtomTypeSet> AtomTypeSetMap;

class AtomType {
 public:
  AtomType(const std::string& version, const std::string& atom_type_name,
           const std::string& element_type_name);
  bool operator<(const AtomType& at) const;
  const Atom& GetElementType() const;
  const std::string& GetAtomTypeName() const;
  const std::string& GetVersion() const;

 private:
  const std::string version;
  std::string atom_type_name;
  Atom element_type;
};

class AtomHash {
 public:
  size_t operator()(const Atom& atom) const;
};

class AtomTypeSet {
 public:
  AtomTypeSet();
  AtomTypeSet(const std::string& version, const std::string& set_name,
              const size_t index);

  bool operator==(const AtomTypeSet& at) const;
  bool operator<(const AtomTypeSet& at) const;
  const AtomSet& GetAtomSet() const;
  const std::string& GetName() const;
  const std::unordered_set<std::string>& GetAtomTypeNames() const {
    return atom_type_names;
  }
  void AddAtom(const Atom& atom);
  void AddAtomType(const std::string& name);
  void SetName(const std::string& name);
  void SetVersion(const std::string& v);
  size_t GetIndex() const;
  const std::string& GetVersion() const;

 private:
  std::string version;
  std::string atom_type_set_name;
  AtomSet contained_atoms;
  std::unordered_set<std::string> atom_type_names;
  size_t index;
};

void CreateAtomTypeSets(const std::list<std::string>& atom_types_file_names,
                        AtomTypeSetMap& atom_type_set);
void CheckAtomTypeXML(const XmlElement_ptr atomtype);
void CheckAtomTypeSetXML(const XmlElement_ptr atomtype_set_ptr);
void AddSimpleAtomTypes(AtomTypeSetMap& atom_type_set,
                        const XmlElement_ptr simpleatomtypes,
                        const std::string& version, size_t& idx);
void AddAtomTypeSets(AtomTypeSetMap& atom_type_set,
                     const XmlElement_ptr atomtype_set_ptrs,
                     const std::string& version, size_t& idx);

std::ostream& operator<<(std::ostream& stream, const AtomTypeSet& type_set);

}  // namespace topology_builder

}  // namespace combi_ff

#endif
