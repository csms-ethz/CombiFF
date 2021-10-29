#include "TopologicalProperty.h"
#include "ContainerOperators.h"

namespace combi_ff {

namespace topology_builder {


/*************************************
TopologicalPropertyBase implementation
*************************************/

/*
Constructors
*/
TopologicalPropertyBase::TopologicalPropertyBase(const std::string& property_name,
                                                 const size_t num_involved_atoms)
  : property_name(property_name),
    num_involved_atoms(num_involved_atoms) {
  involved_atoms.reserve(num_involved_atoms);
}

TopologicalPropertyBase::TopologicalPropertyBase(const std::string& property_name,
                                                 const size_t num_involved_atoms,
                                                 const std::vector<size_t>& involved_atoms,
                                                 const std::unordered_map<std::string, std::string>& attributes)
  : property_name(property_name),
    num_involved_atoms(num_involved_atoms),
    involved_atoms(involved_atoms),
    attributes(attributes) {}

/*
Getter functions
*/

const std::string& TopologicalPropertyBase::GetPropertyName() const {
  return property_name;
}

const std::string& TopologicalPropertyBase::GetPropertyAbbreviation() const {
  return property_abbreviation;
}

size_t TopologicalPropertyBase::GetNumInvolvedAtoms() const {
  return num_involved_atoms;
}

const std::string TopologicalPropertyBase::GetAttribute(const std::string& name)const {
  if (attributes.find(name) != attributes.end())return attributes.find(name)->second;

  else return "";
}

const std::unordered_map<std::string, std::string>& TopologicalPropertyBase::GetAttributes() const {
  return attributes;
}

const std::vector<size_t>& TopologicalPropertyBase::GetInvolvedAtoms() const {
  return involved_atoms;
}

PropertyType TopologicalPropertyBase::GetType() const {
  return type;
}

/*
Setter functions
*/

void TopologicalPropertyBase::SetPropertyName(const std::string& s) {
  property_name = s;
}

void TopologicalPropertyBase::SetAttributes(const std::unordered_map<std::string, std::string>&
                                            attrib) {
  attributes = attrib;
}

void TopologicalPropertyBase::SetInvolvedAtoms(const std::vector<size_t>& vec) {
  involved_atoms = vec;
}

void TopologicalPropertyBase::SetNumInvolvedAtoms(const size_t i) {
  num_involved_atoms = i;
}

void TopologicalPropertyBase::AddAttribute(const std::string& name, const std::string& value) {
  attributes[name] = value;
}

void TopologicalPropertyBase::AddInvolvedAtom(const size_t idx) {
  involved_atoms.push_back(idx);
}

/*
operator<<
*/

std::ostream& operator<<(std::ostream& stream, const TopologicalPropertyBase& prop) {
  stream << prop.GetPropertyName() << " " <<
         prop.GetNumInvolvedAtoms();// << " " //<< prop.GetAttributes() << " ";
  combi_ff::operator<<(stream, prop.GetInvolvedAtoms());
  return stream;
}

/*********************************
TopologicalProperty implementation
*********************************/

/*
template constructors
*/

template<>
TopologicalProperty<bond>::TopologicalProperty() {
  type = bond;
  property_name = "bond";
  property_abbreviation = "bnd";
  num_involved_atoms = 2;
}

template<>
TopologicalProperty<angle>::TopologicalProperty() {
  type = angle;
  property_name = "angle";
  property_abbreviation = "ang";
  num_involved_atoms = 3;
}

template<>
TopologicalProperty<improper_dihedral>::TopologicalProperty() {
  type = improper_dihedral;
  property_name = "improper_dihedral";
  property_abbreviation = "imp";
  num_involved_atoms = 4;
}

template<>
TopologicalProperty<torsional_dihedral>::TopologicalProperty() {
  type = torsional_dihedral;
  property_name = "torsional_dihedral";
  property_abbreviation = "dih";
  num_involved_atoms = 4;
}

/*
template<>
TopologicalProperty<first_neighbor>::TopologicalProperty() {
  type = first_neighbor;
  property_name = "first_neighbor";
  property_abbreviation = "nbr1";
  num_involved_atoms = 2;
}

template<>
TopologicalProperty<second_neighbor>::TopologicalProperty() {
  type = second_neighbor;
  property_name = "second_neighbor";
  property_abbreviation = "nbr2";
  num_involved_atoms = 3;
}

template<>
TopologicalProperty<third_neighbor>::TopologicalProperty() {
  type = third_neighbor;
  property_name = "third_neighbor";
  property_abbreviation = "nbr3";
  num_involved_atoms = 4;
}
*/
/*
sortInvolvedAtoms implementation
*/

template<>
void TopologicalProperty<improper_dihedral>::SortInvolvedAtoms() {
  std::sort(involved_atoms.begin() + 1, involved_atoms.end());
}

template<>
void TopologicalProperty<angle>::SortInvolvedAtoms() {
  if (involved_atoms.front() > involved_atoms.back())
    std::swap(involved_atoms.front(), involved_atoms.back());
}

/*template<>
void TopologicalProperty<second_neighbor>::sortInvolvedAtoms() {
  if(involved_atoms.front() > involved_atoms.back())
    std::swap(involved_atoms.front(), involved_atoms.back());
}*/

template<>
void TopologicalProperty<torsional_dihedral>::SortInvolvedAtoms() {
  if (involved_atoms[1] > involved_atoms[2]) {
    std::swap(involved_atoms.front(), involved_atoms.back());
    std::swap(involved_atoms[1], involved_atoms[2]);
  }
}

/*template<>
void TopologicalProperty<third_neighbor>::sortInvolvedAtoms() {
  if(involved_atoms[1] > involved_atoms[2]) {
    std::swap(involved_atoms.front(), involved_atoms.back());
    std::swap(involved_atoms[1], involved_atoms[2]);
  }
}*/

template<PropertyType t>
void TopologicalProperty<t>::SortInvolvedAtoms() {
  std::sort(involved_atoms.begin(), involved_atoms.end());
}

/*
creates a new TopologicalPropertyBase pointer for the property_name type
*/
std::shared_ptr<TopologicalPropertyBase> CreateNewTopologicalProperty(const PropertyType type) {
  switch (type) {
  case (bond) :
    return std::shared_ptr<TopologicalPropertyBase>(new TopologicalProperty<bond>());
    break;

  case (angle) :
    return std::shared_ptr<TopologicalPropertyBase>(new TopologicalProperty<angle>());
    break;

  case (improper_dihedral) :
    return std::shared_ptr<TopologicalPropertyBase>(new TopologicalProperty<improper_dihedral>());
    break;

  case (torsional_dihedral) :
    return std::shared_ptr<TopologicalPropertyBase>(new TopologicalProperty<torsional_dihedral>());
    break;

  /*
  case(first_neighbor) :
    return std::shared_ptr<TopologicalPropertyBase>(new TopologicalProperty<first_neighbor>());
    break;

  case(second_neighbor) :
    return std::shared_ptr<TopologicalPropertyBase>(new TopologicalProperty<second_neighbor>());
    break;

  case(third_neighbor) :
    return std::shared_ptr<TopologicalPropertyBase>(new TopologicalProperty<third_neighbor>());
    break;
  */
  default :
    throw std::runtime_error("unknown case " + std::to_string(type) +
                             " for creating new TopologicalProperty\n");
  }
}

} //namespace topology_builder

} //namespace Com