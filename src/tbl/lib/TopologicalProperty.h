// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef TOPOLOGICALPROPERTY_H_
#define TOPOLOGICALPROPERTY_H_

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

namespace combi_ff {

namespace topology_builder {

// if new PropertyType is added, also adapt the vector of properties used to
// iterate in mtb.h
typedef enum {
  bond,
  angle,
  improper_dihedral,
  torsional_dihedral,
  num_topological_property_types
} PropertyType;
typedef std::string TopologicalPropertyTypeID;

class TopologicalPropertyBase {
 public:
  virtual ~TopologicalPropertyBase() {}
  TopologicalPropertyBase() = default;
  TopologicalPropertyBase(const std::string& property_name,
                          const size_t num_involved_atoms);
  TopologicalPropertyBase(
      const std::string& property_name, const size_t num_involved_atoms,
      const std::vector<size_t>& involved_atoms,
      const std::unordered_map<std::string, std::string>& attributes);
  const std::string& GetPropertyName() const;
  const std::string& GetPropertyAbbreviation() const;

  void SetPropertyName(const std::string& s);

  size_t GetNumInvolvedAtoms() const;

  void SetAttributes(
      const std::unordered_map<std::string, std::string>& attrib);

  void AddAttribute(const std::string& name, const std::string& value);

  const std::string GetAttribute(const std::string& name) const;

  void AddInvolvedAtom(const size_t idx);

  void SetInvolvedAtoms(const std::vector<size_t>& vec);

  void SetNumInvolvedAtoms(const size_t i);

  const std::unordered_map<std::string, std::string>& GetAttributes() const;

  const std::vector<size_t>& GetInvolvedAtoms() const;

  PropertyType GetType() const;

  virtual void SortInvolvedAtoms() = 0;

 protected:
  std::string property_name{""};
  std::string property_abbreviation{""};
  size_t num_involved_atoms{0};
  std::vector<size_t> involved_atoms{};
  std::unordered_map<std::string, std::string> attributes{};
  PropertyType type;
};

template <const PropertyType type_>
class TopologicalProperty final : public TopologicalPropertyBase {
 private:
 public:
  ~TopologicalProperty() {}
  TopologicalProperty();
  void SortInvolvedAtoms() override;
};

static const std::unordered_map<TopologicalPropertyTypeID, PropertyType>
    topologicalPropertyMap{
        {"bond", bond},
        {"angle", angle},
        {"improper_dihedral", improper_dihedral},
        {"torsional_dihedral", torsional_dihedral}
        //{"first_neighbor", first_neighbor},
        //{"second_neighbor", second_neighbor},
        //{"third_neighbor", third_neighbor}
    };

std::ostream& operator<<(std::ostream& stream,
                         const TopologicalPropertyBase& prop);
std::shared_ptr<TopologicalPropertyBase> CreateNewTopologicalProperty(
    const PropertyType type);

}  // namespace topology_builder

}  // namespace combi_ff

#endif