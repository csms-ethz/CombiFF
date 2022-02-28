// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef PSEUDOATOM_H
#define PSEUDOATOM_H

#include "Matrix.h"
#include "XmlParser.h"
namespace combi_ff {

namespace enu {

class Pseudoatom {
 public:
  Pseudoatom() = default;
  Pseudoatom(const FragmentMatrix& M, const std::string& code,
             const std::string& version, size_t degree);

  const FragmentMatrix& GetMatrix() const;
  const AtomVector<combi_ff::Atom>& GetAtoms() const;
  const std::string& GetCode() const;
  const std::string& GetVersion() const;
  size_t GetDegree() const;

 protected:
  FragmentMatrix M{0};
  const std::string code{"empty_psa"};
  const std::string version{""};
  size_t degree{0};
};

typedef std::unordered_map<std::string, Pseudoatom> PseudoatomMap;
typedef std::vector<Pseudoatom> PseudoatomVector;

void GetNextPseudoatom(PseudoatomMap& pseudoatoms,
                       const XmlElement_ptr pseudoatom,
                       const std::string& version);
void CreatePseudoatoms(PseudoatomMap& pseudoatoms,
                       const std::list<std::string>& pseudoatomSetFileNames);

}  // namespace enu

}  // namespace combi_ff

#endif
