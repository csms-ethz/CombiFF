#ifndef SUBSTRUCTURE_H
#define SUBSTRUCTURE_H

#include "Matrix.h"
#include "Range.h"
#include "XmlParser.h"

namespace combi_ff {

namespace enu {

class AbstractSubstructure {
 public:
  AbstractSubstructure() = default;

  AbstractSubstructure(const size_t N,
                       const AdjacencyVector& v,
                       const std::string& code,
                       const std::string& version) ;

  AbstractSubstructure(const SymmetricalMatrix<size_t>& M,
                       const std::string& code,
                       const std::string& version) ;

  AbstractSubstructure(const AbstractSubstructure& FF);

  const SymmetricalMatrix<size_t>& GetMatrix() const;
  const std::string& GetCode() const;
  const std::string& GetVersion() const;

 protected:
  SymmetricalMatrix<size_t> M {0};
  const std::string code {""};
  const std::string version {""};
};

typedef std::unordered_map<std::string, AbstractSubstructure> AbstractSubstructureMap;


class SubstructureCollection {
 private:
  std::vector<FragmentMatrix> substructure_matrices {std::vector<FragmentMatrix>(0)};
  bool AND {false}, XOR {false};
  const std::string code {""};
  Range r {0, -1};
 public:
  SubstructureCollection() = default;
  SubstructureCollection(const AbstractSubstructure& F,
                         const Range& r,
                         const AtomVector<combi_ff::Atom>& atoms,
                         const bool AND,
                         const bool XOR) ;
  SubstructureCollection(const AbstractSubstructure& F,
                         const Range& r,
                         const FragmentMatrix& M,
                         const bool AND,
                         const bool XOR);
  SubstructureCollection(const AbstractSubstructure& F,
                         const Range& r,
                         std::vector<AtomVector<combi_ff::Atom>>& atoms,
                         const bool AND,
                         const bool XOR) ;


  const std::string& GetCode() const;
  const Range& GetRange() const;
  const std::vector<FragmentMatrix>& GetSubstructureMatrices() const;
  bool GetXOR() const;
  bool GetAND() const;
};

void CreateSubstructures(AbstractSubstructureMap& abstract_substructures,
                         const std::list<std::string>& substructure_file_names);


void GetNextAbstractSubstructure(AbstractSubstructureMap& abstract_substructures,
                                 const XmlElement_ptr substructure,
                                 const std::string& version) ;
} //namespace enu

} //namespace combi_ff

#endif
