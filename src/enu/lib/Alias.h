#ifndef ALIAS_H
#define ALIAS_H

#include "Atom.h"
#include "StringVector.h"
#include "XmlParser.h"
#include <set>

namespace combi_ff{

namespace enu{


class AtomTypeAlias{
public:
	AtomTypeAlias() = default;
	AtomTypeAlias(const std::string& alias_name, const std::string& version);
	const std::string& GetVersion() const;
	const std::string& GetAliasName() const;
	const std::vector<combi_ff::ElementSymbol>& GetContainedAtomTypesConst() const;
	std::vector<combi_ff::ElementSymbol>& GetContainedAtomTypes();
	void SetVersion(const std::string& v);
	void SetAliasName(const std::string& a);
	void AddAtomType(const combi_ff::ElementSymbol& e);
private:
	std::string alias_name {""};
	std::vector<combi_ff::ElementSymbol> contained_atom_types {std::vector<combi_ff::ElementSymbol>(0)};
	std::string version {""};	
};

typedef std::unordered_map<std::string, AtomTypeAlias> AliasMap;

void GetNextAlias(AliasMap& alias_map, const XmlElement_ptr alias, const std::string& version);
void CreateAliases(AliasMap& alias_map, const std::list<std::string>& alias_file_names) ;
void SortElements(combi_ff::StringVector& elements);
bool CompareElems(std::string e1, std::string e2);


} //namespace enu

} //namespace combi_ff


#endif
 