#ifndef ELEMENTTYPES_H_
#define ELEMENTTYPES_H_

#include <iostream>
#include <unordered_map>

namespace combi_ff{

typedef size_t ElementPriority;
typedef size_t Connectivity;

typedef size_t ElementNumber;

typedef std::string ElementSymbol;
typedef std::string ElementIdentifier;

namespace element_types{

static const double hydrogen_mass(1.0079);

struct ElementProperties {
	ElementSymbol symbol {""};
	Connectivity degree {0};
	ElementNumber element_nr {0};
	ElementPriority priority {0};
	double mass {0};
	bool hydrogen_in_smiles {false};

	ElementProperties() = default;
	ElementProperties(ElementSymbol symbol, Connectivity degree, ElementNumber element_nr, ElementPriority priority, double mass, bool hydrogen_in_smiles) : 
		symbol(symbol), degree(degree), element_nr(element_nr),	priority(priority), mass(mass), hydrogen_in_smiles(hydrogen_in_smiles) {}
} ;


//IMPORTANT: element priority has to be unique for every element type!!
static const std::unordered_map  <ElementIdentifier, ElementProperties> element_property_map{
	{"C",  ElementProperties({"C", 4, 6, 0, 12.011, false})},
	{"N",  ElementProperties({"N", 3, 7, 3, 14.007, false})},
	{"P",  ElementProperties({"P", 3, 15, 5, 30.974, false})},
	{"Pf", ElementProperties({"P", 5, 15, 5, 30.974, true})},
	{"O",  ElementProperties({"O", 2, 8, 2, 15.999, false})},
	{"S",  ElementProperties({"S", 2, 16, 4, 32.060, false})},
	{"Sf", ElementProperties({"S", 4, 16, 4, 32.060, true})},
	{"Ss", ElementProperties({"S", 6, 16, 4, 32.060, true})},
	{"H",  ElementProperties({"H", 1, 1, 1, hydrogen_mass, false})},
	{"F",  ElementProperties({"F", 1, 9, 8, 18.998, false})},
	{"Cl", ElementProperties({"Cl", 1, 17, 6, 35.453, false})},
	{"Br", ElementProperties({"Br", 1, 35, 7, 79.904, false})},
	{"I",  ElementProperties({"I", 1, 53, 9, 126.900, false})},
	{"N+", ElementProperties({"N", 4, 7, 20, 14.007, true})},
	{"B",  ElementProperties({"B", 3, 5, 21, 10.81, false})},
	{"Se",  ElementProperties({"Se", 4, 34, 22, 78.96, true})},
	//atom types for sadra's dihedrals
	/*
	{"Me", elementProperties({1, 100, 100, 0})},
	{"OMX", elementProperties({1, 101, 101, 0})},
	{"KOX", elementProperties({1, 102, 102, 0})},
	{"OH", elementProperties({1, 103, 103, 0})},
	{"KdO", elementProperties({1, 104, 104, 0})},
	{"HX", elementProperties({1, 105, 105, 0})},
	*/
	//end of sadra's atom types
	{"*",  ElementProperties({"*",0, 0, 100, 0, false})}
};
} //namespace element_types

} //namespace combi_ff

#endif
