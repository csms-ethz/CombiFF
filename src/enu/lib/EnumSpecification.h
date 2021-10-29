#ifndef ENUMSPECIFICATION_H_
#define ENUMSPECIFICATION_H_

#include "Atom.h"
#include "Range.h"
#include "StringVector.h"

namespace combi_ff{

namespace enu{

/*
EnumSpecifications:
- usedAtomVectors: vector that contains vectors with used atom types, e.g. for C3H4Cl1 and C4Br1F2, usedAtomVectors = [[C,H,Cl], [C,Br,F]]
- lambdaRangesVec: for each AtomVector in usedAtomVectors, lambdaRangesVec contains a vector of LambdaVectors that define the multiplicity of each of the atoms
- usedFamilies: stores the codes of the families that are to be created
- rangedProperties: contains min and max value for the following properties:
					 -rangedProperties[unsatRange] = number of unsaturations
                     -rangedProperties[bondRange] = number of bonds
                     -rangedProperties[SBrange] = number of single bonds (incl H-bonds)
                     -rangedProperties[DBrange] = number of double bonds
                     -rangedProperties[TBrange] = number of triple bonds
                     -rangedProperties[ringRange] = number of rings
- maxDeg: maximum bond degree betwen two atoms
- stereo: stereoisomers?
- explicitHydrogens: should SMILES contain explicit hydrogens?
*/

struct EnumSpecifications {
    std::vector<AtomVector<combi_ff::Atom>> used_atom_vectors {};
    std::vector<std::list<LambdaVector>> lambda_ranges_vec {};
    StringVector used_families {};
    size_t max_degree {4};
    bool stereo {false};
    bool explicit_hydrogens {false};
    RangeVector ranged_properties {num_ranges, Range(0,-1)};

    EnumSpecifications() = default;
};

} //namespace enu

} //namespace combi_ff

#endif