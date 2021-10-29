#include "EnumSpecification.h"

namespace combi_ff{

namespace enu{

EnumSpecifications::EnumSpecifications() :
    used_atom_vectors(std::vector<AtomVector<combi_ff::Atom>>(0)), lambda_ranges_vec(std::vector<std::list<LambdaVector>>(0)), used_families(StringVector(0)),
    max_degree(4), stereo(false), ranged_properties(RangeVector(num_ranges, range(0, -1))) {}

} //namespace enu

} //namespace combi_ff