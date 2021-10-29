#include "Range.h"
#include <assert.h>

namespace combi_ff{

bool IsInRange(size_t i, Range r) {
    if((int)i >= r.first && ((int)i <= r.second || r.second == -1))
        return true;

    else
        return false;
}

bool AreInRange(const std::vector<size_t>& ranges, const RangeVector& ranged_properties){
	assert(ranges.size() == ranged_properties.size());
	for(size_t i = 1; i < ranges.size(); i++){ //don't check unsaturations
		if(!IsInRange(ranges[i], ranged_properties[i]))
			return false;
	}
	return true;
}

} //namespace combi_ff
