#ifndef READLAMBDAS_H
#define READLAMBDAS_H

#include "Range.h"
#include <string>
#include "LambdaVector.h"
#include <list>

namespace combi_ff {

namespace enu {

bool ReadLess(size_t& j, const std::string& formula, std::vector<size_t>& nums);
bool ReadGreater(size_t& j, const std::string& formula, std::vector<size_t>& nums);
bool ReadNumber(size_t& j, const std::string& formula, std::vector<size_t>& nums);
void AddNumsToLambdaRanges(std::vector<size_t>& nums,
                            std::vector<LambdaVector>& lambda_ranges,
                            bool already_sorted);
void AddNumsToLambdaRanges(std::vector<size_t>& nums,
                            std::list<LambdaVector>& lambda_ranges,
                            bool already_sorted);
void ReadLambdas(std::vector<LambdaVector>& lambda_ranges, size_t& j, const std::string& formula);
void ReadLambdas(std::list<LambdaVector>& lambda_ranges, size_t& j, const std::string& formula);
void ReadRange(const std::string& prop, Range& r);
size_t GetNumber(size_t& j, const std::string formula) ;

} //namespace enu

} //namespace combi_ff

#endif
