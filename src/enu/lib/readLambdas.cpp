// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "readLambdas.h"

#include <sstream>

#include "Matrix.h"
#include "exceptions.h"

namespace combi_ff {

namespace enu {

/********************************************************************************************
READ A SEQUENCE/RANGE OF LAMBDA VALUES FOR THE CURRENT ATOM TYPE AND ADD THEM TO
lambda_ranges
********************************************************************************************/
void ReadLambdas(std::list<LambdaVector>& lambda_ranges, size_t& j,
                 const std::string& formula) {
  // if we're at the end of the formula, or the next element is a letter or a
  // '{', assume lambda = 1 e.g. C{CH3} -> C1{CH3}1 and C1H3Cl -> C1H3Cl1
  if (j >= formula.size() || formula[j] == '{') {
    // the new lambda has to be added to all the LambdaVectors in lambda_ranges
    for (auto&& l : lambda_ranges) l.push_back(1);

    j++;
  }

  // if the atom type is followed by a number, use this number as lambda value
  else if (isdigit(formula[j])) {
    size_t n = GetNumber(j, formula);

    // the new lambda has to be added to all the LambdaVectors in lambda_ranges
    for (auto&& lr = lambda_ranges.begin(); lr != lambda_ranges.end(); ++lr)
      lr->push_back(n);
  }

  // if the atom type is followed by a '[', this means that a range or sequence
  // of lambdas is given
  else if (formula[j] == '[') {
    std::vector<size_t> nums;

    /****************************
    ADD THE FIRST NUMBER OR RANGE
    ****************************/
    if (!ReadLess(++j, formula, nums)) {
      if (!ReadGreater(j, formula, nums)) {
        if (!ReadNumber(j, formula, nums))
          throw combi_ff::input_error(
              "couldn't read lambda, unexpected character " +
              std::string(1, formula[j]) + " in " + formula);
      }
    }

    /********************************************************************
    CONTINUE ADDING NUMBERS AND RANGES UNTIL ENCOUNTERING THE CLOSING ']'
    ********************************************************************/
    while (j < formula.size() && formula[j] != ']') {
      // if current character is a comma, either a >, < or a digit has to follow
      if (formula[j] == ',') {
        if (!ReadLess(++j, formula, nums)) {
          if (!ReadGreater(j, formula, nums)) {
            if (!ReadNumber(j, formula, nums))
              throw combi_ff::input_error(
                  "couldn't read lambda, unexpected character " +
                  std::string(1, formula[j]) + " in " + formula);
          }
        }
      }

      // if current character is a -, a digit has to follow, and the range goes
      // from num.back() up to newly read number
      else if (formula[j] == '-') {
        size_t n = GetNumber(++j, formula);
        size_t numPrev = nums.back();

        if (n < numPrev) std::swap(n, numPrev);

        for (size_t i = numPrev; i <= n; i++) nums.push_back(i);
      }

      // unexpected character
      else
        throw combi_ff::input_error("expected a '-' or ',' in the formula " +
                                    formula + ", but encountered " +
                                    formula[j]);
    }

    // add all the unique numbers in num to the lambda_ranges
    AddNumsToLambdaRanges(nums, lambda_ranges, false);
    j++;

  } else if (formula[j] == '*') {
    std::vector<size_t> nums;
    nums.reserve(100);

    for (size_t i = 0; i < 100; i++) nums.push_back(i);

    AddNumsToLambdaRanges(nums, lambda_ranges, true);
    j++;

  } else
    throw combi_ff::input_error("unexpected character " +
                                std::string(1, formula[j]) + " in " + formula);
}

bool ReadLess(size_t& j, const std::string& formula,
              std::vector<size_t>& nums) {
  if (formula[j] == '<') {
    // check if <=
    bool leq(false);

    if (formula[++j] == '=') {
      leq = true;
      j++;
    }

    size_t n = GetNumber(j, formula);

    // add all numbers < n to nums
    for (size_t i = 0; i < n; i++) nums.push_back(i);

    // if <= add n to nums as well
    if (leq) nums.push_back(n);

    return true;

  } else
    return false;
}

bool ReadGreater(size_t& j, const std::string& formula,
                 std::vector<size_t>& nums) {
  if (formula[j] == '>') {
    // check if >=
    bool geq(false);

    if (formula[++j] == '=') {
      geq = true;
      j++;
    }

    size_t n = GetNumber(j, formula);

    // if >= add n to nums
    if (geq) nums.push_back(n);

    // add all numbers greater than n+1 to nums (and smaller than 100)
    for (size_t i = n + 1; i < 100; i++) nums.push_back(i);

    return true;

  } else
    return false;
}

bool ReadNumber(size_t& j, const std::string& formula,
                std::vector<size_t>& nums) {
  if (isdigit(formula[j])) {
    nums.push_back(GetNumber(j, formula));
    return true;

  } else
    return false;
}

void AddNumsToLambdaRanges(std::vector<size_t>& nums,
                           std::vector<LambdaVector>& lambda_ranges,
                           bool already_sorted) {
  if (!already_sorted) {
    std::sort(nums.begin(), nums.end());
    auto&& it = std::unique(nums.begin(), nums.end());
    nums.resize(std::distance(nums.begin(), it));
  }

  std::vector<LambdaVector> lambda_ranges_tmp(0);
  lambda_ranges_tmp.reserve(lambda_ranges.size() * nums.size());

  for (size_t i = 0; i < lambda_ranges.size(); i++) {
    for (size_t n : nums) {
      lambda_ranges_tmp.push_back(lambda_ranges[i]);
      lambda_ranges_tmp.back().push_back(n);
    }
  }

  lambda_ranges = lambda_ranges_tmp;
}

void AddNumsToLambdaRanges(std::vector<size_t>& nums,
                           std::list<LambdaVector>& lambda_ranges,
                           bool already_sorted) {
  if (!already_sorted) {
    std::sort(nums.begin(), nums.end());
    auto&& it = std::unique(nums.begin(), nums.end());
    nums.resize(std::distance(nums.begin(), it));
  }

  std::list<LambdaVector> lambda_ranges_tmp(0);

  for (auto&& lr = lambda_ranges.begin(); lr != lambda_ranges.end(); ++lr) {
    for (size_t n : nums) {
      lambda_ranges_tmp.push_back(*lr);
      lambda_ranges_tmp.back().push_back(n);
    }
  }

  lambda_ranges = lambda_ranges_tmp;
}

/***********
READ A RANGE
***********/
void ReadRange(const std::string& prop_, Range& r) {
  std::istringstream s(prop_);
  std::string prop(""), tmp;

  while (s >> tmp) prop += tmp;

  size_t j = 0;

  if (std::isdigit(prop[j])) {
    size_t num = GetNumber(j, prop);
    r = Range({num, num});

  } else if (prop[j] == '[') {
    size_t a = GetNumber(++j, prop);

    if (prop[j] == ']')
      r = Range(a, a);

    else {
      size_t b = GetNumber(++j, prop);
      r = Range(std::min(a, b), std::max(a, b));
    }

  } else if (prop[j] == '*') {
    r = Range({0, -1});

  } else
    throw combi_ff::input_error("unexpected character " +
                                std::string(1, prop[j]) + " in range " + prop);
}

size_t GetNumber(size_t& j, const std::string& formula) {
  std::string n = "";

  while (j < formula.size() && isdigit(formula[j])) n += formula[j++];

  return std::stoul(n);
}

}  // namespace enu

}  // namespace combi_ff