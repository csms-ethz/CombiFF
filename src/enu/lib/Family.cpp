// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "Family.h"

#include <sstream>

#include "ContainerOperators.h"
#include "EnumSpecification.h"
#include "exceptions.h"
#include "printInfo.h"
#include "readLambdas.h"
#include "version.h"

namespace combi_ff {

namespace enu {

/*******************************
IMPLEMENTATION OF Family METHODS
*******************************/

/***********
CONSTRUCTORS
***********/
Family::Family(const std::string& code, const std::string& version,
               const std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
               const std::vector<std::list<LambdaVector>>& lambda_ranges_vec,
               const size_t max_degree, const RangeVector& ranged_properties,
               const std::vector<SubstructureCollection>& substructures,
               const bool has_pseudoatom)
    : code(code),
      version(version),
      used_atom_vectors(used_atom_vectors),
      lambda_ranges_vec(lambda_ranges_vec),
      max_degree(max_degree),
      ranged_properties(ranged_properties),
      substructures(substructures),
      has_pseudoatom(has_pseudoatom) {}

/*************
GetTER METHODS
*************/
const std::string& Family::GetCode() const { return code; }

const std::string& Family::GetVersion() const { return version; }

const std::vector<AtomVector<combi_ff::Atom>>& Family::GetUsedAtomVectors()
    const {
  return used_atom_vectors;
}

const std::vector<std::list<LambdaVector>>& Family::GetLambdaRangesVec() const {
  return lambda_ranges_vec;
}

size_t Family::GetMaxDeg() const { return max_degree; }

const RangeVector& Family::GetRangedProperties() const {
  return ranged_properties;
}

const std::vector<SubstructureCollection>& Family::GetSubstructures() const {
  return substructures;
}

bool Family::GetHasPseudoatom() const { return has_pseudoatom; }

/***********
PRINT METHOD
***********/

void Family::Print() const {
  const std::string indent1(6, ' ');
  const std::string indent2(indent1 + std::string(2, ' '));
  const std::string indent3(indent2 + std::string(2, ' '));
  std::cout << code << '\n';
  std::cout << indent1 << "->used_atom_vectors:\n";

  for (auto&& at : used_atom_vectors) std::cout << indent2 << at << '\n';

  std::cout << indent1 << "->Restrictions:\n";
  PrintRestrictions(indent2, max_degree, ranged_properties);

  if (substructures.size()) {
    std::cout << indent1 << "->used Substructures:\n";
    PrintSubstructures(indent2, indent3, substructures);
  }
}

/***************************************************
IMPLEMENTATION OF FUNCTIONS THAT HANDLE Family CLASS
***************************************************/

/****************************************************************************************************
CREATE THE VECTOR OF families FROM THE familyDefinitionFile. TO THIS END, FIRST
CREATE TEMPLATE FAMILIES FROM THE DEFINITONS IN THE familyDefinitionFile, AND
THEN GO THROUGH used_families TO ADD THE REQUIRED FAMILIES
****************************************************************************************************/

void CreateFamilies(FamilyVector& families,
                    std::list<std::string>& family_definition_file_names,
                    const StringVector& used_families,
                    const AbstractSubstructureMap& abstr_substructures,
                    const AliasMap& alias_map,
                    const PseudoatomMap& pseudoatoms) {
  if (!family_definition_file_names.size()) {
    if (used_families.size()) {
      std::cout << used_families << '\n';
      throw std::runtime_error(
          "specified the above families to enumerate, but no family Set file "
          "is open\n");

    } else
      throw combi_ff::input_warning("no family Set file\n");
  }

  // only create families if they are used
  if (used_families.size()) {
    std::cout << "creating families\n";

    for (auto&& family_definition_file_name : family_definition_file_names) {
      XmlParserIn parser(family_definition_file_name, XmlParserIn::read_all);
      const XmlElement& families_root = parser.GetTree().GetRoot();
      families_root.CheckTagName("family_definitions");
      families_root.CheckNumberOfChildren_atLeast(1);
      families_root.CheckAttribute("version");
      const std::string& version =
          families_root.attributes.find("version")->second;

      if (version != combi_ff::current_version)
        std::cout << "?Warning: currently running combi_ff version "
                  << combi_ff::current_version << " but family file "
                  << family_definition_file_name << " is version " << version
                  << "\n";

      bool all_fam(false);

      if (std::find(used_families.begin(), used_families.end(), "all") !=
          used_families.end()) {
        // used_families.erase(ff);
        all_fam = true;
      }

      // std::istringstream line(read(familyDefinitionFile));
      // std::string blockName;//("NEWFAMILY");
      // line >> blockName;
      // control vector for code and acronym to check that each code and each
      // acronym occurs only once
      StringVector codes(0), acros(0);

      // while(!familyDefinitionFile.eof())
      for (auto&& family : families_root.children)
        GetNextFamily(version, families, family, used_families, alias_map,
                      pseudoatoms, abstr_substructures, all_fam, codes, acros);
    }
  }

  for (auto&& f : used_families) {
    bool found(false);

    if (f != "all") {
      for (auto&& fam : families) {
        if (fam.GetCode() == f) {
          found = true;
          break;
        }
      }

      if (!found)
        throw combi_ff::input_error("specified family '" + f +
                                    "' in input, but couldn't find this family "
                                    "in the family library file(s)");
    }
  }

  PrintFamilies(families);
}

/******************************************************************************************************************
CREATE A NEW FAMILY BY READING A FAMILY BLOCK IN familyDefinitionFile BETWEEN
THE NEWFAMILY TO THE ENDFAMILY KEYWORDS, BUT IF A FAMILY ACRO/CODE IS NOT
ENCOUNTERED IN THE used_families ACRO/CODE, THE FAMILY CAN BE SKIPPED
******************************************************************************************************************/
void GetNextFamily(const std::string& version, FamilyVector& families,
                   const XmlElement_ptr family,
                   const StringVector& used_families, const AliasMap& alias_map,
                   const PseudoatomMap& pseudoatoms,
                   const AbstractSubstructureMap& abstr_substructures,
                   const bool all_fam, StringVector& codes,
                   StringVector& acros) {
  // declare necessary member variables for a Family
  std::vector<std::list<LambdaVector>> lambda_ranges_vec(0);
  std::vector<LambdaVector> lambda_ranges(0);
  std::vector<AtomVector<combi_ff::Atom>> used_atom_vectors(
      1, AtomVector<combi_ff::Atom>(0));
  RangeVector ranged_properties(num_ranges, {0, -1});
  std::vector<SubstructureCollection> substructures(0);
  PseudoatomMap pseudoatoms_local(0);
  size_t max_degree(5);
  bool has_pseudoatom(false);
  family->CheckTagName("family_definition");
  family->CheckAttribute("code");
  family->CheckAttributeSize(1);
  family->CheckNumberOfChildren_atLeast(1);
  std::string code = family->attributes.find("code")->second;

  if (std::find(codes.begin(), codes.end(), code) != codes.end())
    throw combi_ff::input_error(
        "family with code " + code +
        " occurs more than once in the family Set file");

  if (!all_fam && std::find(used_families.begin(), used_families.end(), code) ==
                      used_families.end()) {
    // skip this family if it's not in the vector of used families
    return;
  }

  codes.push_back(code);
  auto&& family_property_it = family->children.begin();
  const XmlElement_ptr formula = *family_property_it;
  formula->CheckTagName("formula");
  formula->CheckNumberOfChildren_equal(0);
  formula->CheckAttributeSize(0);
  formula->CheckValue();
  GetFormula(formula->value, alias_map, pseudoatoms, used_atom_vectors,
             lambda_ranges_vec, has_pseudoatom, version);
  // while(blockName != "ENDFAMILY") {
  family_property_it++;

  for (; family_property_it != family->children.end(); ++family_property_it) {
    const XmlElement_ptr family_property = *family_property_it;
    const std::string& tag_name = family_property->tag;

    // read MAXDEG block
    if (tag_name == "max_bond_degree") {
      family_property->CheckAttributeSize(0);
      family_property->CheckValue();

      for (auto&& c : family_property->value) {
        if (!std::isdigit(c))
          throw input_error(
              "maximum_bond_degree has to be an integer value, but " +
              std::string(1, c) + " in " + family_property->value +
              "is not a digit\n");
      }

      max_degree = std::stoi(family_property->value);

    } else if (possible_ranges.find(tag_name) != possible_ranges.end())
      GetRestriction(family_property->value,
                     ranged_properties[possible_ranges.find(tag_name)->second]);

    // read SUBSTRUC block
    else if (tag_name == "substructures")
      GetSubstructures(family_property, abstr_substructures, substructures,
                       alias_map, version);

    // unknown input
    else
      throw combi_ff::input_error(
          "in family file, unrecognized option " + tag_name +
          ". Please have a look at the manual for valid input options");
  }

  if (all_fam || std::find(used_families.begin(), used_families.end(), code) !=
                     used_families.end()) {
    std::cout << "    -> added family " << code << "\n";
    families.push_back(Family(code, version, used_atom_vectors,
                              lambda_ranges_vec, max_degree, ranged_properties,
                              substructures, has_pseudoatom));
  }
}

/******************************************************************************************
READ THE FORMULA BLOCK OF THE familyDefinitionFile AND ADD THE CORRESPONDING
ATOMS AND LAMBDAS
 -> note that unlike in the simple processIO.cpp input, the atom types can also
be alias_map here. Additionally, different alias atom types can be defined as
XOR to each other
    -> note that if you want to use AND, just use MAPn, to use n of the same
atom types and if you want to use OR, just use MAPn Mapm ..., to use n of the
same and m of the same or a different type
 -> each AtomVector has a vector of LambdaVectors that it corresponds to
 -> every time, a alias is encountered, the number of used_atom_vectors and
lambda_ranges_vec is increased by a factor corresponding to the size of the
alias, since they're all stored explicitly
 -> xor_idx: has for every atom type. if an atom type doesn't have an XOR
indicator ^, the entry is -1. if two atom types have the same entry, they *have*
to be different
******************************************************************************************/
void GetFormula(const std::string& formula_, const AliasMap& alias_map,
                const PseudoatomMap& pseudoatoms,
                std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
                std::vector<std::list<LambdaVector>>& lambda_ranges_vec,
                bool& has_pseudoatom, const std::string& version) {
  if (lambda_ranges_vec.size())
    throw combi_ff::input_error(
        "please only define one FORMULA entry per family");

  std::list<LambdaVector> lambda_ranges(0);
  lambda_ranges.push_back(LambdaVector(0));
  std::vector<int> xor_idx(0);
  bool has_xor(false);
  size_t j = 0;
  // remove spaces and line breaks
  std::istringstream s(formula_);
  std::string formula(""), tmp;

  while (s >> tmp) formula += tmp;

  // go through FORMULA line, atom type by atom type
  while (j < formula.size()) {
    // case that the current atom type has an XOR
    if (formula[j] == '^') {
      has_xor = true;

      // if there's a number after ^, this is the xor_idx to be used
      if (j + 1 < formula.size() && isdigit(formula[++j]))
        xor_idx.push_back((int)GetNumber(j, formula));

      // unnumbered XORs are thrown toGether under the index 0
      else
        xor_idx.push_back(0);

      // add the atom type and lambda value
      AddAtom(formula, j, alias_map, used_atom_vectors, lambda_ranges, version);
    }

    // case that the current atom type is a united atom
    else if (formula[j] == '{') {
      // not an XOR
      xor_idx.push_back(-1);
      // add the atom type and lambda value
      AddUnitedAtom(formula, j, used_atom_vectors, lambda_ranges);

    } else if (formula[j] == '\'') {
      // not an XOR
      xor_idx.push_back(-1);
      has_pseudoatom = true;
      // add the atom type and lambda value
      AddPseudoatom(formula, j, used_atom_vectors, lambda_ranges, pseudoatoms,
                    version);

    } else if (isalpha(formula[j])) {
      // not an XOR
      xor_idx.push_back(-1);
      // add the atom type and lambda value
      AddAtom(formula, j, alias_map, used_atom_vectors, lambda_ranges, version);

    } else
      throw combi_ff::input_error(
          "format error for \"FORMULA\" in Family: unexpected character " +
          std::string(1, formula[j]) + " in " + formula);
  }

  CleanUpAtomVectorsFam(has_xor, xor_idx, lambda_ranges, used_atom_vectors,
                        lambda_ranges_vec);
}

/**************************************************************************************************************
READ THE RESTRICTION AFTER THE NUNSAT, NBONDS, NSB, NDB, NTB, OR NRINGS KEYWORDS
AND CALL readRange TO ADD THEM
**************************************************************************************************************/
void GetRestriction(const std::string& range_value, Range& r) {
  ReadRange(range_value, r);
}

void GetSubstructures(const XmlElement_ptr substructures_xml,
                      const AbstractSubstructureMap& abstr_substructures,
                      std::vector<SubstructureCollection>& substructures,
                      const AliasMap& alias_map, const std::string& version) {
  substructures_xml->CheckAttributeSize(0);
  substructures_xml->CheckNumberOfChildren_atLeast(1);

  for (auto&& substructure : substructures_xml->children) {
    substructure->CheckTagName("substructure");
    substructure->CheckNumberOfChildren_equal(1);
    substructure->CheckAttribute("substructure_code");
    substructure->CheckAttribute("amount");
    // create an pointer to the AbstractSubstructure with the corresponding code
    // AbstractSubstructure const* prototype = NULL;
    const std::string& code =
        substructure->attributes.find("substructure_code")->second;
    auto prototype_substr = abstr_substructures.find(code);

    if (prototype_substr == abstr_substructures.end())
      throw combi_ff::input_error("didn't find Substructure " + code);

    if (prototype_substr->second.GetVersion() != version)
      std::cout << "?Warning: using substructure " << code << " with version "
                << prototype_substr->second.GetVersion()
                << " but current family is version " << version << "\n";

    const std::string& range_value =
        substructure->attributes.find("amount")->second;
    Range r({0, -1});
    ReadRange(range_value, r);
    auto&& atoms = substructure->GetFirstChild();
    atoms->CheckTagName("atoms");
    atoms->CheckAttributeSize(0);
    atoms->CheckNumberOfChildren_atLeast(1);
    std::vector<int> xor_idx(0);
    std::vector<int> and_idx(0);
    std::string atom_types_in_substructure("");

    for (auto&& atom : atoms->children) {
      atom->CheckTagName("atom");
      atom->CheckAttributeSize(1);
      atom->CheckNoValue();
      atom->CheckAttribute("type");
      atom_types_in_substructure += atom->attributes.find("type")->second + " ";
    }

    std::vector<AtomVector<combi_ff::Atom>> atom_vectors(1);
    size_t j = 0;
    // read the atom names for the current fragment
    AddAtomTypes(atom_types_in_substructure, j, xor_idx, and_idx, atom_vectors,
                 alias_map, version);

    for (auto&& a : atom_vectors) {
      if (a.size() != prototype_substr->second.GetMatrix().GetN())
        throw combi_ff::input_error("atom vector not given correctly for " +
                                    code);
    }

    // check if fragment is supposed to be an AND or XOR
    bool AND(false), XOR(false);

    if (substructure->attributes.find("repetition") !=
        substructure->attributes.end()) {
      const std::string& repetitionType =
          substructure->attributes.find("repetition")->second;

      if (repetitionType == "and")
        AND = true;

      else if (repetitionType == "xor")
        XOR = true;

      else if (repetitionType != "or")
        throw input_error("unknown repetition type " + repetitionType +
                          " for substructure " + code +
                          " (known types are and, xor, or)\n");
    }

    // remove the atom vectors that don't fulfill the AND/XOR restrictions from
    // atom_vectors
    CleanUpAtomVectorsFrag(atom_vectors, xor_idx, and_idx);
    // Add the Fragment whose input parameters were just read to the fragments
    // vector
    substructures.push_back(SubstructureCollection(prototype_substr->second, r,
                                                   atom_vectors, AND, XOR));
  }
}

/*******************************************************
READ AN ATOM NAME AND ADD THE ATOM AND ITS LAMBDA VALUES
*******************************************************/
void AddAtom(const std::string& formula, size_t& j, const AliasMap& alias_map,
             std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
             std::list<LambdaVector>& lambda_ranges,
             const std::string& version) {
  // read the atom name or alias code into nam
  std::string ele("");

  while (j < formula.size() && isalpha(formula[j])) ele += formula[j++];

  // read the lambdas after the type name
  ReadLambdas(lambda_ranges, j, formula);
  bool isAlias(false);
  // check if the current atom type is an alias code
  // auto it = std::find_if( alias_map.begin(), alias_map.end(), [](const alias&
  // element){ return alias.first == ele;} );
  auto alias = alias_map.find(ele);

  if (alias != alias_map.end()) {
    isAlias = true;

    if (alias->second.GetVersion() != version)
      std::cout << "?Warning: using atom type alias " << ele << " with version "
                << alias->second.GetVersion()
                << " but current family is version " << version << "\n";

    std::vector<AtomVector<combi_ff::Atom>> atom_vectors_tmp(0);

    // each atom vector is added as many times as there are elements in the
    // alias, one copy for each of the element in the alias, where the last
    // element is an alias
    for (size_t ii = 0; ii < used_atom_vectors.size(); ii++) {
      for (auto&& alias_atom : alias->second.GetContainedAtomTypesConst()) {
        // add a copy of the previously used used_atom_vectors[ii], and append
        // it with alias_atom
        atom_vectors_tmp.push_back(used_atom_vectors[ii]);

        // case that the alias is a united atom
        if (alias_atom[0] == '{') {
          std::string ele("");
          size_t jj = 1;

          while (jj < alias_atom.size() && alias_atom[jj] != '}' &&
                 isalpha(alias_atom[jj]))
            ele += alias_atom[jj++];

          // add an Atom of type ele (minus the last 'H' char) to the currently
          // last vector in tempAtomsVecs
          atom_vectors_tmp.back().push_back(ele.substr(0, ele.size() - 1));
          // Get number of H atoms in united atom. Default is 1
          size_t n(1);

          if (isdigit(alias_atom[jj])) n = GetNumber(jj, alias_atom);

          // Set the fixed number of H atoms in the Atom that was just added
          atom_vectors_tmp.back().back().SetNumFixedHydrogens(n);
        }

        // if not a unitd atom, directly add the atom to the currently last
        // vector in tempAtomsVecs
        else
          atom_vectors_tmp.back().push_back(alias_atom);
      }
    }

    used_atom_vectors = atom_vectors_tmp;
  }

  // if the current element is not an alias, simply add it to all of the vectors
  // in used_atom_vectors
  if (!isAlias) {
    for (auto&& uA : used_atom_vectors) uA.push_back(ele);
  }
}

/*************************************************************
READ A UNITED ATOM NAME AND ADD THE ATOM AND ITS LAMBDA VALUES
*************************************************************/
void AddUnitedAtom(const std::string& formula, size_t& j,
                   std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
                   std::list<LambdaVector>& lambda_ranges) {
  std::string nam("");

  while (j + 1 < formula.size() && formula[++j] != '}' && isalpha(formula[j]))
    nam += formula[j];

  if (nam.back() != 'H')
    throw combi_ff::input_error(
        "united atom notation only works for hydrogen atom, e.g. {CH3}");

  // shorten nam, s.t. the 'H' is not part of it anymore, but only the atom type
  // name of the first atom
  nam = nam.substr(0, nam.size() - 1);

  // add the atom type
  for (auto&& uA : used_atom_vectors) uA.push_back(nam);

  // determine the number of H atoms in the united atom. Default is 1
  size_t n(1);

  if (isdigit(formula[j])) n = GetNumber(j, formula);

  // Set the fixed number of H atoms of the current atom type
  for (auto&& uA : used_atom_vectors) uA.back().SetNumFixedHydrogens(n);

  if (formula[j++] != '}')
    throw combi_ff::input_error("expected \'}\' in formula " + formula +
                                " but encountered " +
                                std::string(1, formula[j - 1]));

  // read the lambda values of the current atom type
  ReadLambdas(lambda_ranges, j, formula);
}

/*************************************************************
READ A PSEUDOATOM NAME AND ADD THE ATOM AND ITS LAMBDA VALUES
*************************************************************/
void AddPseudoatom(const std::string& formula, size_t& j,
                   std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
                   std::list<LambdaVector>& lambda_ranges,
                   const PseudoatomMap& pseudoatoms,
                   const std::string& version) {
  std::string nam("\'");

  while (j + 1 < formula.size() &&
         formula[++j] != '\'' /* && isalpha(formula[j])*/)
    nam += formula[j];

  nam += '\'';
  auto psa = pseudoatoms.find(nam);

  if (psa == pseudoatoms.end())
    throw combi_ff::input_error("couldn't find Pseudoatom " + nam +
                                " in pseudoatoms");

  if (psa->second.GetVersion() != version)
    std::cout << "?Warning: using pseudoatom " << nam << " with version "
              << psa->second.GetVersion() << " but current family is version "
              << version << "\n";

  // add the atom type
  for (auto&& uA : used_atom_vectors) {
    uA.push_back(Atom(nam, psa->second.GetDegree()));
    uA.back().SetPseudoatom(true);
  }

  // read the lambda values of the current atom type
  ReadLambdas(lambda_ranges, ++j, formula);
}

/************************************************************************************************************
REMOVE ANY DUPLICATES (both in AtomVector and LambdaVectors) AND VIOLATIONS OF
XOR CONDITIONS FROM THE used_atom_vectors
************************************************************************************************************/

void CleanUpAtomVectorsFam(
    bool has_xor, std::vector<int>& xor_idx,
    std::list<LambdaVector>& lambda_ranges,
    std::vector<AtomVector<combi_ff::Atom>>& used_atom_vectors,
    std::vector<std::list<LambdaVector>>& lambda_ranges_vec) {
  // create used_atom_vectors.size() copies of xor_idx in a vector xor_idx_vec,
  // s.t. each vector in used_atom_vectors has its own xor_idx vector. This
  // assures that when the vectors are sorted in the next step, the entries of
  // the xor_idx  vector still correspond to the correct atom types for all of
  // the vectors
  std::vector<std::vector<int>> xor_idx_vec(used_atom_vectors.size(), xor_idx);
  // same goes for lambda_ranges, except here they have to be added to the
  // currently existing lambda_ranges_vec;
  std::vector<std::list<LambdaVector>> lambda_ranges_tmp(
      std::vector<std::list<LambdaVector>>(used_atom_vectors.size(),
                                           lambda_ranges));
  lambda_ranges_vec.insert(lambda_ranges_vec.end(), lambda_ranges_tmp.begin(),
                           lambda_ranges_tmp.end());

  // sort the vectors in used_atom_vectors, toGether with their corresponding
  // vector of LambdaVectors and xor_idx vec
  for (size_t ii = 0; ii < used_atom_vectors.size(); ii++)
    SortAtoms(used_atom_vectors[ii], lambda_ranges_vec[ii], xor_idx_vec[ii]);

  // if there is an XOR in the current FORMULA, erase every entry that doesn't
  // fit the restriction, i.e. every entry, in which two atoms with the same
  // xor_idx value are identical
  if (has_xor) {
    // go through all the atom_vectors in used_atom_vectors
    for (int i = 0; i < (int)used_atom_vectors.size(); i++) {
      bool brk(false);

      // go through all the atom pairs (used_atom_vectors[i][ii],
      // used_atom_vectors[i][jj])
      for (size_t ii = 0; ii < used_atom_vectors.size(); ii++) {
        for (size_t jj = ii + 1; jj < used_atom_vectors[i].size(); jj++) {
          // if two atoms in the current AtomVector are identical, and their
          // xor_idx entries are identical and not equal to -1, the current
          // AtomVector violates the XOR restrictions and needs ot be erased
          // from used_atom_vectors. The corresponding entries also need to be
          // erased from lambda_ranges_vec and xor_idx_vec
          if (used_atom_vectors[i][ii] == used_atom_vectors[i][jj] &&
              xor_idx_vec[i][ii] != -1 &&
              xor_idx_vec[i][ii] == xor_idx_vec[i][jj]) {
            used_atom_vectors.erase(used_atom_vectors.begin() + i);
            lambda_ranges_vec.erase(lambda_ranges_vec.begin() + i);
            xor_idx_vec.erase(xor_idx_vec.begin() + i);
            i--;
            brk = true;
            break;
          }
        }

        // don't keep looping over the current atom pairs, if the corresponding
        // AtomVector is already erased
        if (brk) break;
      }
    }
  }

  // go through all the atom_vectors in used_atom_vectors, and combine the
  // entries that correspond to an identical atom type
  for (size_t i = 0; i < used_atom_vectors.size(); i++) {
    // go through all the pairs (used_atom_vectors[i][ii],
    // used_atom_vectors[i][jj])
    for (size_t ii = 0; ii < used_atom_vectors[i].size(); ii++) {
      for (int jj = (int)ii + 1; jj < (int)used_atom_vectors[i].size(); jj++) {
        // if the atom entry ii is identical to the atom entry jj, combine the
        // two by adding their lambda values and erasing the jj entries in the
        // lambda_ranges_vec and the used_atom_vectors
        if (used_atom_vectors[i][ii] == used_atom_vectors[i][jj]) {
          // for(size_t k = 0; k < lambda_ranges_vec[i].size(); k++) {
          for (auto&& lr : lambda_ranges_vec[i]) {
            // lambda_ranges_vec[i][k][ii] += lambda_ranges_vec[i][k][jj];
            lr[ii] += lr[jj];
            // lambda_ranges_vec[i][k].erase(lambda_ranges_vec[i][k].begin() +
            // jj);
            lr.erase(lr.begin() + jj);
          }

          used_atom_vectors[i].erase(used_atom_vectors[i].begin() + jj);
          jj--;
        }
      }
    }
  }

  // erase the atoms in the atom_vectors that have a lambda value of 0. If the
  // thereby created AtomVector is already in used_atom_vectors, simply erase it
  // and add the additional lambda values, and otherwise, just add the newly
  // created AtomVector with its lambdaValues
  size_t s = used_atom_vectors.size();

  for (size_t i = 0; i < s; i++) {
    bool erase(false);

    for (auto&& lr = lambda_ranges_vec[i].begin();
         lambda_ranges_vec[i].size() && lr != lambda_ranges_vec[i].end();
         lr++) {
      AtomVector<combi_ff::Atom> curr_atom_vec = used_atom_vectors[i];

      if (erase) lr--;

      bool zeroEntry(false);

      // search for zero entries in all entries of lambda_ranges_vec[i][jj] and
      // erase the corresponding atom types
      for (int k = 0; k < (int)lr->size(); k++) {
        if ((*lr)[k] == 0) {
          zeroEntry = true;
          curr_atom_vec.erase(curr_atom_vec.begin() + k);
          lr->erase(lr->begin() + k);
          k--;
        }
      }

      erase = (false);

      // if current AtomVector was changed, check if it has to be added to
      // used_atom_vectors, or if it's already in there
      if (zeroEntry) {
        auto lr2 = lr;
        LambdaVector lambda_vec_tmp = *lr;
        lr++;
        lambda_ranges_vec[i].erase(lr2);
        erase = true;
        // lr--;
        bool found(false);

        // std::cout << &lambda_vec_tmp << '\n';
        // check if curr_atom_vec is found in used_atom_vectors and add the
        // newly created LambdaVector
        for (size_t iii = 0; iii < used_atom_vectors.size(); iii++) {
          if (used_atom_vectors[iii] == curr_atom_vec) {
            found = true;
            // check if the corresponding LambdaVector that was just changed is
            // already in the AtomVector's lamdaRanges if not, add it to
            // lambda_ranges_vec[iii]
            lambda_ranges_vec[iii].push_back(lambda_vec_tmp);
            break;
          }
        }

        // if curr_atom_vec is not in used_atom_vectors, add it toGether with
        // its lambda vector
        if (!found) {
          used_atom_vectors.push_back(curr_atom_vec);
          std::list<LambdaVector> ll(1, lambda_vec_tmp);
          lambda_ranges_vec.push_back(ll);
        }
      }
    }
  }

  // remove empty atom_vectors
  for (size_t i = 0; i < used_atom_vectors.size(); i++) {
    if (!used_atom_vectors[i].size()) {
      used_atom_vectors.erase(used_atom_vectors.begin() + i);
      lambda_ranges_vec.erase(lambda_ranges_vec.begin() + i);
      i--;
    }
  }

  // erase all atom_vectors that appear more than once in used_atom_vectors. Add
  // the lambda_ranges of the deleted AtomVector to the ones of the one that's
  // not deleted
  int resize = 0;

  for (int ii = 0; ii < ((int)used_atom_vectors.size() - resize); ii++) {
    for (int jj = ii + 1; jj < ((int)used_atom_vectors.size() - resize); jj++) {
      if (used_atom_vectors[ii] == used_atom_vectors[jj]) {
        lambda_ranges_vec[ii].insert(lambda_ranges_vec[ii].end(),
                                     lambda_ranges_vec[jj].begin(),
                                     lambda_ranges_vec[jj].end());
        std::swap(*(used_atom_vectors.begin() + jj),
                  *(used_atom_vectors.end() - 1 - resize));
        std::swap(*(lambda_ranges_vec.begin() + jj),
                  *(lambda_ranges_vec.end() - 1 - resize));
        resize++;
        jj--;
      }
    }
  }

  for (size_t i = 0; i < used_atom_vectors.size(); i++) {
    for (auto&& lr : lambda_ranges_vec[i]) {
      if (lr.size() != used_atom_vectors[i].size())
        throw std::runtime_error(
            "size of lr is not equal to size of used_atom_vectors");
    }
  }

  used_atom_vectors.resize(used_atom_vectors.size() - resize);
  lambda_ranges_vec.resize(lambda_ranges_vec.size() - resize);

  // make sure that all entries in the lambda_ranges are unique for all the
  // used_atom_vectors
  for (size_t i = 0; i < used_atom_vectors.size(); i++) {
    lambda_ranges_vec[i].sort();
    auto&& it =
        std::unique(lambda_ranges_vec[i].begin(), lambda_ranges_vec[i].end());
    lambda_ranges_vec[i].resize(
        std::distance(lambda_ranges_vec[i].begin(), it));
  }

  for (size_t ii = 0; ii < used_atom_vectors.size(); ii++) {
    for (size_t jj = ii + 1; jj < used_atom_vectors.size(); jj++) {
      if (used_atom_vectors[ii].size() > used_atom_vectors[jj].size()) {
        /*for(auto && lambda : lambda_ranges_vec[ii])
            assert(lambda.size() == used_atom_vectors[ii].size());

        for(auto && lambda : lambda_ranges_vec[jj])
            assert(lambda.size() == used_atom_vectors[jj].size());*/
        used_atom_vectors[ii].swap(used_atom_vectors[jj]);
        std::swap(lambda_ranges_vec[ii], lambda_ranges_vec[jj]);
        /*for(auto && lambda : lambda_ranges_vec[ii])
            assert(lambda.size() == used_atom_vectors[ii].size());

        for(auto && lambda : lambda_ranges_vec[jj])
            assert(lambda.size() == used_atom_vectors[jj].size());*/

      } else if (used_atom_vectors[ii].size() == used_atom_vectors[jj].size() &&
                 used_atom_vectors[ii] > used_atom_vectors[jj]) {
        /*for(auto && lambda : lambda_ranges_vec[ii])
            assert(lambda.size() == used_atom_vectors[ii].size());

        for(auto && lambda : lambda_ranges_vec[jj])
            assert(lambda.size() == used_atom_vectors[jj].size());*/
        used_atom_vectors[ii].swap(used_atom_vectors[jj]);
        std::swap(lambda_ranges_vec[ii], lambda_ranges_vec[jj]);
        /*for(auto && lambda : lambda_ranges_vec[ii]) {
            assert(lambda.size() == used_atom_vectors[ii].size());;
        }

        for(auto && lambda : lambda_ranges_vec[jj])
            assert(lambda.size() == used_atom_vectors[jj].size());*/
      }
    }
  }
}

/******************************************************************************
READ THE ATOM TYPES GIVEN IN THE INPUT, IN ORDER TO DEFINE THE CURRENT FRAGMENT
 ->atom_vectors contains all the atom vectors that are allowed for the current
fragment
 ->if none of the elements are alias_map, it contains just one atom vector
 ->otherwise, it will explicitly contain every possible combination of the
regular and alias atom types
******************************************************************************/

void AddAtomTypes(std::string& frg, size_t& j, std::vector<int>& xor_idx,
                  std::vector<int>& and_idx,
                  std::vector<AtomVector<combi_ff::Atom>>& atom_vectors,
                  const AliasMap& alias_map, const std::string& version) {
  std::string name("");

  // go through all the characters between the brackets [...] in frg
  while (j < frg.size()) {
    // case that atom name starts with a letter, or is an asterisk => add an
    // entry of atom type * to all atom_vectors
    if (frg[j] == '*') {
      for (auto&& a : atom_vectors) a.push_back(Atom("*"));

      // current atom is neither an XOR, nor an AND, indicated by the index -1
      xor_idx.push_back(-1);
      and_idx.push_back(-1);
    }

    // case that the atom type is a simple atom, or an alias
    else if (isalpha(frg[j])) {
      while (frg[j] != ' ') name += frg[j++];

      and_idx.push_back(-1);
      xor_idx.push_back(-1);
      // current atom is neither an XOR, nor an AND, indicated by the index -1
      AddAtomType(name, alias_map, atom_vectors, version);
    }

    // case that current atom is a united atom
    else if (frg[j] == '{') {
      // read the united atom name
      while (frg[++j] != '}' && isalpha(frg[j])) name += frg[j];

      // make sure that the united atom contains hydrogen
      if (name.back() != 'H')
        throw combi_ff::input_error(
            "united atom notation only works for hydrogen at the moment, e.g. "
            "{CH3}");

      // read number of hydrogen atoms. default is 1
      size_t n(1);

      if (isdigit(frg[j])) n = GetNumber(j, frg);

      // add the atom type toGether with the fixed number of hydrogen atoms
      for (auto&& a : atom_vectors) {
        a.push_back(name.substr(0, name.size() - 1));
        a.back().SetNumFixedHydrogens(n);
      }

      // current atom is neither an XOR, nor an AND, as indicated by index -1
      xor_idx.push_back(-1);
      and_idx.push_back(-1);

      if (frg[j++] != '}')
        throw combi_ff::input_error("expected \'}\' in " + frg +
                                    " at position " + std::to_string(j - 1) +
                                    " but encountered " +
                                    std::string(1, frg[j - 1]));

      // reSet name string
      name = "";
    }

    // case that current atom is an XOR
    else if (frg[j] == '^') {
      // read XOR idx. default is 0
      int n(0);

      if (isdigit(frg[++j])) n = (int)GetNumber(j, frg);

      xor_idx.push_back(n);

      // case that current atom is also an AND
      if (frg[j] == '&') {
        n = 0;

        if (isdigit(frg[++j])) n = (int)GetNumber(j, frg);

        and_idx.push_back(n);
      }

      // otherwise, not an AND, indicated by index -1
      else
        and_idx.push_back(-1);

      // read the atom name and add the atom to atom_vectors
      while (isalpha(frg[j])) name += frg[j++];

      AddAtomType(name, alias_map, atom_vectors, version);
    }

    // case that current atom is an AND
    else if (frg[j] == '&') {
      // read AND idx. default is 0
      int n(0);

      if (isdigit(frg[++j])) n = (int)GetNumber(j, frg);

      and_idx.push_back(n);

      // case that current atom is also an XOR
      if (frg[j] == '^') {
        n = 0;

        if (isdigit(frg[++j])) n = (int)GetNumber(j, frg);

        xor_idx.push_back(n);
      }

      // otherwise not an XOR, indicated by index -1
      else
        xor_idx.push_back(-1);

      // read the atom name and add the atom to atom_vectors
      while (isalpha(frg[j])) name += frg[j++];

      AddAtomType(name, alias_map, atom_vectors, version);
    }

    j++;
  }
}

/*******************************************************
USED TO DETERMINE IF AN ATOM WITH NAME name IS AN ALIAS
 ->if it's an alias, add all possible combinations to the atom_vectors
 ->if not, simply add the atom to the back of each AtomVector in atom_vectors
*******************************************************/
void AddAtomType(std::string& name, const AliasMap& alias_map,
                 std::vector<AtomVector<combi_ff::Atom>>& atom_vectors,
                 const std::string& version) {
  bool isAlias(false);
  // go through alias_map
  auto alias = alias_map.find(name);

  if (alias != alias_map.end()) {
    if (alias->second.GetVersion() != version)
      std::cout << "?Warning: using atom type alias " << name
                << " with version " << alias->second.GetVersion()
                << " but current family is version " << version << "\n";

    isAlias = true;
    size_t n = atom_vectors.size();
    std::vector<AtomVector<combi_ff::Atom>> atom_vectors_tmp;
    atom_vectors_tmp.reserve(n *
                             alias->second.GetContainedAtomTypesConst().size());

    // add to all currently existing atom_vectors ...addatom
    for (size_t i = 0; i < n; i++) {
      //... all possible alias_map
      for (auto map : alias->second.GetContainedAtomTypesConst()) {
        size_t j = 0;

        // case that an element alias is a united atom
        if (map[j] == '{') {
          std::string nam("");

          while (map[++j] != '}' && isalpha(map[j])) nam += map[j];

          if (nam.back() != 'H')
            throw combi_ff::input_error(
                "united atom notation only works for hydrogen atm, e.g. {CH3}");

          // read number of hydrogen atoms. default is 1
          size_t n(1);

          if (isdigit(map[j])) n = GetNumber(j, map);

          // add a copy of the i-th AtomVector, and append the current atom to
          // it
          atom_vectors_tmp.push_back(atom_vectors[i]);
          atom_vectors_tmp.back().push_back(
              Atom(nam.substr(0, nam.size() - 1)));
          atom_vectors_tmp.back().back().SetNumFixedHydrogens(n);

          if (map[j] != '}')
            throw combi_ff::input_error("expected \'}\' in " + map +
                                        " at position " + std::to_string(j) +
                                        " but encountered " +
                                        std::string(1, map[j]));
        }

        // otherwise directly add the atom defined by the current alias
        else {
          atom_vectors_tmp.push_back(atom_vectors[i]);
          atom_vectors_tmp.back().push_back(Atom(map));
        }
      }
    }

    // Set the atom_vectors to the atom_vectors_tmp
    atom_vectors = atom_vectors_tmp;
  }

  // if the atom is not an element alias, just add it to all the atom_vectors
  if (!isAlias) {
    for (auto&& a : atom_vectors) a.push_back(Atom(name));
  }

  // reSet the name string
  name = "";
}

/*************************************************
CLEAN UP THE ATOM VECTORS OF A FRAGMENT COLLECTION
*************************************************/
void CleanUpAtomVectorsFrag(
    std::vector<AtomVector<combi_ff::Atom>>& atom_vectors,
    std::vector<int>& xor_idx, std::vector<int>& and_idx) {
  // test all atom_vectors for violations of XOR or AND within their atoms
  for (int i = 0; i < (int)atom_vectors.size(); i++) {
    auto&& a = atom_vectors[i];
    bool del(false);

    // test all pairs of atoms for violations of XOR or AND
    for (size_t ii = 0; ii < a.size(); ii++) {
      // save XOR and AND idx, of first atom of the pair. save the atom name as
      // well. just for convenience, since it's reused
      int currXOR(xor_idx[ii]), currAND(and_idx[ii]);
      std::string curr_atom_type(a[ii].GetElementSymbol());

      for (size_t jj = ii + 1; jj < a.size(); jj++) {
        // if atom ii is an XOR
        if (currXOR != -1) {
          // if atom ii and jj have the same XOR index, but are of the same atom
          // type => violation!
          if (currXOR == xor_idx[jj] &&
              curr_atom_type == a[jj].GetElementSymbol())
            del = true;
        }

        // if current AtomVector was not deleted, and atom ii is an AND
        if (!del && currAND != -1) {
          // if atom ii and atom jj have the same AND index, but don't have the
          // same atom type => violation!
          if (currAND == and_idx[jj] &&
              curr_atom_type != a[jj].GetElementSymbol())
            del = true;
        }

        // if a violation was found, erase the current AtomVector
        if (del) {
          atom_vectors.erase(atom_vectors.begin() + i);
          i--;
          break;
        }
      }

      if (del) break;
    }
  }
}

/***********************
PRINT FAMILY INFORMATION
***********************/
void PrintFamilies(FamilyVector& families) {
  if (families.size()) {
    std::cout << "  -Families:\n";

    for (auto&& f : families) {
      std::cout << "    -";
      f.Print();
    }

  } else
    std::cout << "  -no Families\n";

  std::cout << "******************************************************\n";
}

}  // namespace enu

}  // namespace combi_ff