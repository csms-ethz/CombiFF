// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "InputOutput.h"

#include <stdio.h>

#include <algorithm>

#include "ContainerOperators.h"
#include "XmlParser.h"
#include "exceptions.h"
#include "printInfo.h"
#include "readLambdas.h"

namespace combi_ff {

namespace enu {

const std::string IOFileProperties::output_file_name_default = "isomerlist.xml";

// constructor for IOFileProperties
IOFileProperties::IOFileProperties()
    : rnd(rand()),
      input_file_names(FileNameVector(num_input_files)),
      output_dir(""),
      input_file_categories(StringVector(num_input_files)),
      input_file_name_building_blocks(StringVector(num_input_files)) {
  for (size_t f = 0; f < num_input_files; f++) {
    auto type =
        possible_input_file_properties.find(static_cast<InputFileType>(f));

    if (type == possible_input_file_properties.end())
      throw std::logic_error("input file type " + std::to_string(f) +
                             " not found\n");

    input_file_categories[type->first] = type->second.category;
    input_file_name_building_blocks[type->first] = type->second.name_block;
  }
}

// constructor for InputOutput
InputOutput::InputOutput(const StringVector& arguments) : arguments(arguments) {
  // read arguments to Get input and output file names
  ReadInputArguments();
}

// destructor for InputOutput
InputOutput::~InputOutput() {
  remove(io_file_properties.file_name_tmp.c_str());
}

// read the input options and use them to define the enum specifications and
// IOFiles
void InputOutput::ReadInputArguments() {
  if (!arguments.size()) {
    PrintInputOptions();
    throw combi_ff::input_error("no arguments found.");

  } else if (arguments[0].find("help") != std::string::npos) {
    PrintInputOptions();
    throw combi_ff::help_exception();
  }

  // go through the arguments vector
  for (size_t i = 0; i < arguments.size(); i++) GetNextInputOption(i);
}

// read the keywords in the arguments vector and perform the corresponding
// action
void InputOutput::GetNextInputOption(size_t& i) {
  // arg is the current argument
  std::string arg = arguments[i];

  if (arg.front() != '-')
    throw combi_ff::input_error(
        "expected a keyword starting with \'-\', but got \'" + arg +
        "\'. To view input options, run the program with command line argument "
        "\'help\'");

  if (arg == "-input")
    AddInputFile(i);

  else if (arg == "-formula")
    AddAtoms(i);

  else if (arg == "-families")
    AddFamilies(i);

  else if (arg == "-output")
    AddOutputFile(i);

  else if (arg == "-output_directory")
    AddOutputDir(i);

  else if (arg == "-max_bond_degree")
    AddMaxDegree(i);

  else if (possible_ranges.find(arg.substr(1, arg.size() - 1)) !=
           possible_ranges.end())
    AddRestriction(i,
                   possible_ranges.find(arg.substr(1, arg.size() - 1))->second);

  else if (arg == "-stereo")
    AddStereo();

  else if (possible_input_files.find(arg) != possible_input_files.end())
    ReadFileNames(i, possible_input_files.find(arg)->second);

  else
    throw combi_ff::input_error("unrecognized keyword \'" + arg +
                                "\' in input. To view input options, run the "
                                "program with command line argument \'help\'");
}

// read the content of a given input file and add it to arguments
void InputOutput::AddInputFile(size_t& i) {
  std::ifstream arguments_file;

  if (i + 1 > arguments.size() - 1)
    throw combi_ff::input_error("expected argument after -input");

  // the next argument in arguments corresponds to the filename after -input
  std::string arg = arguments[++i];
  arguments_file.open(arg.c_str());
  // read the input file, and add all the arguments to the arguments vector
  std::string next;

  while (arguments_file >> next) {
    // don't continue reading after a comment
    if (next[0] == '#') getline(arguments_file, next);

    // otherwise, add the current string to arguments vector
    else
      arguments.push_back(next);
  }

  arguments_file.close();
}

// read the formula given after the -atoms keyword
//  ->for each atom type, first read the atom type name, and then read the
//  different lambdas
void InputOutput::AddAtoms(size_t& i) {
  // create a local atom vec and a local vector of LambdaVectors
  AtomVector<combi_ff::Atom> used_atoms(0);
  std::list<LambdaVector> lambda_ranges(0);
  lambda_ranges.push_back(LambdaVector(0));
  // read in the whole formula given on the current line
  std::string formula("");

  while (i < arguments.size() - 1 && arguments[i + 1][0] != '-')
    formula += arguments[++i];

  size_t j = 0;

  // go through formula to consecutively add atom types and ranges
  while (j < formula.size()) {
    // if an atom name starts with '{', it's a united atom
    if (formula[j] == '{') AddUnitedAtom(j, formula, used_atoms);

    // if an atom name starts with a letter, read the name
    else if (isalpha(formula[j]))
      AddAtom(j, formula, used_atoms);

    else
      throw combi_ff::input_error(
          "format error in formula: expected an alphabetical character at "
          "position " +
          std::to_string(j) + " in " + formula);

    // read the lambdas corresponding to the current atom type
    ReadLambdas(lambda_ranges, j, formula);
  }

  // make sure that all the atoms in the final used_atoms vector are unique. if
  // necessary, combine two entries in used_atoms/lambda_ranges
  for (j = 0; j < used_atoms.size(); j++) {
    for (int k = (int)j + 1; k < (int)used_atoms.size(); k++) {
      if (used_atoms[j] == used_atoms[k]) {
        // std::cout << "combining atom " << j << " and " << k << " as they have
        // the same type\n";
        for (auto&& l : lambda_ranges) {
          l[j] += l[k];
          l.erase(l.begin() + k);
        }

        used_atoms.erase(used_atoms.begin() + k);
        k--;
      }
    }
  }

  // sort Atoms in decreasing order of degrees, and alphabetically, and make
  // sure that all the resulting vectors in lambda_ranges are unique
  SortAtoms(used_atoms, lambda_ranges);
  lambda_ranges.sort();
  auto&& it = std::unique(lambda_ranges.begin(), lambda_ranges.end());
  lambda_ranges.resize(std::distance(lambda_ranges.begin(), it));
  // add used_atoms and lambda_ranges to the usedAtomVectors and the
  // lambda_rangesVec
  enum_spec.used_atom_vectors.push_back(used_atoms);
  enum_spec.lambda_ranges_vec.push_back(lambda_ranges);
}

// adds files of inputType to tempInputFile and opens the tempInputFile in
// inputFile
void InputOutput::ReadFileNames(size_t& i, const InputFileType inputType) {
  const std::string& name =
      io_file_properties.input_file_name_building_blocks[inputType];
  std::string line;
  std::string keyword = arguments[i];
  std::string arg = arguments[++i];

  // read the next argument
  do {
    std::cout << "using " << arg << " as " << name << " input.\n";
    io_file_properties.input_file_names[inputType].push_back(arg);

    if (++i == arguments.size()) break;

    arg = arguments[i];
  } while (arg[0] != '-');

  i--;
}

// read the family names after the -families keyword
void InputOutput::AddFamilies(size_t& i) {
  while (i + 1 < arguments.size() && arguments[i + 1][0] != '-')
    enum_spec.used_families.push_back(arguments[++i]);

  std::sort(enum_spec.used_families.begin(), enum_spec.used_families.end());
  auto&& it = std::unique(enum_spec.used_families.begin(),
                          enum_spec.used_families.end());

  if (it != enum_spec.used_families.end()) {
    while (it != enum_spec.used_families.end()) {
      std::cout << "!Error:listed " << *it << " more than once\n";
      ++it;
    }

    throw combi_ff::input_error(
        "Please make sure not to list a family twice in the -families "
        "argument");
  }
}

// read in the output file name for the direct enumeration after the -output
// keyword
void InputOutput::AddOutputFile(size_t& i) {
  if (i + 1 >= arguments.size())
    throw combi_ff::input_error("expected argument after -output");

  if (io_file_properties.output_file_name == "" ||
      io_file_properties.output_file_name ==
          IOFileProperties::output_file_name_default)
    io_file_properties.output_file_name = arguments[++i];

  else
    throw combi_ff::input_error("Please specify only one output file name");
}
// read in the output directory after the -output_directory keyword
void InputOutput::AddOutputDir(size_t& i) {
  if (i + 1 >= arguments.size())
    throw combi_ff::input_error("expected argument after -output_directory");

  io_file_properties.output_dir = arguments[++i];

  if (io_file_properties.output_dir.back() != '/')
    io_file_properties.output_dir += '/';
}

// read in the maximum bond degree parameter after the -maxDeg keyword
void InputOutput::AddMaxDegree(size_t& i) {
  if (i + 1 >= arguments.size())
    throw combi_ff::input_error("expected integer argument after -maxDeg");

  try {
    enum_spec.max_degree = std::stoi(arguments[++i]);

  } catch (const std::invalid_argument& ia) {
    throw combi_ff::input_error(
        "expected integer argument after -maxDeg, but encountered " +
        arguments[i] + " . " + ia.what());
  }
}
// read in the range of a rangedProperty
void InputOutput::AddRestriction(size_t& i,
                                 const size_t position_in_range_vec) {
  Range& r = enum_spec.ranged_properties[position_in_range_vec];
  std::string prop("");

  while (i < arguments.size() - 1 && arguments[i + 1][0] != '-')
    prop += arguments[++i];

  ReadRange(prop, r);
}
// Set stereo to true if the -stereo keyword is used
void InputOutput::AddStereo() {
  std::cout << "enumerating molecules including stereoisomerism\n";
  enum_spec.stereo = true;
}

// read the atom type of an atom and add it to used_atoms
void InputOutput::AddAtom(size_t& j, std::string& formula,
                          AtomVector<combi_ff::Atom>& used_atoms) {
  std::string nam;

  while (j < formula.size() && isalpha(formula[j]) && formula[j] != '[')
    nam += formula[j++];

  used_atoms.push_back(Atom(nam));
}
// read the atom type of a united atom, Set the number of fixed hydrogens and
// add the united atom to used_atoms
void InputOutput::AddUnitedAtom(size_t& j, std::string& formula,
                                AtomVector<combi_ff::Atom>& used_atoms) {
  std::string nam;

  // read the atom name
  while (++j < formula.size() && formula[j] != '}' && isalpha(formula[j]))
    nam += formula[j];

  if (nam.back() != 'H')
    throw combi_ff::input_error(
        "Error: united atom notation only works for hydrogen, use e.g. {CH3}");

  // add an atom with typeName "nam" to the used_atoms vector, leaving out the
  // last letter corresponding to the hydrogen atom
  used_atoms.push_back(Atom(nam.substr(0, nam.size() - 1)));

  // read the number of hydrogen atoms in the united atom
  if (isdigit(formula[j])) {
    size_t n = GetNumber(j, formula);
    used_atoms.back().SetNumFixedHydrogens(n);
  }

  // if there's no number, it's assumed that there's 1 hydrogen
  else
    used_atoms.back().SetNumFixedHydrogens(1);

  // after the name and number have been read, the united atom formula needs to
  // be closed with a '}' advance j to the next position after the '}'
  if (formula[j++] != '}')
    throw combi_ff::input_error("expected \'}\' at end of united atom");
}

/*GetTER FUNCTIONS*/
std::list<std::string>& InputOutput::GetInputFileNamesAt(size_t idx) {
  return io_file_properties.input_file_names[idx];
}
const std::string& InputOutput::GetOutputFileName() const {
  return io_file_properties.output_file_name;
}
const IOFileProperties& InputOutput::GetIOFilProps() const {
  return io_file_properties;
}
const EnumSpecifications& InputOutput::GetEnumSpec() const { return enum_spec; }

/**********************************************
PRINT THE INPUT OUTPUT OPTIONS THAT CAN BE USED
**********************************************/
void InputOutput::PrintInputOptions() {
  std::cout
      << "input options are:\n\n"
      << "-input: used to specify an input file that can contain one or "
         "several of the following arguments.\n"
      << "   - Note that these arguments could also be used directly on the "
         "command line\n"
      << "   - lines starting with a # are considered comments and are "
         "ignored\n\n"
      << "-output_directory: used to specify in which directory the '.fie' "
         "files should be saved\n\n"
      << "-output: used to specify the file name of the '.fie' file that is "
         "created for the formulas given by -atoms\n\n"
      << "-atoms: used to specify a molecular formula. The following formats "
         "can be used to specify the amounts and types of the atoms:\n"
      << "   - Xn, where X is an atom type and n is a number\n"
      << "   - X[n-m], where X is an atom type and n-m is a range of nubmers\n"
      << "   - X[n,m,l,...,], where n,m,l,... are specific numbers to be used\n"
      << "   - the notation of [n-m] and [n,m,...] can be combined at will, "
         "e.g. [1-3,4,10-8, 7]\n"
      << "   - X*, in order to specify any amount of atom type X (in practice "
         "at the moment simply corresponds to range of 0-100)\n"
      << "   - X>n or X>=n, where X is an atom type and n is a number, in "
         "order to specify that there should at least n+1, or n, atoms of type "
         "X\n"
      << "   - X<n or X<=n, where X is an atom type and n is a number, in "
         "order to specify that there should be at most n-1, or n, atoms of "
         "type X\n"
      << "   - {XHn} can be used instead of X for all of the above options in "
         "order to specify a united atom, where X is an atom type and n is a "
         "number,\n     specifying the amount of hydrogen atoms in the united "
         "atom\n"
      << "   - the -atoms keyword can be used several times, in order to "
         "specify different molecules\n"
      << "   - example: -atoms C10H12Cl>1 F[2,3,5-4] -atoms C10 {CH3}4 H*\n\n"
      << "-max_bond_degree: used to specify maximum bond degree between two "
         "atoms\n"
      << "   - specify by using -maxDeg n, where n is the desired maximum "
         "degree number\n"
      << "   - defaults to 4, if not specified\n\n"
      << "-unsaturations: used to specify the desired number of unsaturations. "
         "The following formats can be used:\n"
      << "   - -unsaturations n, in order to specify that there should be "
         "exactly n unsaturations\n"
      << "   - -unsaturations [n-m], in order to specify that there should be "
         "between n and m unsaturations\n"
      << "   - -unsaturations <n, or <=n, in order to specify that there "
         "should be at most n-1, or n, unsaturations\n"
      << "   - -unsaturations >n, or >=n, in order to specify that there "
         "should be at least n+1, or n, unsaturations\n"
      << "   - defaults to unrestricted, if not specified\n\n"
      << "-total_bonds: used to specify the desired number of bonds\n"
      << "   - same formats as for -unsaturations can be used\n\n"
      << "-single_bonds: used to specify the desired number of single bonds "
         "(incl. H-bonds)\n"
      << "   - same formats as for -unsaturations can be used\n\n"
      << "-double_bonds: used to specify the desired number of double bonds\n"
      << "   - same formats as for -unsaturations can be used\n\n"
      << "-triple_bonds: used to specify the desired number of triple bonds\n"
      << "   - same formats as for -unsaturations can be used\n\n"
      << "-quadruple_bonds: used to specify the desired number of quadruple "
         "bonds\n"
      << "   - same formats as for -unsaturations can be used\n\n"
      << "-cycles-: used to specify the desired number of rings\n"
      << "   - same formats as for -unsaturations can be used\n\n"
      << "-substructure_files: used to specify the input file(s) in which the "
         "Substructures are defined. Several files can be listed in a row\n\n"
      << "-family_files: used to specify the input file(s) in which the Family "
         "Library is defined. Several files can be listed in a row\n\n"
      << "-element_alias_files: used to specify the input file(s) in which the "
         "elementAliases are defined. Several files can be listed in a row\n\n"
      << "-pseudoatom_files: used to specify the input file(s) in which the "
         "pseudoatoms are defined. Several files can be listed in a row\n\n"
      << "-families: used to specify the desired families, as defined in the "
         "Family Library file(s)\n"
      << "   - specified by giving the CODE of the family, as defined in the "
         "Family Library file(s)\n"
      << "   - several families can be listed after each other, or the -family "
         "keyword can be used several times in order to specify multiple "
         "familes\n\n"
      << "-stereo: used to specify that stereoisomerism should also be "
         "considered in the enumeration\n\n";
}

}  // namespace enu

}  // namespace combi_ff