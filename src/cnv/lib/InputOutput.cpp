// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "InputOutput.h"

#include <iostream>

#include "exceptions.h"

namespace combi_ff {

namespace cnv {

/**********
CONSTRUCTOR
**********/
InputOutput::InputOutput(
    std::list<std::pair<std::string, std::string>>& input_list,
    std::vector<bool>& print_options, InputOption& input,
    std::string& output_file_name,
    combi_ff::StringVector& family_enumeration_file_names, int& argc,
    char* argv[])
    : input_list(input_list),
      print_options(print_options),
      input(input),
      output_file_name(output_file_name),
      family_enumeration_file_names(family_enumeration_file_names),
      input_file_names(combi_ff::StringVector(0)) {
  // add command line arguments to arguments vector
  arguments = (combi_ff::StringVector(0));

  for (int ii = 1; ii < argc; ii++) arguments.push_back(std::string(argv[ii]));
}

/*******************************************************************************
READ THE GIVEN OR DEFAULT TBFragmentFiles, familyLibraryFiles, and moleculeFiles
 use them to define the required properties
*******************************************************************************/
void InputOutput::ReadArguments() {
  if (arguments.size() == 0) {
    InputOutput::PrintInputOptions();
    throw combi_ff::input_error("no arguments found");

  } else if (arguments.front().find("help") != std::string::npos) {
    InputOutput::PrintInputOptions();
    throw combi_ff::help_exception();
  }

  std::string arg;

  // go through the arguments vector
  // for(size_t i = 0; i < arguments.size(); i++) {
  for (size_t i = 0; i < arguments.size(); i++) {
    // arg is the current argument
    arg = arguments[i];

    /*if (arg == "-fie") {
      if (!family_enumeration_file_names.size()) {
        while (i + 1 < arguments.size())
          family_enumeration_file_names.push_back(arguments[++i]);
      }

    } else */
    if (arg.front() == '-') {
      size_t found_option(0);

      if (arg == "-h") {
        InputOutput::PrintInputOptions();
        throw combi_ff::help_exception();
      }

      size_t find = arg.find("I");

      if (find != (size_t)std::string::npos) {
        found_option++;
        InputOutput::ReadInputOption(i + find, arguments);
        arg.erase(arg.begin() + find);
      }

      find = arg.find("O");

      if (find != (size_t)std::string::npos) {
        found_option++;
        InputOutput::ReadOutputOption(i + find, arguments);
        arg.erase(arg.begin() + find);
      }

      find = arg.find("i");

      if (find != (size_t)std::string::npos) {
        found_option++;

        if (i + find >= arguments.size())
          throw combi_ff::input_error("expected argument for -i option");

        GetInputFromFile(i + find);
        arg.erase(arg.begin() + find);
      }

      find = arg.find("o");

      if (find != (size_t)std::string::npos) {
        found_option++;

        if (i + find >= arguments.size())
          throw combi_ff::input_error("expected argument for -o option");

        output_file_name = arguments[i + find];
        arg.erase(arg.begin() + find);
      }

      /*find = arg.find("s");

      if (find != (size_t)std::string::npos) {
        found_option++;

        if (i + find >= arguments.size())
          throw combi_ff::input_error("expected argument for -s option");

        InputOutput::ReadFieFileNamesFromSetupFile(arguments[i + find]);
        arg.erase(arg.begin() + find);
      }*/

      i += found_option;

      if (arg != "-")
        throw combi_ff::input_error("uncrecognized option " +
                                    arg.substr(1, arg.size() - 1));

      else if (!found_option)
        throw combi_ff::input_error("uncrecognized option " + arg);

    } else {
      while (i < arguments.size() && arguments[i].front() != '-')
        input_list.push_back(std::pair<std::string, std::string>(
            "command_line", arguments[i++]));

      if (i < arguments.size() && arguments[i].front() == '-') i--;
    }
  }
}

/*************************************************************
READ THE CONTENT OF A GIVEN INPUT FILE AND ADD IT TO arguments
*************************************************************/
void InputOutput::GetInputFromFile(const size_t pos) {
  std::ifstream arguments_file;
  // the next argument in arguments corresponds to the filename after @input
  // std::string arg = arguments[++i];
  input_file_names = combi_ff::StringVector(1, arguments[pos]);

  for (int i = 0; i < (int)input_file_names.size(); i++) {
    for (int j = 0; j < (int)input_file_names[i].size(); j++) {
      if ((input_file_names[i][j] == ' ' || input_file_names[i][j] == ',') &&
          input_file_names[i].size() > 1) {
        input_file_names.push_back(input_file_names[i].substr(
            j + 1, input_file_names[i].size() - (j + 1)));
        input_file_names[i] = input_file_names[i].substr(0, j);
        j--;
      }
    }

    if (!input_file_names[i].size()) {
      input_file_names.erase(input_file_names.begin() + i);
      i--;
    }
  }

  for (auto&& fileName : input_file_names) {
    arguments_file.open(fileName.c_str());
    // read the input file, and add all the arguments to the arguments vector
    std::string next;

    while (getline(arguments_file, next)) {
      // don't continue reading after a comment
      if (next.front() == '#' || next.front() == '%')
        getline(arguments_file, next);

      // otherwise, add the current string to arguments vector
      else {
        if (next.back() == 13)  // catch carriage return
          next.pop_back();

        input_list.push_back(
            std::pair<std::string, std::string>(fileName, next));
      }
    }

    arguments_file.close();
  }
}

void InputOutput::ReadInputOption(const size_t pos,
                                  const combi_ff::StringVector& arguments) {
  if (pos >= arguments.size())
    throw combi_ff::input_error("expected argument after -I");

  std::string arg = arguments[pos];

  if (arg == "frm")
    input = formula;

  else if (arg == "smi")
    input = smiles;

  // else if (arg == "nam")
  //   input = name;

  // else if (arg == "fmi")
  //   input = family_enumeration;

  else if (arg == "mat")
    input = matrix;

  else
    throw combi_ff::input_error("unrecognized option for -I");
}

const std::string InputOutput::IncompatibilityMessage(
    const std::string& input, const std::string& output) const {
  return "input option \'" + input + "\' not compatible with output option \'" +
         output;
}

void InputOutput::ReadOutputOption(const size_t pos,
                                   const combi_ff::StringVector& arguments) {
  if (pos >= arguments.size())
    throw combi_ff::input_error("expected argument after -I");

  std::string arg = arguments[pos];
  bool found_option(false);
  size_t find = arg.find("smi");

  if (find != (size_t)std::string::npos) {
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "smi"));

    /*else if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "smi"));*/

    found_option = true;
    print_options[cnv::print_canon_smiles] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 3);
  }

  find = arg.find("frm");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "frm"));*/

    found_option = true;
    print_options[cnv::print_canon_formula] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 3);
  }

  /*find = arg.find("nam");

  if (find != (size_t)std::string::npos) {
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "nam"));

    else if (input == smiles)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("smi", "nam"));

    else if (input == family_enumeration)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("fmi", "nam"));

    else if (input == matrix)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("mat", "nam"));

    found_option = true;
    print_options[cnv::print_canon_name] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 3);
  }*/

  find = arg.find("mass");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "mass"));*/

    found_option = true;
    print_options[cnv::print_mass] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 4);
  }

  /*find = arg.find("fmi");

  if (find != (size_t)std::string::npos) {
    if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "fmi"));

    else if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "fmi"));

    else if (input == family_enumeration)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("fmi", "fmi"));

    found_option = true;
    print_options[cnv::print_family_enumeration] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 3);
  }*/

  find = arg.find("num_sat");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "num_sat"));*/

    found_option = true;
    print_options[cnv::print_n_unsaturations] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 7);
  }

  find = arg.find("num_bnd");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "num_mul"));

    else */
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "num_bnd"));

    found_option = true;
    print_options[cnv::print_n_bonds] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 7);
  }

  find = arg.find("num_sin");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "num_mul"));

    else */
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "num_sin"));

    found_option = true;
    print_options[cnv::print_n_single_bonds] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 7);
  }

  find = arg.find("num_mul");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "num_mul"));

    else */
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "num_mul"));

    found_option = true;
    print_options[cnv::print_n_multiple_bonds] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 7);
  }

  find = arg.find("num_dbl");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "num_dbl"));

    else */
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "num_dbl"));

    found_option = true;
    print_options[cnv::print_n_double_bonds] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 7);
  }

  /*find = arg.find("num_arm");

  if (find != (size_t)std::string::npos) {

    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "num_arm"));

    found_option = true;
    print_options[cnv::print_n_aromatic_bonds] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 7);
  }*/

  find = arg.find("num_tri");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "num_tri"));

    else */
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "num_tri"));

    found_option = true;
    print_options[cnv::print_n_triple_bonds] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 7);
  }

  find = arg.find("num_cyc");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "num_cyc"));

    else */
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "num_cyc"));

    found_option = true;
    print_options[cnv::print_n_cycles] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 7);
  }

  find = arg.find("atmV");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "atmV"));

    else */
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "atmV"));

    found_option = true;
    print_options[cnv::print_canon_atom_vector] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 4);
  }

  find = arg.find("stack");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "stack"));

    else */
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "stack"));

    found_option = true;
    print_options[cnv::print_stack] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 5);
  }

  find = arg.find("mat");

  if (find != (size_t)std::string::npos) {
    /*if (input == name)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("nam", "mat"));

    else */
    if (input == formula)
      throw combi_ff::input_error(
          InputOutput::IncompatibilityMessage("frm", "mat"));

    found_option = true;
    print_options[cnv::print_matrix] = true;
    arg.erase(arg.begin() + find, arg.begin() + find + 3);
  }

  while (arg.find(" ") != std::string::npos)
    arg.erase(arg.begin() + arg.find(" "));

  while (arg.find(",") != std::string::npos)
    arg.erase(arg.begin() + arg.find(","));

  if (arg != "")
    throw combi_ff::input_error("unrecognized output option " + arg);

  else if (!found_option)
    throw combi_ff::input_error("expected an argument for -O");
}

void InputOutput::ReadFieFileNamesFromSetupFile(
    const std::string& setup_file_name) {
  std::ifstream setup_file(setup_file_name);

  if (!setup_file.is_open())
    throw combi_ff::input_error("manual setup file " + setup_file_name +
                                " not open");

  std::string next("");

  while (setup_file >> next) {
    if (next.front() == '#')
      std::getline(setup_file, next);

    else
      family_enumeration_file_names.push_back(next);
  }
}

const combi_ff::StringVector& InputOutput::GetInputFileNames() const {
  return input_file_names;
}

/**********************************************
PRINT THE INPUT OUTPUT OPTIONS THAT CAN BE USED
**********************************************/
void InputOutput::PrintInputOptions() {
  std::cout
      << "input options are:\n\n"
      << " -I : used to define the format of the input. takes only one "
         "argument. possible arguments are:\n"
      << "      - smi:   smiles strings\n"
      << "      - frm:   molecular formulas\n"
      //<< "      - nam:   molecule names\n"
      //<< "      - fmi:   family molecular identifiers (i.e. codes from the fie
      //"
      //   "files)\n"
      << "      - mat:   adjacency matrix\n\n"
      << " -O : used to define the desired output. takes at least one "
         "argument.\n"
      << "      arguments can be written all together without any whitespaces, "
         "or spearated by commas, or within single quotes and separated by "
         "commas and/or whitespaces\n"
      << "      possible arguments are:\n"
      << "      - smi:     print the canonical smiles string\n"
      << "      - frm:     print the canonical formula\n"
      << "      - mass:    print the molecular mass\n"
      //<< "      - fmi:     print the code of the corresponding canonical "
      //   "smiles string from the fie files\n"
      << "      - num_sat: print number of unsaturations\n"
      << "      - num_bnd: print number of bonds\n"
      << "      - num_sin: print number of single bonds\n"
      << "      - num_mul: print number of multiple bonds\n"
      << "      - num_dbl: print number of double bonds\n"
      //<< "      - num_arm: print number of aromatic bonds\n"
      << "      - num_tri: print number of triple bonds\n"
      << "      - num_cyc: print number of cycles\n"
      << "      - stack:   print the stack of the canonical adjacency matrix\n"
      << "      - atmV:    print the atom vector of the canonical matrix\n"
      << "      - mat:     print the canonical adjacency matrix\n\n"
      << " -i : used to define an input file. several files can be given "
         "separated by commas, or within single quotes separated by commas "
         "and/or whitespaces\n"
      << "      lists of smiles strings can be separated by commas, periods, "
         "and/or whitespaces.\n"
      << "      lists of formulas, names, or family identifiers can be "
         "separated by commas and/or whitespaces\n"
      << "      Note: if no input file is given, the list of arguments (or a "
         "single argument) can be given directly on the command line.\n"
      << " -o : used to define an output file\n\n"
      //<< " -s : used to overwrite the default option usr/default.cnv for the "
      //   "Setup file\n\n"
      << " -h : display help message\n\n";
}

}  // namespace cnv

}  // namespace combi_ff