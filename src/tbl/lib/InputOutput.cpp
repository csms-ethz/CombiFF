// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "InputOutput.h"

#include <cassert>
#include <iostream>
#include <sstream>

#include "exceptions.h"

namespace combi_ff {

namespace topology_builder {

// constructor for IOFileProperties
IOFileProperties::IOFileProperties()
    : fie_families(StringVector(0)),
      input_file_names(FileNameVector(num_input_files)),
      output_dir(""),
      output_dir_molecule_decompositions(""),
      output_dir_molecules_with_macros(""),
      output_dir_mtb(""),
      input_file_categories(StringVector(num_input_files)),
      input_file_name_building_blocks(StringVector(num_input_files)) {
  for (size_t f = 0; f < num_input_files; f++) {
    auto type = possible_input_file_props.find(static_cast<InputFileType>(f));

    if (type == possible_input_file_props.end())
      throw std::logic_error("input file type " + std::to_string(f) +
                             " not found\n");

    input_file_categories[type->first] = type->second.category;
    input_file_name_building_blocks[type->first] = type->second.name_block;
  }
}

/**********
InputOutput
**********/
InputOutput::InputOutput(const std::vector<std::string>& arguments)
    : fragments_to_use(std::vector<std::string>(0)), arguments(arguments) {
  i = 0;
  ReadInput();
}

/******************************
open the input and output files
 use the input files to define the paramters given by main
******************************/
void InputOutput::ReadInput() {
  // read arguments and open corresponding input and output files. use input
  // files to define required properties
  ReadArguments();
}

/*******************************************************************************
READ THE GIVEN OR DEFAULT TBFragmentFiles, familyLibraryFiles, and moleculeFiles
 use them to define the required properties
*******************************************************************************/
void InputOutput::ReadArguments() {
  if (arguments.size() == 0) {
    PrintInputOptions();
    throw input_error("no arguments found.\n");

  } else if (arguments[0].find("help") != std::string::npos) {
    PrintInputOptions();
    throw help_exception();
  }

  std::string arg;

  // go through the arguments vector
  // for(size_t i = 0; i < arguments.size(); i++) {
  for (i = 0; i < arguments.size(); i++) {
    // arg is the current argument
    arg = arguments[i];

    // command line option -input . add file contents to arguments vector
    if (arg == "-input") GetInputFile();

    // command line option -fragments
    else if (possible_input_files.find(arg) != possible_input_files.end())
      ReadFileNames(i, possible_input_files.find(arg)->second);

    else if (arg == "-families") {
      while (i + 1 < arguments.size()) {
        arg = arguments[++i];

        if (arg[0] == '-') {
          i--;
          break;

        } else
          io_file_properties.fie_families.push_back(arg);
      }
    }

    // command line optoin -outputDir
    else if (arg == "-output_directory") {
      io_file_properties.output_dir = arguments[++i];

      if (io_file_properties.output_dir.back() != '/')
        io_file_properties.output_dir += '/';

      if (io_file_properties.output_dir_molecule_decompositions.size()) {
        std::cout
            << "?Warning: using output directory for molecule decompositions '"
            << io_file_properties.output_dir_molecule_decompositions
            << "'. ignoring general output directory '"
            << io_file_properties.output_dir << "'\n";
      }

      if (io_file_properties.output_dir_molecules_with_macros.size()) {
        std::cout
            << "?Warning: using output directory for molecules with macros '"
            << io_file_properties.output_dir_molecules_with_macros
            << "'. ignoring general output directory '"
            << io_file_properties.output_dir << "'\n";
      }

      if (io_file_properties.output_dir_mtb.size()) {
        std::cout << "?Warning: using output directory for mtb '"
                  << io_file_properties.output_dir_mtb
                  << "'. ignoring general output directory '"
                  << io_file_properties.output_dir << "'\n";
      }

    } else if (arg == "-output_directory_decompositions") {
      io_file_properties.output_dir_molecule_decompositions = arguments[++i];

      if (io_file_properties.output_dir_molecule_decompositions.back() != '/')
        io_file_properties.output_dir_molecule_decompositions += '/';

      if (io_file_properties.output_dir.size()) {
        std::cout
            << "?Warning: using output directory for molecule decompositions '"
            << io_file_properties.output_dir_molecule_decompositions
            << "'. ignoring general output directory '"
            << io_file_properties.output_dir << "'\n";
      }

    } else if (arg == "-output_directory_molecule_macros") {
      io_file_properties.output_dir_molecules_with_macros = arguments[++i];

      if (io_file_properties.output_dir_molecules_with_macros.back() != '/')
        io_file_properties.output_dir_molecules_with_macros += '/';

      if (io_file_properties.output_dir.size()) {
        std::cout
            << "?Warning: using output directory for molecules with macros '"
            << io_file_properties.output_dir_molecules_with_macros
            << "'. ignoring general output directory '"
            << io_file_properties.output_dir << "'\n";
      }

    } else if (arg == "-output_directory_mtb") {
      io_file_properties.output_dir_mtb = arguments[++i];

      if (io_file_properties.output_dir_mtb.back() != '/')
        io_file_properties.output_dir_mtb += '/';

      if (io_file_properties.output_dir.size()) {
        std::cout << "?Warning: using output directory for mtb '"
                  << io_file_properties.output_dir_mtb
                  << "'. ignoring general output directory '"
                  << io_file_properties.output_dir << "'\n";
      }
    }

    // command line option -fragments_to_use
    else if (arg == "-fragments_to_use") {
      while (i < arguments.size() - 1 && arguments[i + 1][0] != '-')
        fragments_to_use.push_back(arguments[++i]);

    } else if (arg == "-united_atom")
      io_file_properties.united_atom = true;

    else if (arg == "-third_neighbor_exclusions")
      io_file_properties.third_neighbor_exclusions = true;

    else if (arg == "-unique_torsionals")
      io_file_properties.unique_torsionals = true;

    else {
      PrintInputOptions();
      std::cout << std::endl
                << std::endl
                << "*******************" << std::endl
                << std::endl;
      throw input_error("Unrecognized input option " + arg +
                        ". See above for possible options\n");
    }
  }
}

// read names of the input files
void InputOutput::ReadFileNames(size_t& i, const InputFileType input_type) {
  const std::string& name =
      io_file_properties.input_file_name_building_blocks[input_type];
  std::string line;
  std::string keyword = arguments[i];
  std::string arg = arguments[++i];

  // read the next argument
  do {
    std::cout << "using " << arg << " as " << name << " input.\n";
    io_file_properties.input_file_names[input_type].push_back(arg);

    // read the next argument
    if (++i == arguments.size()) break;

    arg = arguments[i];
  } while (arg[0] != '-');

  i--;
}

/*************************************************************
READ THE CONTENT OF A GIVEN INPUT FILE AND ADD IT TO arguments
*************************************************************/
void InputOutput::GetInputFile() {
  std::ifstream arguments_file;
  // the next argument in arguments corresponds to the filename after -input
  std::string arg = arguments[++i];
  arguments_file.open(arg.c_str());

  if (!arguments_file.is_open()) throw std::runtime_error(arg + " not open\n");

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

const IOFileProperties& InputOutput::GetIOFilProps() const {
  return io_file_properties;
}

const std::vector<std::string> InputOutput::GetFragmentsToUse() const {
  return fragments_to_use;
}

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
      << "-fragments: used to specify the input file(s) in which the Topology "
         "Builder Fragments are defined. Several files can be listed in a row\n"
      << "   - defaults to ./fragments/fragments_alkanes, if not specified\n\n"
      << "-fragments_to_use: used to specify the Fragments to use by listing "
         "the fragment codes\n"
      << "   - if not specified, all fragments are used\n\n"
      << "-families: for which families to create the decomposition\n\n"
      << "-united_atoms: should united atoms be used?\n\n"
      << "-third_neighbor_exclusions: should third neighbors be excluded?\n\n"
      << "-unique_torsionals: should only one torsional dihedral per central "
         "bond be added to the topology?\n\n"
      << "-output_directory: used to specify the path where general output "
         "files should be stored\n\n"
      << "-output_directory_decompositions: used to specify the path where the "
         "molecule decomposition files should be stored\n\n"
      << "-output_directory_molecule_macros: used ot specify the path where "
         "the molecules with macros files should be stored\n\n"
      << "-output_directory_mtb: used to specify the path where the mtb files "
         "should be stored\n\n";
}

}  // namespace topology_builder

}  // namespace combi_ff