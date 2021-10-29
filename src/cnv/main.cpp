#include<iostream>
#include<fstream>
#include "SmilesHandler.h"
#include "FormulaHandler.h"
#include "NameHandler.h"
#include "FamilyIdentifierHandler.h"
#include "AdjacencyMatrixHandler.h"
#include <ctype.h>
#include <ostream>

void convertFormula(const std::string& formula, std::ofstream& output_file,
                    std::vector<bool> print_options) ;



int main(int argc, char* argv[]) {
  try {
    std::list<std::pair<std::string, std::string>> input_list(0);
    std::string output_file_name("");
    std::string setup_file_name("usr/default.cnv");
    combi_ff::StringVector fie_file_names(0);
    /*
      print_options[print_canon_smiles] : print canonical smiles
      printOptinos[print_canon_name]   : print canonical name
      print_options[print_canon_formula] : print canonical formula
      print_options[print_mass]   : print mass
      print_options[print_family_enumeration]    : print family identifier
      print_options[print_n_unsaturations] : print number of unsaturations
      print_options[print_n_multiple_bonds]    : print number of mulitple bonds
      print_options[print_n_double_bonds]    : print number of double bonds
      print_options[print_n_quadruple_bonds]    : print number of quadruple bonds
      print_options[print_n_aromatic_bonds]    : print number of aromatic bonds
      print_options[print_n_triple_bonds]    : print number of triple bonds
      print_options[print_n_cycles]   : print number of cycles
      print_options[print_canon_atom_vector]  : print atom vector of canonical matrix
      print_options[print_stack]  : print stack of the canonical adjacency matrix
      print_options[print_matrix]    : print canonical adjacency matrix
    */
    std::vector<bool> print_options(combi_ff::cnv::num_print_options, false);
    combi_ff::cnv::InputOption input(combi_ff::cnv::not_set);
    combi_ff::cnv::InputOutput IO(input_list,
                                  print_options,
                                  input,
                                  output_file_name,
                                  fie_file_names,
                                  argc,
                                  argv);

    /**************************************************************
    OPEN INPUT AND OUTPUT FILES. READ IN ALL THE VARIOUS PARAMETERS
    **************************************************************/
    try {
      IO.ReadArguments();
    }

    //catch incorrect user input
    catch (combi_ff::input_error& e) {
      std::cerr << "!Input Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }

    //case that './cnv help' was used, this is not a failure, but the program should stop
    catch (combi_ff::help_exception& e) {
      return EXIT_SUCCESS;
    }

    //catch any other errors
    catch (std::exception& e) {
      std::cerr << "!Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }

    if (!input_list.size()) {
      std::cerr << "!Error: no input found to convert."
                << " Please give a (list of) smiles, molecular formulas, names, or family identifiers\n";
      return EXIT_FAILURE;
    }

    std::ofstream output_file;

    if (output_file_name != "") {
      output_file.open(output_file_name.c_str());

      if (!output_file.is_open()) {
        std::cerr << "!Error: output file " << output_file_name << " not open\n";
        return EXIT_FAILURE;
      }
    }

    if (!fie_file_names.size() && (input == combi_ff::cnv::family_enumeration ||
                                   print_options[combi_ff::cnv::print_family_enumeration])) {
      std::ifstream setup_file(setup_file_name);
      std::string next;

      while (setup_file >> next) {
        if (next[0] == '#')
          std::getline(setup_file, next);

        else
          fie_file_names.push_back(next);
      }

      setup_file.close();
    }

    const int column_width(20);

    if (input == combi_ff::cnv::not_set) {
      std::cerr << "# no input arguments -I found, choosing smi by default\n";
      input = combi_ff::cnv::smiles;
    }

    std::list<std::string> input_list_no_file_names;

    for (auto && arg : input_list)
      input_list_no_file_names.push_back(arg.second);

    if (input == combi_ff::cnv::smiles)
      combi_ff::cnv::SmilesHandler(output_file,
                                   fie_file_names,
                                   column_width,
                                   print_options,
                                   input_list).Run();

    else if (input == combi_ff::cnv::formula)
      combi_ff::cnv::FormulaHandler(output_file,
                                    fie_file_names,
                                    column_width,
                                    print_options,
                                    input_list).Run();

    else if (input == combi_ff::cnv::name)
      combi_ff::cnv::NameHandler(output_file,
                                 fie_file_names,
                                 column_width,
                                 print_options,
                                 input_list).Run();

    else if (input == combi_ff::cnv::family_enumeration)
      combi_ff::cnv::FamilyIdentifierHandler(output_file,
                                             fie_file_names,
                                             column_width,
                                             print_options,
                                             input_list).Run();

    else if (input == combi_ff::cnv::matrix)
      combi_ff::cnv::AdjacencyMatrixHandler(output_file,
                                            fie_file_names,
                                            column_width,
                                            print_options,
                                            input_list).Run();

    else
      throw combi_ff::input_error("unknown input " + input);

    output_file.close();
  }

//catch incorrect user input
  catch (combi_ff::input_error& e) {
    std::cerr << "!Input Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  //case that './cnv help' was used, this is not a failure, but the program should stop
  catch (combi_ff::help_exception& e) {
    return EXIT_SUCCESS;
  }

  //catch any other errors
  catch (std::exception& e) {
    std::cerr << "!Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
















