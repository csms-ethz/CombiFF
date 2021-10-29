#include "TblRun.h"
#include "exceptions.h"
#include <stdio.h>

int main(int argc, char* argv[]) {
  try {
    //add command line arguments to arguments vector
    if (argc > 1) {
      std::vector<std::string> arguments(argv + 1, argv + argc);
      //run the topology builder
      combi_ff::topology_builder::TblRun(arguments).run();

    } else
      std::cout << "no command line arguments found\n";
  }

  //catch incorrect user input
  catch (combi_ff::input_error& e) {
    std::cout << "!Input Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  //case that './tbl help' was used, this is not a failure, but the program should stop
  catch (combi_ff::help_exception& e) {
    return EXIT_SUCCESS;
  }

  //catch any other errors
  catch (std::exception& e) {
    std::cout << "!Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  //success
  return EXIT_SUCCESS;
}


