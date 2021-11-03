//
// Created by benja on 02.11.2021.
//

#include "../cnv/main.cpp"

#include "pycnv.h"


//pybindsWrapper
void pycnv(char inputFile, char outputFile, char inputFormat, char desiredOutput){
    std::vector<char> argv;
    argv.push_back('-I');
    argv.push_back(inputFormat);
    argv.push_back('-O');
    argv.push_back(desiredOutput);
    argv.push_back('-i');
    argv.push_back(inputFile);
    argv.push_back('-o');
    argv.push_back(outputFile);
    int argc = argv.size();
    char* cargv = &argv[0];

    main(argc, &cargv);
}


