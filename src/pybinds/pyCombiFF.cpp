//
// Created by benja on 02.11.2021.
//

#ifndef PYCOMBIFF_H
#define PYCOMBIFF_H
#include <pybind11/pybind11.h>

#include "./pycnv.h"



namespace py = pybind11;

PYBIND11_MODULE(pycombiff, m) {
    m.doc() = "pyCombiFF FUn text";
    m.attr("version");

    m.def("cnv",
          &pycnv,
          py::return_value_policy::automatic,
          "This function performs cnv.\n"
          "\n"
          "Parameters\n"
          "-----------\n"
          "data: Iterable[Iterable[Number]]\n"
          "\tis the data to be clustered.\n"
          "\n"
          "RETURNS\n"
          "-------\n"
          "bool:\n"
          "\tWas successful?\n"
          "\n",
          py::arg("inputFile"),
          py::arg("outputFile"),
          py::arg("inputFormat")="smi",
          py::arg("desiredOutput")="smi");
}

#endif //PYCOMBIFF_H