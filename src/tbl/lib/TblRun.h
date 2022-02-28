// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef TBLRUN_H_
#define TBLRUN_H_

#include <ctype.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "MatrixTbl.h"

namespace combi_ff {

namespace topology_builder {

class TblRun {
 private:
  const std::vector<std::string>& arguments;

 public:
  TblRun(const std::vector<std::string>& arguments);

  void run();
};

}  // namespace topology_builder

}  // namespace combi_ff

#endif
