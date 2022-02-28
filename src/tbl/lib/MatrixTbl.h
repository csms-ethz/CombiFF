// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef AdjacencyMatrixTBL_H_
#define AdjacencyMatrixTBL_H_

#include "Matrix.h"

namespace combi_ff {
namespace topology_builder {

typedef combi_ff::AdjacencyMatrix<double, combi_ff::CnvAtom> AdjacencyMatrix;
typedef combi_ff::SymmetricalMatrix<double> FragmentMatrixTbl;

}  // namespace topology_builder

}  // namespace combi_ff

#endif