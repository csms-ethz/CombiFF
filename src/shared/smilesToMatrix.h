// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "Matrix.h"

namespace combi_ff {

struct SmilesBlock;

combi_ff::AdjacencyMatrix<double, combi_ff::CnvAtom> ConvertToMatrix(
    const std::string& smiles);

bool GetBasicInformationFromSmiles(
    const std::string& smiles, combi_ff::AtomVector<combi_ff::CnvAtom>& atoms,
    size_t& num_atoms, size_t& num_rings, std::vector<bool>& is_aromatic);

void CreateNameVector(std::vector<combi_ff::SmilesBlock>& names,
                      const std::string& smiles,
                      combi_ff::AtomVector<combi_ff::CnvAtom>& atoms,
                      size_t& num_atoms);

void CreateMatrixFromSmilesBlocks(
    combi_ff::AdjacencyMatrix<double, combi_ff::CnvAtom>& A,
    const std::string& smiles, const std::vector<combi_ff::SmilesBlock>& names,
    const size_t num_rings);

}  // namespace combi_ff
