// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "TblRun.h"

#include "TblFragment.h"
#include "exceptions.h"
#include "moleculeMacros.h"
#include "smilesToMatrix.h"
//#include <filesystem>
#include "MoleculeDecomposer.h"
#include "mtb.h"
#include "replacement.h"
#include "version.h"

namespace combi_ff {

namespace topology_builder {

TblRun::TblRun(const std::vector<std::string>& arguments)
    : arguments(arguments) {}

void TblRun::run() {
  // open input and output files. read in all the various parameters
  topology_builder::InputOutput io(arguments);
  const IOFileProperties& io_file_properties = io.GetIOFilProps();
  // create the map of atom type Sets
  AtomTypeSetMap atom_type_set(0);
  CreateAtomTypeSets(io_file_properties.input_file_names[atomtypes_file],
                     atom_type_set);
  // create the fragments used for the molecule decomposition
  std::vector<TblFragment> tbl_fragments;
  CreateTblFragments(io_file_properties.input_file_names[fragment_file],
                     tbl_fragments, io.GetFragmentsToUse(), atom_type_set);

  for (auto&& filename_familyIsomerEnumeration :
       io_file_properties.fie_families) {
    std::cout
        << "\n\n"
        << "**********************************************************\n"
        << "family isomer enumeration file is "
        << filename_familyIsomerEnumeration << '\n'
        << "**********************************************************\n\n";
    FamilyDecomposer decomposer(filename_familyIsomerEnumeration,
                                io_file_properties, tbl_fragments);
    std::string filename_molecules_with_macros("");
    CreateMoleculesWithMacros(
        decomposer.GetFilename(), filename_molecules_with_macros,
        decomposer.GetFamilyCode(), tbl_fragments, io_file_properties);
    CreateMTB(filename_molecules_with_macros, decomposer.GetFamilyCode(),
              io_file_properties);

    ReplaceMacrosByParameters(filename_molecules_with_macros,
                              decomposer.GetFamilyCode(), io_file_properties);
  }
}

}  // namespace topology_builder

}  // namespace combi_ff