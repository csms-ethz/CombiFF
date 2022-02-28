// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef ENUMRUN_H_
#define ENUMRUN_H_

#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include "StringVector.h"
#include "version.h"

namespace combi_ff {

namespace enu {

class FamilySpecifications;
class InputOutput;
class IOFileProperties;
class Family;

class EnumRun {
 public:
  EnumRun(const StringVector& arguments);

  void Run();

  void EnumerateDirect(const InputOutput& io);

  void EnumerateFamilies(const FamilySpecifications& family_spec,
                         const InputOutput& io);

  void WriteOutputDirect(const IOFileProperties& io_files,
                         const std::string& output_file_name, const bool stereo,
                         const std::chrono::milliseconds time,
                         const size_t num_isomers);

  void WriteOutput(const IOFileProperties& io_files,
                   const std::string& output_file_name, const bool stereo,
                   const std::chrono::milliseconds time,
                   const size_t num_isomers);

  void WriteOutputFamily(const IOFileProperties& io_files, const bool stereo,
                         const std::chrono::milliseconds time,
                         const Family& family, const size_t num_isomers);

  void PrintSummary(const IOFileProperties& io_files,
                    const std::string& output_file_name, const size_t numIso,
                    const std::chrono::milliseconds time) const;

 private:
  const StringVector& arguments;
  std::ofstream output_file;
};

}  // namespace enu

}  // namespace combi_ff

#endif
