#include "EnumRun.h"

#include "ContainerOperators.h"
#include "Enumerator.h"
#include "FamilySpecifications.h"
#include "InputOutput.h"
#include "Time.h"
#include "printInfo.h"

namespace combi_ff {

namespace enu {

// constructor
EnumRun::EnumRun(const StringVector& arguments) : arguments(arguments) {}

// driver function that runs the whole programs
void EnumRun::Run() {
  // initialize random seed for temporary input file names
  srand((uint)time(NULL));
  // declare parameters and vectors used to store the input information and
  // enumSpecifications, and open input/output files
  InputOutput io(arguments);
  // create FamilySpecifications
  FamilySpecifications family_spec(io);
  // construct the molecules corresponding to input given directly on the
  // command line (if any direct input was given)
  EnumerateDirect(io);
  // construct the molecules corresponding to the required families(if any
  // familes are given)
  EnumerateFamilies(family_spec, io);
}

void EnumRun::EnumerateDirect(const InputOutput& io) {
  if (io.GetEnumSpec().used_atom_vectors.size()) {
    // print initial information
    std::cout << "Restrictions:\n";
    PrintRestrictions("  ", io.GetEnumSpec().max_degree,
                      io.GetEnumSpec().ranged_properties);
    // start timer
    Time start(Clock::now());
    // create Enumerator instance and use it to enumerate the isomers of the
    // given family
    Enumerator enumerator(io.GetEnumSpec(), io.GetIOFilProps().file_name_tmp);
    enumerator.EnumerateIsomers();
    // stop timer
    const std::chrono::milliseconds duration(GetDuration(start, Clock::now()));
    const std::string output_file_name =
        io.GetOutputFileDirectory() + io.GetOutputFileName();
    // print final information and close output file;
    WriteOutputDirect(io.GetIOFilProps(), output_file_name,
                      io.GetEnumSpec().stereo, io.GetEnumSpec().count_only,
                      duration, enumerator.GetNumIsomers());
  }
}

void EnumRun::EnumerateFamilies(const FamilySpecifications& family_spec,
                                const InputOutput& io) {
  for (auto&& family : family_spec.GetFamilies()) {
    // start timer
    Time start(Clock::now());
    // create Enumerator instance and use it to enumerate the isomers of the
    // given family
    Enumerator enumerator(family, family_spec.GetPseudoatoms(),
                          io.GetEnumSpec(), io.GetIOFilProps().file_name_tmp);
    enumerator.EnumerateIsomers();
    // stop timer
    const std::chrono::milliseconds duration(GetDuration(start, Clock::now()));
    // print final information and close output file
    WriteOutputFamily(io.GetIOFilProps(), io.GetEnumSpec().stereo,
                      io.GetEnumSpec().count_only, duration, family,
                      enumerator.GetNumIsomers());
  }
}

/*********************
WRITE THE FINAL OUTPUT
*********************/
void EnumRun::WriteOutput(const IOFileProperties& io_files,
                          const std::string& output_file_name,
                          const bool stereo, const bool count_only,
                          const std::chrono::milliseconds time,
                          const size_t num_isomers) {
  PrintSummary(io_files, output_file_name, num_isomers, time);

  if (!count_only) {
    std::ifstream tmp_file(io_files.file_name_tmp);

    if (!tmp_file.is_open())
      throw std::runtime_error(io_files.file_name_tmp + " not open\n");

    /*outputFile << "# " << numIsomers << (stereo ? " stereo" : " constitutional
    ") << "isomers\n"; outputFile << "# enumeration took " << time << '\n';
    outputFile <<
    "#==========================================================================\n"
               << "#" << '\n';*/
    output_file << xml_indent << "<enumeration_type type=\""
                << (stereo ? "stereo" : "constitutional") << "\"/>\n"
                << xml_indent << "<number_of_isomers>" << num_isomers
                << "</number_of_isomers>\n"
                << xml_indent << "<enumeration_time>" << time
                << "</enumeration_time>\n";

    if (num_isomers) {
      output_file << xml_indent << "<isomer_lists>\n"
                  << tmp_file.rdbuf()  // append tempFile to outputFile
                  << xml_indent << "</isomer_lists>\n";

    } else
      output_file << xml_indent << "<isomer_lists/>\n";

    tmp_file.close();
  }
}

/************************************************
WRITE THE FINAL OUTPUT FOR THE DIRECT ENUMERATION
*************************************************/
void EnumRun::WriteOutputDirect(const IOFileProperties& io_files,
                                const std::string& output_file_name,
                                const bool stereo, const bool count_only,
                                const std::chrono::milliseconds time,
                                const size_t num_isomers) {
  if (!count_only) {
    output_file.open(output_file_name);

    if (!output_file.is_open())
      throw std::runtime_error(output_file_name + " is not open\n");

    /*outputFile  <<
       "#==========================================================================\n"
                << "# isomer enumeration file\n"
                << "# created by ENU (version " << combi_ff::current_version <<
       ")\n";*/
    output_file
        << "<?xml version='1.0' encoding='UTF-8'?>\n"
        << "<!DOCTYPE isomer_enumeration SYSTEM \"isomer_enumeration.dtd\">\n"
        << "<isomer_enumeration enu_version=\"" << combi_ff::current_version
        << "\">\n";
    WriteOutput(io_files, output_file_name, stereo, count_only, time,
                num_isomers);
    output_file << "</isomer_enumeration>";
    output_file.close();
  } else {
    const std::string output_file_name_ = "";
    PrintSummary(io_files, output_file_name_, num_isomers, time);
  }
}

/************************************************
WRITE THE FINAL OUTPUT FOR THE FAMILY ENUMERATION
*************************************************/
void EnumRun::WriteOutputFamily(const IOFileProperties& io_files,
                                const bool stereo, const bool count_only,
                                const std::chrono::milliseconds time,
                                const Family& family,
                                const size_t num_isomers) {
  if (!count_only) {
    const std::string output_file_name(
        io_files.output_dir + "family_isomer_enumeration_" + family.GetCode() +
        (stereo ? "_stereo" : "") + ".xml");
    output_file.open(output_file_name);

    if (!output_file.is_open())
      throw std::runtime_error(output_file_name + " is not open\n");

    /*outputFile  <<
       "#==========================================================================\n"
                << "# isomer enumeration file for family " << family.GetCode()
       << " (version " << family.GetVersion() << ")\n"
                << "# created by ENU (version " << combi_ff::current_version <<
       ")\n";*/
    output_file << "<?xml version='1.0' encoding='UTF-8'?>\n"
                << "<!DOCTYPE family_isomer_enumeration SYSTEM "
                   "\"family_isomer_enumeration.dtd\">\n"
                << "<family_isomer_enumeration enu_version=\""
                << combi_ff::current_version << "\" family_code=\""
                << family.GetCode() << "\" family_version=\""
                << family.GetVersion() << "\">\n";
    WriteOutput(io_files, output_file_name, stereo, count_only, time,
                num_isomers);
    output_file << "</family_isomer_enumeration>";
    output_file.close();
  } else {
    const std::string output_file_name_ = "";
    PrintSummary(io_files, output_file_name_, num_isomers, time);
  }
}

/***************************
WRITE SUMMARY OF ENUMERATION
***************************/
void EnumRun::PrintSummary(const IOFileProperties& io_files,
                           const std::string& output_file_name,
                           const size_t num_isomers,
                           const std::chrono::milliseconds time) const {
  std::cout << "\n\n"
            << "********************************************************\n"
            << "successfully finished the current enumeration.\n"
            << "********************************************************\n"
            << "Summary:\n"
            << "********\n"
            << "used family Set file(s):       "
            << io_files.input_file_names[family_file] << "\n"
            << "used substructure Set file(s): "
            << io_files.input_file_names[substructure_file] << "\n"
            << "used alias Set file(s):        "
            << io_files.input_file_names[alias_file] << "\n"
            << "used pseudoatom Set file(s):   "
            << io_files.input_file_names[pseudoatom_file] << "\n";

  if (output_file_name.size())
    std::cout << "wrote output to:               " << output_file_name << "\n";

  std::cout << "********************************************************\n"
            << "total number of isomers:       " << num_isomers << "\n"
            << "********************************************************\n"
            << "enumeration took:              " << time << '\n'
            << "********************************************************\n";
}

}  // namespace enu

}  // namespace combi_ff