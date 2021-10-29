#ifndef SMILESHANDLER_H_
#define SMILESHANDLER_H_

#include "Handler.h"
#include "AdjacencyMatrixCnv.h"

namespace combi_ff {

namespace cnv {

class SmilesHandler : public cnv::Handler {
 public:
  SmilesHandler(std::ofstream& output_file,
                combi_ff::StringVector& fie_file_names,
                const int column_width,
                std::vector<bool>& print_options,
                std::list<std::pair<std::string, std::string>>& input_list)
    : cnv::Handler(output_file, fie_file_names, column_width, print_options, input_list),
      max_size(20) {}


  void Run();
  void PrintFirstLine() ;
  void ConvertSmiles(const std::string& smiles_orig,
                     std::string& smiles_canon,
                     cnv::AdjacencyMatrix& A);
  void PrintOutput(const std::string& smiles_orig,
                   const std::string& smiles_canon,
                   const std::string& fmi,
                   const cnv::AdjacencyMatrix& A) ;
  void FindFamilyIdentifier(const std::string& smiles_canon, std::string& fmi);

 private:
  size_t max_size;

};

} //namespace cnv

} //namespace combi_ff

#endif