#include "NameHandler.h"
#include <iomanip>
#include <algorithm>

namespace combi_ff {

namespace cnv {

void NameHandler::Run() {
  cnv::Handler::Run(print_canon_name);
  NameHandler::PrintFirstLine();

  //determine the different names from the input_list
  for (auto name_it = input_list.begin(); name_it != input_list.end(); ++name_it) {
    for (size_t i = 0; i < name_it->second.size(); i++) {
      if ((name_it->second)[i] == ' ' || (name_it->second)[i] == ',') {
        auto next = ++name_it;
        name_it--;
        input_list.insert(next, {"", name_it->second.substr(i + 1, name_it->second.size() - i)});
        name_it->second = name_it->second.substr(0, i);
      }
    }
  }

  for (auto && name : input_list) {
    if (name.second.size())
      NameHandler::ConvertName(name.second);
  }
}

void NameHandler::ConvertName(const std::string& name) {
  std::string name_canon(name);
  std::transform(name_canon.begin(), name_canon.end(), name_canon.begin(), tolower);

  for (size_t i = 0; i < name_canon.size(); i++) {
    if (name_canon[i] == ' ') {
      name_canon.erase(name_canon.begin() + i);
      i--;
    }
  }

  NameHandler::PrintOutput(name, name_canon);
}

void NameHandler::PrintOutput(const std::string& name_orig, const std::string& name_canon) {
  std::ostream* out;

  if (output_file.is_open())
    out = &output_file;

  else
    out = &std::cout;

  *out << std::setw(column_width) << std::left << name_orig << " ";
  *out << std::setw(column_width) << std::left << name_canon << " ";
  *out << '\n';
}


void NameHandler::PrintFirstLine() {
  Handler::PrintFirstLine("# originalName ");
}

} //namespace cnv

} //namespace combi_ff