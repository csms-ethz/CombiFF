#include "mtb.h"
#include "version.h"
#include "TopologicalProperty.h"
#include "InputOutput.h"

namespace combi_ff {

namespace topology_builder {

void CreateMTB(const std::string& filename_molecules_with_macros,
               const std::string& family_code,
               const IOFileProperties& io_file_properties) {
  std::cout << "\ncreating mtb\n\n";
  size_t counter(0);
  XmlParserIn molecules_with_macros_in(filename_molecules_with_macros,
                                       XmlParserIn::read_element_by_element);
  const std::string filename_mtb((io_file_properties.output_dir_mtb.size() ?
                                  io_file_properties.output_dir_mtb : io_file_properties.output_dir)
                                 + family_code + ".mtb");
  molecules_with_macros_in.GetRoot().CheckAttribute("family_code");
  molecules_with_macros_in.GetRoot().CheckAttribute("family_version");
  molecules_with_macros_in.GetRoot().CheckAttribute("enu_version");
  molecules_with_macros_in.GetRoot().CheckAttribute("tbl_version");
  std::ofstream mtb_file(filename_mtb);
  mtb_file << "TITLE\n"
          << "topology file for family " << molecules_with_macros_in.GetRoot().attributes["family_code"]
          << " (version " << molecules_with_macros_in.GetRoot().attributes["family_version"] << ")"
          << " created by enu version " << molecules_with_macros_in.GetRoot().attributes["enu_version"]
          << " and tbl version " << molecules_with_macros_in.GetRoot().attributes["tbl_version"]
          << ". This file was created by tbl version " << combi_ff::current_version << "\n"
          << "END\n"
          << "FORCEFIELD\n"
          << combi_ff::topology_builder::GROMOS_force_field << '\n'
          << "END\n"
          << "MAKETOPVERSION\n"
          << combi_ff::topology_builder::maketop_version << "\n"
          << "END\n";

  while (molecules_with_macros_in.ReadUntilElement("molecule_with_macros")) {
    auto&& molecule_with_macros = molecules_with_macros_in.GetRoot().GetLastChild();

    if (molecule_with_macros->children.size() < 1) {
      std::cout << "no decomposition found for this molecule\n";
      continue;
    }

    std::cout << "molecule nr. " << ++counter << ": "
              << molecule_with_macros->attributes["smiles"] << '\n';
    mtb_file   << "MTBUILDBLSOLUTE\n"
              << "# smiles: " << molecule_with_macros->attributes["smiles"] << "\n"
              << "# RNME\n"
              << molecule_with_macros->attributes["isomer_id"] << "\n"
              << "# number of atoms, number of preceding exclusion\n"
              << "# NMAT,NLIN\n";
    auto&& molecule_with_macros_property_it = molecule_with_macros->children.begin();
    auto&& atoms = *molecule_with_macros_property_it;
    atoms->CheckTagName("atoms");
    atoms->CheckAttribute("amount");
    mtb_file << std::setw(5)
            << std::right
            << atoms->attributes["amount"]
            << " "
            << std::setw(4)
            << std::right
            << 0
            << '\n'
            << "# atoms\n"
            << "#ATOM ANM  IACM MASS        CGM ICGM MAE MSAE\n";
    size_t i = 0;
    std::unordered_map<std::string, size_t> atoms_ANM_to_ATOM;

    for (auto && atom : atoms->children)
      atoms_ANM_to_ATOM[atom->attributes["atom_id"]] = ++i;

    for (auto && atom : atoms->children) {
      atom->CheckAttribute("atom_id");
      const std::string& atomtype = atom->attributes["atom_type"];
      mtb_file << std::setw(5)
              << std::right
              << atoms_ANM_to_ATOM[atom->attributes["atom_id"]]
              << " " << atom->attributes["atom_id"]
              << " typ_" << atomtype
              << " mass_" << atomtype
              << " " << "chg_" << atomtype
              << " " << 1;
      std::vector<size_t> exclusions{};

      if (atom->children.size()) {
        auto&& excluded_atoms = atom->GetFirstChild();

        for (auto && ea : excluded_atoms->children) {
          ea->CheckTagName("excluded_atom");

          if (atoms_ANM_to_ATOM[ea->value] > atoms_ANM_to_ATOM[atom->attributes["atom_id"]])
            exclusions.push_back(atoms_ANM_to_ATOM[ea->value]);
        }
      }

      mtb_file << " " << exclusions.size();

      for (auto && ea : exclusions)
        mtb_file << " " << ea;

      mtb_file << '\n';
    }

    molecule_with_macros_property_it++;

    for (auto prop : property_types_in_mtb) {
      if (molecule_with_macros_property_it == molecule_with_macros->children.end()) {
        mtb_file << "#" << std::setw(4) << std::right << parameter_count_comments[prop] << std::endl;
        mtb_file << "0\n";
        break;
      }

      auto&& cur_property = *molecule_with_macros_property_it;
      mtb_file << "#" << std::setw(4) << std::right << parameter_count_comments[prop] << std::endl;

      if (topologicalPropertyMap.find(cur_property->tag.substr(0,
                                                              cur_property->tag.size() - 1))->second == prop) {
        cur_property->CheckAttribute("amount");
        mtb_file << std::setw(5) << std::right << cur_property->attributes["amount"] << '\n';
        std::shared_ptr<TopologicalPropertyBase> new_property = CreateNewTopologicalProperty(prop);
        mtb_file << "#";
        mtb_file << std::setw(4) << std::right << parameter_comments[0];

        for (size_t i = 1; i < new_property->GetNumInvolvedAtoms(); i++)
          mtb_file << std::setw(5) << std::right << parameter_comments[i];

        mtb_file << std::setw(6) << std::right << "MCB\n";

        for (auto && property_child : cur_property->children) {
          property_child->CheckTagName(new_property->GetPropertyName());
          property_child->CheckNumberOfChildren_equal(1);
          property_child->GetFirstChild()->CheckNumberOfChildren_equal(new_property->GetNumInvolvedAtoms());

          for (auto && involved_atom : property_child->GetFirstChild()->children)
            mtb_file << " " << std::setw(4) << std::right << atoms_ANM_to_ATOM[involved_atom->value];

          mtb_file << "    " << property_child->attributes["parameter"] << "\n";
        }

        molecule_with_macros_property_it++;

      } else
        mtb_file << "    0\n";
    }

    mtb_file << "# NEX\n"
            << "0\n"
            << "END\n";
  }
}

} //namespace topology_builder

} //namespace combi_ff