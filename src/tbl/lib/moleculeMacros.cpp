// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "moleculeMacros.h"

#include <set>

#include "InputOutput.h"
#include "Linkage.h"
#include "TblAtom.h"
#include "TblFragment.h"
#include "version.h"

namespace combi_ff {

namespace topology_builder {

/***************************
AtomInMolecule implemenation
***************************/

/*
constructor
*/
AtomInMolecule::AtomInMolecule(
    size_t idx, const std::string& par, const std::string& id,
    const std::vector<std::string>& connected_link_atoms)
    : idx(idx),
      neighbors(0),
      atom_id(id),
      parameter(par),
      connected_link_atoms(connected_link_atoms) {}

// use the XML file molecule_decompositions_
void CreateMoleculesWithMacros(
    const std::string& filename_molecule_decomposition,
    std::string& filename_molecules_with_macros, const std::string& family_code,
    const std::vector<TblFragment>& tbl_fragments,
    const IOFileProperties& io_file_properties) {
  std::cout << "\ncreating molecules with macros\n\n";
  XmlParserIn molecule_decompositions_in(filename_molecule_decomposition,
                                         XmlParserIn::read_element_by_element);
  filename_molecules_with_macros =
      ((io_file_properties.output_dir_molecules_with_macros.size()
            ? io_file_properties.output_dir_molecules_with_macros
            : io_file_properties.output_dir) +
       "molecules_with_macros_" + family_code + ".xml");
  size_t counter(0);
  XmlParserOut molecules_macros(filename_molecules_with_macros,
                                "molecules_with_macros");
  molecules_macros.GetRoot().attributes["tbl_version"] =
      combi_ff::current_version;
  molecule_decompositions_in.GetRoot().CheckAttribute("family_code");
  molecules_macros.GetRoot().attributes["tbl_version"] =
      combi_ff::current_version;
  molecules_macros.GetRoot().attributes["family_code"] =
      molecule_decompositions_in.GetRoot().attributes["family_code"];
  molecules_macros.GetRoot().attributes["family_version"] =
      molecule_decompositions_in.GetRoot().attributes["family_version"];
  molecules_macros.GetRoot().attributes["enu_version"] =
      molecule_decompositions_in.GetRoot().attributes["enu_version"];
  molecules_macros.WriteHead();

  // for(auto && molecule_decomposition :
  // molecule_decompositions_in.GetRoot().children) {
  while (
      molecule_decompositions_in.ReadUntilElement("molecule_decomposition")) {
    auto&& molecule_decomposition =
        molecule_decompositions_in.GetRoot().GetLastChild();
    std::cout << "molecule nr. " << ++counter << ": "
              << molecule_decomposition->attributes["smiles"] << '\n';

    if (molecule_decomposition->children.size() != 2) {
      std::cout << "no decomposition found for this molecule\n";
      continue;
    }

    molecules_macros.GetRoot().AddElement("molecule_with_macros");
    XmlElement& molecule_macros = *(molecules_macros.GetRoot().GetLastChild());
    CreateMoleculeWithMacros(molecule_decomposition, molecule_macros,
                             tbl_fragments, io_file_properties.united_atom,
                             io_file_properties.third_neighbor_exclusions,
                             io_file_properties.unique_torsionals);
    molecules_macros.WriteAndRemoveLastChild();
  }

  molecules_macros.WriteTail();
}

struct less {
  bool operator()(const std::pair<std::string, const AtomTypeSet*>& a,
                  const std::pair<std::string, const AtomTypeSet*>& b) const {
    return *(a.second) < *(b.second);
  }
};

void CreateMoleculeWithMacros(XmlElement_ptr molecule_decomposition,
                              XmlElement& molecule_macros,
                              const std::vector<TblFragment>& tbl_fragments,
                              const bool united_atom,
                              const bool third_neighbor_exclusions,
                              const bool unique_torsionals) {
  molecule_macros.SetAttribute(*molecule_decomposition, "smiles");
  molecule_macros.SetAttribute(*molecule_decomposition, "molecule_name");
  molecule_macros.SetAttribute(*molecule_decomposition, "isomer_id");
  // molecule_macros.AddElement("atoms");
  // auto&& atoms = molecule_macros.GetLastChild();
  molecule_decomposition->CheckNumberOfChildren_equal(2);
  auto&& fragments = molecule_decomposition->GetFirstChild();
  fragments->CheckTagName("fragments");
  std::unordered_map<std::string, const TblFragment*> used_fragments;
  // pair: 1st element is atom_id, 2nd element is atomtype
  std::unordered_map<std::string, std::pair<std::string, std::string>>
      core_atom_types;
  std::vector<std::vector<topo>> topological_properties(
      num_topological_property_types);
  std::vector<std::pair<std::string, const AtomTypeSet*>> atms(0);

  for (auto&& fragment : fragments->children) {
    fragment->CheckTagName("fragment");
    fragment->CheckAttribute("id");
    std::string fragment_id = fragment->attributes["id"];
    std::string fragCode = fragment->value;
    bool found = false;

    if (used_fragments.find(fragment_id) != used_fragments.end())
      throw std::runtime_error("fragment ID " + fragment_id +
                               " used more than once in molecule " +
                               molecule_macros.attributes["smiles"]);

    for (auto&& frag : tbl_fragments) {
      if (frag.GetCode() == fragCode) {
        used_fragments[fragment_id] = &frag;

        if (frag.GetVersion() != combi_ff::current_version)
          std::cout << "?Warning: currently running combi_ff version "
                    << combi_ff::current_version << " but fragment " << fragCode
                    << " has version " << frag.GetVersion() << "\n";

        for (auto&& TBatom : frag.GetTblAtoms()) {
          if (TBatom.GetLinkageType() == core) {
            atms.push_back({fragment_id + "_" + TBatom.GetAtomID(),
                            &(TBatom.GetAtomTypes())});
          }

          core_atom_types[fragment_id + "_" + TBatom.GetAtomID()] =
              std::pair<std::string, std::string>(
                  fragment_id + "_" + TBatom.GetAtomID(),
                  frag.GetCoreAtom(TBatom.GetAtomID())
                      ->GetAtomTypes()
                      .GetName());
        }

        for (auto&& prop : frag.GetTopologicalProperties()) {
          topo newtopo;
          newtopo.property = prop;
          newtopo.frag = &frag;
          newtopo.fragment_id = fragment_id;
          topological_properties
              [topologicalPropertyMap.find(prop->GetPropertyName())->second]
                  .push_back(newtopo);
        }

        found = true;
        break;
      }
    }

    if (!found)
      throw std::runtime_error("fragment " + fragment_id +
                               " not found in tbl_fragments\n");
  }

  auto&& linkages = molecule_decomposition->GetLastChild();
  linkages->CheckTagName("linkages");
  std::unordered_map<std::string, std::vector<std::string>> connected_fragments;

  for (auto&& linkage : linkages->children) {
    linkage->CheckTagName("linkage");
    linkage->CheckNumberOfChildren_equal(2);
    std::vector<std::string> involved_atoms;
    std::vector<const TblFragment*> frags;
    std::vector<const std::string*> fragment_ids;
    involved_atoms.reserve(2);
    frags.reserve(2);
    fragment_ids.reserve(2);

    for (auto&& involved_fragment : linkage->children) {
      involved_fragment->CheckTagName("involved_fragment");
      involved_fragment->GetFirstChild()->CheckTagName("fragment_id");
      involved_fragment->GetLastChild()->CheckTagName("linksite");
      involved_atoms.push_back(involved_fragment->GetLastChild()->value);
      frags.push_back(
          used_fragments[involved_fragment->GetFirstChild()->value]);
      fragment_ids.push_back(&involved_fragment->GetFirstChild()->value);
    }

    // pair: 1st element is atom_id, 2nd element is atomtype of core atom
    core_atom_types[*fragment_ids[0] + "_" + involved_atoms[0]].first =
        *fragment_ids[1] + "_" +
        frags[1]->GetCoreAtom(involved_atoms[1])->GetAtomID();
    core_atom_types[*fragment_ids[0] + "_" + involved_atoms[0]].second =
        frags[1]->GetCoreAtom(involved_atoms[1])->GetAtomTypes().GetName();
    core_atom_types[*fragment_ids[1] + "_" + involved_atoms[1]].first =
        *fragment_ids[0] + "_" +
        frags[0]->GetCoreAtom(involved_atoms[0])->GetAtomID();
    core_atom_types[*fragment_ids[1] + "_" + involved_atoms[1]].second =
        frags[0]->GetCoreAtom(involved_atoms[0])->GetAtomTypes().GetName();
    connected_fragments[*fragment_ids[0] + "_" +
                        frags[0]->GetCoreAtom(involved_atoms[0])->GetAtomID()]
        .push_back(frags[1]->GetCode() + "." + involved_atoms[1]);
    connected_fragments[*fragment_ids[1] + "_" +
                        frags[1]->GetCoreAtom(involved_atoms[1])->GetAtomID()]
        .push_back(frags[0]->GetCode() + "." + involved_atoms[0]);
  }

  std::cout << "connected fragments:\n";

  for (auto&& cc : connected_fragments) {
    // std::sort(cc.second.begin(), cc.second.end()); // don't sort, as
    // fragments arent necessarily symmetric
    std::cout << cc.first << " " << cc.second[0] << std::endl;
  }

  std::sort(atms.begin(), atms.end(), less());
  std::vector<AtomInMolecule> indexed_atoms;
  indexed_atoms.reserve(atms.size());
  std::unordered_map<std::string, AtomInMolecule> atoms_in_molecule;
  molecule_macros.AddElement("atoms");
  auto&& atoms = molecule_macros.GetLastChild();
  auto&& at = atms.front();
  atoms->AddElement("atom");
  size_t idx = 1;
  atoms->GetLastChild()->attributes["atom_id"] =
      (at.second->GetAtomSet().begin()->GetElementSymbol()) +
      std::to_string(idx);
  atoms->GetLastChild()->attributes["atom_type"] = at.second->GetName();

  if (connected_fragments[at.first].size()) {
    atoms->GetLastChild()->attributes["atom_type"] += "/";

    for (auto&& tt : connected_fragments[at.first])
      atoms->GetLastChild()->attributes["atom_type"] += tt + ",";

    atoms->GetLastChild()->attributes["atom_type"].pop_back();
  }

  atoms_in_molecule.emplace(
      at.first,
      AtomInMolecule(0, at.second->GetName(),
                     (at.second->GetAtomSet().begin()->GetElementSymbol()) +
                         std::to_string(idx),
                     connected_fragments[at.first]));
  indexed_atoms.push_back(
      AtomInMolecule(0, at.second->GetName(),
                     (at.second->GetAtomSet().begin()->GetElementSymbol()) +
                         std::to_string(idx),
                     connected_fragments[at.first]));
  indexed_atoms.back().atomtype = core_atom_types[at.first].second;
  idx++;

  for (size_t i = 1; i < atms.size(); i++) {
    auto&& at = atms[i];
    atoms->AddElement("atom");
    // if (atms[i - 1].second->GetAtomSet().begin()->GetElementSymbol() !=
    // at.second->GetAtomSet().begin()->GetElementSymbol())
    //   idx = 1;
    atoms->GetLastChild()->attributes["atom_id"] =
        (at.second->GetAtomSet().begin()->GetElementSymbol()) +
        std::to_string(idx);
    atoms->GetLastChild()->attributes["atom_type"] = at.second->GetName();

    if (connected_fragments[at.first].size()) {
      atoms->GetLastChild()->attributes["atom_type"] += "/";

      for (auto&& tt : connected_fragments[at.first])
        atoms->GetLastChild()->attributes["atom_type"] += tt + ",";

      atoms->GetLastChild()->attributes["atom_type"].pop_back();
    }

    atoms_in_molecule.emplace(
        at.first,
        AtomInMolecule(i, at.second->GetName(),
                       (at.second->GetAtomSet().begin()->GetElementSymbol()) +
                           std::to_string(idx),
                       connected_fragments[at.first]));
    indexed_atoms.push_back(
        AtomInMolecule(i, at.second->GetName(),
                       (at.second->GetAtomSet().begin()->GetElementSymbol()) +
                           std::to_string(idx),
                       connected_fragments[at.first]));
    indexed_atoms.back().atomtype = core_atom_types[at.first].second;
    idx++;
  }

  atoms->attributes["amount"] = std::to_string(atoms->children.size());
  std::vector<std::list<std::shared_ptr<TopologicalPropertyBase>>>
      topological_properties_in_molecule(num_topological_property_types);

  for (auto&& prop_type : topological_properties) {
    if (prop_type.size()) {
      PropertyType type =
          topologicalPropertyMap
              .find(prop_type.front().property->GetPropertyName())
              ->second;

      for (auto&& prop : prop_type) {
        topological_properties_in_molecule[type].push_back(
            CreateNewTopologicalProperty(type));
        auto&& cur_property = topological_properties_in_molecule[type].back();

        if (type == bond)
          cur_property->AddAttribute("degree",
                                     prop.property->GetAttribute("degree"));

        for (auto&& idx : prop.property->GetInvolvedAtoms()) {
          auto&& atom = prop.frag->GetTblAtom(idx);
          cur_property->AddInvolvedAtom(
              atoms_in_molecule[core_atom_types[prop.fragment_id + "_" +
                                                atom.GetAtomID()]
                                    .first]
                  .idx);
        }

        cur_property->SortInvolvedAtoms();
        std::vector<std::string> involved_atoms_codes(0);

        for (auto&& idx : cur_property->GetInvolvedAtoms()) {
          for (auto&& idx_ : prop.property->GetInvolvedAtoms()) {
            auto&& atom = prop.frag->GetTblAtom(idx_);

            if (atoms_in_molecule[core_atom_types[prop.fragment_id + "_" +
                                                  atom.GetAtomID()]
                                      .first]
                    .idx == idx) {
              involved_atoms_codes.push_back(
                  core_atom_types[prop.fragment_id + "_" + atom.GetAtomID()]
                      .first);
              break;
            }
          }
        }

        std::string parameter_name =
            (prop.property->GetAttribute("parameter").size()
                 ? prop.property->GetAttribute("parameter") + "_"
                 : "");
        parameter_name += GetParameterNameFromFragment(
            *prop.frag, prop.fragment_id, *prop.property, core_atom_types,
            involved_atoms_codes, connected_fragments);
        cur_property->AddAttribute("parameter", parameter_name);
      }

      // for (auto && pp : topological_properties_in_molecule[type])
      //   pp->SortInvolvedAtoms();

      if (topological_properties_in_molecule[type].size() > 1) {
        topological_properties_in_molecule[type].sort();
        auto&& it = topological_properties_in_molecule[type].begin();
        auto it_next = it;
        it_next++;

        for (; it_next != topological_properties_in_molecule[type].end();
             it++) {
          if ((*it)->GetInvolvedAtoms() == (*it_next)->GetInvolvedAtoms()) {
            topological_properties_in_molecule[type].erase(it_next);
            it_next = it;
            it--;
          }

          it_next++;
        }
      }
    }
  }

  // determine dihedrals
  topology_builder::AdjacencyMatrix A(atoms_in_molecule.size());

  for (auto&& bd : topological_properties_in_molecule[bond]) {
    double deg;
    const std::string& bond_degree = bd->GetAttribute("degree");

    if (bond_degree == "single")
      deg = 1;

    else if (bond_degree == "double")
      deg = 2;

    else if (bond_degree == "triple")
      deg = 3;

    else if (bond_degree == "quadruple")
      deg = 4;

    else if (bond_degree == "aromatic")
      deg = 1.5;

    else
      throw std::runtime_error("unknown bond degree " + bond_degree + "\n");

    A.SetElement(bd->GetInvolvedAtoms()[0], bd->GetInvolvedAtoms()[1], deg);
    indexed_atoms[bd->GetInvolvedAtoms()[0]].neighbors.push_back(
        bd->GetInvolvedAtoms()[1]);
    indexed_atoms[bd->GetInvolvedAtoms()[1]].neighbors.push_back(
        bd->GetInvolvedAtoms()[0]);
    // auto&& involved_atoms = bd->GetFirstChild();
    // auto&& involved_atom1 = involved_atoms->GetFirstChild()->value;
    // auto&& involved_atom2 = involved_atoms->GetLastChild()->value;
    // size_t deg;
    // if(bd->attributes["par"]);
    // A.SetElement(AtomInMolecule[involved_atom1].idx, AtomInMolecule);
  }

  /*for (auto && bd : topological_properties_in_molecule[bond]) {
    topological_properties_in_molecule[first_neighbor].push_back(CreateNewTopologicalProperty(first_neighbor));
    auto&& cur_first_neighbor =
  topological_properties_in_molecule[first_neighbor].back();
    cur_first_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[0]].idx);
    cur_first_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[1]].idx);
    cur_first_neighbor->SortInvolvedAtoms();
    std::string parameter_name = GetParameterNameAutomatically(A,
  cur_first_neighbor, indexed_atoms, core_atom_types);
    cur_first_neighbor->AddAttribute("parameter", parameter_name);
  }*/

  for (auto&& bd : topological_properties_in_molecule[bond]) {
    if (indexed_atoms[bd->GetInvolvedAtoms()[0]].neighbors.size() > 1) {
      // angle(s) found
      for (auto&& nbr_idx :
           indexed_atoms[bd->GetInvolvedAtoms()[0]].neighbors) {
        if (nbr_idx != indexed_atoms[bd->GetInvolvedAtoms()[1]].idx) {
          topological_properties_in_molecule[angle].push_back(
              CreateNewTopologicalProperty(angle));
          auto&& cur_angle = topological_properties_in_molecule[angle].back();
          cur_angle->AddInvolvedAtom(
              indexed_atoms[bd->GetInvolvedAtoms()[1]].idx);
          cur_angle->AddInvolvedAtom(
              indexed_atoms[bd->GetInvolvedAtoms()[0]].idx);
          cur_angle->AddInvolvedAtom(nbr_idx);
          cur_angle->SortInvolvedAtoms();
          std::string parameter_name = GetParameterNameAutomatically(
              A, cur_angle, indexed_atoms, core_atom_types);
          cur_angle->AddAttribute("parameter", parameter_name);
          /*
          topological_properties_in_molecule[second_neighbor].push_back(CreateNewTopologicalProperty(second_neighbor));
          auto&& cur_second_neighbor =
          topological_properties_in_molecule[second_neighbor].back();
          cur_second_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[1]].idx);
          cur_second_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[0]].idx);
          cur_second_neighbor->AddInvolvedAtom(nbr_idx);
          cur_second_neighbor->SortInvolvedAtoms();
          parameter_name = GetParameterNameAutomatically(A, cur_second_neighbor,
          indexed_atoms, core_atom_types);
          cur_second_neighbor->AddAttribute("parameter", parameter_name);
          */
        }
      }
    }

    if (indexed_atoms[bd->GetInvolvedAtoms()[1]].neighbors.size() > 1) {
      // angle(s) found
      for (auto&& nbr_idx :
           indexed_atoms[bd->GetInvolvedAtoms()[1]].neighbors) {
        if (nbr_idx != indexed_atoms[bd->GetInvolvedAtoms()[0]].idx) {
          topological_properties_in_molecule[angle].push_back(
              CreateNewTopologicalProperty(angle));
          auto&& cur_angle = topological_properties_in_molecule[angle].back();
          cur_angle->AddInvolvedAtom(
              indexed_atoms[bd->GetInvolvedAtoms()[0]].idx);
          cur_angle->AddInvolvedAtom(
              indexed_atoms[bd->GetInvolvedAtoms()[1]].idx);
          cur_angle->AddInvolvedAtom(nbr_idx);
          cur_angle->SortInvolvedAtoms();
          std::string parameter_name = GetParameterNameAutomatically(
              A, cur_angle, indexed_atoms, core_atom_types);
          cur_angle->AddAttribute("parameter", parameter_name);
          /*
          topological_properties_in_molecule[second_neighbor].push_back(CreateNewTopologicalProperty(second_neighbor));
          auto&& cur_second_neighbor =
          topological_properties_in_molecule[second_neighbor].back();
          cur_second_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[0]].idx);
          cur_second_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[1]].idx);
          cur_second_neighbor->AddInvolvedAtom(nbr_idx);
          cur_second_neighbor->SortInvolvedAtoms();
          parameter_name = GetParameterNameAutomatically(A, cur_second_neighbor,
          indexed_atoms, core_atom_types);
          cur_second_neighbor->AddAttribute("parameter", parameter_name);
          */
        }
      }
    }

    if (indexed_atoms[bd->GetInvolvedAtoms()[0]].neighbors.size() > 1 &&
        indexed_atoms[bd->GetInvolvedAtoms()[1]].neighbors.size() > 1) {
      // torsional found
      // topological_properties_in_molecule[torsional_dihedral].push_back(CreateNewTopologicalProperty(torsional_dihedral));
      // topological_properties_in_molecule[third_neighbor].push_back(CreateNewTopologicalProperty(third_neighbor));
      // auto&& cur_torsional_dihedral =
      // topological_properties_in_molecule[torsional_dihedral].back();
      //  auto&& cur_third_neighbor =
      //  topological_properties_in_molecule[third_neighbor].back();
      for (const auto nbr1 :
           indexed_atoms[bd->GetInvolvedAtoms()[0]].neighbors) {
        for (const auto nbr2 :
             indexed_atoms[bd->GetInvolvedAtoms()[1]].neighbors) {
          if (nbr1 != indexed_atoms[bd->GetInvolvedAtoms()[1]].idx &&
              nbr2 != indexed_atoms[bd->GetInvolvedAtoms()[0]].idx) {
            topological_properties_in_molecule[torsional_dihedral].push_back(
                CreateNewTopologicalProperty(torsional_dihedral));
            auto&& cur_torsional_dihedral =
                topological_properties_in_molecule[torsional_dihedral].back();
            cur_torsional_dihedral->AddInvolvedAtom(nbr1);
            cur_torsional_dihedral->AddInvolvedAtom(
                indexed_atoms[bd->GetInvolvedAtoms()[0]].idx);
            cur_torsional_dihedral->AddInvolvedAtom(
                indexed_atoms[bd->GetInvolvedAtoms()[1]].idx);
            cur_torsional_dihedral->AddInvolvedAtom(nbr2);
            cur_torsional_dihedral->SortInvolvedAtoms();
            // cur_third_neighbor->SortInvolvedAtoms();
            std::string parameter_name = GetParameterNameAutomatically(
                A, cur_torsional_dihedral, indexed_atoms, core_atom_types);
            cur_torsional_dihedral->AddAttribute("parameter", parameter_name);
          }
        }
      }

      /*
      if (indexed_atoms[bd->GetInvolvedAtoms()[0]].neighbors[0] ==
      indexed_atoms[bd->GetInvolvedAtoms()[1]].idx){
        cur_torsional_dihedral->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[0]].neighbors[1]);
        //cur_third_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[0]].neighbors[1]);
      }

      else{
        cur_torsional_dihedral->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[0]].neighbors[0]);
        //cur_third_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[0]].neighbors[0]);
      }

      cur_torsional_dihedral->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[0]].idx);
      cur_torsional_dihedral->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[1]].idx);
      //cur_third_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[0]].idx);
      //cur_third_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[1]].idx);

      if (indexed_atoms[bd->GetInvolvedAtoms()[1]].neighbors[0] ==
      indexed_atoms[bd->GetInvolvedAtoms()[0]].idx){
        cur_torsional_dihedral->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[1]].neighbors[1]);
       //
      cur_third_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[1]].neighbors[1]);
      }

      else{
        cur_torsional_dihedral->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[1]].neighbors[0]);
       //
      cur_third_neighbor->AddInvolvedAtom(indexed_atoms[bd->GetInvolvedAtoms()[1]].neighbors[0]);
      }*/
      // parameter_name = GetParameterNameAutomatically(A, cur_third_neighbor,
      // indexed_atoms, core_atom_types);
      //  cur_third_neighbor->AddAttribute("parameter", parameter_name);
    }
  }

  /*
  remove bonds with united hydrogens
  */
  if (united_atom) {
    auto&& it = topological_properties_in_molecule[bond].begin();

    for (; it != topological_properties_in_molecule[bond].end(); it++) {
      for (auto&& at : (*it)->GetInvolvedAtoms()) {
        if (indexed_atoms[at].atomtype.find("united") != std::string::npos) {
          auto it_erase = it;
          it--;
          topological_properties_in_molecule[bond].erase(it_erase);
          break;
        }
      }

      if (!topological_properties_in_molecule[bond].size()) break;
    }
  }

  /*
  remove duplicate topological properties (e.g. angle a-b-c = c-b-a, after
  sorting)
  */
  for (int property = angle; property < num_topological_property_types;
       property++) {
    if (topological_properties_in_molecule[property].size() > 1) {
      topological_properties_in_molecule[property].sort();
      auto&& it = topological_properties_in_molecule[property].begin();
      auto it_next = it;
      it_next++;

      for (; it_next != topological_properties_in_molecule[property].end();
           it++) {
        if ((*it)->GetInvolvedAtoms() == (*it_next)->GetInvolvedAtoms()) {
          topological_properties_in_molecule[property].erase(it_next);
          it_next = it;
          it--;
        }

        it_next++;
      }
    }

    if (united_atom) {
      auto&& it = topological_properties_in_molecule[property].begin();

      for (; it != topological_properties_in_molecule[property].end(); it++) {
        for (auto&& at : (*it)->GetInvolvedAtoms()) {
          if (indexed_atoms[at].atomtype.find("united") != std::string::npos) {
            auto it_erase = it;
            it--;
            topological_properties_in_molecule[property].erase(it_erase);
            break;
          }
        }

        if (!topological_properties_in_molecule[property].size()) break;
      }
    }
  }

  if (unique_torsionals) {
    if (topological_properties_in_molecule[torsional_dihedral].size() > 1) {
      for (auto it =
               topological_properties_in_molecule[torsional_dihedral].begin();
           it != topological_properties_in_molecule[torsional_dihedral].end();
           it++) {
        auto it_next = it;
        it_next++;

        for (; it_next !=
               topological_properties_in_molecule[torsional_dihedral].end();
             it_next++) {
          if ((*it)->GetInvolvedAtoms()[1] ==
                  (*it_next)->GetInvolvedAtoms()[1] &&
              (*it)->GetInvolvedAtoms()[2] ==
                  (*it_next)->GetInvolvedAtoms()[2]) {
            auto it_prev = it_next;
            it_prev--;
            topological_properties_in_molecule[torsional_dihedral].erase(
                it_next);
            it_next = it_prev;
          }
        }
      }
    }
  }

  for (auto&& prop_type : topological_properties_in_molecule) {
    if (prop_type.size()) {
      molecule_macros.AddElement(prop_type.front()->GetPropertyName() + "s");
      auto&& cur_element = molecule_macros.GetLastChild();
      cur_element->attributes["amount"] = std::to_string(prop_type.size());

      for (auto&& prop : prop_type) {
        cur_element->AddElement(prop->GetPropertyName());
        auto&& cur_property = cur_element->GetLastChild();
        cur_property->attributes["parameter"] = prop->GetAttribute("parameter");
        cur_property->AddElement("involved_atoms");
        auto&& involved_atoms = cur_property->GetLastChild();

        for (auto&& atom : prop->GetInvolvedAtoms()) {
          involved_atoms->AddElement("involved_atom");
          involved_atoms->GetLastChild()->value = indexed_atoms[atom].atom_id;
        }
      }
    }
  }

  if (united_atom) {
    for (auto&& at_it = atoms->children.begin(); at_it != atoms->children.end();
         at_it++) {  // atoms
      if ((*at_it)->attributes["atom_type"].find("united") !=
          std::string::npos) {
        auto erase_it = at_it;
        at_it--;
        molecule_macros.GetFirstChild()->RemoveChild(erase_it);
      }
    }

    // reSet amount of atoms
    atoms->attributes["amount"] = std::to_string(atoms->children.size());
  }

  if (!third_neighbor_exclusions) {
    for (auto&& at : atoms->children) {
      const std::string atom_id = at->attributes["atom_id"];
      std::vector<const AtomInMolecule*> excluded_atoms{};

      for (auto&& bd : topological_properties_in_molecule[bond]) {
        if (indexed_atoms[bd->GetInvolvedAtoms().front()].atom_id == atom_id ||
            indexed_atoms[bd->GetInvolvedAtoms().back()].atom_id == atom_id) {
          excluded_atoms.push_back(&indexed_atoms[bd->GetInvolvedAtoms()[0]]);
          excluded_atoms.push_back(&indexed_atoms[bd->GetInvolvedAtoms()[1]]);
        }
      }

      for (auto&& ang : topological_properties_in_molecule[angle]) {
        if (indexed_atoms[ang->GetInvolvedAtoms().front()].atom_id == atom_id ||
            indexed_atoms[ang->GetInvolvedAtoms().back()].atom_id == atom_id ||
            indexed_atoms[ang->GetInvolvedAtoms()[1]].atom_id == atom_id) {
          excluded_atoms.push_back(&indexed_atoms[ang->GetInvolvedAtoms()[0]]);
          excluded_atoms.push_back(&indexed_atoms[ang->GetInvolvedAtoms()[1]]);
          excluded_atoms.push_back(&indexed_atoms[ang->GetInvolvedAtoms()[2]]);
        }
      }

      std::sort(excluded_atoms.begin(), excluded_atoms.end());
      excluded_atoms.erase(
          std::unique(excluded_atoms.begin(), excluded_atoms.end()),
          excluded_atoms.end());

      if (excluded_atoms.size()) {
        at->AddElement("excluded_atoms");

        for (auto&& ea : excluded_atoms) {
          at->GetFirstChild()->AddElement("excluded_atom");
          at->GetFirstChild()->GetLastChild()->value = ea->atom_id;
        }
      }
    }

  } else {
    for (auto&& at : atoms->children) {
      const std::string atom_id = at->attributes["atom_id"];
      std::vector<const AtomInMolecule*> excluded_atoms;

      for (auto&& bd : topological_properties_in_molecule[bond]) {
        if (indexed_atoms[bd->GetInvolvedAtoms().front()].atom_id == atom_id ||
            indexed_atoms[bd->GetInvolvedAtoms().back()].atom_id == atom_id) {
          excluded_atoms.push_back(&indexed_atoms[bd->GetInvolvedAtoms()[0]]);
          excluded_atoms.push_back(&indexed_atoms[bd->GetInvolvedAtoms()[1]]);
        }
      }

      for (auto&& ang : topological_properties_in_molecule[angle]) {
        if (indexed_atoms[ang->GetInvolvedAtoms().front()].atom_id == atom_id ||
            indexed_atoms[ang->GetInvolvedAtoms().back()].atom_id == atom_id ||
            indexed_atoms[ang->GetInvolvedAtoms()[1]].atom_id == atom_id) {
          excluded_atoms.push_back(&indexed_atoms[ang->GetInvolvedAtoms()[0]]);
          excluded_atoms.push_back(&indexed_atoms[ang->GetInvolvedAtoms()[1]]);
          excluded_atoms.push_back(&indexed_atoms[ang->GetInvolvedAtoms()[2]]);
        }
      }

      for (auto&& tors :
           topological_properties_in_molecule[torsional_dihedral]) {
        if (indexed_atoms[tors->GetInvolvedAtoms().front()].atom_id ==
                atom_id ||
            indexed_atoms[tors->GetInvolvedAtoms().back()].atom_id == atom_id ||
            indexed_atoms[tors->GetInvolvedAtoms()[1]].atom_id == atom_id ||
            indexed_atoms[tors->GetInvolvedAtoms()[2]].atom_id == atom_id) {
          excluded_atoms.push_back(&indexed_atoms[tors->GetInvolvedAtoms()[0]]);
          excluded_atoms.push_back(&indexed_atoms[tors->GetInvolvedAtoms()[1]]);
          excluded_atoms.push_back(&indexed_atoms[tors->GetInvolvedAtoms()[2]]);
          excluded_atoms.push_back(&indexed_atoms[tors->GetInvolvedAtoms()[3]]);
        }
      }

      std::sort(excluded_atoms.begin(), excluded_atoms.end());
      excluded_atoms.erase(
          std::unique(excluded_atoms.begin(), excluded_atoms.end()),
          excluded_atoms.end());

      if (excluded_atoms.size()) {
        at->AddElement("excluded_atoms");

        for (auto&& ea : excluded_atoms) {
          at->GetFirstChild()->AddElement("excluded_atom");
          at->GetFirstChild()->GetLastChild()->value = ea->atom_id;
        }
      }
    }
  }
}

bool operator<(const std::pair<std::string, const AtomTypeSet*>& at1,
               const std::pair<std::string, const AtomTypeSet*>& at2) {
  const Atom& atom1 = *(at1.second->GetAtomSet().begin());
  const Atom& atom2 = *(at2.second->GetAtomSet().begin());

  if (atom1 < atom2)
    return true;

  else
    return false;
}

bool operator<(const std::shared_ptr<TopologicalPropertyBase>& b1,
               const std::shared_ptr<TopologicalPropertyBase>& b2) {
  for (size_t i = 0; i < b1->GetNumInvolvedAtoms(); i++) {
    if (b1->GetInvolvedAtoms()[i] < b2->GetInvolvedAtoms()[i])
      return true;

    else if (b1->GetInvolvedAtoms()[i] > b2->GetInvolvedAtoms()[i])
      return false;
  }

  return false;
}

const std::string GetParameterNameFromFragment(
    const TblFragment& frag, const std::string& fragment_id,
    const TopologicalPropertyBase& prop,
    const std::unordered_map<std::string, std::pair<std::string, std::string>>&
        core_atom_types,
    const std::vector<std::string>& involved_atoms_codes,
    const std::unordered_map<std::string, std::vector<std::string>>&
        connected_fragments) {
  std::string parameter_name = prop.GetPropertyAbbreviation();

  if (prop.GetInvolvedAtoms().size()) {
    parameter_name +=
        "_" + core_atom_types
                  .find(fragment_id + "_" +
                        frag.GetTblAtom(prop.GetInvolvedAtoms()[0]).GetAtomID())
                  ->second.second;
    PropertyType type =
        topologicalPropertyMap.find(prop.GetPropertyName())->second;

    if (type != improper_dihedral) {
      for (size_t i = 0; i < prop.GetInvolvedAtoms().size() - 1; i++) {
        const size_t idx1 = prop.GetInvolvedAtoms()[i];
        const size_t idx2 = prop.GetInvolvedAtoms()[i + 1];

        if (frag.GetMatrix().GetElement(idx1, idx2) == 1)
          parameter_name += "-";

        else if (frag.GetMatrix().GetElement(idx1, idx2) == 1.5)
          parameter_name += ":";

        else if (frag.GetMatrix().GetElement(idx1, idx2) == 2)
          parameter_name += "=";

        else if (frag.GetMatrix().GetElement(idx1, idx2) == 3)
          parameter_name += "#";

        else
          parameter_name += ",";

        parameter_name +=
            core_atom_types
                .find(fragment_id + "_" + frag.GetTblAtom(idx2).GetAtomID())
                ->second.second;
      }

    } else {
      for (size_t i = 0; i < prop.GetInvolvedAtoms().size() - 1; i++) {
        const size_t idx1 = prop.GetInvolvedAtoms()[0];
        const size_t idx2 = prop.GetInvolvedAtoms()[i + 1];

        if (frag.GetMatrix().GetElement(idx1, idx2) == 1)
          parameter_name += ",-";

        else if (frag.GetMatrix().GetElement(idx1, idx2) == 1.5)
          parameter_name += ",:";

        else if (frag.GetMatrix().GetElement(idx1, idx2) == 2)
          parameter_name += ",=";

        else if (frag.GetMatrix().GetElement(idx1, idx2) == 3)
          parameter_name += ",#";

        else
          parameter_name += ",";

        parameter_name +=
            core_atom_types
                .find(fragment_id + "_" + frag.GetTblAtom(idx2).GetAtomID())
                ->second.second;
      }
    }

    parameter_name += "/";

    for (auto&& ia : involved_atoms_codes) {
      if (connected_fragments.find(ia) != connected_fragments.end()) {
        for (auto&& cf : connected_fragments.find(ia)->second)
          parameter_name += cf + ",";

        if (connected_fragments.find(ia)->second.size())
          parameter_name.back() = '%';
      }
    }

    if (involved_atoms_codes.size()) parameter_name.pop_back();
  }

  return parameter_name;
}

// create the parameter name for automatically determined properties (currently:
// angles and torsionals)
const std::string GetParameterNameAutomatically(const topology_builder::AdjacencyMatrix& A,
                                                std::shared_ptr<TopologicalPropertyBase>& cur_property,
                                                const std::vector<AtomInMolecule>& indexed_atoms,
                                                const std::unordered_map<std::string, std::pair<std::string, std::string>>& core_atom_types/*,
                                                const std::vector<std::string>& involved_atoms_codes,
                                                const std::unordered_map<std::string, std::vector<std::string>>& connected_fragments*/) {
  std::string parameter_name = cur_property->GetPropertyAbbreviation() + "_";
  parameter_name +=
      indexed_atoms[cur_property->GetInvolvedAtoms()[0]].parameter;

  for (size_t i = 0; i < cur_property->GetInvolvedAtoms().size() - 1; i++) {
    const size_t idx1 = cur_property->GetInvolvedAtoms()[i];
    const size_t idx2 = cur_property->GetInvolvedAtoms()[i + 1];

    if (A.GetElement(idx1, idx2) == 1)
      parameter_name += "-";

    else if (A.GetElement(idx1, idx2) == 1.5)
      parameter_name += ":";

    else if (A.GetElement(idx1, idx2) == 2)
      parameter_name += "=";

    else if (A.GetElement(idx1, idx2) == 3)
      parameter_name += "#";

    else
      parameter_name += ",";

    parameter_name += indexed_atoms[idx2].parameter;
  }

  parameter_name += "/";

  for (auto&& atmidx : cur_property->GetInvolvedAtoms()) {
    for (auto&& linkatoms : indexed_atoms[atmidx].connected_link_atoms)
      parameter_name += linkatoms + ",";

    parameter_name.back() = '%';
  }

  parameter_name.pop_back();
  return parameter_name;
}

}  // namespace topology_builder

}  // namespace combi_ff