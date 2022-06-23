// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "MoleculeDecomposer.h"

#include <iomanip>

#include "TblFragment.h"
#include "smilesToMatrix.h"
#include "version.h"

namespace combi_ff {

namespace topology_builder {

/*
!!!!!!!!!!!!!!!!!!!!HalfLink.fragIdx is one higher than the C++ index of the
fragment in the fragments vector!!!!!!!! i.e. for frag[0] the fragIdx is one!!!!
*/

FamilyDecomposer::FamilyDecomposer(
    const std::string& filename_family_isomer_enumeration,
    const IOFileProperties& io_file_properties,
    const std::vector<TblFragment>& tbl_fragments)
    : filename_family_isomer_enumeration(filename_family_isomer_enumeration),
      io_file_properties(&io_file_properties),
      tbl_fragments(&tbl_fragments) {
  CreateFamilyMoleculeDecompositions();
}

void FamilyDecomposer::CreateFamilyMoleculeDecompositions() {
  XmlParserIn family_isomer_enumeration_parser(
      filename_family_isomer_enumeration, XmlParserIn::read_element_by_element);
  family_isomer_enumeration_parser.ReadUntilElement("enumeration_time");
  auto&& family_isomer_enumeration_root =
      family_isomer_enumeration_parser.GetRoot();
  CheckXMLFormat(family_isomer_enumeration_root);
  filename_molecule_decomposition =
      (io_file_properties->output_dir_molecule_decompositions.size()
           ? io_file_properties->output_dir_molecule_decompositions
           : io_file_properties->output_dir) +
      "molecule_decompositions_" + family_code + ".xml";
  XmlParserOut molecule_decompositions_parser(filename_molecule_decomposition,
                                              "molecule_decompositions");
  XmlElement* molecule_decompositions_root =
      molecule_decompositions_parser.GetTree().GetRootPointer();
  molecule_decompositions_root->attributes["tbl_version"] =
      combi_ff::current_version;
  molecule_decompositions_root->SetAttribute(family_isomer_enumeration_root,
                                             "family_code");
  molecule_decompositions_root->SetAttribute(family_isomer_enumeration_root,
                                             "family_version");
  molecule_decompositions_root->SetAttribute(family_isomer_enumeration_root,
                                             "enu_version");
  molecule_decompositions_parser.WriteHead();
  size_t num_isomers = GetNumIsomers(family_isomer_enumeration_root);

  if (num_isomers) {
    size_t counter(0);

    while (family_isomer_enumeration_parser.ReadUntilElement("isomer_list")) {
      // while(family_isomer_enumeration_parser.ReadUntilElement("isomer")) {
      for (auto&& isomer_list :
           family_isomer_enumeration_parser.GetCurNode()->children) {
        for (auto&& isomer : isomer_list->children) {
          std::cout << isomer->tag << std::endl;
          auto&& isomer_child_ptr = isomer->children.begin();
          (*isomer_child_ptr)->CheckTagName("constitutional_SMILES");
          (*isomer_child_ptr)->CheckValue();
          const std::string& smiles = (*isomer_child_ptr)->value;
          std::cout << "molecule nr. " << ++counter << ": " << smiles << '\n';
          const topology_builder::AdjacencyMatrix A(ConvertToMatrix(smiles));
          molecule_decompositions_root->AddElement("molecule_decomposition");
          XmlElement& molecule_decomposition =
              *(molecule_decompositions_root->GetLastChild());
          molecule_decomposition.attributes["smiles"] = smiles;
          molecule_decomposition.SetAttribute(*isomer, "isomer_id");
          MoleculeDecomposer(smiles, *tbl_fragments, molecule_decomposition);
          molecule_decompositions_parser.WriteAndRemoveLastChild();
        }
      }
    }
  }

  molecule_decompositions_parser.WriteTail();
}

void FamilyDecomposer::CheckXMLFormat(
    XmlElement& family_isomer_enumeration_root) {
  family_isomer_enumeration_root.CheckAttributeSize(3);
  family_isomer_enumeration_root.CheckAttribute("enu_version");
  family_isomer_enumeration_root.CheckAttribute("family_code");
  family_isomer_enumeration_root.CheckAttribute("family_version");
  family_code = family_isomer_enumeration_root.attributes["family_code"];
  family_isomer_enumeration_root.CheckNumberOfChildren_equal(3);
}

size_t FamilyDecomposer::GetNumIsomers(
    const XmlElement& family_isomer_enumeration_root) {
  auto&& enumeration_info_ptr = family_isomer_enumeration_root.children.begin();
  (*enumeration_info_ptr)->CheckTagName("enumeration_type");
  (*enumeration_info_ptr)->CheckAttribute("type");
  const std::string enu_type = (*enumeration_info_ptr)->attributes["type"];
  enumeration_info_ptr++;
  (*enumeration_info_ptr)->CheckTagName("number_of_isomers");
  const std::string num_isomers = (*enumeration_info_ptr)->value;
  enumeration_info_ptr++;
  (*enumeration_info_ptr)->CheckTagName("enumeration_time");
  const std::string enutime = (*enumeration_info_ptr)->value;
  return std::stoi(num_isomers);
}

const std::string& FamilyDecomposer::GetFilename() const {
  return filename_molecule_decomposition;
}

const std::string& FamilyDecomposer::GetFamilyCode() const {
  return family_code;
}

MoleculeDecomposer::MoleculeDecomposer(
    const std::string& smiles, const std::vector<TblFragment>& tbl_fragments,
    XmlElement& molecule_decomposition)
    : A(ConvertToMatrix(smiles).ExtendToFullMatrix()),
      N(A.GetN()),
      tbl_fragments(tbl_fragments),
      is_assigned_core_atom(N, false),
      is_assigned_core_atom_all_true(N, true),
      is_assigned_link_atom(N, 0),
      core_atoms(N, NULL),
      link_atoms(N, std::vector<const TblAtom*>(0)),
      fragments(0),
      half_links(0),
      links(0),
      all_assigned(false) {
  CreateMoleculeDecomposition(molecule_decomposition);
}

void MoleculeDecomposer::CreateMoleculeDecomposition(
    XmlElement& molecule_decomposition) {
  for (size_t i = 0; i < tbl_fragments.size(); i++) {
    const TblFragment& frag = (tbl_fragments)[i];

    while (FindBondLinkingMatch(frag)) {
      if (is_assigned_core_atom == is_assigned_core_atom_all_true) {
        all_assigned = true;
        break;
      }
    }

    if (all_assigned) break;
  }

  if (!CheckValidity(molecule_decomposition)) return;

  // sort the links by their fragment indices
  for (size_t i = 0; i < links.size(); i++) {
    for (size_t j = i + 1; j < links.size(); j++) {
      if (links[i].GetHalfLink1()->GetFragmentIndex() >
          links[j].GetHalfLink1()->GetFragmentIndex())
        std::swap(links[i], links[j]);

      else if (links[i].GetHalfLink1()->GetFragmentIndex() ==
                   links[j].GetHalfLink1()->GetFragmentIndex() &&
               links[i].GetHalfLink2()->GetFragmentIndex() >
                   links[j].GetHalfLink2()->GetFragmentIndex())
        std::swap(links[i], links[j]);
    }
  }

  for (auto&& frg : fragments) std::cout << frg->GetCode() << " ";

  std::cout << std::endl;
  molecule_decomposition.AddElement("fragments");
  auto&& fragments_element = *(molecule_decomposition.GetLastChild());
  size_t i = 0;

  for (auto&& frag : fragments) {
    fragments_element.AddElement("fragment");
    fragments_element.GetLastChild()->attributes["id"] =
        "f" + std::to_string(++i);
    fragments_element.GetLastChild()->value = frag->GetCode();

    if (frag->GetVersion() != combi_ff::current_version)
      std::cout << "?Warning: currently running combi_ff version "
                << combi_ff::current_version << " but fragment "
                << frag->GetCode() << " has version " << frag->GetVersion()
                << "\n";
  }

  molecule_decomposition.AddElement("linkages");
  auto&& linkages_element = *(molecule_decomposition.GetLastChild());

  for (auto&& link : links) {
    linkages_element.AddElement("linkage");
    auto&& linkage_element = *(linkages_element.GetLastChild());
    linkage_element.AddElement("involved_fragment");
    auto&& involved_fragment_element1 = *(linkage_element.GetLastChild());
    involved_fragment_element1.AddElement("fragment_id");
    involved_fragment_element1.GetLastChild()->value =
        "f" + std::to_string(link.GetHalfLink1()->GetFragmentIndex());
    involved_fragment_element1.AddElement("linksite");
    involved_fragment_element1.GetLastChild()->value =
        link.GetHalfLink1()->GetLinkatomName();
    linkage_element.AddElement("involved_fragment");
    auto&& involved_fragment_element2 = *(linkage_element.GetLastChild());
    involved_fragment_element2.AddElement("fragment_id");
    involved_fragment_element2.GetLastChild()->value =
        "f" + std::to_string(link.GetHalfLink2()->GetFragmentIndex());
    involved_fragment_element2.AddElement("linksite");
    involved_fragment_element2.GetLastChild()->value =
        link.GetHalfLink2()->GetLinkatomName();
  }
}

bool MoleculeDecomposer::FindBondLinkingMatch(const TblFragment& frag) {
  const size_t m = A.GetN();
  const size_t n = frag.GetMatrix().GetN();

  if (n > m) return false;

  const AtomVector<combi_ff::CnvAtom> atoms = A.GetAtomVector();
  const FragmentMatrixTbl& fragment_matrix = frag.GetMatrix();
  const std::vector<TblAtom>& fragment_tbl_atoms = frag.GetTblAtoms();
  std::vector<const AtomSet*> fragment_atom_set(0);

  for (auto&& at : fragment_tbl_atoms)
    fragment_atom_set.push_back(&(at.GetAtomTypes().GetAtomSet()));

  std::vector<const std::unordered_set<std::string>*> fragment_atom_names_sets(
      0);

  for (auto&& at : fragment_tbl_atoms)
    fragment_atom_names_sets.push_back(&(at.GetAtomTypes().GetAtomTypeNames()));

  ComparisonMatrix M(n, m);

  for (size_t i = 0; i < n; i++) {
    bool potential_match(false);

    if (fragment_tbl_atoms[i].IsCoreAtom()) {
      // std::cout << "case core\n";
      for (size_t j = 0; j < m; j++) {
        if (fragment_atom_set[i]->find(atoms[j].GetElementSymbol()) !=
            fragment_atom_set[i]->end()) {
          if (!(is_assigned_core_atom)[j]) {
            if (!is_assigned_link_atom[j]) {
              potential_match = true;
              M.SetElement(i, j, true);

            } else {
              bool incompatible_type(false);

              for (auto&& linkatom : link_atoms[j]) {
                if (linkatom->GetAtomTypes().GetAtomTypeNames().find(
                        fragment_tbl_atoms[i].GetAtomTypes().GetName()) ==
                    linkatom[j].GetAtomTypes().GetAtomTypeNames().end()) {
                  incompatible_type = true;
                  break;
                }
              }

              if (!incompatible_type) {
                potential_match = true;
                M.SetElement(i, j, true);
              }
            }
          }
        }
      }

    } else {
      // std::cout << "case link\n";
      for (size_t j = 0; j < m; j++) {
        if (fragment_atom_set[i]->find(atoms[j].GetElementSymbol()) !=
                fragment_atom_set[i]->end() ||
            fragment_atom_set[i]->find(atoms[j].GetUnitedAtomSymbol()) !=
                fragment_atom_set[i]->end()) {
          if (!is_assigned_core_atom[j] ||
              fragment_atom_names_sets[i]->find(
                  core_atoms[j]->GetAtomTypes().GetName()) !=
                  fragment_atom_names_sets[i]->end()) {
            potential_match = true;
            M.SetElement(i, j, true);
          }
        }
      }
    }

    if (!potential_match) {
      // std::cout << "found no match for " <<
      // fragment_tbl_atoms[i].GetAtomTypes() << std::endl;
      return false;
    }
  }

  // M.print();
  int k = 0;
  bool match = UllmannMatch(M, k, n, m, fragment_matrix, fragment_tbl_atoms);

  if (match) {
    // std::cout << "found match for " << frag.GetCode() << '\n';
    // std::cout << "  " << fragAtoms << '\n';
    size_t j;
    bool unlinked = false;

    // i are the indices of the atoms in the fragment
    for (size_t i = 0; i < n; i++) {
      j = GetMatchIndex(M, i);

      if (fragment_tbl_atoms[i].IsCoreAtom()) {
        if (!(is_assigned_core_atom)[j]) unlinked = true;
      }
    }

    if (unlinked) {
      fragments.push_back(&frag);

      // i are the indices of the atoms in the fragment
      for (size_t i = 0; i < n; i++) {
        // j is the index of the atom in the adjacency matrix that is matched to
        // atom i in the fragment
        j = GetMatchIndex(M, i);

        // if atom i in the fragment is a core atom, simply assign it as a core
        // atom
        if (fragment_tbl_atoms[i].IsCoreAtom()) {
          // std::cout << "  case core\n";
          (is_assigned_core_atom)[j] = true;
          core_atoms[j] = &fragment_tbl_atoms[i];

        } else {
          (is_assigned_link_atom)[j]++;
          link_atoms[j].push_back(&fragment_tbl_atoms[i]);

          // check if atom j in the adjacency matrix is already an assigned core
          // atom
          if ((is_assigned_core_atom)[j]) {
            // ii is the index of the core atom in the fragment that is
            // connected to the link atom
            size_t ii = frag.GetCoreIndex(i);
            // k is the index of the atom in the adjacency matrix that is
            // matched to the core atom ii in the fragment
            size_t k = GetMatchIndex(M, ii);
            // core is the index of the atom in the adjacency matrix that is
            // matched to link atom i in the fragment
            size_t core = GetMatchIndex(M, i);

            // search in the existing half_links for the matching halflink
            for (auto&& halflink = half_links.begin();
                 halflink != half_links.end(); halflink++) {
              if (halflink->GetAtomIndex() == k &&
                  halflink->GetCoreIndex() == core) {
                std::string name_link_atom_2 =
                    fragment_tbl_atoms[i].GetAtomID();
                std::string name_link_atom_1 = halflink->GetLinkatomName();
                links.push_back(Link(
                    HalfLink(halflink->GetAtomIndex(), halflink->GetCoreIndex(),
                             halflink->GetFragmentIndex(), name_link_atom_1),
                    HalfLink(k, core, fragments.size(), name_link_atom_2)));
                half_links.erase(halflink);
                (is_assigned_link_atom)[k]++;
                link_atoms[k].push_back(&fragment_tbl_atoms[ii]);
                break;
              }
            }

          } else
            half_links.push_back(
                HalfLink(j, GetMatchIndex(M, frag.GetCoreIndex(i)),
                         fragments.size(), fragment_tbl_atoms[i].GetAtomID()));
        }
      }

      return true;
    }

    return true;
  }

  return false;
}

bool MoleculeDecomposer::CheckValidity(XmlElement& molecule_decomposition) {
  if (!A.hasCycle() && fragments.size() != links.size() + 1) {
    molecule_decomposition.value =
        "\n" + std::string(molecule_decomposition.distance * 2 + 2, ' ') +
        "<!-- unable to find a valid decomposition for this molecule with the "
        "given fragments because\n" +
        std::string(molecule_decomposition.distance * 2 + 2, ' ') +
        " number of fragments (" + std::to_string(fragments.size()) +
        ") is not equal to number of links plus one " + "(" +
        std::to_string(links.size()) +
        " + 1 = " + std::to_string(links.size() + 1) + ")\n" +
        std::string(molecule_decomposition.distance * 2 + 2, ' ') +
        "fragments are:\n" +
        std::string(molecule_decomposition.distance * 2 + 2, ' ');

    for (auto&& frag : fragments)
      molecule_decomposition.value += frag->GetCode() + " ";

    molecule_decomposition.value +=
        "\n" + std::string(molecule_decomposition.distance * 2 + 2, ' ') +
        "linkages are:\n" +
        std::string(molecule_decomposition.distance * 2 + 2, ' ');

    for (auto&& link : links)
      molecule_decomposition.value +=
          fragments[link.GetHalfLink1()->GetFragmentIndex() - 1]->GetCode() +
          "[" + link.GetHalfLink1()->GetLinkatomName() + ":" +
          link.GetHalfLink2()->GetLinkatomName() + "]" +
          fragments[link.GetHalfLink2()->GetFragmentIndex() - 1]->GetCode() +
          " ";

    molecule_decomposition.value +=
        "-->\n" + std::string(molecule_decomposition.distance * 2, ' ');
    return false;
  }

  for (size_t i = 0; i < A.GetN(); i++) {
    if (!(is_assigned_core_atom[i] ||
          (A.GetAtomVector())[i].GetUnitedAtomSymbol() == "H")) {
      molecule_decomposition.value =
          "\n" + std::string(molecule_decomposition.distance * 2 + 2, ' ') +
          "<!-- unable to find a valid decomposition for this molecule with "
          "the given fragments because\n" +
          std::string(molecule_decomposition.distance * 2 + 2, ' ') + " atom " +
          std::to_string(i) + " (" +
          A.GetAtomVector()[i].GetUnitedAtomSymbol() +
          ") is not an assigned core atom\n" +
          std::string(molecule_decomposition.distance * 2 + 2, ' ') +
          "fragments are:\n" +
          std::string(molecule_decomposition.distance * 2 + 2, ' ');

      for (auto&& frag : fragments)
        molecule_decomposition.value += frag->GetCode() + " ";

      molecule_decomposition.value +=
          "\n" + std::string(molecule_decomposition.distance * 2 + 2, ' ') +
          "linkages are:\n" +
          std::string(molecule_decomposition.distance * 2 + 2, ' ');

      for (auto&& link : links)
        molecule_decomposition.value +=
            fragments[link.GetHalfLink1()->GetFragmentIndex() - 1]->GetCode() +
            "[" + link.GetHalfLink1()->GetLinkatomName() + ":" +
            link.GetHalfLink2()->GetLinkatomName() + "]" +
            fragments[link.GetHalfLink2()->GetFragmentIndex() - 1]->GetCode() +
            " ";

      molecule_decomposition.value +=
          "-->\n" + std::string(molecule_decomposition.distance * 2, ' ');
      return false;
    }
  }

  if (half_links.size() != 0) {
    molecule_decomposition.value =
        "\n" + std::string(molecule_decomposition.distance * 2 + 2, ' ') +
        "<!-- unable to find a valid decomposition for this molecule with the "
        "given fragments because\n" +
        std::string(molecule_decomposition.distance * 2 + 2, ' ') +
        " not all link atoms are linked\n" +
        std::string(molecule_decomposition.distance * 2 + 2, ' ') +
        "fragments are:\n" +
        std::string(molecule_decomposition.distance * 2 + 2, ' ');

    for (auto&& frag : fragments)
      molecule_decomposition.value += frag->GetCode() + " ";

    molecule_decomposition.value +=
        "\n" + std::string(molecule_decomposition.distance * 2 + 2, ' ') +
        "linkages are:\n" +
        std::string(molecule_decomposition.distance * 2 + 2, ' ');

    for (auto&& link : links)
      molecule_decomposition.value +=
          fragments[link.GetHalfLink1()->GetFragmentIndex() - 1]->GetCode() +
          "[" + link.GetHalfLink1()->GetLinkatomName() + ":" +
          link.GetHalfLink2()->GetLinkatomName() + "]" +
          fragments[link.GetHalfLink2()->GetFragmentIndex() - 1]->GetCode() +
          " ";

    molecule_decomposition.value +=
        "-->\n" + std::string(molecule_decomposition.distance * 2, ' ');
    return false;
  }

  return true;
}

size_t MoleculeDecomposer::GetMatchIndex(const ComparisonMatrix& M,
                                         const size_t i) {
  for (size_t j = 0; j < M.GetM(); j++) {
    if (M.GetElement(i, j)) return j;
  }

  throw std::runtime_error(
      "couldn't find match index for the above ComparisonMatrix for i=" +
      std::to_string(i) + "\n");
}

bool MoleculeDecomposer::UllmannMatch(
    ComparisonMatrix& M, int k, const size_t n, const size_t m,
    const FragmentMatrixTbl& fragment_matrix,
    const std::vector<TblAtom>& fragment_tbl_atoms) {
  if (k == (int)n) {
    return true;
  }

  ComparisonMatrix M_save(n, m);

  for (size_t l = 0; l < m; l++) {
    if (M.GetElement(k, l)) {
      M_save = M;

      for (size_t j = 0; j < m; j++) M.SetElement(k, j, false);

      for (size_t i = 0; i < n; i++) M.SetElement(i, l, false);

      M.SetElement(k, l, true);

      if (Refine(M, k, n, m, fragment_matrix, fragment_tbl_atoms)) {
        if (UllmannMatch(M, k + 1, n, m, fragment_matrix, fragment_tbl_atoms))
          return true;
      }

      M = M_save;
    }
  }

  return false;
}

bool MoleculeDecomposer::Refine(
    combi_ff::ComparisonMatrix& M, int k, const size_t n, const size_t m,
    const FragmentMatrixTbl& fragment_matrix,
    const std::vector<TblAtom>& fragment_tbl_atoms) {
  bool changed(true), valid(false), found(false);

  while (changed) {
    changed = false;

    for (size_t i = k + 1; i < n; i++) {
      for (size_t j = 0; j < m; j++) {
        if (M.GetElement(i, j)) {
          valid = true;

          for (auto&& i2 : fragment_tbl_atoms[i].GetNeighbours()) {
            found = false;

            for (auto&& j2 : A.GetAtomVector()[j].GetNeighbours()) {
              if (M.GetElement(i2, j2)) {
                if (fragment_matrix.GetElement(i, i2) == A.GetElement(j, j2)) {
                  found = true;
                  break;
                }
              }
            }

            if (!found) {
              valid = false;
              break;
            }
          }

          if (!valid) {
            M.SetElement(i, j, false);
            changed = true;
            bool nonzero = false;

            for (size_t h = 0; h < m; h++) {
              if (M.GetElement(i, h)) nonzero = true;
            }

            if (!nonzero) return false;
          }
        }
      }
    }
  }

  return true;
}

}  // namespace topology_builder

}  // namespace combi_ff
