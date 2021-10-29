#include "Enumerator.h"
#include "printInfo.h"
#include "HydrogenDistribution.h"
#include <math.h>
#include "matchingEnu.h"
#include "StereoGenerator.h"
#include "Family.h"
#include "EnumSpecification.h"
#include "XmlParser.h"

namespace combi_ff {

namespace enu {

Enumerator::Enumerator(const enu::Family& family,
                       const enu::PseudoatomMap& pseudoatoms,
                       const bool stereo,
                       const std::string& output_file_name) :
  code(family.GetCode()),
  lambda_ranges_vec(family.GetLambdaRangesVec()),
  used_atom_vectors(family.GetUsedAtomVectors()),
  ranged_properties(family.GetRangedProperties()),
  substructures(&family.GetSubstructures()),
  pseudoatoms(&pseudoatoms),
  max_degree(family.GetMaxDeg()),
  stereo(stereo),
  num_isomers(0),
  num_isomers_for_listing(0),
  has_pseudoatom(family.GetHasPseudoatom()),
  formula_idx(0),
  ranges(num_ranges, 0) {
  std::cout << "******************************************************\n"
            << "finding isomers for family " << family.GetCode() << "...\n\n";
  output_file.open(output_file_name);

  if (!output_file.is_open())
    throw std::runtime_error("couldn't open " + output_file_name);
}

Enumerator::Enumerator(const enu::EnumSpecifications& enum_spec,
                       const std::string& output_file_name) :
  code("mol"), lambda_ranges_vec(enum_spec.lambda_ranges_vec),
  used_atom_vectors(enum_spec.used_atom_vectors),
  ranged_properties(enum_spec.ranged_properties),
  substructures(NULL),
  pseudoatoms(NULL),
  max_degree(enum_spec.max_degree),
  stereo(enum_spec.stereo),
  deg_unsaturations(0),
  num_isomers(0),
  num_isomers_for_listing(0),
  has_pseudoatom(false),
  formula_idx(0),
  ranges(num_ranges, 0) {
  std::cout << "******************************************************\n"
            << "finding isomers for formula given directly in input file...\n\n";
  output_file.open(output_file_name);

  if (!output_file.is_open())
    throw std::runtime_error("couldn't open " + output_file_name);
}


/*************
GetTER METHODS
*************/
size_t Enumerator::GetNumIsomers() const {
  return num_isomers;
}


/***************************************************************
FUNCTION TO PREPARE FOR MOLECULE ENUMERATION FOR A GIVEN used_atom_vectors
****************************************************************/
void Enumerator::EnumerateIsomers() {
  /************************************************************************
  FOR A GIVEN AtomVector atom_types, LOOP OVER ALL THE DIFFERENT LAMBDA VECTORS
  *************************************************************************/
  for (size_t k = 0; k < used_atom_vectors.size(); k++) {
    const combi_ff::AtomVector<combi_ff::Atom>& atom_types = used_atom_vectors[k];
    const std::list<combi_ff::LambdaVector>& lambda_ranges = lambda_ranges_vec[k];

    if (!atom_types.size())
      throw std::runtime_error("atom_types vector is empty in Enumerator::enumerateIsomers()");

    for (auto && lambda : lambda_ranges) {
      if (lambda.size() != atom_types.size())
        throw std::runtime_error("lambdaRanges vector doesn't have same size as aotmTypes vector in Enumerator::enumerateIsomers()");

      size_t sum_degree(0);
      double degree_unsaturations_(0);

      for (size_t i = 0; i < lambda.size(); i++) {
        degree_unsaturations_ += ((double)atom_types[i].GetDegree() - 2.) *
                                 (double)lambda[i];
        sum_degree += (atom_types[i].GetDegree()) * lambda[i];
      }

      ranges[range_unsaturations] = deg_unsaturations
                                    = (int)floor((degree_unsaturations_ / 2.) + 1);

      /***********************************************************
      IF sum_degree and degree_unsaturations MEET THE REQUIREMENTS; PRODUCE THE ISOMERS
      ************************************************************/
      if (!(sum_degree % 2) &&
          IsInRange(deg_unsaturations, ranged_properties[range_unsaturations])) {
        enumerator_arguments.formula_ptr = std::unique_ptr<const std::string>
                                           (new std::string(Enumerator::CreateCanonicalFormula(lambda, atom_types)));
        //print information on the new class of molecules,
        std::cout << "formula nr. "
                  << ++formula_idx << ": "
                  << *enumerator_arguments.formula_ptr
                  << '\n';
        //determine if the local maxDeg can be further restricted due to a small DoU, as
        // -> if DoU is zero, the max deg can be restricted to 1, if it's one the max deg can be restricted to 2 etc.
        // -> thus, select the minimum of maxDeg and DoU + 1 for max_degree_local
        size_t max_degree_local = std::min((int)max_degree, deg_unsaturations + 1);
        num_isomers_for_formula = 0;
        std::cout << " -> found " << num_isomers_for_formula << std::flush << '\r';
        //find all molecules that fit the restrictions
        output_file << std::string(xml_indent_size * 2,
                                   ' ') << "<isomer_list formula=\"" << *enumerator_arguments.formula_ptr << "\"";
        const auto cur_pos = output_file.tellp();
        output_file << ">\n";
        Enumerator::EnumerateFormula(atom_types, lambda, max_degree_local, stereo);

        if (num_isomers_for_formula)
          output_file << std::string(xml_indent_size * 2, ' ') << "</isomer_list>\n";

        else {
          output_file.seekp(cur_pos);
          output_file << "/>\n";
        }

        num_isomers += num_isomers_for_formula;
        std::cout << " -> found " << num_isomers_for_formula
                  << " isomers for the formula " << *enumerator_arguments.formula_ptr <<
                  " with the given restrictions.\n"
                  << " -> total number of isomers in family so far is " << num_isomers << "\n";
      }
    }
  }

  output_file.close();
}


/****************************************************************
FUNCTION THAT CONTROLS THE ISOMER ENUMERATION FOR A GIVEN FORMULA
*****************************************************************/
void Enumerator::EnumerateFormula(const AtomVector<combi_ff::Atom>& atom_types,
                                  const LambdaVector& lambda, const size_t max_degree,
                                  const bool stereo) {
  AtomVector<combi_ff::Atom> atoms(0), non_H_atoms(0);
  atoms.reserve(std::accumulate(lambda.begin(), lambda.end(), 0));
  non_H_atoms.reserve(std::accumulate(lambda.begin(), lambda.end() - 1, 0));

  for (size_t i = 0; i < lambda.size() - 1; i++) {
    combi_ff::AtomVector<combi_ff::Atom> a(lambda[i], atom_types[i]);
    atoms.insert(atoms.end(), a.begin(), a.end());
    non_H_atoms.insert(non_H_atoms.end(), a.begin(), a.end());
  }

  combi_ff::AtomVector<combi_ff::Atom> a(lambda.back(), atom_types.back());
  atoms.insert(atoms.end(), a.begin(), a.end());

  if (atom_types.back().GetElementSymbol() != "H")
    non_H_atoms.insert(non_H_atoms.end(), a.begin(), a.end());

  //don't continue if there are only H atoms in the atom vector
  if (!non_H_atoms.size())
    return;

  //calculate total number of atoms N, number of non-hydrogen atoms Nhat, and number of hydrogen atoms NH
  const size_t N(atoms.size());
  const size_t N_hat(non_H_atoms.size());
  const size_t N_hyd(N - N_hat);
  //define connectivity vector, which describes how the hydrogen atoms are distributed among the non-hydrogen atoms
  //i.e. if atom with index 2 is connected to 3 hydrogen atoms, Hcon[2] is equal to 3
  enu::HydrogenDistribution hydrogen_dist(lambda, non_H_atoms, (int)N_hyd);
  //count the number of implicit hydrogen atoms and the number of implicit atoms from pseudoatoms
  enumerator_arguments.n_pseudoatoms = 0;
  size_t n_fixed_H(0);
  Enumerator::GetNImplicitAtoms(non_H_atoms, n_fixed_H);
  size_t N_full(N + enumerator_arguments.n_pseudoatoms + n_fixed_H);
  /*************************************************************
  LOOP OVER ALL THE POSSIBLE DISTRIBUTIONS OF THE HYDROGEN ATOMS
  *************************************************************/

  while (hydrogen_dist.GetNextDistribution()) {
    const combi_ff::ConnectivityVector H_con = hydrogen_dist.GetHcon();
    //make local copy of non_H_atoms, where the implicit hydrogens are added to the atoms
    combi_ff::AtomVector<combi_ff::Atom> atoms_hat(non_H_atoms);

    for (size_t i = 0; i < N_hat; i++)
      if (H_con[i])
        atoms_hat[i].SetNumHydrogens(H_con[i]);

    //sort atomsHat and create the new lambda vector lambdaHat
    combi_ff::LambdaVector lambda_hat = Sort(atoms_hat);
    //prepare an atom vector that contains the fixed hydrogens from the united atoms, in order to extend the
    //found matrices in the next step
    enumerator_arguments.full_matrix_BU_ptr = std::unique_ptr<AdjacencyMatrix>
                                              (new AdjacencyMatrix(N_full));
    enumerator_arguments.extended_matrix_BU_ptr = std::unique_ptr<AdjacencyMatrix>
                                                  (new AdjacencyMatrix(N_hat + enumerator_arguments.n_pseudoatoms));
    Enumerator::GetFullAndExtendedMatrixSkeletons(N_full, atoms_hat);
    //define the adjacency matrix A, and the three supplementary matrices M, L, C
    enumerator_arguments.A = std::shared_ptr<AdjacencyMatrix>(new AdjacencyMatrix(
                                                                N_hat, lambda_hat, atoms_hat));
    enumerator_arguments.max_fill_alg_ptr = std::unique_ptr<MaxFillAlg>
                                            (new MaxFillAlg(
                                               *enumerator_arguments.A, stereo, max_degree));

    if ((*enumerator_arguments.A).GetN() == 1) {
      if (atoms_hat[0].GetDegree() != 0)
        return;

      if (Enumerator::CheckMolecule())
        WriteMolecule();

    } else {
      while (enumerator_arguments.max_fill_alg_ptr->GetNextCanonicalMatrix()) {
        if (Enumerator::CheckMolecule())
          WriteMolecule();
      }
    }
  }
}


/**************************************************************************************
COUNT THE NUMER OF IMPLICIT HYDROGENS AND THE NUMBER OF IMPLICIT ATOMS FROM PSEUDOATOMS
***************************************************************************************/
void Enumerator::GetNImplicitAtoms(const AtomVector<combi_ff::Atom>& non_H_atoms,
                                   size_t& n_fixed_H) {
  for (size_t k = 0; k < non_H_atoms.size(); k++) {
    if (non_H_atoms[k].IsPseudoatom()) {
      auto psa = pseudoatoms->find(non_H_atoms[k].GetElementSymbol());
      n_fixed_H += psa->second.GetAtoms()[0].GetNumFixedHydrogens();

      for (size_t i = 1; i < psa->second.GetAtoms().size(); i++) {
        if (psa->second.GetAtoms()[i].GetElementSymbol() != "H") {
          enumerator_arguments.n_pseudoatoms++;

          if (psa->second.GetAtoms()[i].GetNumFixedHydrogens())
            n_fixed_H += psa->second.GetAtoms()[i].GetNumFixedHydrogens();

        } else
          n_fixed_H++;
      }

    } else if (
      non_H_atoms[k].GetNumFixedHydrogens()) //note: pseudoatoms can't have a fixed NH, if this is needed, another larger pseudoatom should be created
      n_fixed_H += non_H_atoms[k].GetNumFixedHydrogens();
  }
}

/***************************************************************************************************************
CREATE SKELETONS FOR THE fullMatrix (WITH EXPLICIT HYDROGENS) AND THE extendedMatrix (WITH EXPLICIT PSEUDOATOMS)
***************************************************************************************************************/
void Enumerator::GetFullAndExtendedMatrixSkeletons(const size_t N_full,
                                                   const AtomVector<combi_ff::Atom>& atoms_hat) {
  AtomVector<combi_ff::Atom> atoms_full(atoms_hat);
  AtomVector<combi_ff::Atom> atoms_extended(atoms_hat);
  atoms_full.reserve(N_full);

  for (size_t k = 0; k < atoms_full.size(); k++) {
    if (atoms_full[k].IsPseudoatom()) {
      auto psa = pseudoatoms->find(atoms_full[k].GetElementSymbol());

      for (size_t i = 1; i < psa->second.GetAtoms().size(); i++) {
        atoms_full.push_back(psa->second.GetAtoms()[i]);
        atoms_extended.push_back(psa->second.GetAtoms()[i]);

        if (psa->second.GetAtoms()[i].GetNumFixedHydrogens()) {
          atoms_full.back().SetNumFixedHydrogens(0);
          atoms_full.back().SetUnitedAtomSymbol(atoms_full.back().GetElementSymbol());
        }
      }

      const size_t nh = atoms_full[k].GetNumHydrogens();
      atoms_extended[k] = psa->second.GetAtoms()[0];
      atoms_extended[k].SetNumHydrogens(nh);
      atoms_full[k] = psa->second.GetAtoms()[0];
      atoms_full[k].SetNumHydrogens(nh);
    }

    if (atoms_full[k].GetNumFixedHydrogens()) {
      atoms_full[k].SetNumFixedHydrogens(0);
      atoms_full[k].SetUnitedAtomSymbol(atoms_full[k].GetElementSymbol());

    } else if (atoms_full[k].GetNumHydrogens()) {
      atoms_full[k].SetNumHydrogens(0);
      atoms_full[k].SetUnitedAtomSymbol(atoms_full[k].GetElementSymbol());
    }
  }

  for (size_t k = atoms_full.size(); k < N_full; k++)
    atoms_full.push_back(Atom("H"));

  enumerator_arguments.full_matrix_BU_ptr->SetAtomVector(atoms_full);
  enumerator_arguments.extended_matrix_BU_ptr->SetAtomVector(atoms_extended);
  size_t pseudoatom_idx(atoms_hat.size());

  if (has_pseudoatom) {
    for (size_t k = 0; k < atoms_hat.size(); k++) {
      if (atoms_hat[k].IsPseudoatom()) {
        auto psa = pseudoatoms->find(atoms_hat[k].GetElementSymbol());

        for (size_t jj = 1; jj < psa->second.GetMatrix().GetN(); jj++) {
          enumerator_arguments.full_matrix_BU_ptr->SetElement(k, jj + pseudoatom_idx - 1,
                                                              psa->second.GetMatrix().GetElement(0, jj));
          enumerator_arguments.extended_matrix_BU_ptr->SetElement(k,
                                                                  jj + pseudoatom_idx - 1,
                                                                  psa->second.GetMatrix().GetElement(0, jj));
        }

        for (size_t ii = 1; ii < psa->second.GetMatrix().GetN(); ii++) {
          for (size_t jj = ii + 1; jj < psa->second.GetMatrix().GetN(); jj++) {
            //if(psa->GetAtoms()[jj].GetElementSymbol() != "H") {
            enumerator_arguments.full_matrix_BU_ptr->SetElement(ii + pseudoatom_idx - 1,
                                                                jj + pseudoatom_idx - 1, psa->second.GetMatrix().GetElement(ii, jj));
            enumerator_arguments.extended_matrix_BU_ptr->SetElement(ii + pseudoatom_idx - 1,
                                                                    jj + pseudoatom_idx - 1, psa->second.GetMatrix().GetElement(ii, jj));
            //}
          }
        }

        pseudoatom_idx += psa->second.GetMatrix().GetN() - 1;
      }
    }
  }
}


/*
test, whether the found molecule fits the criteria of the restrictions
*/
bool Enumerator::CheckMolecule() {
  if (substructures != NULL && substructures->size()) {
    enumerator_arguments.full_matrix_ptr = std::unique_ptr<AdjacencyMatrix>
                                           (new AdjacencyMatrix(*enumerator_arguments.full_matrix_BU_ptr));
    GetFullMatrix();
  }

  //Get number of single, double, triple bonds, and number of rings in the molecule
  std::fill(ranges.begin() + 1, ranges.end(), 0);
  enumerator_arguments.A->GetNumMultipleBonds(ranges);
  ranges[range_rings] = deg_unsaturations - (ranges[range_double_bonds] + 2 *
                                             ranges[range_triple_bonds] + 3 * ranges[range_quadruple_bonds]);

  //only accept the matrix, if the bond/ring restrictions are fulfilled
  if (AreInRange(ranges, ranged_properties) && (substructures == NULL ||
                                                FindFragMatches(*enumerator_arguments.full_matrix_ptr, *substructures)))
    return true;

  else
    return false;
}

/*
generate the SMILES (and stereo SMILES if specified) and write them to the output_file
*/
void Enumerator::WriteMolecule() {
  AdjacencyMatrix* mat(NULL);
  RepresentationSystem* u0(NULL);
  RepresentationSystem* u_automorph(NULL);
  RepresentationSystem* u_id(NULL);

  if (has_pseudoatom) {
    enumerator_arguments.extended_matrix_ptr = std::shared_ptr<AdjacencyMatrix>
                                               (new AdjacencyMatrix(*enumerator_arguments.extended_matrix_BU_ptr));
    GetExtendedMatrix(u0, u_automorph, u_id);
    mat = &*enumerator_arguments.extended_matrix_ptr;

  } else {
    mat = &*enumerator_arguments.A;
    u0 = &enumerator_arguments.max_fill_alg_ptr->GetU0_non_const();
    u_automorph = &enumerator_arguments.max_fill_alg_ptr->GetUAutomorph_non_const();
    u_id = &enumerator_arguments.max_fill_alg_ptr->GetUId_non_const();
  }

  //test if this molecule is aromatic and has already been found
  bool canonical(true);

  if (ranges[range_double_bonds] >= 3 && ranges[range_rings] >= 1 &&
      FindBenzMatch(*mat, canonical, *u0) && !canonical)
    return;

  SmilesGeneratorEnu smiles_gen(*mat);
  smiles_gen.GenerateSmiles();
  num_isomers_for_formula++; //another molecule found
  num_isomers_for_listing++;
  PrintConstitutionalIsomer(smiles_gen.GetSmiles());

  if (stereo
      && (ranges[range_double_bonds] || mat->GetAtomVector().front().GetDegree() >= 4
          || (mat->GetAtomVector().front().GetDegree() == 3 &&
              mat->GetAtomVector().front().GetNumTotalHydrogens() == 1)))
    EnumerateStereoSmiles(ranges[range_rings], *u_automorph, *u_id, smiles_gen);

  else
    PrintState(num_isomers_for_formula);

  PrintClosingIsomerTag();

  if (has_pseudoatom) {
    delete u0;
    delete u_automorph;
    delete u_id;
  }
}

void Enumerator::PrintConstitutionalIsomer(const std::string& smiles) {
  output_file << std::string(xml_indent_size * 3,
                             ' ') << "<isomer isomer_id=\"" << code << "_" << std::setw(
                9) << std::setfill('0') << std::right <<
              std::to_string(num_isomers_for_listing) << "\">\n"
              /* << std::string(xml_indent_size * 3, ' ') << "<formula>" << *enumerator_arguments.formula << "</formula>\n"*/
              << std::string(xml_indent_size * 4,
                             ' ') << "<constitutional_SMILES>" << smiles << "</constitutional_SMILES>\n";
}

void Enumerator::PrintClosingIsomerTag() {
  output_file << std::string(xml_indent_size * 3, ' ') << "</isomer>\n";
}

/*****************************************************************************************************
FILL THE MATRIX full_matrix_ptr FROM THE ADJACENCY MATRIX A, S.T. ALL HYDROGEN ATOMS ARE GIVEN EXPLICITLY
*****************************************************************************************************/
void Enumerator::GetFullMatrix() {
  size_t N = enumerator_arguments.A->GetN();
  size_t H_index(N + enumerator_arguments.n_pseudoatoms);

  for (size_t i = 0; i < N; i++) {
    for (size_t j = i + 1; j < N; j++)
      enumerator_arguments.full_matrix_ptr->SetElement(i, j,
                                                       enumerator_arguments.A->GetElement(i, j));

    size_t deg = enumerator_arguments.full_matrix_ptr->AccumulateRow(i);

    for (size_t j = deg + 1;
         j <= enumerator_arguments.full_matrix_ptr->GetAtomVector()[i].GetDegree(); j++)
      enumerator_arguments.full_matrix_ptr->SetElement(i, H_index++, 1);
  }

  for (size_t i = N; i < N + enumerator_arguments.n_pseudoatoms; i++) {
    size_t deg = enumerator_arguments.full_matrix_ptr->AccumulateRow(i);

    for (size_t j = deg + 1;
         j <= enumerator_arguments.full_matrix_ptr->GetAtomVector()[i].GetDegree()/*A.GetAtomVector()[i].GetNumHydrogens()*/;
         j++)
      enumerator_arguments.full_matrix_ptr->SetElement(i, H_index++, 1);
  }
}


/************************************************
CREATE ADJACENCY MATRIX WITH EXPLICIT PSEUDOATOMS
*************************************************/
void Enumerator::GetExtendedMatrix(RepresentationSystem*& u0,
                                   RepresentationSystem*& u_automorph, RepresentationSystem*& u_id) {
  size_t N = enumerator_arguments.A->GetN();

  for (size_t i = 0; i < N; i++) {
    for (size_t j = i + 1; j < N; j++)
      enumerator_arguments.extended_matrix_ptr->SetElement(i, j,
                                                           enumerator_arguments.A->GetElement(i, j));
  }

  enumerator_arguments.extended_matrix_ptr->SortAtomVector();
  LambdaVector lambda_extended =
    enumerator_arguments.extended_matrix_ptr->GetLambda();
  size_t N_extended = enumerator_arguments.extended_matrix_ptr->GetN();
  std::vector<size_t> num_perms_extended(N_extended);
  int idx(0);

  for (size_t i = 0; i < lambda_extended.size(); i++) {
    for (size_t j = 0; j < lambda_extended[i]; j++)
      num_perms_extended[idx++] = (lambda_extended[i] - j);
  }

  u0 = (new RepresentationSystem(N_extended - 1, StabilizerVector(1)));
  u_id = (new RepresentationSystem(N_extended - 1, StabilizerVector(1)));
  u_automorph = (new RepresentationSystem(N_extended - 1, StabilizerVector(1)));

  for (size_t i = 0; i < u0->size(); i++) {
    (*u0)[i].resize(num_perms_extended[i]);
    (*u_automorph)[i].reserve(num_perms_extended[i]);
    (*u_id)[i][0] = (Permutations(1, {i, i}));
    (*u_automorph)[i][0] = (Permutations(1, {i, i}));
    (*u0)[i][0] = Permutations(1, {i, i});

    for (size_t j = 1; j < num_perms_extended[i]; j++) {
      (*u0)[i][j] = (Permutations(1, {i, i + j}));
    }
  }

  enumerator_arguments.extended_matrix_ptr->MakeCanonical(*u0);

  if (stereo) {
    //Get automorphism group of extended matrix
    PermutationIterator perm_it(*u0);

    while (perm_it.GetNextPermutation()) {
      size_t jpos;

      if (enumerator_arguments.extended_matrix_ptr->EqualTo(*
                                                            (perm_it.GetPermutedIndices()),
                                                            jpos)) {
        (*u_automorph)[perm_it.GetSmallestDiffIndex()].push_back(*
                                                                 (perm_it.GetCombinedPermutation()));
        perm_it.SetCurrentIndexToSmallestDiffIndex();

      } else {
        if (perm_it.GetCurrentIndex() > jpos)
          perm_it.SetCurrentIndex(jpos);
      }
    }
  }
}



/********************************************************************************************
USE SmilesGeneratorEnu CLASS TO GENERATE STEREO CONFORMATIONS AND PRINT THE CORRESPONDING OUTPUT
*********************************************************************************************/
void Enumerator::EnumerateStereoSmiles(const size_t n_rings,
                                       const RepresentationSystem& u_automorph,
                                       const RepresentationSystem& u_id,
                                       const SmilesGeneratorEnu& smiles_gen) {
  std::shared_ptr<AdjacencyMatrix> mat;

  if (has_pseudoatom)
    mat = enumerator_arguments.extended_matrix_ptr;

  else
    mat = enumerator_arguments.A;

  //mat->print();
  //std::cout << smilesGen.GetSmiles() << std::endl;
  //std::cout << smilesGen.GetVisitedIndices() << std::endl;
  StereoGenerator stereo_gen(*mat, n_rings, u_automorph, u_id, smiles_gen);
  stereo_gen.GenerateStereoSmiles();
  std::vector<std::tuple<std::string, int, std::pair<int, int>>> stereo_smiles_vec =
    stereo_gen.GetStereoSmiles();

  if (stereo_smiles_vec.size()) {
    output_file << std::string(xml_indent_size * 4, ' ') << "<stereo_isomers>\n"
                << std::string(xml_indent_size * 5,
                               ' ') << "<num_stereoisomers>" << stereo_smiles_vec.size() <<
                "</num_stereoisomers>\n";

    for (size_t n = 0; n < stereo_smiles_vec.size(); n++) {
      PrintState(num_isomers_for_formula + (size_t)n);
      const std::string stereo_smiles = std::get<0>(stereo_smiles_vec[n]);
      const int enantiomer = std::get<1>(stereo_smiles_vec[n]);
      const int n_tet_centers = std::get<2>(stereo_smiles_vec[n]).first;
      const int n_ct_centers = std::get<2>(stereo_smiles_vec[n]).second;
      output_file  << std::string(xml_indent_size * 5,
                                  ' ') << "<stereo_isomer stereo_id=\""
                   << code << "_" << std::setw(9) << std::setfill('0') << std::right <<
                   std::to_string(num_isomers_for_listing)
                   << "_" << std::setw(4) << std::to_string(n + 1) << "\">\n"
                   << std::string(xml_indent_size * 6,
                                  ' ') << "<stereo_SMILES>" << stereo_smiles << "</stereo_SMILES>\n"
                   << std::string(xml_indent_size * 6,
                                  ' ') <<  "<num_tetrahedral_stereocenters>" << n_tet_centers <<
                   "</num_tetrahedral_stereocenters>\n"
                   << std::string(xml_indent_size * 6,
                                  ' ') << "<num_cis_trans_stereocenters>" << n_ct_centers <<
                   "</num_cis_trans_stereocenters>\n";

      if (enantiomer != -1)
        output_file << std::string(xml_indent_size * 6, ' ') << "<enantiomer>"
                    << code << "_" << std::setw(9) << std::setfill('0') << std::right <<
                    std::to_string(num_isomers_for_listing)
                    << "_" << std::setw(4) << std::to_string(enantiomer + 1) << "</enantiomer>\n";

      output_file << std::string(xml_indent_size * 5, ' ') << "</stereo_isomer>\n";
      //output_file << formula << " " << smiles << " ";
      /*std::string smiles = std::get<0>(stereoSMILES[n]);
      int enantiomer = std::get<1>(stereoSMILES[n]);
      int numStereoCentersLocal = std::get<2>(stereoSMILES[n]).first;
      int numCTBondsLocal = std::get<2>(stereoSMILES[n]).second;*/
      /*  if(enantiomer == -1)
          output_file << "M[" << numStereoCentersLocal << "," << numCTBondsLocal << "]" << '\n';

        else
          output_file << "C[" << numStereoCentersLocal << "," << numCTBondsLocal << "]" << "(" << std::setw(4) << std::setfill('0') << std::right << enantiomer
                     << ")" << '\n';*/
    }

    output_file << std::string(xml_indent_size * 4, ' ') << "</stereo_isomers>\n";
    num_isomers_for_formula += stereo_smiles_vec.size() - 1;

  } else {
    //printIsomerToFile(smilesGen.GetSmiles(), "A");
    PrintState(num_isomers_for_formula);
  }
}

/*******************
PRINT CONSOLE OUTPUT
********************/
void Enumerator::PrintState(const size_t n_cur_isomers) {
  if (!(n_cur_isomers % 100)) {
    //std::cout << ". " << std::flush;
    //if(!(n_cur_isomers % 10000)) {
    std::cout << " -> found " << n_cur_isomers << std::flush << '\r';
    /*std::cout << '\n'
          << std::setw(3) << std::setfill('0') << n_cur_isomers/1000000000 << "'"
              << std::setw(3) << std::setfill('0') << (n_cur_isomers/1000000)%1000 << "'"
              << std::setw(3) << std::setfill('0') << (n_cur_isomers/1000)%1000 << "'"
              << std::setw(3) << std::setfill('0') << (n_cur_isomers%1000) << " ";*/
    //}
  }
}

/********************************************************************
CREATES A CANONICAL FORMULA STRING FROM A GIVEN AtomVector AND LambdaVector
BY SORTING BY ATOM PRIORITY AND EXTRACTING IMPLICIT HYDROGENS
*********************************************************************/
std::string Enumerator::CreateCanonicalFormula(const LambdaVector& l,
                                               const AtomVector<combi_ff::Atom>& used_atoms) {
  AtomVector<combi_ff::Atom> atoms_full(0);

  //extract implicit atoms pseudoatoms ot an explicit atom vector atoms_full
  for (size_t i = 0; i < l.size(); i++) {
    if (used_atoms[i].IsPseudoatom()) {
      //const Pseudoatom* psa = NULL;
      auto psa = pseudoatoms->find(used_atoms[i].GetElementSymbol());

      for (auto && a : psa->second.GetAtoms()) {
        AtomVector<combi_ff::Atom> atmp(l[i], a);
        atoms_full.insert(atoms_full.end(), atmp.begin(), atmp.end());
      }

    } else {
      AtomVector<combi_ff::Atom> atmp(l[i], used_atoms[i]);
      atoms_full.insert(atoms_full.end(), atmp.begin(), atmp.end());
    }
  }

  return CreateCanonicalFormulaFromAtomVector<combi_ff::Atom>(atoms_full);
}

} //namespace enu

} //namespace combi_ff