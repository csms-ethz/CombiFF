#include "smilesToMatrix.h"
#include "SmilesGenerator.h"
#include <algorithm>

namespace combi_ff {


combi_ff::AdjacencyMatrix<double, combi_ff::CnvAtom> ConvertToMatrix(
  const std::string& smiles) {
  combi_ff::AtomVector<combi_ff::CnvAtom> atoms(0);
  size_t num_atoms(0);
  size_t num_rings(0);

  if (smiles.size() == 1 && !isalpha(smiles[0]))
    throw std::runtime_error("can't convert a smiles string that only consists of one symbol and doesn't contain an atom: "
                             + smiles);

  std::vector<bool> is_aromatic(0);

  if (!GetBasicInformationFromSmiles(smiles, atoms, num_atoms, num_rings,
                                     is_aromatic))
    return combi_ff::AdjacencyMatrix<double, combi_ff::CnvAtom>(0);

  std::vector<combi_ff::SmilesBlock> names(num_atoms);
  CreateNameVector(names, smiles, atoms, num_atoms);
  combi_ff::AdjacencyMatrix<double, combi_ff::CnvAtom> A(atoms.size());
  A.SetAtomVector(atoms);
  A.SetIsAromaticCarbon(is_aromatic);
  CreateMatrixFromSmilesBlocks(A, smiles, names, num_rings);
  return A;
}

bool  GetBasicInformationFromSmiles(const std::string& smiles,
                                    combi_ff::AtomVector<combi_ff::CnvAtom>& atoms, size_t& num_atoms,
                                    size_t& num_rings,
                                    std::vector<bool>& is_aromatic) {
  size_t n_triple(0), n_quadruple(0), degree(0);
  double n_double(0);

  for (size_t i = 0; i < smiles.size(); i++) {
    if (smiles[i] == '=')
      n_double += 1;

    else if (smiles[i] == '#')
      n_triple++;

    else if (smiles[i] == '$')
      n_quadruple++;

    else if (smiles[i] == '%') {
      //assert(i + 2 < smiles.size());
      num_rings++;
      i += 2;

    } else if (std::isdigit(smiles[i]))
      num_rings++;

    else if (std::isupper(smiles[i])) {
      if (i > 1 && (smiles[i - 2] == '[' || smiles[i - 1] == '@')) {
        continue; // don't add hydrogens e.g. in [SH2] or [C@H]
      }

      try {
        if (i + 1 < smiles.size() && (smiles[i + 1] == 'l' || smiles[i + 1] == 'r' ||
                                      smiles[i + 1] == 'g' || smiles[i + 1] == 'e' ||
                                      smiles[i + 1] == 'd' || smiles[i + 1] == 't' ||
                                      smiles[i + 1] == 'u' || smiles[i + 1] == 'a')) {
          std::string s = smiles.substr(i, 2);
          //assert(s == "Br" || s == "Cl");
          atoms.push_back(combi_ff::CnvAtom(s));
          degree += atoms.back().GetDegree();
          i++;

        } else {
          atoms.push_back(combi_ff::CnvAtom(std::string(1, smiles[i])));
          degree += atoms.back().GetDegree();
        }
      }

      //catch incorrect user input -> unknown atom types
      catch (combi_ff::input_error& e) {
        std::cerr << "!Input Error: " << e.what() << std::endl;
        return false;
      }

      is_aromatic.push_back(false);
      num_atoms++;

    } else if (std::islower(smiles[i])) {
      num_atoms++;
      std::string s("");
      s += (smiles[i]);

      if (i + 1 < smiles.size() && smiles.substr(i, 2) == "se") {
        s += smiles[i + 1];
        i++;
      }

      std::transform(s.begin(), s.begin() + 1, s.begin(), toupper);
      is_aromatic.push_back(true);

      try {
        atoms.push_back(combi_ff::CnvAtom(s));
      }

      //catch incorrect user input -> unknown atom types
      catch (combi_ff::input_error& e) {
        std::cerr << "!Input Error: " << e.what() << std::endl;
        return false;
      }

      if (i > 1 && (smiles[i - 2] == '[' || smiles[i - 1] == '@')) {
        continue; // don't add hydrogens e.g. in [SH2] or [C@H]
      }

      n_double += 0.5;
      degree += atoms.back().GetDegree();
    }
  }
  
  //assert(nR % 2 == 0);
  /*num_rings /= 2;
  int n_unsaturations = (int)num_rings + (int)n_double + 2 * (int)n_triple + 3 *
                        (int)n_quadruple;

  if (n_unsaturations - 1 + (int)atoms.size() - 0.5 * (double)degree > 0)
    std::cerr <<
              "?Warning: nUnsat larger than zero in smilesToMatrix transformation of " <<
              smiles << " (unsat= " << n_unsaturations << " #atoms = " <<
              atoms.size() << " deg = " << degree << ")" << '\n';
  */
  return true;
}

void CreateNameVector(std::vector<combi_ff::SmilesBlock>& names,
                      const std::string& smiles, combi_ff::AtomVector<combi_ff::CnvAtom>& atoms,
                      size_t& num_atoms) {
  size_t curr_atom(0);

  for (size_t i = 0; i < smiles.size(); i++) {
    if (smiles[i] == '=' && smiles[i + 1] != '%' && !std::isdigit(smiles[i + 1]))
      names[curr_atom].bond_type = "=";

    else if (smiles[i] == '#' && smiles[i + 1] != '%' &&
             !std::isdigit(smiles[i + 1]))
      names[curr_atom].bond_type = "#";

    else if (smiles[i] == '$' && smiles[i + 1] != '%' &&
             !std::isdigit(smiles[i + 1]))
      names[curr_atom].bond_type = "$";

    else if (smiles[i] == '%') {
      //assert(i + 2 < smiles.size());
      if (smiles[i - 1] == '=')
        names[curr_atom - 1].ring_indices.emplace_back("=%",
                                                       std::stoi(smiles.substr(i + 1, 2)));

      else if (smiles[i - 1] == '#')
        names[curr_atom - 1].ring_indices.emplace_back("#%" ,
                                                       std::stoi(smiles.substr(i + 1, 2)));

      else if (smiles[i - 1] == '$')
        names[curr_atom - 1].ring_indices.emplace_back("$%" ,
                                                       std::stoi(smiles.substr(i + 1, 2)));

      else
        names[curr_atom - 1].ring_indices.emplace_back("%",
                                                       std::stoi(smiles.substr(i + 1, 2)));

      i += 2;

    } else if (std::isdigit(smiles[i])) {
      if (smiles[i - 1] == '=')
        names[curr_atom - 1].ring_indices.emplace_back("=" , std::stoi(smiles.substr(i,
                                                                                     1)));

      else if (smiles[i - 1] == '#')
        names[curr_atom - 1].ring_indices.emplace_back("#" , std::stoi(smiles.substr(i,
                                                                                     1)));

      else if (smiles[i - 1] == '$')
        names[curr_atom - 1].ring_indices.emplace_back("$" , std::stoi(smiles.substr(i,
                                                                                     1)));

      else
        names[curr_atom - 1].ring_indices.emplace_back("", std::stoi(smiles.substr(i,
                                                                                   1)));

    } else if (std::isupper(smiles[i])) {
      std::string s("");

      if (i + 1 < smiles.size() && (smiles[i + 1] == 'l' || smiles[i + 1] == 'r')) {
        s = smiles.substr(i, 2);
        //assert(s == "Br" || s == "Cl");
        i++;

      } else
        s += smiles[i];

      names[curr_atom].element_name = s;
      curr_atom++;

    } else if (std::islower(smiles[i])) {
      std::string s("");
      s += (smiles[i]);
      names[curr_atom].element_name = s;
      curr_atom++;

    } else if (smiles[i] == '/' || smiles[i] == '\\')
      std::cerr << "?Warning: removing stereo information\n";

    else if (smiles[i] == '(') {
      if (curr_atom >= names.size())
        throw combi_ff::input_warning("encountered unexpected symbol " + std::string(1,
                                                                                     smiles[i]) + " at position " + std::to_string(i) + " in " + smiles);

      names[curr_atom].opening_braces += '(';

    } else if (smiles[i] == ')')
      names[curr_atom - 1].closing_braces += ')';

    else if (smiles[i] == '[') {
      if (++i == smiles.size())
        throw combi_ff::input_warning("unclosed [ encountered in " + smiles);

      if (std::isupper(smiles[i])) {
        //names[curr_atom].element_name = '[';
        std::string s("");

        if ((i + 1 < smiles.size()) && (smiles[i + 1] == 'l' || smiles[i + 1] == 'r' ||
                                        smiles[i + 1] == 'e')) {
          s = smiles.substr(i, 2);
          //assert(s == "Br" || s == "Cl");
          i++;

        } else
          s += smiles[i];

        names[curr_atom].element_name += s;

        if (i + 1 < smiles.size() && smiles[i + 1] != '@')
          names[curr_atom].element_name.insert(0, "[");

        curr_atom++;

      } else if (std::islower(smiles[i])) {
        names[curr_atom].element_name = smiles[i];

        if (i + 1 < smiles.size() && smiles[i + 1] != '@')
          names[curr_atom].element_name.insert(0, "[");

        curr_atom++;

      } else {
        throw combi_ff::input_warning("encountered unexpected symbol " + std::string(1,
                                                                                     smiles[i]) + " at position " + std::to_string(i) + " in " + smiles);
      }

      if (++i == smiles.size())
        throw combi_ff::input_warning("unclosed [ encountered in " + smiles);

      if (smiles[i] == '@') {
        std::cerr << "?Warning: removing stereo information\n";
        i++;

        if (i < smiles.size() && smiles[i] == '@')
          i++;

      } else if (smiles[i] == '+' || smiles[i] == '-') {
        if (names[curr_atom - 1].element_name.front() != '[')
          names[curr_atom - 1].element_name.insert(0, "[");

        names[curr_atom - 1].element_name += smiles[i];
        names[curr_atom - 1].formal_charge += smiles[i];
        atoms[curr_atom - 1].SetFormalCharge(atoms[curr_atom - 1].GetFormalCharge() +
                                             smiles[i]);
        i++;

        if (i < smiles.size() && std::isdigit(smiles[i])) {
          names[curr_atom - 1].element_name += smiles[i];
          names[curr_atom - 1].formal_charge += smiles[i];
          atoms[curr_atom - 1].SetFormalCharge(atoms[curr_atom - 1].GetFormalCharge() +
                                               smiles[i]);
          i++;
        }
      }

      if (smiles[i] == 'H') {
        if (++i == smiles.size())
          throw combi_ff::input_warning("unclosed [ encountered in " + smiles);

        if (names[curr_atom - 1].element_name.front() == '[') {
          size_t num_h = 1;

          if (std::isdigit(smiles[i])) {
            num_h = std::stoi(std::string(1, smiles[i]));

            if (++i == smiles.size())
              throw combi_ff::input_warning("unclosed [ encountered in " + smiles);
          }

          names[curr_atom - 1].element_name += "H" + std::to_string(num_h);
          atoms[curr_atom - 1].SetDegree(atoms[curr_atom - 1].GetDegree() + num_h);
          atoms[curr_atom - 1].SetNumHydrogens(num_h);
          atoms[curr_atom - 1].SetHydrogensInSmiles(true);
        }

        if (smiles[i] != ']')
          throw combi_ff::input_warning("encountered " + std::string(1,
                                                                     smiles[i]) + " at position " + std::to_string(i) + " in " + smiles +
                                        ". but expected ]") ;

        else if (names[curr_atom - 1].element_name.front() == '[')
          names[curr_atom - 1].element_name += ']';

      } else if (smiles[i] == ']') {
        if (names[curr_atom - 1].element_name.front() == '[')
          names[curr_atom - 1].element_name += ']';

      } else {
        throw combi_ff::input_warning("encountered unexpected symbol " + std::string(1,
                                                                                     smiles[i]) + " at position " + std::to_string(i) + " in " + smiles);
      }

    } else if ((smiles[i] == '=' || smiles[i] == '#') && (isdigit(smiles[i + 1]) ||
                                                          smiles[i + 1] == '%')) {
      // nothing is to be done, this is just s.t. there is no error recognized in the final else clause!
      // DO NOT REMOVE
    } else if (smiles[i] == '-') {
      // ignore explicit single bonds
    } else if (smiles[i] == '@' || smiles[i] == '\\' || smiles[i] == '/')
      std::cerr << "?Warning: removing stereo information\n";

    else
      throw std::runtime_error("encountered unexpected symbol " + std::string(1,
                                                                              smiles[i]) + " at position " + std::to_string(i) + " in " + smiles);
  }

  std::string reconstructedSmiles("");

  for (auto && nn : names) {
    //assert(nn.opening_braces.size() <= 1 && nn.closing_braces.size() <= 1);
    reconstructedSmiles += nn.opening_braces + nn.bond_type + nn.element_name;

    for (auto && r : nn.ring_indices)
      reconstructedSmiles += r.first + std::to_string(r.second);

    reconstructedSmiles += nn.closing_braces;
  }

  // check if reconstructed smiles is equal to original one (ignore explicit single bonds)
  //std::string smiles_backup = smiles;
  //smiles_backup.erase(std::remove(smiles_backup.begin(), smiles_backup.end(),
  //                                '-'), smiles_backup.end());

  if (reconstructedSmiles != smiles)
    std::cerr << ("?Warning: reconstructed smiles: \'" + reconstructedSmiles +
                  "\' is different from original smiles \'" + smiles + "\'\n");

//remove redundant braces
  std::stack<size_t> braceStack;

  for (size_t ii = 0; ii < names.size(); ii++) {
    for (size_t jj = 0; jj < names[ii].opening_braces.size(); jj++)
      braceStack.push(ii);

    if (names[ii].closing_braces.size() > 1) {
      for (size_t jj = 1 ; jj < names[ii].closing_braces.size(); jj++) {
        names[braceStack.top()].opening_braces =
          names[braceStack.top()].opening_braces.substr(0,
                                                        names[braceStack.top()].opening_braces.size() - 1);
        braceStack.pop();
      }

      names[ii].closing_braces = names[ii].closing_braces[0];

    } else if (names[ii].closing_braces.size() == 1)
      braceStack.pop();
  }

  /*for(auto && n : names)
    n.print();

  std::cout << std::endl;*/
}

void CreateMatrixFromSmilesBlocks(
  combi_ff::AdjacencyMatrix<double, combi_ff::CnvAtom>& A,
  const std::string& smiles, const std::vector<combi_ff::SmilesBlock>& names,
  const size_t nR) {
  std::vector<std::pair<int, int>> ring_indices(nR, std::pair<int, int>(-1, -1));
  std::vector<int> ringClosureBondDeg(nR, -1);
  std::stack<size_t> atmIdx;

  if (names[0].ring_indices.size()) {
    size_t r;

    for (size_t j = 0; j < names[0].ring_indices.size(); j++) {
      auto&& closure = names[0].ring_indices[j];
      r = closure.second - 1 ;
      ring_indices[r] = std::pair<int, int> (0, j);
    }
  }

  if (names[0].opening_braces != "") {
    //std::cout << "adding " << 0 << std::endl;
    atmIdx.push(0);
  }

  for (size_t i = 1; i < names.size(); i++) {
    auto&& name = names[i];

    if (name.opening_braces == "" && names[i - 1].closing_braces == "") {
      if (A.GetElement(i - 1, i))
        throw std::runtime_error("already Set value of " + std::to_string(
                                   i) + " and " + std::to_string(i - 1));

      if ((std::islower(name.element_name.front()) ||
           (name.element_name.front() == '[' && std::islower(name.element_name[1])))
          && (std::islower(names[i - 1].element_name.front()) ||
              (names[i - 1].element_name.front() == '[' &&
               std::islower(names[i - 1].element_name[1]))))
        A.SetElement(i - 1, i, 1.5);

      else if (name.bond_type == "")
        A.SetElement(i - 1, i, 1);

      else if (name.bond_type == "=")
        A.SetElement(i - 1, i, 2);

      else if (name.bond_type == "#")
        A.SetElement(i - 1, i, 3);

      else if (name.bond_type == "$")
        A.SetElement(i - 1, i, 4);

    } else if (name.opening_braces == "" && names[i - 1].closing_braces != "") {
      size_t nbr = atmIdx.top();

      if (A.GetElement(i, nbr))
        throw std::runtime_error("already Set value of " + std::to_string(
                                   i) + " and " + std::to_string(nbr));

      if (std::islower(name.element_name.front()) &&
          std::islower(names[nbr].element_name.front()))
        A.SetElement(i, nbr, 1.5);

      else if (name.bond_type == "")
        A.SetElement(i, nbr, 1);

      else if (name.bond_type == "=")
        A.SetElement(i, nbr, 2);

      else if (name.bond_type == "#")
        A.SetElement(i, nbr, 3);

      else if (name.bond_type == "$")
        A.SetElement(i, nbr, 4);

      //encountered the final branch of an atom (because previous closing braces are ')')
      if (!atmIdx.size())
        throw std::runtime_error("cannot pop atmIdx");

      else
        atmIdx.pop();

    } else if (name.opening_braces != "" && names[i - 1].closing_braces == "") {
      //std::cout << "adding " << i - 1 << std::endl;
      atmIdx.push(i - 1);

      if (A.GetElement(i - 1, i))
        throw std::runtime_error("already Set value of " + std::to_string(
                                   i) + " and " + std::to_string(i - 1));

      if (std::islower(name.element_name.front()) &&
          std::islower(names[i - 1].element_name.front()))
        A.SetElement(i - 1, i, 1.5);

      else if (name.bond_type == "")
        A.SetElement(i - 1, i, 1);

      else if (name.bond_type == "=")
        A.SetElement(i - 1, i, 2);

      else if (name.bond_type == "#")
        A.SetElement(i - 1, i, 3);

      else if (name.bond_type == "$")
        A.SetElement(i - 1, i, 4);

    } else if (name.opening_braces != "" && names[i - 1].closing_braces != "") {
      size_t nbr = atmIdx.top();
      //std::cout << "connecting to " << atmIdx.top() << std::endl;

      if (A.GetElement(i - 1, i))
        throw std::runtime_error("already Set value of " + std::to_string(
                                   i) + " and " + std::to_string(nbr));

      if (std::islower(name.element_name.front()) &&
          std::islower(names[nbr].element_name.front()))
        A.SetElement(i, nbr, 1.5);

      else if (name.bond_type == "")
        A.SetElement(i, nbr, 1);

      else if (name.bond_type == "=")
        A.SetElement(i, nbr, 2);

      else if (name.bond_type == "#")
        A.SetElement(i, nbr, 3);

      else if (name.bond_type == "$")
        A.SetElement(i, nbr, 4);
    }

    if (name.ring_indices.size()) {
      size_t r;

      for (size_t j = 0; j < name.ring_indices.size(); j++) {
        auto&& closure = name.ring_indices[j];
        r = closure.second - 1;

        if (ring_indices[r].first != -1) {
          if (std::islower(names[ring_indices[r].first].element_name.front()) &&
              std::islower(name.element_name.front()))
            A.SetElement(ring_indices[r].first, i, 1.5);

          else if (closure.first[0] == '=' ||
                   names[ring_indices[r].first].ring_indices[ring_indices[r].second].first[0] ==
                   '=')
            A.SetElement(ring_indices[r].first, i, 2);

          else if (closure.first[0] == '#' ||
                   names[ring_indices[r].first].ring_indices[ring_indices[r].second].first[0] ==
                   '#')
            A.SetElement(ring_indices[r].first, i, 3);

          else if (closure.first[0] == '$' ||
                   names[ring_indices[r].first].ring_indices[ring_indices[r].second].first[0] ==
                   '$')
            A.SetElement(ring_indices[r].first, i, 4);

          else
            A.SetElement(ring_indices[r].first, i, 1);

          ring_indices[r] = std::pair<int, int>(-1, -1);

        } else
          ring_indices[r] = std::pair<int, int> (i, j);
      }
    }
  }

  for (size_t i = 0; i < A.GetN(); i++) {
    for (size_t j = i + 1; j < A.GetN(); j++) {
      if (A.GetElement(i, j) == 1.5 && !A.AreInSameCycle(i, j))
        A.SetElement(i, j, 1);
    }
  }

  if (names.back().closing_braces != "")
    atmIdx.pop();

  if (atmIdx.size()) {
    while (atmIdx.size()) {
      std::cerr << atmIdx.top() << " ";
      atmIdx.pop();
    }

    std::cerr << std::endl;
    throw combi_ff::input_warning("braces mismatch in SMILES");
  }

  for (size_t i = 0; i < A.GetAtomVector().size(); i++) {
    double rowSum = A.AccumulateRow(i);

    if (int(rowSum * 10) % 10) {
      //std::cerr << "?Warning: bond degree of atom " << std::to_string(
        //          i) << " in " << smiles << " is not an integer (ok for fused aromatic rings)\n";
    }

    //assert(int(rowSum * 10) % 10 == 0);
    //assert(rowSum <= (A.GetAtomVector())[i].GetDegree());
    if (rowSum > (double)A.GetAtomVector()[i].GetDegree() &&
        !A.GetAtom(i).GetFormalCharge().size()) {
      //A.GetAtom(i).SetDegree((size_t)
      //                       rowSum); // e.g. for S, where different valences are possible
    }

    if ((double)A.GetAtomVector()[i].GetDegree() >
        rowSum && !A.GetAtom(
          i).GetHydrogenInSmiles()/*&& A.GetAtomVector()[i].GetElementSymbol() == "C"*/) // add implicit hydrogens
      (A.GetAtom(i)).SetNumHydrogens(size_t((double)A.GetAtomVector()[i].GetDegree() -
                                            rowSum));
  }
}

} //namespace combi_ff
