#include "HydrogenDistribution.h"
#include "Matrix.h"

namespace combi_ff {

namespace enu {

HydrogenDistribution::HydrogenDistribution(const LambdaVector& lambda,
                                           const AtomVector<combi_ff::Atom>& non_H_atoms, const int N_hyd)
  : i(-1),
    non_H_atoms(non_H_atoms),
    N_hyd(N_hyd),
    N_hat(non_H_atoms.size()),
    canonical(false) {
  //define lambda_bar, the LambdaVector of the non-hydrogen atoms
  //if there are any hydrogen atoms, lambda_bar is lambda without the last entry, otherwise it's just lambda
  lambda_bar = LambdaVector(N_hyd ? LambdaVector(lambda.begin(), lambda.end() - 1) : lambda);
  H_con = ConnectivityVector(N_hat, 0);
  //define degree_vec_bar, the degree vector of the non-hydrogen atoms, containing the degree of each of the atoms in non_H_atom
  degree_vec_bar = std::vector<int>(N_hat, 0);

  for (size_t i = 0; i < N_hat; i++)
    degree_vec_bar[i] = (int)non_H_atoms[i].GetDegree();

  //create supplementary vector M (analogous to matrix M of the maximum filling algorithm, but vector instead of matrix)
  //M[i] stores the maximum possible number of hydrogen atoms that atom i can be connected to
  //this is the minimum of the atom's bond degree minus one (we want a connected graph) and the total number of hydrogen atoms
  M = std::vector<int>(N_hat, 0);

  if (N_hat > 1) {
    for (size_t j = 0; j < N_hat; j++) {
      if (!non_H_atoms[j].GetHasFixedHydrogens())
        M[j] = std::min(degree_vec_bar[j] - 1, N_hyd);
    }

  } else {
    if (non_H_atoms[0].GetHasFixedHydrogens())
      return;

    else if (degree_vec_bar[0] == N_hyd)
      M[0] = degree_vec_bar[0];

    else
      return;
  }
}


/************************************************************************
Get THE LEXICOGRAPHICALLY NEXT SMALLER VALID HYDROGEN CONNECTIVITY VECTOR
************************************************************************/
bool HydrogenDistribution::GetNextDistribution() {
  canonical = false;

  //case that this is the first iteration
  if (i == -1)
    ForwardStep();

  //continue backstepping until there's either a canonical decomposition, or no other decomposition is possible
  while (!canonical && i > 0)
    BackwardStep();

  if (canonical)
    return true;

  else
    return false;
}

/*********************************************************************************************
FORWARD STEP ANALOGOUS TO THE MAXIMUM FILLING ALGORITHM, BUT WITH A VECTOR INSTEAD OF A MATRIX
*********************************************************************************************/
void HydrogenDistribution::ForwardStep() {
  //if i reaches end of H_con vector
  if (++i == (int)N_hat) {
    //if the sum over H_con is equal to the number of H atoms -> valid distribution -> test for canonicity
    if (std::accumulate(H_con.begin(), H_con.end(), 0) == N_hyd)
      CanonicityTest();

    //otherwise, backstep
    else
      BackwardStep();
  }

  //otherwise determine highest possible entry of H_con[i] and continue with forward step
  //-> highest possible value is determined as the minimum of the value of M[i] and the sum over H_con, i.e. the already distributed H atoms
  else {
    H_con[i] = std::min(M[i], N_hyd - std::accumulate(H_con.begin(), H_con.end(), 0));
    ForwardStep();
  }
}

/**********************************************************************************************
BACKWARD STEP ANALOGOUS TO THE MAXIMUM FILLING ALGORITHM, BUT WITH A VECTOR INSTEAD OF A MATRIX
**********************************************************************************************/
void HydrogenDistribution::BackwardStep() {
  //if i is back at the becinning of H_con vector, there are no more canonical decompositions
  if (i == 0)
    return;

  //check if the next lower H_con[i] can be decreased by one
  if (H_con[--i] > 0) {
    H_con[i]--;

    //if this makes H_con canonical up to i, continue with forward step
    if (LexOrderUpTo(i))
      ForwardStep();

    //otherwise perform another backward step on the same H_con[i] (increase i by one, since it will be decreased at beginning of backward step)
    else {
      i++;
      BackwardStep();
    }
  }

  //otherwise, perform a backward step for the next lower entry of H_con
  else
    BackwardStep();
}

/***************************************************
CHECK IF THE FIRST max ENTRIES IN H_con ARE CANONICAL
***************************************************/
bool HydrogenDistribution::LexOrderUpTo(const size_t max) {
  //jmin and jmax are used as the borders in H_con for a certain atom type
  size_t jmin(0), jmax;

  //go through all atom types
  for (size_t i = 0; i < lambda_bar.size(); i++) {
    //jmax is the highest index of the current atom type, described by lambda_bar[i]
    //jmin is the lowest index of the current atom type, describbed by lambda_bar[i]
    //between jmin and jmax, the elements of H_con have to be sorted in decreasing order, otherwise not canonical
    jmax = std::accumulate(lambda_bar.begin(), lambda_bar.begin() + 1 + i, 0);

    for (size_t j = jmin; j < jmax - 1 && j < max; j++) {
      if (H_con[j] < H_con[j + 1])
        return false;
    }

    //update jmin
    jmin = jmax;

    if (jmin > max)
      break;
  }

  return true;
}

/***********************
CANONICITY TEST FOR H_con
***********************/
void HydrogenDistribution::CanonicityTest() {
  //degree_vec_hat is the degree vector of atomsHat, when applying the hydrogen distribution of H_con
  DegreeVector degree_vec_hat(N_hat);

  for (size_t j = 0; j < N_hat; j++)
    degree_vec_hat[j] = degree_vec_bar[j] - H_con[j];

  //in order to establish canonicity, confirm that
  // - the sum over all hydrogen atoms distributed through H_con is equal to the total number of hydrogen atoms
  // - the sum of the degree of the united atoms is even (since nbonds = 1/2 sum (degrees) and we can't have half bonds)
  // - the values of the entries in H_con are decreasing within each lambda in lambda_bar
  if (std::accumulate(H_con.begin(), H_con.end(), 0) == N_hyd
      && std::accumulate(degree_vec_hat.begin(), degree_vec_hat.end(), 0) % 2 == 0
      && LexOrder())
    canonical = true;

  return;
}


/************************************************************
TESTS IF H_con IS LEXICOGRAPHICALLY ORDERED FOR ALL ATOM TYPES
************************************************************/
bool HydrogenDistribution::LexOrder() {
  //jmin and jmax are used as the borders in H_con for a certain atom type
  size_t jmin(0), jmax;

  for (size_t i = 0; i < lambda_bar.size(); i++) {
    //jmax is the highest index of the current atom type, described by lambda_bar[i]
    //jmin is the lowest index of the current atom type, describbed by lambda_bar[i]
    //between jmin and jmax, the elements of H_con have to be sorted in decreasing order, otherwise not canonical
    jmax = std::accumulate(lambda_bar.begin(), lambda_bar.begin() + 1 + i, 0);

    for (size_t j = jmin; j < jmax - 1; j++) {
      if (H_con[j] < H_con[j + 1])
        return false;
    }

    //update jmin;
    jmin = jmax;
  }

  return true;
}

const ConnectivityVector& HydrogenDistribution::GetHcon() {
  return H_con;
}

} //namespace enu

} //namespace combi_ff