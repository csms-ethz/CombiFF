// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef HYDROGENDISTRIBUTION_H
#define HYDROGENDISTRIBUTION_H

#include "Atom.h"
#include "LambdaVector.h"

namespace combi_ff {

namespace enu {

class HydrogenDistribution {
 public:
  HydrogenDistribution(const LambdaVector& lambda,
                       const AtomVector<combi_ff::Atom>& non_H_atoms,
                       const int NH);

  bool GetNextDistribution();
  void CanonicityTest();
  void BackwardStep();
  void ForwardStep();

  bool LexOrder();
  bool LexOrderUpTo(const size_t max);
  const ConnectivityVector& GetHcon();

 private:
  int i;
  ConnectivityVector H_con;
  std::vector<int> M;
  LambdaVector lambda_bar;
  const AtomVector<combi_ff::Atom>& non_H_atoms;
  std::vector<int> degree_vec_bar;
  const int N_hyd;
  const size_t N_hat;
  bool canonical;
};

}  // namespace enu

}  // namespace combi_ff

#endif
