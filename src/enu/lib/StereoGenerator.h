// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef STEREOGENERATOR_H_
#define STEREOGENERATOR_H_

#include "AdjacencyMatrixEnu.h"
#include "SmilesGenerator.h"

namespace combi_ff {

namespace enu {

typedef std::vector<int> Config;

// enum {clockwise, counterclockwise, notTetStereo = -1} tetrahedralOrientation;
// enum {trans, cis, notETStereo = -1} CTdirection;

struct CisTransPair {
  size_t idx1;
  size_t idx2;
  size_t nbr1;
  size_t nbr2;

  CisTransPair(size_t idx1, size_t idx2, size_t nbr1, size_t nbr2)
      : idx1(idx1), idx2(idx2), nbr1(nbr1), nbr2(nbr2) {}

  CisTransPair() : idx1(0), idx2(0), nbr1(0), nbr2(0) {}
};

struct Configuration {
  Config valid_config_tetra;
  Config valid_config_ct;

  Configuration() : valid_config_tetra(Config(0)), valid_config_ct(Config(0)) {}
  Configuration(const Config& valid_config_tetra, const Config& valid_config_ct)
      : valid_config_tetra(valid_config_tetra),
        valid_config_ct(valid_config_ct) {}
};

class StereoGenerator {
 public:
  StereoGenerator(const enu::AdjacencyMatrix& A, const size_t n_rings,
                  const RepresentationSystem& u_automorph,
                  const RepresentationSystem& u_id,
                  const SmilesGeneratorEnu& smiles_gen);

  void GenerateStereoSmiles();
  void FindPotentialStereoCenters();
  void TestIfPotentialTetrahedralStereoCenter(const Atom& a, const size_t idx);
  void TestIfPotentialCisTransStereoCenter(const Atom& a, const size_t idx);
  void FindTrueStereoCenters();
  void FindPotentialParaStereoCenters();
  void WriteStereoSmiles(std::vector<Configuration>& validConfigurations);

  size_t GetNumStereoCenters() const;
  size_t GetNumCTBonds() const;

  void FindShortestPath(std::vector<size_t>& minimalPath, size_t i,
                        size_t j) const;
  bool IsSmaller(Config& perm, Config& orig) const;

  void FindValidConfigurations(std::vector<Configuration>& validConfigurations);
  void FindValidTrueCTConfigurations(
      std::vector<Config>& validConfigurationsCT);
  void FindValidTrueTetrahedralConfigurations(
      std::vector<Config>& validConfigurationsTetra);
  void AddValidTrueConfigurations(
      std::vector<Configuration>& validConfigurations,
      std::vector<Config>& validConfigurationsTetra,
      std::vector<Config>& validConfigurationsCT);
  void FindValidParaConfigurations(
      std::vector<Configuration>& validConfigurations);

  bool NeighborsArePermuted(const Atom& a,
                            const std::vector<size_t>& permutedIndices);

  std::vector<std::tuple<std::string, int, std::pair<int, int>>>&
  GetStereoSmiles();

  void NeighborOrder(const size_t idx, const size_t idxPerm,
                     const std::vector<size_t>& permutedIndices,
                     std::vector<size_t>& neighborsOrigOrder,
                     std::vector<size_t>& neighborsPermOrder) const;

 private:
  const enu::AdjacencyMatrix& A;
  const size_t N;
  const RepresentationSystem& u_automorph;
  const RepresentationSystem& u_id;
  size_t num_stereo_centers;
  size_t num_ct_bonds;
  const std::vector<size_t>& visited_indices;

  const std::vector<SmilesBlock>& smiles_blocks;
  const std::vector<size_t>& coming_from;
  const std::vector<std::vector<size_t>>& going_to;
  const std::vector<std::vector<size_t>>& ring_connections;
  std::vector<bool> is_in_cycle;
  std::vector<size_t> potential_stereo_idx;
  std::vector<bool> true_stereo;
  std::vector<bool> potential_cis_trans_stereo;
  std::vector<size_t> potential_cis_trans_stereo_idx;
  std::vector<size_t> true_cis_trans_stereo_idx;
  std::vector<int> double_bond_partner;
  std::vector<int> single_bond_partner;
  std::vector<size_t> true_stereo_idx;
  std::vector<size_t> potential_para_stereo_idx;
  std::vector<size_t> potential_para_cis_trans_stereo_idx;
  std::vector<CisTransPair> true_cis_trans_pairs;
  std::vector<CisTransPair> para_cis_trans_pairs;

  std::vector<int> position_in_true_stereo_idx;
  std::vector<bool> true_cis_trans_stereo;
  std::vector<size_t> visited_idx;

  // store stereo smiles, whether the stereoisomer is an enantiomer, and the
  // number of stereocenters and stereo CT bonds
  std::vector<std::tuple<std::string, int, std::pair<int, int>>> stereo_smiles;
};

Config operator++(Config& b);

}  // namespace enu

}  // namespace combi_ff

#endif
