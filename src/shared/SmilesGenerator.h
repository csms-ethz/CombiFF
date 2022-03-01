// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#ifndef SMILESGENERATOR_H
#define SMILESGENERATOR_H

#include <math.h>

#include "Matrix.h"

namespace combi_ff {

struct SmilesBlock {
  std::string opening_braces{""};
  std::string closing_braces{""};
  std::string element_name{""};
  std::string formal_charge{""};
  std::string bond_type{""};
  // ring_indices: first element is bond degree, 2nd element is ring index (e.g.
  // C=1 -> ("=",1))
  std::vector<std::pair<std::string, size_t>> ring_indices;

  size_t size() {
    return opening_braces.size() + closing_braces.size() + element_name.size() +
           formal_charge.size() + bond_type.size() + ring_indices.size();
  }
  void Print() const {
    std::cout << opening_braces << element_name << formal_charge << bond_type;

    for (auto&& ri : ring_indices)
      std::cout << "(" << ri.first << "," << ri.second << ")";

    std::cout << closing_braces;
  }
};

template <typename T, typename AtomClass>
class SmilesGenerator {
 public:
  SmilesGenerator(const combi_ff::AdjacencyMatrix<T, AtomClass>& A);
  void GenerateSmiles();
  void GenerateSmilesNonCanon();
  std::vector<size_t> Partition();
  size_t GetInv(const size_t idx);
  std::string CreateSmiles();
  std::vector<size_t> Sort(std::vector<size_t>& curr_partition);
  void Tiebreak(std::vector<size_t>& curr_partition);
  void Refine();
  void SortRingNodes(
      std::vector<std::tuple<size_t, size_t, size_t>>& ring_nodes,
      std::vector<size_t>& node_visited_at_position);
  void TestValidity(const std::string& smiles, const size_t& nArom,
                    const size_t& ring);
  const std::string GetSmiles() const;
  const std::vector<size_t>& GetVisitedIndices() const;
  const std::vector<size_t>& GetComingFrom() const;
  const std::vector<std::vector<size_t>>& GetGoingTo() const;
  const std::vector<std::vector<size_t>>& GetRingConnections() const;
  const std::vector<SmilesBlock>& GetSmilesBlocks() const;
  static bool IsLarger(const size_t i, const size_t j) {
    if (i > j)
      return true;

    else
      return false;
  }

 private:
  const combi_ff::AdjacencyMatrix<T, AtomClass>& A;
  const size_t N;
  const std::vector<bool>& is_aromatic_carbon;
  std::vector<std::vector<size_t>> partition_classes;
  std::vector<size_t> indices;
  std::vector<SmilesBlock> smiles_blocks;
  std::vector<size_t> visited_indices;
  std::vector<size_t> coming_from;
  std::vector<std::vector<size_t>> going_to;
  std::vector<std::vector<size_t>> ring_connections;
  std::string smiles;

  std::vector<size_t> keys;
  std::vector<size_t> partitions_to_sort;
  std::vector<size_t> changed_indices;
  std::vector<size_t> affected_partition_classes;
  std::vector<std::vector<size_t>> sorted_neighbors;

  std::vector<bool> changed;
};

typedef SmilesGenerator<size_t, combi_ff::Atom> SmilesGeneratorEnu;
typedef SmilesGenerator<double, combi_ff::CnvAtom> SmilesGeneratorCnv;

template <typename T, typename AtomClass>
SmilesGenerator<T, AtomClass>::SmilesGenerator(
    const combi_ff::AdjacencyMatrix<T, AtomClass>& A)
    : A(A),
      N(A.GetN()),
      is_aromatic_carbon(A.GetIsAromaticCarbonConst()),
      partition_classes(std::vector<std::vector<size_t>>(N)),
      smiles_blocks(std::vector<SmilesBlock>(N)),
      visited_indices(std::vector<size_t>(0)),
      coming_from(std::vector<size_t>(N, 0)),
      going_to(std::vector<std::vector<size_t>>(N, std::vector<size_t>(0))),
      ring_connections(
          std::vector<std::vector<size_t>>(N, std::vector<size_t>(0))),
      smiles(std::string("")) {
  partitions_to_sort.reserve(N);
  keys.reserve(N);
  changed_indices.reserve(N);
  affected_partition_classes.reserve(N);
  affected_partition_classes.reserve(N);
  changed_indices.reserve(N);
  changed = std::vector<bool>(N);
}

template <typename T, typename AtomClass>
void SmilesGenerator<T, AtomClass>::GenerateSmiles() {
  indices = Partition();
  bool done = false;

  for (auto&& p : partition_classes) p.reserve(N);

  while (!done) {
    Refine();
    done = true;

    for (auto&& p : partition_classes) {
      if (p.size() != 1) {
        done = false;
        break;
      }
    }
  }

  for (size_t i = 0; i < N; i++) indices[i] = partition_classes[i].front();

  smiles = CreateSmiles();
}

template <typename T, typename AtomClass>
void SmilesGenerator<T, AtomClass>::GenerateSmilesNonCanon() {
  indices = std::vector<size_t>(N);
  std::iota(indices.begin(), indices.end(), 0);

  smiles = CreateSmiles();
}

template <typename T, typename AtomClass>
std::vector<size_t> SmilesGenerator<T, AtomClass>::Partition() {
  std::vector<size_t> invariants(N);

  for (size_t i = 0; i < N; i++) invariants[i] = (GetInv(i));

  std::vector<size_t> keys(invariants);
  std::sort(keys.begin(), keys.end());
  std::vector<size_t>::iterator it = std::unique(keys.begin(), keys.end());
  keys.resize(std::distance(keys.begin(), it));
  std::vector<size_t> indices(N, 1);
  int ni;
  int idx = 0;

  for (size_t i = 0; i < keys.size(); i++) {
    ni = 0;

    for (size_t j = 0; j < N; j++) {
      if (invariants[j] == keys[i]) {
        ni++;
        indices[j] = idx;
      }
    }

    idx = idx + ni;
  }

  for (size_t i = 0; i < partition_classes.size(); i++)
    partition_classes[indices[i]].push_back(i);

  for (size_t i = 0; i < partition_classes.size(); i++) {
  }

  return indices;
}
template <typename T, typename AtomClass>
size_t SmilesGenerator<T, AtomClass>::GetInv(const size_t idx) {
  size_t inv(0);
  const combi_ff::Atom& atom = A.GetAtomVector()[idx];
  inv += 10000000 * atom.GetNumConnections();

  for (auto&& i : atom.GetNeighbours()) {
    if (A.GetElement(idx, i)) inv += 100000;
  }

  inv += 1000 * atom.GetElementNumber();
  inv += (atom.GetNumHydrogens() + atom.GetNumFixedHydrogens());
  return inv;
}
template <typename T, typename AtomClass>
std::string SmilesGenerator<T, AtomClass>::CreateSmiles() {
  std::string smiles;
  std::stack<size_t> stack;
  size_t nArom(0);
  // indices: indices[i] = j means that atom j has priority i, aka that at
  // priority i is atom j

  for (size_t i = 0; i < indices.size(); i++) {
    combi_ff::ElementSymbol symbol =
        A.GetAtomVector()[indices[i]].GetElementSymbol();

    if (A.GetAtomVector()[indices[i]].GetFormalCharge().size())
      smiles_blocks[indices[i]].formal_charge =
          A.GetAtomVector()[indices[i]].GetFormalCharge();

    if (is_aromatic_carbon[indices[i]]) {
      nArom++;
      smiles_blocks[indices[i]].element_name = symbol;
      std::transform(smiles_blocks[indices[i]].element_name.begin(),
                     smiles_blocks[indices[i]].element_name.end(),
                     smiles_blocks[indices[i]].element_name.begin(), tolower);

    } else
      smiles_blocks[indices[i]].element_name = symbol;

    if (A.GetAtomVector()[indices[i]].GetHydrogenInSmiles()) {
      if (A.GetAtomVector()[indices[i]].GetNumTotalHydrogens()) {
        smiles_blocks[indices[i]].element_name.insert(0, "[");
        size_t nH = A.GetAtomVector()[indices[i]].GetNumTotalHydrogens();
        smiles_blocks[indices[i]].element_name +=
            "H" + std::to_string(nH) + "]";

      } else {
        smiles_blocks[indices[i]].element_name.insert(0, "[");
        smiles_blocks[indices[i]].element_name += "]";
      }
    }

    smiles_blocks[indices[i]].ring_indices =
        std::vector<std::pair<std::string, size_t>>(0);
  }

  // std::cout << &smiles_blocks << std::endl;
  stack.push(indices[0]);
  std::vector<bool> written(indices.size(), false);
  size_t ring = 0;
  std::vector<size_t> node_visited_at_position(indices.size(), -1);
  // std::vector<std::vector<size_t>> going_to(indices.size(),
  // std::vector<size_t>(0));
  std::vector<std::tuple<size_t, size_t, size_t>> ring_nodes;
  size_t currId(0);
  size_t idx = 0;
  visited_indices.resize(N);

  while (!stack.empty()) {
    size_t i = stack.top();
    stack.pop();

    if (!written[i]) {
      // std::cout << "current node " << i << std::endl;
      visited_indices[idx++] = i;
      node_visited_at_position[i] = currId++;
      written[i] = true;
      bool found = false;
      SmilesBlock& name = smiles_blocks[i];

      if (i != indices[0]) going_to[coming_from[i]].push_back(i);

      for (int j = (int)indices.size() - 1; j >= 0; j--) {
        if (A.GetElement(i, indices[j])) {
          // std::cout << "  neighbour " << indices.at(j) << std::endl;
          // openIdx[i].insert(openIdx[i].begin(), j);
          if (!written[indices[j]]) {
            if (found) {
              // std::cout << i << " " << indices.at(j) << std::endl;
              // sif(smiles_blocks[indices[j]].opening_braces == "")
              smiles_blocks[indices[j]].opening_braces = "(";
            }

            found = true;
            // std::cout << "visiting " << indices.at(j) << std::endl;
            stack.push(indices[j]);

            if (A.GetElement(i, indices[j]) == 2 &&
                !(is_aromatic_carbon[indices[j]] && is_aromatic_carbon[i]))
              smiles_blocks[indices[j]].bond_type = "=";

            else if (A.GetElement(i, indices[j]) == 3)
              smiles_blocks[indices[j]].bond_type = "#";

            else if (A.GetElement(i, indices[j]) == 4)
              smiles_blocks[indices[j]].bond_type = "$";

            // std::cout << " . coming from at " << indices[j] << " is " << i <<
            // std::endl;
            coming_from[indices[j]] = i;
            // std::cout << &coming_from << std::endl;

          } else if (coming_from[i] != indices[j]) {
            ring++;

            if (ring > 99)
              throw std::runtime_error(
                  "can't handle more than 99 rings in a SMILES string");

            ring_nodes.push_back(
                std::tuple<size_t, size_t, size_t>(indices[j], i, ring));

            if (A.GetElement(i, indices[j]) == 1 ||
                (double)A.GetElement(i, indices[j]) == 1.5) {
              smiles_blocks[indices[j]].ring_indices.push_back(
                  std::pair<std::string, size_t>("", ring));
              name.ring_indices.push_back(
                  std::pair<std::string, size_t>("", ring));
              name.bond_type = "";

            } else if (A.GetElement(i, indices[j]) == 2) {
              if (is_aromatic_carbon[indices[j]])
                smiles_blocks[indices[j]].ring_indices.push_back(
                    std::pair<std::string, size_t>("", ring));

              else
                smiles_blocks[indices[j]].ring_indices.push_back(
                    std::pair<std::string, size_t>("=", ring));

              name.bond_type = "";
              name.ring_indices.push_back(
                  std::pair<std::string, size_t>("", ring));

            } else if (A.GetElement(i, indices[j]) == 3) {
              smiles_blocks[indices[j]].ring_indices.push_back(
                  std::pair<std::string, size_t>(
                      "#", ring));  // += '#' + std::to_string(ring);
              name.bond_type = "";
              name.ring_indices.push_back(std::pair<std::string, size_t>(
                  "",
                  ring));  /// += std::to_string(ring);

            } else if (A.GetElement(i, indices[j]) == 4) {
              smiles_blocks[indices[j]].ring_indices.push_back(
                  std::pair<std::string, size_t>(
                      "$", ring));  // += '$' + std::to_string(ring);
              name.bond_type = "";
              name.ring_indices.push_back(std::pair<std::string, size_t>(
                  "",
                  ring));  // += std::to_string(ring);
            }

            if (A.GetElement(coming_from[i], i) == 2 &&
                (!is_aromatic_carbon[i] || !is_aromatic_carbon[coming_from[i]]))
              name.bond_type = "=";

            else if (A.GetElement(coming_from[i], i) == 3)
              name.bond_type = "#";

            else if (A.GetElement(coming_from[i], i) == 4)
              name.bond_type = "$";
          }
        }
      }

      if (!found) name.closing_braces += ")";

      // smiles_blocks[i] = name;
    }
  }

  smiles_blocks[visited_indices.back()].closing_braces = "";

  if (ring > 0) {
    std::vector<size_t> chronological_ring_indices(ring + 1, 0);
    size_t idx(1);

    for (size_t i = 0; i < visited_indices.size() - 1; i++) {
      if (going_to[visited_indices[i]].size() >= 1)
        smiles_blocks[going_to[visited_indices[i]].back()].opening_braces = "";
    }

    SortRingNodes(ring_nodes, node_visited_at_position);

    for (auto&& rn : ring_nodes)
      chronological_ring_indices[std::get<2>(rn)] = idx++;

    for (auto&& n : smiles_blocks) {
      for (auto&& nn : n.ring_indices)
        nn.second = (chronological_ring_indices[nn.second]);

      for (size_t i = 0; i < n.ring_indices.size(); i++) {
        for (size_t j = i + 1; j < n.ring_indices.size(); j++) {
          if (n.ring_indices[i].second > n.ring_indices[j].second) {
            auto temp = n.ring_indices[i];
            n.ring_indices[i] = n.ring_indices[j];
            n.ring_indices[j] = temp;
          }
        }
      }
    }
  }

  for (auto&& rn : ring_nodes) {
    ring_connections[std::get<0>(rn)].push_back(std::get<1>(rn));
    ring_connections[std::get<1>(rn)].push_back(std::get<0>(rn));
  }

  for (auto&& v : visited_indices) {
    auto&& n = smiles_blocks[v];
    smiles += n.opening_braces + n.bond_type;

    if (n.formal_charge.size())
      smiles += "[" + n.element_name + n.formal_charge + "]";

    else
      smiles += n.element_name;

    // Hack for explicit hydrogen
    //      for(size_t ii = 0; ii < A.GetAtomVector()[v].GetNumTotalHydrogens();
    //      ii++){
    //        smiles += "(H)";
    //      }
    // end hack

    for (auto&& ri : n.ring_indices) {
      smiles += ri.first;

      if (ri.second <= 9)
        smiles += std::to_string(ri.second);

      else
        smiles += "%" + std::to_string(ri.second);
    }

    smiles += n.closing_braces;
  }

  // test the validity of the smiles string . can be commented out at some point
  // testValidity(smiles, nArom, ring);
  return smiles;
}
template <typename T, typename AtomClass>
std::vector<size_t> SmilesGenerator<T, AtomClass>::Sort(
    std::vector<size_t>& curr_partition) {
  sorted_neighbors.resize(curr_partition.size());
  const combi_ff::AtomVector<AtomClass>& atoms = A.GetAtomVector();
  size_t idx;

  for (size_t i = 0; i < curr_partition.size(); i++) {
    idx = curr_partition[i];
    sorted_neighbors[i].resize(atoms[idx].GetNeighbours().size());
    const auto& curr_neighbor_vector = atoms[idx].GetNeighbours();

    // for(auto && j : atoms[idx].GetNeighbours())
    auto neighbor_it = curr_neighbor_vector.begin();
    for (size_t j = 0; j < sorted_neighbors[i].size(); j++) {
      sorted_neighbors[i][j] = (indices[*neighbor_it++]);
    }

    std::sort(sorted_neighbors[i].begin(), sorted_neighbors[i].end(), IsLarger);
  }

  size_t n = curr_partition.size();
  changed_indices.resize(N);
  std::fill(changed.begin(), changed.end(), false);
  idx = 0;

  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      if (sorted_neighbors[i] < sorted_neighbors[j]) {
        indices[curr_partition[j]]++;

        if (!changed[curr_partition[j]]) {
          changed[curr_partition[j]] = true;
          changed_indices[idx++] = (curr_partition[j]);
        }

      } else if (sorted_neighbors[j] < sorted_neighbors[i]) {
        indices[curr_partition[i]]++;

        if (!changed[curr_partition[i]]) {
          changed[curr_partition[i]] = true;
          changed_indices[idx++] = (curr_partition[i]);
        }
      }
    }
  }

  changed_indices.resize(idx);
  return changed_indices;
}
template <typename T, typename AtomClass>
void SmilesGenerator<T, AtomClass>::Tiebreak(
    std::vector<size_t>& curr_partition) {
  indices[curr_partition.back()] += (curr_partition.size() - 1);
  partition_classes[indices[curr_partition.back()]].push_back(
      curr_partition.back());
  curr_partition.pop_back();
}
template <typename T, typename AtomClass>
void SmilesGenerator<T, AtomClass>::Refine() {
  keys.assign(indices.begin(), indices.end());
  partitions_to_sort.resize(0);
  std::sort(keys.begin(), keys.end());
  std::vector<size_t>::iterator it = std::unique(keys.begin(), keys.end());
  keys.resize(std::distance(keys.begin(), it));

  for (size_t i = 0; i < N; i++) {
    if (partition_classes[i].size() > 1) partitions_to_sort.push_back(i);
  }

  auto&& atoms = A.GetAtomVector();

  while (!partitions_to_sort.empty()) {
    std::vector<size_t>& curr_partition =
        partition_classes[partitions_to_sort.back()];
    partitions_to_sort.pop_back();
    // sort the current partitions and store which ones were changed in
    // changed_indices
    Sort(curr_partition);
    affected_partition_classes.resize(0);

    for (auto&& idx : changed_indices) {
      for (auto&& j : atoms[idx].GetNeighbours())
        affected_partition_classes.push_back(indices[j]);
    }

    std::sort(affected_partition_classes.begin(),
              affected_partition_classes.end());
    auto&& it = std::unique(affected_partition_classes.begin(),
                            affected_partition_classes.end());
    affected_partition_classes.resize(
        std::distance(affected_partition_classes.begin(), it));

    for (size_t j = 0; j < partition_classes.size(); j++)
      partition_classes[j].resize(0);

    for (size_t j = 0; j < indices.size(); j++)
      partition_classes[indices[j]].push_back(j);

    for (size_t j = 0; j < affected_partition_classes.size(); j++) {
      if (partition_classes[affected_partition_classes[j]].size() > 1)
        partitions_to_sort.push_back(affected_partition_classes[j]);
    }

    std::sort(partitions_to_sort.begin(), partitions_to_sort.end());
    it = std::unique(partitions_to_sort.begin(), partitions_to_sort.end());
    partitions_to_sort.resize(std::distance(partitions_to_sort.begin(), it));
  }

  for (int i = (int)partition_classes.size() - 1; i >= 0; i--) {
    std::vector<size_t>& curr_partition = partition_classes[i];

    if (curr_partition.size() > 1) {
      Tiebreak(curr_partition);
      break;
    }
  }
}
template <typename T, typename AtomClass>
void SmilesGenerator<T, AtomClass>::SortRingNodes(
    std::vector<std::tuple<size_t, size_t, size_t>>& ring_nodes,
    std::vector<size_t>& node_visited_at_position) {
  for (size_t i = 0; i < ring_nodes.size(); i++) {
    for (size_t j = i + 1; j < ring_nodes.size(); j++) {
      if (node_visited_at_position[std::get<0>(ring_nodes[i])] >
          node_visited_at_position[std::get<0>(ring_nodes[j])]) {
        std::tuple<size_t, size_t, size_t> temp = ring_nodes[i];
        ring_nodes[i] = ring_nodes[j];
        ring_nodes[j] = temp;
        // std::swap(ring_nodes.begin() + i, ring_nodes.begin() +  j);

      } else if (node_visited_at_position[std::get<0>(ring_nodes[i])] ==
                 node_visited_at_position[std::get<0>(ring_nodes[j])]) {
        if (node_visited_at_position[std::get<1>(ring_nodes[i])] >
            node_visited_at_position[std::get<1>(ring_nodes[j])]) {
          std::tuple<size_t, size_t, size_t> temp = ring_nodes[i];
          ring_nodes[i] = ring_nodes[j];
          ring_nodes[j] = temp;
        }
      }
    }
  }
}
template <typename T, typename AtomClass>
void SmilesGenerator<T, AtomClass>::TestValidity(const std::string& smiles,
                                                 const size_t& nArom,
                                                 const size_t& ring) {
  {
    std::stack<size_t> braceOpen;
    size_t nD(nArom / 2), nT(0), nQ(0);

    for (auto&& v : visited_indices) {
      if (smiles_blocks[v].opening_braces != "") braceOpen.push(1);

      if (smiles_blocks[v].closing_braces != "") {
        if (braceOpen.empty()) {
          std::cout << "\nbrace mismatch in SMILES: " << smiles << "\n";
          abort();

        } else
          braceOpen.pop();
      }

      if (smiles_blocks[v].bond_type != "") {
        if (smiles_blocks[v].bond_type == "=")
          nD++;

        else if (smiles_blocks[v].bond_type == "#")
          nT++;

        else
          nQ++;
      }

      for (auto&& r : smiles_blocks[v].ring_indices) {
        if (r.first != "") {
          if (r.first == "=")
            nD++;

          else if (r.first == "#")
            nT++;

          else
            nQ++;
        }
      }
    }

    /*size_t nSingle(0), nDouble(0), nTriple(0), nQuadruple(0);
    A.GetNumMultipleBonds(nSingle, nDouble, nTriple, nQuadruple);

    if(nD != nDouble || nT != nTriple || nQ != nQuadruple || ring != nRings) {
        std::cout << "\nwrong number of double/triple/quadruple bonds in SMILES:
    " << smiles << "\n"
                  << "number of double bonds in Matrix:     " << nDouble << "\n"
                  << "number of double bonds in smiles:     " << nD << "\n"
                  << "number of triple bonds in Matrix:     " << nTriple << "\n"
                  << "number of double bonds in smiles:     " << nT << "\n"
                  << "number of quadruple bonds in Matrix:  " << nQuadruple <<
    "\n"
                  << "number of quadruple bonds in Matrix:  " << nQ << "\n"
                  << "number of rings in Matrix:            " << nRings << "\n"
                  << "number of rings in smiles:            " << ring << "\n"
                  << "Adjacency Matrix:\n";
        A.print();
        throw std::runtime_error("smiles " + smiles + " not valid");
    }

    if(!braceOpen.empty())
        throw std::runtime_error("\nbrace mismatch in SMILES: " + smiles);*/
  }
}
template <typename T, typename AtomClass>
const std::string SmilesGenerator<T, AtomClass>::GetSmiles() const {
  return smiles;
}
template <typename T, typename AtomClass>
const std::vector<size_t>& SmilesGenerator<T, AtomClass>::GetVisitedIndices()
    const {
  return visited_indices;
}
template <typename T, typename AtomClass>
const std::vector<size_t>& SmilesGenerator<T, AtomClass>::GetComingFrom()
    const {
  return coming_from;
}
template <typename T, typename AtomClass>
const std::vector<std::vector<size_t>>&
SmilesGenerator<T, AtomClass>::GetGoingTo() const {
  return going_to;
}
template <typename T, typename AtomClass>
const std::vector<std::vector<size_t>>&
SmilesGenerator<T, AtomClass>::GetRingConnections() const {
  return ring_connections;
}
template <typename T, typename AtomClass>
const std::vector<SmilesBlock>& SmilesGenerator<T, AtomClass>::GetSmilesBlocks()
    const {
  return smiles_blocks;
}

}  // namespace combi_ff

#endif
