#ifndef LINK_H
#define LINK_H

#include <cstdlib>
#include <string>

namespace combi_ff {

namespace topology_builder {

/*
!!!!!!!!!!!!!!!!!!!!fragment_index is one higher than the C++ index of the fragment in fragments!!!!!!!!
i.e. for frag[0] the fragment_index is one!!!!
*/

class HalfLink {

 public:

  HalfLink() = default;
  HalfLink(size_t atom_index, size_t core_index, size_t fragment_index, std::string linkatom_name);

  size_t GetAtomIndex() const;
  size_t GetCoreIndex() const;
  size_t GetFragmentIndex() const;
  std::string GetLinkatomName() const;

 private:
  size_t atom_index {0};
  size_t core_index {0};
  size_t fragment_index {0};
  std::string linkatom_name {""};
};

class Link {
 public:
  Link() = default;
  Link(HalfLink half_link_1, HalfLink half_link_2);

  const HalfLink* GetHalfLink1() const;
  const HalfLink* GetHalfLink2() const;

 private:
  HalfLink half_link_1 {HalfLink()};
  HalfLink half_link_2 {HalfLink()};
};

} //namespace topology_builder

} //namespace combi_ff


#endif
