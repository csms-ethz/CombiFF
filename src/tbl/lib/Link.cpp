// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include "Link.h"

namespace combi_ff {

namespace topology_builder {

HalfLink::HalfLink(size_t atom_index, size_t core_index, size_t fragment_index,
                   const std::string& linkatom_name)
    : atom_index(atom_index),
      core_index(core_index),
      fragment_index(fragment_index),
      linkatom_name(linkatom_name) {}

size_t HalfLink::GetAtomIndex() const { return atom_index; }

size_t HalfLink::GetCoreIndex() const { return core_index; }

size_t HalfLink::GetFragmentIndex() const { return fragment_index; }

std::string HalfLink::GetLinkatomName() const { return linkatom_name; }

Link::Link(const HalfLink& half_link_1, const HalfLink& half_link_2)
    : half_link_1(half_link_1), half_link_2(half_link_2) {}

const HalfLink* Link::GetHalfLink1() const { return &half_link_1; }
const HalfLink* Link::GetHalfLink2() const { return &half_link_2; }

}  // namespace topology_builder

}  // namespace combi_ff
