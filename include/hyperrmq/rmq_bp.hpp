#pragma once

#include "hyperrmq/rmm_tree.hpp"

namespace hyperrmq {

template <uint32_t c = 8, uint32_t b = 1024>
struct RMQBP {
    uint32_t num_of_nodes;
    RMMTree<c, b> rmm_tree;

    RMQBP();
    explicit RMQBP(const std::vector<int32_t>& values);

    uint32_t query(uint32_t i, uint32_t j) const;

    uint64_t evaluate_memory_consumption() const;
};

template <uint32_t c, uint32_t b>
RMQBP<c, b>::RMQBP() : num_of_nodes(0) {}

template <uint32_t c, uint32_t b>
RMQBP<c, b>::RMQBP(const std::vector<int32_t>& values)
    : num_of_nodes(values.size()), rmm_tree(cartesian_tree_bp(values)) {}

template <uint32_t c, uint32_t b>
uint32_t RMQBP<c, b>::query(uint32_t i, uint32_t j) const {
    assert(0 <= i && i <= j && j < num_of_nodes);
    return rmm_tree.rank1(
        rmm_tree.rmq(rmm_tree.select1(i), rmm_tree.select1(j) + 1));
}

template <uint32_t c, uint32_t b>
uint64_t RMQBP<c, b>::evaluate_memory_consumption() const {
    return rmm_tree.evaluate_memory_consumption();
}

}  // namespace hyperrmq
