#pragma once

#include "hyperrmq/rmm_tree.hpp"

namespace hyperrmq {

template <uint64_t c = 8, uint64_t b = 1024>
struct RMQBP {
    uint64_t num_of_nodes;
    RMMTree<c, b> rmm_tree;

    RMQBP();
    explicit RMQBP(const std::vector<int64_t>& values);

    uint64_t query(uint64_t i, uint64_t j) const;

    uint64_t evaluate_memory_consumption() const;
};

template <uint64_t c, uint64_t b>
RMQBP<c, b>::RMQBP() : num_of_nodes(0) {}

template <uint64_t c, uint64_t b>
RMQBP<c, b>::RMQBP(const std::vector<int64_t>& values)
    : num_of_nodes(values.size()), rmm_tree(cartesian_tree_bp(values)) {}

template <uint64_t c, uint64_t b>
uint64_t RMQBP<c, b>::query(uint64_t i, uint64_t j) const {
    assert(0 <= i && i <= j && j < num_of_nodes);
    return rmm_tree.rank1(
        rmm_tree.rmq(rmm_tree.select1(i), rmm_tree.select1(j) + 1));
}

template <uint64_t c, uint64_t b>
uint64_t RMQBP<c, b>::evaluate_memory_consumption() const {
    return rmm_tree.evaluate_memory_consumption();
}

}  // namespace hyperrmq
