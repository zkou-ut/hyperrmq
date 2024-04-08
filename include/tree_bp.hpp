#pragma once

#include <cstdint>
#include <vector>

#include "bit_array.hpp"

namespace average_case_optimal_rmq {

struct TreeBP {
    uint64_t n;
    BitArray bp;

    TreeBP();
    explicit TreeBP(uint64_t n, const BitArray& bp);
    explicit TreeBP(uint64_t n, BitArray&& bp);
    explicit TreeBP(const std::string& bp_str);

    bool operator==(const TreeBP& rhs) const;
    bool operator!=(const TreeBP& rhs) const;

    bool operator<(const TreeBP& rhs) const;

    std::string to_string() const;

    uint32_t naive_fwdsearch(uint32_t j, int32_t d) const;
    int32_t naive_bwdsearch(uint32_t j, int32_t d) const;
    int32_t naive_minexcess(uint32_t j, int32_t k) const;

    uint32_t naive_open(uint32_t index) const;
    uint32_t naive_close(uint32_t index) const;

    uint32_t naive_lca(uint32_t u_inorder, uint32_t v_inorder) const;

    uint64_t evaluate_memory_consumption() const;
};

TreeBP cartesian_tree_bp(const std::vector<int32_t>& values);

}  // namespace average_case_optimal_rmq
