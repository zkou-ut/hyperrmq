#pragma once

#include <cstdint>
#include <vector>

#include "hyperrmq/bit_array.hpp"

namespace hyperrmq {

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

    uint64_t naive_fwdsearch(uint64_t j, int64_t d) const;
    int64_t naive_bwdsearch(uint64_t j, int64_t d) const;
    int64_t naive_minexcess(uint64_t j, int64_t k) const;

    uint64_t naive_open(uint64_t index) const;
    uint64_t naive_close(uint64_t index) const;

    uint64_t naive_lca(uint64_t u_inorder, uint64_t v_inorder) const;

    uint64_t evaluate_memory_consumption() const;
};

template <typename T>
TreeBP cartesian_tree_bp(const std::vector<T>& values);

}  // namespace hyperrmq
