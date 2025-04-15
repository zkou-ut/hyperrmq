#pragma once

#include <cstdint>
#include <vector>

#include "hyperrmq/binary_tree.hpp"

namespace hyperrmq {

struct CartesianTree {
    int64_t n;
    BinaryTree binary_tree;

    explicit CartesianTree(const std::vector<int64_t> &values);
};

}  // namespace hyperrmq
