#pragma once

#include <cstdint>
#include <vector>

#include "hyperrmq/binary_tree.hpp"

namespace hyperrmq {

struct CartesianTree {
    int32_t n;
    BinaryTree binary_tree;

    explicit CartesianTree(const std::vector<int32_t> &values);
};

}  // namespace hyperrmq
