#pragma once

#include <cstdint>
#include <vector>

#include "binary_tree.hpp"

namespace average_case_optimal_rmq {

struct CartesianTree {
    int32_t n;
    BinaryTree binary_tree;

    explicit CartesianTree(const std::vector<int32_t> &values);
};

}  // namespace average_case_optimal_rmq
