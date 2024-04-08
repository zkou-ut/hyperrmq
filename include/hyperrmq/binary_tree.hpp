#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace average_case_optimal_rmq {

struct BinaryTree {
    int32_t None = -1;

    int32_t n;
    std::vector<int32_t> parent;
    std::vector<int32_t> left;
    std::vector<int32_t> right;
    int32_t root;

    BinaryTree();

    explicit BinaryTree(const std::vector<int32_t> parent_inorder);

    std::string to_BP_string(std::string left_delim = "(",
                             std::string middle_delim = ")",
                             std::string right_delim = "");

    bool operator==(const BinaryTree &other) const {
        return (this->n == other.n) && (this->parent == other.parent);
    }
};

}  // namespace average_case_optimal_rmq
