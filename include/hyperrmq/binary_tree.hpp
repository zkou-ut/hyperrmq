#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace hyperrmq {

struct BinaryTree {
    int64_t None = -1;

    int64_t n;
    std::vector<int64_t> parent;
    std::vector<int64_t> left;
    std::vector<int64_t> right;
    int64_t root;

    BinaryTree();

    explicit BinaryTree(const std::vector<int64_t> parent_inorder);

    std::string to_BP_string(std::string left_delim = "(",
                             std::string middle_delim = ")",
                             std::string right_delim = "");

    bool operator==(const BinaryTree &other) const {
        return (this->n == other.n) && (this->parent == other.parent);
    }
};

}  // namespace hyperrmq
