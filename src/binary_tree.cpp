#include "hyperrmq/binary_tree.hpp"

namespace hyperrmq {

BinaryTree::BinaryTree() : n(0), root(None) {}

BinaryTree::BinaryTree(const std::vector<int32_t> parent_inorder)
    : parent(parent_inorder), root(None) {
    n = parent.size();
    left.resize(n, None);
    right.resize(n, None);
    for (int i = 0; i < n; i++) {
        if (parent[i] == None) {
            root = i;
        } else {
            if (i < parent[i]) {
                left[parent[i]] = i;
            } else {
                right[parent[i]] = i;
            }
        }
    }
}

std::string BinaryTree::to_BP_string(std::string left_delim,
                                     std::string middle_delim,
                                     std::string right_delim) {
    auto bp = [&](auto self, int32_t v) -> std::string {
        std::string result = "";
        result += left_delim;
        if (left[v] != None) result += self(self, left[v]);
        result += middle_delim;
        if (right[v] != None) result += self(self, right[v]);
        result += right_delim;
        return result;
    };
    return root == None ? "" : bp(bp, root);
}

}  // namespace hyperrmq
