#include "hyperrmq/binary_tree.hpp"

#include <cstdint>

#include "gtest/gtest.h"

namespace average_case_optimal_rmq {
namespace {

TEST(BinaryTreeTest, ZeroNodes) {
    BinaryTree binary_tree;
    ASSERT_EQ(binary_tree.left, std::vector<int32_t>({}));
    ASSERT_EQ(binary_tree.right, std::vector<int32_t>({}));
    ASSERT_EQ(binary_tree.to_BP_string("(", ")", ""), "");
}

TEST(BinaryTreeTest, OneNode) {
    BinaryTree binary_tree({-1});
    ASSERT_EQ(binary_tree.left, std::vector<int32_t>({-1}));
    ASSERT_EQ(binary_tree.right, std::vector<int32_t>({-1}));
    ASSERT_EQ(binary_tree.to_BP_string("(", ")", ""), "()");
}

TEST(BinaryTreeTest, TwoNodesLeft) {
    BinaryTree binary_tree({1, -1});
    ASSERT_EQ(binary_tree.left, std::vector<int32_t>({-1, 0}));
    ASSERT_EQ(binary_tree.right, std::vector<int32_t>({-1, -1}));
    ASSERT_EQ(binary_tree.to_BP_string("(", ")", ""), "(())");
}

TEST(BinaryTreeTest, TwoNodesRight) {
    BinaryTree binary_tree({-1, 0});
    ASSERT_EQ(binary_tree.left, std::vector<int32_t>({-1, -1}));
    ASSERT_EQ(binary_tree.right, std::vector<int32_t>({1, -1}));
    ASSERT_EQ(binary_tree.to_BP_string("(", ")", ""), "()()");
}

TEST(BinaryTreeTest, FiveNodes) {
    BinaryTree binary_tree({1, 3, 1, -1, 3});
    ASSERT_EQ(binary_tree.left, std::vector<int32_t>({-1, 0, -1, 1, -1}));
    ASSERT_EQ(binary_tree.right, std::vector<int32_t>({-1, 2, -1, 4, -1}));
    ASSERT_EQ(binary_tree.to_BP_string("(", ")", ""), "((())())()");
}

TEST(BinaryTreeTest, Equal) {
    BinaryTree binary_tree_1({1, 3, 1, -1, 3});
    BinaryTree binary_tree_2({1, 3, 1, -1, 3});
    BinaryTree binary_tree_3({1, 3, 1, 4, -1});
    BinaryTree binary_tree_4({1, 3, 1, -1});
    ASSERT_TRUE(binary_tree_1 == binary_tree_2);
    ASSERT_FALSE(binary_tree_1 == binary_tree_3);
    ASSERT_FALSE(binary_tree_1 == binary_tree_4);
}

}  // namespace
}  // namespace average_case_optimal_rmq
