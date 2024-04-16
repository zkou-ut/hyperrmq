
#include "hyperrmq/cartesian_tree.hpp"

#include "gtest/gtest.h"

namespace hyperrmq {
namespace {

TEST(CartesianTreeTest, FiveNodes) {
    CartesianTree cartesian_tree({3, 2, 4, 1, 5});
    BinaryTree binary_tree({1, 3, 1, -1, 3});
    ASSERT_EQ(cartesian_tree.binary_tree, binary_tree);
}

TEST(CartesianTreeTest, SameValues) {
    CartesianTree cartesian_tree({3, 3, 3});
    BinaryTree binary_tree({-1, 0, 1});
    ASSERT_EQ(cartesian_tree.binary_tree, binary_tree);
}

TEST(CartesianTreeTest, TenNodes) {
    CartesianTree cartesian_tree({1, 7, 3, 2, 4, 5, 0, 9, 8, 6});
    BinaryTree binary_tree({6, 2, 3, 0, 3, 4, -1, 8, 9, 6});
    ASSERT_EQ(cartesian_tree.binary_tree, binary_tree);
}

}  // namespace
}  // namespace hyperrmq
