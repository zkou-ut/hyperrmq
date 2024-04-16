#include "hyperrmq/hypersuccinct_binary_tree.hpp"

#include <random>

#include "gtest/gtest.h"
#include "hyperrmq/binary_tree.hpp"
#include "hyperrmq/cartesian_tree.hpp"
#include "hyperrmq/rmq_bp.hpp"
#include "hyperrmq/runs.hpp"

namespace hyperrmq {
namespace {

TEST(HypersuccinctBinaryTreeTest, ScanBPByIncSmall) {
    TreeBP tree("(())(())()");
    HypersuccinctBinaryTree<> hs(tree, 2);

    auto node = HypersuccinctBinaryTree<>::Node(&hs, 0, 0, false);
    for (int i = 0; i < tree.bp.size(); i++) {
        ASSERT_TRUE(node.is_valid());
        ASSERT_EQ(node.access(), tree.bp.get(i));
        node.inc();
    }
    ASSERT_FALSE(node.is_valid());
}

TEST(HypersuccinctBinaryTreeTest, ScanBPByIncLarge) {
    const int n = 1e3;
    std::mt19937 engine(0);
    std::vector<int32_t> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    HypersuccinctBinaryTree<> hs(tree, 5);

    auto node = HypersuccinctBinaryTree<>::Node(&hs, 0, 0, false);
    for (int i = 0; i < tree.bp.size(); i++) {
        ASSERT_TRUE(node.is_valid());
        ASSERT_EQ(node.access(), tree.bp.get(i));
        node.inc();
    }
    ASSERT_FALSE(node.is_valid());
}

TEST(HypersuccinctBinaryTreeTest, ScanBPByDecSmall) {
    TreeBP tree("(())(())()");
    HypersuccinctBinaryTree<> hs(tree, 2);

    auto node = hs.inorder_to_node(hs.num_of_nodes - 1);
    for (int i = tree.bp.size() - 1; i >= 0; i--) {
        ASSERT_TRUE(node.is_valid());
        ASSERT_EQ(node.access(), tree.bp.get(i));
        node.dec();
    }
    ASSERT_FALSE(node.is_valid());
}

TEST(HypersuccinctBinaryTreeTest, ScanBPByDecLarge) {
    const int n = 1e3;
    std::mt19937 engine(0);
    std::vector<int32_t> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    HypersuccinctBinaryTree<> hs(tree, 5);

    auto node = hs.inorder_to_node(hs.num_of_nodes - 1);
    for (int i = tree.bp.size() - 1; i >= 0; i--) {
        ASSERT_TRUE(node.is_valid());
        ASSERT_EQ(node.access(), tree.bp.get(i));
        node.dec();
    }
    ASSERT_FALSE(node.is_valid());
}

void tree_queries_stress_test(const BinaryTree &expected,
                              const HypersuccinctBinaryTree<> &actual) {
    ASSERT_EQ(expected.root, actual.node_to_inorder(actual.root()));

    for (int i = 0; i < expected.n; i++) {
        ASSERT_EQ(
            expected.parent[i],
            actual.node_to_inorder(actual.parent(actual.inorder_to_node(i))));
        ASSERT_EQ(expected.left[i], actual.node_to_inorder(actual.left_child(
                                        actual.inorder_to_node(i))));
        ASSERT_EQ(expected.right[i], actual.node_to_inorder(actual.right_child(
                                         actual.inorder_to_node(i))));

        ASSERT_EQ(expected.left[i] == -1 && expected.right[i] == -1,
                  actual.is_leaf(actual.inorder_to_node(i)));

        int expected_child_label =
            (i == expected.root                       ? -1
             : i == expected.left[expected.parent[i]] ? 0
                                                      : 1);
        ASSERT_EQ(expected_child_label,
                  actual.child_label(actual.inorder_to_node(i)));

        int expected_leftmost_desc = i;
        while (expected.left[expected_leftmost_desc] != -1) {
            expected_leftmost_desc = expected.left[expected_leftmost_desc];
        }
        ASSERT_EQ(expected_leftmost_desc,
                  actual.node_to_inorder(
                      actual.leftmost_desc(actual.inorder_to_node(i))));

        int expected_rightmost_desc = i;
        while (expected.right[expected_rightmost_desc] != -1) {
            expected_rightmost_desc = expected.right[expected_rightmost_desc];
        }
        ASSERT_EQ(expected_rightmost_desc,
                  actual.node_to_inorder(
                      actual.rightmost_desc(actual.inorder_to_node(i))));

        int expected_subtree_size =
            expected_rightmost_desc - expected_leftmost_desc + 1;
        ASSERT_EQ(expected_subtree_size,
                  actual.subtree_size(actual.inorder_to_node(i)));

        if (expected.parent[i] != -1) {
            ASSERT_TRUE(
                actual.is_ancestor(actual.inorder_to_node(expected.parent[i]),
                                   actual.inorder_to_node(i)));
            ASSERT_FALSE(
                actual.is_ancestor(actual.inorder_to_node(i),
                                   actual.inorder_to_node(expected.parent[i])));
        }
    }
};

void is_ancestor_lca_stress_test(const RMQBP<> &rmq,
                                 const HypersuccinctBinaryTree<> &actual,
                                 int q = 1000, int seed = 0) {
    std::mt19937 engine(seed);

    auto lca = [&](int u, int v) {
        return rmq.query(std::min(u, v), std::max(u, v));
    };
    for (int i = 0; i < q; i++) {
        auto u = engine() % actual.num_of_nodes;
        auto v = engine() % actual.num_of_nodes;
        auto node_u = actual.inorder_to_node(u);
        auto node_v = actual.inorder_to_node(v);
        ASSERT_EQ(lca(u, v) == u, actual.is_ancestor(node_u, node_v));
        ASSERT_EQ(lca(u, v),
                  actual.node_to_inorder(actual.lca(node_u, node_v)));
    }
};

TEST(HypersuccinctBinaryTreeTest, OneNode) {
    BinaryTree expected({-1});
    HypersuccinctBinaryTree<> actual(TreeBP(expected.to_BP_string()), 2);

    tree_queries_stress_test(expected, actual);
}

TEST(HypersuccinctBinaryTreeTest, LeftRootPath) {
    BinaryTree expected({-1, 0, 1, 2, 3, 4});
    HypersuccinctBinaryTree<> actual(TreeBP(expected.to_BP_string()), 2);

    tree_queries_stress_test(expected, actual);
}

TEST(HypersuccinctBinaryTreeTest, RightRootPath) {
    BinaryTree expected({1, 2, 3, 4, 5, -1});
    HypersuccinctBinaryTree<> actual(TreeBP(expected.to_BP_string()), 2);

    tree_queries_stress_test(expected, actual);
}

TEST(HypersuccinctBinaryTreeTest, FiveNodes) {
    BinaryTree expected({1, 3, 1, -1, 3});
    HypersuccinctBinaryTree<> actual(TreeBP(expected.to_BP_string()), 2);

    tree_queries_stress_test(expected, actual);
}

TEST(HypersuccinctBinaryTreeTest, CartesianStressTest) {
    std::mt19937 engine(0);
    for (int n = 100; n <= 300; n += 100) {
        std::vector<int32_t> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), engine);

        CartesianTree ct(perm);
        BinaryTree expected = ct.binary_tree;

        HypersuccinctBinaryTree<> actual(cartesian_tree_bp(perm), 7);

        tree_queries_stress_test(expected, actual);
        is_ancestor_lca_stress_test(RMQBP<>(perm), actual);
    }
}

TEST(HypersuccinctBinaryTreeTest, IncreasingRuns) {
    const int n = 200;
    std::vector<int32_t> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    for (int s = 0; s < 5; s++) {
        random_roughly_fixed_incresing_runs(perm, 1 << s, 0);

        CartesianTree ct(perm);
        BinaryTree expected = ct.binary_tree;

        HypersuccinctBinaryTree<> actual(cartesian_tree_bp(perm), 7);

        tree_queries_stress_test(expected, actual);
        is_ancestor_lca_stress_test(RMQBP<>(perm), actual);
    }
}

}  // namespace
}  // namespace hyperrmq
