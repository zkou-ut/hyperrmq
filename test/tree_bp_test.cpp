#include "hyperrmq/tree_bp.hpp"

#include <random>
#include <set>
#include <stack>

#include "gtest/gtest.h"
#include "hyperrmq/cartesian_tree.hpp"

namespace hyperrmq {
namespace {

TEST(TreeBPTest, TenNodes) {
    std::vector<int> perm = {1, 7, 3, 2, 4, 5, 0, 9, 8, 6};
    CartesianTree ct(perm);
    TreeBP ctbp = cartesian_tree_bp(perm);
    ASSERT_EQ(ct.binary_tree.to_BP_string(), ctbp.to_string());
}

TEST(TreeBPTest, EqualityAndInequality) {
    TreeBP ctbp1 = cartesian_tree_bp({1, 0, 2});
    TreeBP ctbp2 = cartesian_tree_bp({2, 0, 1});
    TreeBP ctbp3 = cartesian_tree_bp({0, 1, 2});

    ASSERT_TRUE(ctbp1 == ctbp2);
    ASSERT_FALSE(ctbp1 != ctbp2);

    ASSERT_FALSE(ctbp1 == ctbp3);
    ASSERT_TRUE(ctbp1 != ctbp3);

    std::set<TreeBP> s = {ctbp1, ctbp2, ctbp3};
    ASSERT_EQ(s.size(), 2);
}

TEST(TreeBPTest, SetFromPermutations) {
    const std::vector<int> catalan = {
        1,    1,    2,     5,     14,     42,     132,     429,
        1430, 4862, 16796, 58786, 208012, 742900, 2674440, 9694845};

    for (int n = 1; n <= 7; n++) {
        std::vector<int> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        std::set<TreeBP> trees;
        do {
            trees.insert(cartesian_tree_bp(perm));
        } while (std::next_permutation(perm.begin(), perm.end()));
        ASSERT_EQ(trees.size(), catalan[n]);
    }
}

TEST(TreeBPTest, BuildFromBP) {
    std::vector<int> perm = {1, 7, 3, 2, 4, 5, 0, 9, 8, 6};
    TreeBP ctbp = cartesian_tree_bp(perm);
    TreeBP another_ctbp(ctbp.to_string());
    ASSERT_EQ(ctbp, another_ctbp);
}

TEST(TreeBPTest, StressTestPermutation) {
    const int n = 10000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    for (int test_case = 0; test_case < 10; test_case++) {
        shuffle(perm.begin(), perm.end(), engine);
        CartesianTree ct(perm);
        TreeBP ctbp = cartesian_tree_bp(perm);
        ASSERT_EQ(ct.binary_tree.to_BP_string(), ctbp.to_string());
    }
}

TEST(TreeBPTest, StressTestSimilarValues) {
    const int n = 10000;
    std::mt19937 engine(0);
    for (int test_case = 0; test_case < 10; test_case++) {
        std::vector<int> values(n);
        for (int i = 0; i < n; i++) {
            values[i] = engine() % 100;
        }
        CartesianTree ct(values);
        TreeBP ctbp = cartesian_tree_bp(values);
        ASSERT_EQ(ct.binary_tree.to_BP_string(), ctbp.to_string());
    }
}

TEST(TreeBPTest, NaiveFwdSearch) {
    TreeBP tree("(()(()()))((()))");
    ASSERT_EQ(tree.naive_fwdsearch(0, -1), 17);
    ASSERT_EQ(tree.naive_fwdsearch(1, -1), 10);
    ASSERT_EQ(tree.naive_fwdsearch(2, -1), 3);
    ASSERT_EQ(tree.naive_fwdsearch(3, -1), 10);
    ASSERT_EQ(tree.naive_fwdsearch(4, -1), 9);
    ASSERT_EQ(tree.naive_fwdsearch(5, -1), 6);
    ASSERT_EQ(tree.naive_fwdsearch(6, -1), 9);
    ASSERT_EQ(tree.naive_fwdsearch(7, -1), 8);
    ASSERT_EQ(tree.naive_fwdsearch(8, -1), 9);
    ASSERT_EQ(tree.naive_fwdsearch(9, -1), 10);
    ASSERT_EQ(tree.naive_fwdsearch(10, -1), 17);
    ASSERT_EQ(tree.naive_fwdsearch(11, -1), 16);
    ASSERT_EQ(tree.naive_fwdsearch(12, -1), 15);
    ASSERT_EQ(tree.naive_fwdsearch(13, -1), 14);
    ASSERT_EQ(tree.naive_fwdsearch(14, -1), 15);
    ASSERT_EQ(tree.naive_fwdsearch(15, -1), 16);
    ASSERT_EQ(tree.naive_fwdsearch(2, -2), 10);
    ASSERT_EQ(tree.naive_fwdsearch(2, -3), 17);
    ASSERT_EQ(tree.naive_fwdsearch(5, -3), 10);
    ASSERT_EQ(tree.naive_fwdsearch(5, -4), 17);
}

TEST(TreeBPTest, NaiveBwdSearch) {
    TreeBP tree("(()(()()))((()))");

    ASSERT_EQ(tree.naive_bwdsearch(1, -1), 0);
    ASSERT_EQ(tree.naive_bwdsearch(2, -1), 1);
    ASSERT_EQ(tree.naive_bwdsearch(3, -1), 0);
    ASSERT_EQ(tree.naive_bwdsearch(4, -1), 3);
    ASSERT_EQ(tree.naive_bwdsearch(5, -1), 4);
    ASSERT_EQ(tree.naive_bwdsearch(6, -1), 3);
    ASSERT_EQ(tree.naive_bwdsearch(7, -1), 6);
    ASSERT_EQ(tree.naive_bwdsearch(8, -1), 3);
    ASSERT_EQ(tree.naive_bwdsearch(9, -1), 0);
    ASSERT_EQ(tree.naive_bwdsearch(10, -1), -1);
    ASSERT_EQ(tree.naive_bwdsearch(11, -1), 10);
    ASSERT_EQ(tree.naive_bwdsearch(12, -1), 11);
    ASSERT_EQ(tree.naive_bwdsearch(13, -1), 12);
    ASSERT_EQ(tree.naive_bwdsearch(14, -1), 11);
    ASSERT_EQ(tree.naive_bwdsearch(15, -1), 10);
    ASSERT_EQ(tree.naive_bwdsearch(16, -1), -1);

    ASSERT_EQ(tree.naive_bwdsearch(2, -2), 0);
    ASSERT_EQ(tree.naive_bwdsearch(2, -3), -1);

    ASSERT_EQ(tree.naive_bwdsearch(5, -3), 0);
    ASSERT_EQ(tree.naive_bwdsearch(5, -4), -1);

    ASSERT_EQ(tree.naive_bwdsearch(13, -3), 10);
    ASSERT_EQ(tree.naive_bwdsearch(13, -4), -1);
}

void exhaustive_lca_rmq_test(const std::vector<int> &values) {
    auto naive_rmq = [&](int first, int last) -> int {
        if (first > last) {
            std::swap(first, last);
        }
        auto min_val = values[first];
        auto min_idx = first;
        for (int i = first + 1; i < last + 1; i++) {
            if (min_val > values[i]) {
                min_val = values[i];
                min_idx = i;
            }
        }
        return min_idx;
    };

    auto tree = cartesian_tree_bp(values);

    for (int i = 0; i < values.size(); i++) {
        for (int j = 0; j < values.size(); j++) {
            ASSERT_EQ(tree.naive_lca(i, j), naive_rmq(i, j));
        }
    }
}

TEST(TreeBPTest, OpenCloseStressTest) {
    const int n = 100000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);

    std::stack<int> st;
    std::vector<std::pair<int, int>> paren_pairs;
    for (int i = 0; i < tree.bp.size(); i++) {
        if (tree.bp.get(i) == 0) {
            st.push(i);
        } else {
            paren_pairs.push_back({st.top(), i});
            st.pop();
        }
    }

    for (auto &&[l, r] : paren_pairs) {
        ASSERT_EQ(tree.naive_open(r), l);
        ASSERT_EQ(tree.naive_close(l), r);
    }
}

TEST(TreeCoveringTest, NaiveLCASmallTest) {
    std::vector<int> perm = {3, 1, 4, 1, 5};
    exhaustive_lca_rmq_test(perm);
}

TEST(TreeCoveringTest, NaiveLCAStressTest) {
    std::mt19937 engine(0);
    for (int test_case = 0; test_case < 10; test_case++) {
        int n = 1 << test_case;
        std::vector<int> values(n);
        for (int i = 0; i < n; i++) {
            values[i] = engine() % 10;
        }
        exhaustive_lca_rmq_test(values);
    }
}

}  // namespace
}  // namespace hyperrmq
