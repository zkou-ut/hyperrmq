#include "rmm_tree.hpp"

#include <random>
#include <stack>

#include "gtest/gtest.h"
#include "tree_bp.hpp"

namespace average_case_optimal_rmq {
namespace {

TEST(RmmTreeTest, IdxLeafConversion) {
    RMMTree<8, 16> rmm(cartesian_tree_bp(std::vector<int>(80)));
    ASSERT_EQ(rmm.num_of_blocks, 10);

    ASSERT_EQ(rmm.leaf_to_idx(16), 0);
    ASSERT_EQ(rmm.leaf_to_idx(17), 1);
    ASSERT_EQ(rmm.leaf_to_idx(18), 2);
    ASSERT_EQ(rmm.leaf_to_idx(19), 3);
    ASSERT_EQ(rmm.leaf_to_idx(10), 4);
    ASSERT_EQ(rmm.leaf_to_idx(11), 5);
    ASSERT_EQ(rmm.leaf_to_idx(12), 6);
    ASSERT_EQ(rmm.leaf_to_idx(13), 7);
    ASSERT_EQ(rmm.leaf_to_idx(14), 8);
    ASSERT_EQ(rmm.leaf_to_idx(15), 9);

    ASSERT_EQ(rmm.idx_to_leaf(0), 16);
    ASSERT_EQ(rmm.idx_to_leaf(1), 17);
    ASSERT_EQ(rmm.idx_to_leaf(2), 18);
    ASSERT_EQ(rmm.idx_to_leaf(3), 19);
    ASSERT_EQ(rmm.idx_to_leaf(4), 10);
    ASSERT_EQ(rmm.idx_to_leaf(5), 11);
    ASSERT_EQ(rmm.idx_to_leaf(6), 12);
    ASSERT_EQ(rmm.idx_to_leaf(7), 13);
    ASSERT_EQ(rmm.idx_to_leaf(8), 14);
    ASSERT_EQ(rmm.idx_to_leaf(9), 15);
}

TEST(RmmTreeTest, IdxLeafIdentity) {
    for (int node = 100; node <= 1000; node += 100) {
        RMMTree<8, 16> rmm(cartesian_tree_bp(std::vector<int>(node)));
        for (int leaf = 0; leaf < rmm.num_of_blocks; leaf++) {
            int num = rmm.idx_to_leaf(leaf);
            int leaf2 = rmm.leaf_to_idx(num);
            ASSERT_EQ(leaf, leaf2);
        }
    }
}

TEST(RmmTreeTest, IsLeftmost) {
    RMMTree<8, 16> rmm(cartesian_tree_bp(std::vector<int>(80)));
    ASSERT_EQ(rmm.num_of_blocks, 10);

    std::set<int> leftmost = {1, 2, 4, 8, 16};
    for (int v = 1; v < rmm.num_of_blocks * 2; v++) {
        ASSERT_EQ(rmm.is_leftmost(v), bool(leftmost.count(v)));
    }
}

TEST(RmmTreeTest, IsRightmost) {
    RMMTree<8, 16> rmm(cartesian_tree_bp(std::vector<int>(80)));
    ASSERT_EQ(rmm.num_of_blocks, 10);

    std::set<int> rightmost = {1, 3, 7, 15};
    for (int v = 1; v < rmm.num_of_blocks * 2; v++) {
        ASSERT_EQ(rmm.is_rightmost(v), bool(rightmost.count(v)));
    }
}

TEST(RmmTreeTest, VertexToInterval) {
    RMMTree<8, 16> rmm(cartesian_tree_bp(std::vector<int>(80)));
    ASSERT_EQ(rmm.num_of_blocks, 10);

    std::vector<std::pair<uint32_t, uint32_t>> results_from_one = {
        {0, 10}, {0, 6},  {6, 10}, {0, 4}, {4, 6},
        {6, 8},  {8, 10}, {0, 2},  {2, 4}};

    for (int v = 0; v < 9; v++) {
        auto expected = results_from_one[v];
        ASSERT_EQ(rmm.vertex_to_interval(v + 1), expected);
    }

    for (uint32_t i = 0; i < 10; i++) {
        ASSERT_EQ(rmm.vertex_to_interval(rmm.idx_to_leaf(i)),
                  std::make_pair(i, i + 1));
    }
}

TEST(RmmTreeTest, NodeExcess) {
    const int n = 10000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 16> rmm(tree);

    std::vector<int32_t> actual(rmm.num_of_blocks * 2);
    for (int i = 0; i < rmm.num_of_blocks; i++) {
        actual[rmm.idx_to_leaf(i)] =
            rmm.excess_sample[i + 1] - rmm.excess_sample[i];
    }
    for (int i = rmm.num_of_blocks - 1; i >= 1; i--) {
        actual[i] = actual[2 * i] + actual[2 * i + 1];
    }

    for (int v = 1; v < 2 * rmm.num_of_blocks; v++) {
        ASSERT_EQ(rmm.node_excess(v), actual[v]);
    }
}

TEST(RmmTreeTest, ExcessPrefixSum) {
    const int n = 10000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 16> rmm(tree);
    ASSERT_EQ(rmm.excess_prefix_sum(0), 0);
    uint32_t expected = 0;
    for (int i = 0; i < tree.bp.size(); i++) {
        expected += (tree.bp.get(i) == 0) ? 1 : -1;
        ASSERT_EQ(rmm.excess_prefix_sum(i + 1), expected);
    }
}

TEST(RmmTreeTest, Rank0) {
    const int n = 10000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 16> rmm(tree);
    ASSERT_EQ(rmm.rank0(0), 0);
    uint32_t expected = 0;
    for (int i = 0; i < tree.bp.size(); i++) {
        expected += (tree.bp.get(i) == 0);
        ASSERT_EQ(rmm.rank0(i + 1), expected);
    }
}

TEST(RmmTreeTest, Rank1) {
    const int n = 10000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 16> rmm(tree);
    ASSERT_EQ(rmm.rank1(0), 0);
    uint32_t expected = 0;
    for (int i = 0; i < tree.bp.size(); i++) {
        expected += (tree.bp.get(i) == 1);
        ASSERT_EQ(rmm.rank1(i + 1), expected);
    }
}

TEST(RmmTreeTest, SelectRankIdentity) {
    const int n = 10000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 16> rmm(tree);
    for (int i = 0; i < tree.bp.size(); i++) {
        if (tree.bp.get(i) == 0) {
            ASSERT_EQ(i, rmm.select0(rmm.rank0(i)));
        } else {
            ASSERT_EQ(i, rmm.select1(rmm.rank1(i)));
        }
    }
}

TEST(RmmTreeTest, MinExcessInequality) {
    const int n = 10000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 16> rmm(tree);
    for (int i = 1; i < rmm.num_of_blocks * 2; i++) {
        ASSERT_GE(rmm.node_excess(i), rmm.min_tree[i]);
    }
}

TEST(RmmTreeTest, FwdBlock) {
    TreeBP tree("(()(()()))((()))");
    RMMTree<8, 16> rmm(tree);
    ASSERT_EQ(rmm.fwdblock(0, -1), std::make_pair(0, 16u));
    ASSERT_EQ(rmm.fwdblock(1, -1), std::make_pair(-1, 10u));
    ASSERT_EQ(rmm.fwdblock(2, -1), std::make_pair(-1, 3u));
    ASSERT_EQ(rmm.fwdblock(3, -1), std::make_pair(-1, 10u));
    ASSERT_EQ(rmm.fwdblock(4, -1), std::make_pair(-1, 9u));
    ASSERT_EQ(rmm.fwdblock(5, -1), std::make_pair(-1, 6u));
    ASSERT_EQ(rmm.fwdblock(6, -1), std::make_pair(-1, 9u));
    ASSERT_EQ(rmm.fwdblock(7, -1), std::make_pair(-1, 8u));
    ASSERT_EQ(rmm.fwdblock(8, -1), std::make_pair(-1, 9u));
    ASSERT_EQ(rmm.fwdblock(9, -1), std::make_pair(-1, 10u));
    ASSERT_EQ(rmm.fwdblock(10, -1), std::make_pair(0, 16u));
    ASSERT_EQ(rmm.fwdblock(11, -1), std::make_pair(-1, 16u));
    ASSERT_EQ(rmm.fwdblock(12, -1), std::make_pair(-1, 15u));
    ASSERT_EQ(rmm.fwdblock(13, -1), std::make_pair(-1, 14u));
    ASSERT_EQ(rmm.fwdblock(14, -1), std::make_pair(-1, 15u));
    ASSERT_EQ(rmm.fwdblock(15, -1), std::make_pair(-1, 16u));

    ASSERT_EQ(rmm.fwdblock(2, -2), std::make_pair(-2, 10u));
    ASSERT_EQ(rmm.fwdblock(2, -3), std::make_pair(-2, 16u));

    ASSERT_EQ(rmm.fwdblock(5, -3), std::make_pair(-3, 10u));
    ASSERT_EQ(rmm.fwdblock(5, -4), std::make_pair(-3, 16u));
}

TEST(RmmTreeTest, FwdSearchSmall) {
    TreeBP tree("(()(()()))((()))");
    RMMTree<8, 8> rmm(tree);
    ASSERT_EQ(rmm.fwdsearch(0, -1), 17);
    ASSERT_EQ(rmm.fwdsearch(1, -1), 10);
    ASSERT_EQ(rmm.fwdsearch(2, -1), 3);
    ASSERT_EQ(rmm.fwdsearch(3, -1), 10);
    ASSERT_EQ(rmm.fwdsearch(4, -1), 9);
    ASSERT_EQ(rmm.fwdsearch(5, -1), 6);
    ASSERT_EQ(rmm.fwdsearch(6, -1), 9);
    ASSERT_EQ(rmm.fwdsearch(7, -1), 8);
    ASSERT_EQ(rmm.fwdsearch(8, -1), 9);
    ASSERT_EQ(rmm.fwdsearch(9, -1), 10);
    ASSERT_EQ(rmm.fwdsearch(10, -1), 17);
    ASSERT_EQ(rmm.fwdsearch(11, -1), 16);
    ASSERT_EQ(rmm.fwdsearch(12, -1), 15);
    ASSERT_EQ(rmm.fwdsearch(13, -1), 14);
    ASSERT_EQ(rmm.fwdsearch(14, -1), 15);
    ASSERT_EQ(rmm.fwdsearch(15, -1), 16);

    ASSERT_EQ(rmm.fwdsearch(2, -2), 10);
    ASSERT_EQ(rmm.fwdsearch(2, -3), 17);

    ASSERT_EQ(rmm.fwdsearch(5, -3), 10);
    ASSERT_EQ(rmm.fwdsearch(5, -4), 17);
}

TEST(RmmTreeTest, FwdSearchVerySmall) {
    TreeBP tree("()");
    RMMTree<8, 16> rmm(tree);
    ASSERT_EQ(rmm.fwdsearch(0, -1), 3);
    ASSERT_EQ(rmm.fwdsearch(1, -1), 2);
    ASSERT_EQ(rmm.fwdsearch(1, -2), 3);
}

TEST(RmmTreeTest, FwdSearchStressTest) {
    const int n = 1000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 16> rmm(tree);

    auto naive_fwdsearch = [&](int i, int d) -> int {
        int j = i, d_prime = 0;
        while (j < rmm.num_of_bits) {
            d_prime += 1 - 2 * tree.bp.get(j);
            j++;
            if (d_prime == d) {
                return j;
            }
        }
        return rmm.num_of_bits + 1;
    };

    for (int i = 0; i < rmm.num_of_bits; i++) {
        int d = -1;
        while (true) {
            auto expected = naive_fwdsearch(i, d);
            auto actual = rmm.fwdsearch(i, d);
            ASSERT_EQ(expected, actual);
            if (expected == rmm.num_of_bits + 1) {
                break;
            }
            d--;
        }
    }
}

TEST(RmmTreeTest, CloseStressTest) {
    const int n = 100000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 256> rmm(tree);

    std::stack<int> st;
    std::map<int, int> close_result;
    for (int i = 0; i < tree.bp.size(); i++) {
        if (tree.bp.get(i) == 0) {
            st.push(i);
        } else {
            close_result[st.top()] = i;
            st.pop();
        }
    }

    for (auto &&[l, r] : close_result) {
        ASSERT_EQ(rmm.close(l), r);
    }
}

TEST(RmmTreeTest, BwdBlock) {
    TreeBP tree("(()(()()))((()))");
    RMMTree<8, 16> rmm(tree);

    ASSERT_EQ(rmm.bwdblock(1, -1), std::make_pair(-1, 0u));
    ASSERT_EQ(rmm.bwdblock(2, -1), std::make_pair(-1, 1u));
    ASSERT_EQ(rmm.bwdblock(3, -1), std::make_pair(-1, 0u));
    ASSERT_EQ(rmm.bwdblock(4, -1), std::make_pair(-1, 3u));
    ASSERT_EQ(rmm.bwdblock(5, -1), std::make_pair(-1, 4u));
    ASSERT_EQ(rmm.bwdblock(6, -1), std::make_pair(-1, 3u));
    ASSERT_EQ(rmm.bwdblock(7, -1), std::make_pair(-1, 6u));
    ASSERT_EQ(rmm.bwdblock(8, -1), std::make_pair(-1, 3u));
    ASSERT_EQ(rmm.bwdblock(9, -1), std::make_pair(-1, 0u));
    ASSERT_EQ(rmm.bwdblock(10, -1), std::make_pair(0, 0u));
    ASSERT_EQ(rmm.bwdblock(11, -1), std::make_pair(-1, 10u));
    ASSERT_EQ(rmm.bwdblock(12, -1), std::make_pair(-1, 11u));
    ASSERT_EQ(rmm.bwdblock(13, -1), std::make_pair(-1, 12u));
    ASSERT_EQ(rmm.bwdblock(14, -1), std::make_pair(-1, 11u));
    ASSERT_EQ(rmm.bwdblock(15, -1), std::make_pair(-1, 10u));
    ASSERT_EQ(rmm.bwdblock(16, -1), std::make_pair(0, 0u));

    ASSERT_EQ(rmm.bwdblock(2, -2), std::make_pair(-2, 0u));
    ASSERT_EQ(rmm.bwdblock(2, -3), std::make_pair(-2, 0u));

    ASSERT_EQ(rmm.bwdblock(5, -3), std::make_pair(-3, 0u));
    ASSERT_EQ(rmm.bwdblock(5, -4), std::make_pair(-3, 0u));

    ASSERT_EQ(rmm.bwdblock(13, -3), std::make_pair(-3, 10u));
    ASSERT_EQ(rmm.bwdblock(13, -4), std::make_pair(-3, 0u));
}

TEST(RmmTreeTest, BwdSearch) {
    TreeBP tree("(()(()()))((()))");
    RMMTree<8, 8> rmm(tree);

    ASSERT_EQ(rmm.bwdsearch(1, -1), 0);
    ASSERT_EQ(rmm.bwdsearch(2, -1), 1);
    ASSERT_EQ(rmm.bwdsearch(3, -1), 0);
    ASSERT_EQ(rmm.bwdsearch(4, -1), 3);
    ASSERT_EQ(rmm.bwdsearch(5, -1), 4);
    ASSERT_EQ(rmm.bwdsearch(6, -1), 3);
    ASSERT_EQ(rmm.bwdsearch(7, -1), 6);
    ASSERT_EQ(rmm.bwdsearch(8, -1), 3);
    ASSERT_EQ(rmm.bwdsearch(9, -1), 0);
    ASSERT_EQ(rmm.bwdsearch(10, -1), -1);
    ASSERT_EQ(rmm.bwdsearch(11, -1), 10);
    ASSERT_EQ(rmm.bwdsearch(12, -1), 11);
    ASSERT_EQ(rmm.bwdsearch(13, -1), 12);
    ASSERT_EQ(rmm.bwdsearch(14, -1), 11);
    ASSERT_EQ(rmm.bwdsearch(15, -1), 10);
    ASSERT_EQ(rmm.bwdsearch(16, -1), -1);

    ASSERT_EQ(rmm.bwdsearch(2, -2), 0);
    ASSERT_EQ(rmm.bwdsearch(2, -3), -1);

    ASSERT_EQ(rmm.bwdsearch(5, -3), 0);
    ASSERT_EQ(rmm.bwdsearch(5, -4), -1);

    ASSERT_EQ(rmm.bwdsearch(13, -3), 10);
    ASSERT_EQ(rmm.bwdsearch(13, -4), -1);
}

TEST(RmmTreeTest, BwdSearchStressTest) {
    const int n = 1000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 16> rmm(tree);

    auto naive_bwdsearch = [&](int i, int d) -> int {
        int j = i - 1, d_prime = 0;
        while (j >= 0) {
            d_prime -= 1 - 2 * tree.bp.get(j);
            if (d_prime == d) {
                return j;
            }
            j--;
        }
        return -1;
    };

    for (int i = 1; i <= rmm.num_of_bits; i++) {
        int d = -1;
        while (true) {
            auto expected = naive_bwdsearch(i, d);
            auto actual = rmm.bwdsearch(i, d);
            ASSERT_EQ(expected, actual);
            if (expected == -1) {
                break;
            }
            d--;
        }
    }
}

TEST(RmmTreeTest, OpenStressTest) {
    const int n = 100000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 256> rmm(tree);

    std::stack<int> st;
    std::map<int, int> open_result;
    for (int i = 0; i < tree.bp.size(); i++) {
        if (tree.bp.get(i) == 0) {
            st.push(i);
        } else {
            open_result[i] = st.top();
            st.pop();
        }
    }

    for (auto &&[r, l] : open_result) {
        ASSERT_EQ(rmm.open(r), l);
    }
}

TEST(RmmTreeTest, MinBlockSmall) {
    TreeBP tree("(()(()()))((()))");
    RMMTree<8, 16> rmm(tree);
    auto naive_minblock = [&](int i, int j) -> std::pair<int, int> {
        int d = 0, m = 0;
        for (int k = i; k < j; k++) {
            d += 1 - 2 * tree.bp.get(k);
            m = std::min(m, d);
        }
        return {m, d};
    };

    for (int i = 0; i < 17; i++) {
        for (int j = i; j < 17; j++) {
            ASSERT_EQ(rmm.minblock(i, j), naive_minblock(i, j));
        }
    }
}

TEST(RmmTreeTest, MinBlockStressTest) {
    const int n = 1000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 1024> rmm(tree);

    for (int i = 0; i < n; i++) {
        int d = 0, m_expected = 0;
        ASSERT_EQ(rmm.minblock(i, i), std::make_pair(m_expected, d));
        for (int j = i; j < n; j++) {
            d += 1 - 2 * tree.bp.get(j);
            m_expected = std::min(m_expected, d);
            ASSERT_EQ(rmm.minblock(i, j + 1), std::make_pair(m_expected, d));
        }
    }
}

TEST(RmmTreeTest, MinExcessStressTest) {
    const int n = 1000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 16> rmm(tree);

    for (int i = 0; i < n; i++) {
        int d = 0, m_expected = 0;
        ASSERT_EQ(rmm.minexcess(i, i), m_expected);
        for (int j = i; j < n; j++) {
            d += 1 - 2 * tree.bp.get(j);
            m_expected = std::min(m_expected, d);
            ASSERT_EQ(rmm.minexcess(i, j + 1), m_expected);
        }
    }
}

TEST(RmmTreeTest, RmqStressTest) {
    const int n = 1000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree<8, 16> rmm(tree);

    for (int i = 0; i < n; i++) {
        int d = 0, m = 0;
        uint32_t expected = i;
        for (int j = i; j < n; j++) {
            d += 1 - 2 * tree.bp.get(j);
            if (m > d) {
                m = d;
                expected = j;
            }
            ASSERT_EQ(rmm.rmq(i, j + 1), expected);
        }
    }
}

}  // namespace
}  // namespace average_case_optimal_rmq
