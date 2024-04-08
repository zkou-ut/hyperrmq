#include "microtree_array.hpp"

#include <random>

#include "gtest/gtest.h"

namespace average_case_optimal_rmq {
namespace {

TEST(MicrotreeArrayTest, EditableSmall) {
    EditableMicrotreeArray mts(3, 4);
    std::vector<std::string> bps = {"()", "()()", "(())", "(())(())"};

    for (int i = 0; i < bps.size(); i++) {
        for (int j = 0; j < bps[i].size(); j++) {
            mts.set(i, j, bps[i][j] == ')');
        }
    }

    for (int i = 0; i < bps.size(); i++) {
        for (int j = 0; j < bps[i].size(); j++) {
            ASSERT_EQ(mts.get(i, j), bps[i][j] == ')');
        }
    }

    for (int i = 0; i < bps.size(); i++) {
        ASSERT_EQ(mts[i], TreeBP(bps[i]));
    }

    std::vector<TreeBP> trees;
    for (auto &&bp : bps) {
        trees.push_back(TreeBP(bp));
    }
    EditableMicrotreeArray mts2 = trees;
    EditableMicrotreeArray mts3(4, trees.size());
    for (int i = 0; i < trees.size(); i++) {
        mts3.set(i, trees[i]);
    }

    ASSERT_EQ(mts, mts2);
    ASSERT_TRUE(mts == mts2);
    ASSERT_FALSE(mts != mts2);

    ASSERT_EQ(mts, mts3);
    ASSERT_TRUE(mts == mts3);
    ASSERT_FALSE(mts != mts3);

    mts.set(0, mts2[1]);
    ASSERT_NE(mts, mts2);
    ASSERT_FALSE(mts == mts2);
    ASSERT_TRUE(mts != mts2);

    mts.set(0, mts2[0]);
    ASSERT_EQ(mts, mts2);
    ASSERT_TRUE(mts == mts2);
    ASSERT_FALSE(mts != mts2);
}

std::string random_tree_bp(int node_count, std::mt19937 engine) {
    std::string bp;
    int open = 0, close = 0;
    while (open + close < 2 * node_count) {
        bool use_open;
        if (open == close) {
            use_open = true;
        } else if (open == node_count) {
            use_open = false;
        } else {
            use_open = (engine() % 2);
        }
        open += use_open;
        close += !use_open;
        bp.push_back(use_open ? '(' : ')');
    }
    return bp;
};

TEST(MicrotreeArrayTest, EditableStressTestBitWise) {
    int tree_count = 100, block_size = 100;
    std::vector<std::string> bps(tree_count);

    std::mt19937 engine(0);

    for (int i = 0; i < tree_count; i++) {
        bps[i] = random_tree_bp(2 * block_size - 1, engine);
    }

    EditableMicrotreeArray mts(block_size, tree_count);

    for (int i = 0; i < bps.size(); i++) {
        for (int j = 0; j < bps[i].size(); j++) {
            mts.set(i, j, bps[i][j] == ')');
        }
    }

    for (int i = 0; i < bps.size(); i++) {
        for (int j = 0; j < bps[i].size(); j++) {
            ASSERT_EQ(mts.get(i, j), bps[i][j] == ')');
        }
    }

    for (int i = 0; i < bps.size(); i++) {
        ASSERT_EQ(mts[i], TreeBP(bps[i]));
    }
}

TEST(MicrotreeArrayTest, EditableStressTest) {
    int tree_count = 10, block_size = 10;
    int query_count = 100;

    std::mt19937 engine(0);

    std::vector<TreeBP> expected(tree_count);
    EditableMicrotreeArray actual(block_size, tree_count);
    ASSERT_EQ(EditableMicrotreeArray(expected), actual);
    for (int q = 0; q < query_count; q++) {
        TreeBP tree(random_tree_bp(engine() % (2 * block_size), engine));
        int i = engine() % tree_count;
        expected[i] = tree;
        actual.set(i, tree);
        ASSERT_EQ(EditableMicrotreeArray(expected), actual);
    }
}

template <typename CompressedMicrotreeSplitRankArray>
void test_microtree_split_rank_array_helper() {
    static_assert(std::is_base_of_v<CompressedMicrotreeSplitRankArrayInterface,
                                    CompressedMicrotreeSplitRankArray>);
    std::mt19937 engine(0);

    for (int t = 0; t < 100; t++) {
        int n = engine() % 100 + 1;
        int B = engine() % 10 + 1;
        std::vector<std::pair<TreeBP, uint32_t>> vec(n);
        for (int i = 0; i < n; i++) {
            vec[i].first =
                TreeBP(random_tree_bp(engine() % (2 * B - 1) + 1, engine));
            vec[i].second = engine() % (vec[i].first.n + 1);
        }
        MicrotreeSplitRankArray array(vec);
        CompressedMicrotreeSplitRankArray compressed(array);
        ASSERT_EQ(vec.size(), array.size());
        ASSERT_EQ(array.size(), compressed.size());
        for (int i = 0; i < n; i++) {
            ASSERT_EQ(vec[i], array[i]);
            ASSERT_EQ(array[i], compressed[i]);
        }
    }
};

TEST(MicrotreeArrayTest, MicrotreeSplitRankArrayHuffmanNaiveStressTest) {
    test_microtree_split_rank_array_helper<
        CompressedMicrotreeSplitRankArrayHuffmanNaive<>>();
}

TEST(MicrotreeArrayTest, MicrotreeSplitRankArrayHuffmanStressTest) {
    test_microtree_split_rank_array_helper<
        CompressedMicrotreeSplitRankArrayHuffman<>>();
}

TEST(MicrotreeArrayTest, MicrotreeSplitRankArrayArithmeticDepthStressTest) {
    test_microtree_split_rank_array_helper<
        CompressedMicrotreeSplitRankArrayArithmetic<true>>();
}

TEST(MicrotreeArrayTest, MicrotreeSplitRankArrayArithmeticBreadthStressTest) {
    test_microtree_split_rank_array_helper<
        CompressedMicrotreeSplitRankArrayArithmetic<false>>();
}

TEST(MicrotreeArrayTest, MicrotreeSplitRankArrayAllArithmeticStressTest) {
    test_microtree_split_rank_array_helper<
        CompressedMicrotreeSplitRankArrayAllArithmetic>();
}

}  // namespace
}  // namespace average_case_optimal_rmq
