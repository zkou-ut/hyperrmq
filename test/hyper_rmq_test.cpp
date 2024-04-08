#include "hyper_rmq.hpp"

#include <random>
#include <stack>

#include "gtest/gtest.h"
#include "rmq_bp.hpp"
#include "runs.hpp"
#include "tree_bp.hpp"
#include "tree_covering.hpp"

namespace average_case_optimal_rmq {
namespace {

TEST(HyperRMQTest, RestoreBPFromChunk) {
    const int n = 1000, B = 7;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    HyperRMQ<128> hyper_rmq(perm, B);
    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree rmm_tree(tree);

    auto [upsilon, tree_and_split_rank] = tree_covering(B, tree);
    BitArray chunk_concat(2 * n);
    int idx = 0;
    for (int i = 0; i < upsilon.bp.size(); i++) {
        auto chunk = hyper_rmq.get_chunk(i);
        chunk_concat.write_interval(idx, chunk);
        idx += chunk.size();
    }
    ASSERT_EQ(idx, 2 * n);

    ASSERT_EQ(tree.bp, chunk_concat);
}

TEST(HyperRMQTest, InorderTauConversionTest) {
    const int n = 1000, B = 7;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    HyperRMQ<64> hyper_rmq(perm, B);
    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree rmm_tree(tree);

    auto [upsilon, microtrees] = legacy_tree_covering(B, tree);
    std::vector<BitArray> chunks;
    chunks.reserve(upsilon.bp.size());
    for (int i = 0; i < upsilon.bp.size(); i++) {
        chunks.push_back(hyper_rmq.get_chunk(i));
    }
    int count = 0;
    for (int chunk_idx = 0; chunk_idx < chunks.size(); chunk_idx++) {
        auto chunk = chunks[chunk_idx];
        int close = 0;
        for (int bit_idx = 0; bit_idx < chunk.size(); bit_idx++) {
            auto bit = chunk.get(bit_idx);
            if (bit) {
                std::pair<uint32_t, uint32_t> tau_expected = {chunk_idx, close};
                std::pair<uint32_t, uint32_t> tau_actual =
                    hyper_rmq.inorderselect(count);
                ASSERT_EQ(tau_expected, tau_actual);
            }
            auto inorder_actual = hyper_rmq.inorder({chunk_idx, close});
            auto inorder_expected = count;
            ASSERT_EQ(inorder_expected, inorder_actual);

            count += bit;
            close += bit;
        }
    }
}

TEST(HyperRMQTest, FewNodes) {
    std::mt19937 engine(0);
    const int B = 3;
    for (int n = 1; n < 20; n++) {
        std::vector<int> perm(n);
        std::generate(perm.begin(), perm.end(), engine);

        RMQBP<> rmq_bp(perm);
        HyperRMQ<32, CompressedMicrotreeSplitRankArrayHuffmanNaive<32>>
            rmq_huffman_naive(perm, B);
        HyperRMQ<32> rmq_huffman(perm, B);
        HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<true>>
            rmq_depth(perm, B);
        HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<false>>
            rmq_breadth(perm, B);
        HyperRMQ<32, CompressedMicrotreeSplitRankArrayAllArithmetic>
            rmq_all_arith(perm, B);

        for (int l = 0; l < n; l++) {
            for (int r = l; r < n; r++) {
                auto expected = rmq_bp.query(l, r);
                ASSERT_EQ(expected, rmq_huffman_naive.query(l, r));
                ASSERT_EQ(expected, rmq_huffman.query(l, r));
                ASSERT_EQ(expected, rmq_depth.query(l, r));
                ASSERT_EQ(expected, rmq_breadth.query(l, r));
                ASSERT_EQ(expected, rmq_all_arith.query(l, r));
            }
        }
    }
}

TEST(HyperRMQTest, QuerySeventyTest) {
    const int n = 70, B = 6;
    std::vector<int32_t> perm = {
        6,  7,  5,  9,  8,  4,  12, 11, 13, 14, 15, 10, 3,  18, 17, 16, 21, 22,
        20, 24, 25, 23, 19, 27, 26, 30, 29, 28, 33, 32, 31, 2,  35, 34, 36, 1,
        39, 42, 41, 40, 43, 38, 46, 45, 44, 49, 48, 51, 50, 52, 47, 55, 54, 56,
        53, 37, 58, 59, 57, 63, 62, 65, 64, 61, 67, 66, 60, 70, 69, 68};

    RMQBP<> rmq_bp(perm);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayHuffmanNaive<32>>
        rmq_huffman_naive(perm, B);
    HyperRMQ<32> rmq_huffman(perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<true>> rmq_depth(
        perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<false>>
        rmq_breadth(perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayAllArithmetic> rmq_all_arith(
        perm, B);

    for (int l = 0; l < n; l++) {
        for (int r = l; r < n; r++) {
            auto expected = rmq_bp.query(l, r);
            ASSERT_EQ(expected, rmq_huffman_naive.query(l, r));
            ASSERT_EQ(expected, rmq_huffman.query(l, r));
            ASSERT_EQ(expected, rmq_depth.query(l, r));
            ASSERT_EQ(expected, rmq_breadth.query(l, r));
            ASSERT_EQ(expected, rmq_all_arith.query(l, r));
        }
    }
}

TEST(HyperRMQTest, QueryStressTestSmall) {
    const int n = 100, B = 6;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    RMQBP<> rmq_bp(perm);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayHuffmanNaive<32>>
        rmq_huffman_naive(perm, B);
    HyperRMQ<32> rmq_huffman(perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<true>> rmq_depth(
        perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<false>>
        rmq_breadth(perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayAllArithmetic> rmq_all_arith(
        perm, B);

    for (int l = 0; l < n; l++) {
        for (int r = l; r < n; r++) {
            auto expected = rmq_bp.query(l, r);
            ASSERT_EQ(expected, rmq_huffman_naive.query(l, r));
            ASSERT_EQ(expected, rmq_huffman.query(l, r));
            ASSERT_EQ(expected, rmq_depth.query(l, r));
            ASSERT_EQ(expected, rmq_breadth.query(l, r));
            ASSERT_EQ(expected, rmq_all_arith.query(l, r));
        }
    }
}

TEST(HyperRMQTest, QueryStressTestMedium) {
    const int n = 10000, B = 20, query_count = 3000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    RMQBP<16, 1024> rmq_bp(perm);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayHuffmanNaive<32>>
        rmq_huffman_naive(perm, B);
    HyperRMQ<32> rmq_huffman(perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<true>> rmq_depth(
        perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<false>>
        rmq_breadth(perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayAllArithmetic> rmq_all_arith(
        perm, B);

    for (int q = 0; q < query_count; q++) {
        int l = engine() % n;
        int r = engine() % n;
        if (l > r) {
            std::swap(l, r);
        }

        auto expected = rmq_bp.query(l, r);
        ASSERT_EQ(expected, rmq_huffman_naive.query(l, r));
        ASSERT_EQ(expected, rmq_huffman.query(l, r));
        ASSERT_EQ(expected, rmq_depth.query(l, r));
        ASSERT_EQ(expected, rmq_breadth.query(l, r));
        ASSERT_EQ(expected, rmq_all_arith.query(l, r));
    }
}

TEST(HyperRMQTest, QueryStressTestVariousWidth) {
    const int n = 10000, B_huffman = 5, B_arith = 50, query_count = 300;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    RMQBP<16, 1024> rmq_bp(perm);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayHuffmanNaive<32>>
        rmq_huffman_naive(perm, B_huffman);
    HyperRMQ<32> rmq_huffman(perm, B_huffman);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<true>> rmq_depth(
        perm, B_arith);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<false>>
        rmq_breadth(perm, B_arith);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayAllArithmetic> rmq_all_arith(
        perm, B_arith);

    for (int w = 1; w <= n; w <<= 1) {
        for (int q = 0; q < query_count; q++) {
            int l = engine() % (n - w + 1);
            int r = l + w - 1;
            if (l > r) {
                std::swap(l, r);
            }

            auto expected = rmq_bp.query(l, r);
            ASSERT_EQ(expected, rmq_huffman_naive.query(l, r));
            ASSERT_EQ(expected, rmq_huffman.query(l, r));
            ASSERT_EQ(expected, rmq_depth.query(l, r));
            ASSERT_EQ(expected, rmq_breadth.query(l, r));
            ASSERT_EQ(expected, rmq_all_arith.query(l, r));
        }
    }
}

TEST(HyperRMQTest, QueryStressTestLargeB) {
    const int n = 1000000, B = 100000, query_count = 10;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    RMQBP<16, 1024> rmq_bp(perm);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayHuffmanNaive<32>>
        rmq_huffman_naive(perm, B);
    HyperRMQ<32> rmq_huffman(perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<true>> rmq_depth(
        perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<false>>
        rmq_breadth(perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayAllArithmetic> rmq_all_arith(
        perm, B);

    for (int q = 0; q < query_count; q++) {
        int l = engine() % n;
        int r = engine() % n;
        if (l > r) {
            std::swap(l, r);
        }

        auto expected = rmq_bp.query(l, r);
        ASSERT_EQ(expected, rmq_huffman_naive.query(l, r));
        ASSERT_EQ(expected, rmq_huffman.query(l, r));
        ASSERT_EQ(expected, rmq_depth.query(l, r));
        ASSERT_EQ(expected, rmq_breadth.query(l, r));
        ASSERT_EQ(expected, rmq_all_arith.query(l, r));
    }
}

TEST(HyperRMQTest, QueryStressTestMediumIncreasingRuns) {
    const int n = 10000, B = 20, query_count = 3000;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    random_roughly_fixed_incresing_runs(perm, perm.size() / 20);

    RMQBP<> rmq_bp(perm);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayHuffmanNaive<32>>
        rmq_huffman_naive(perm, B);
    HyperRMQ<32> rmq_huffman(perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<true>> rmq_depth(
        perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayArithmetic<false>>
        rmq_breadth(perm, B);
    HyperRMQ<32, CompressedMicrotreeSplitRankArrayAllArithmetic> rmq_all_arith(
        perm, B);

    for (int q = 0; q < query_count; q++) {
        int l = engine() % n;
        int r = engine() % n;
        if (l > r) {
            std::swap(l, r);
        }

        auto expected = rmq_bp.query(l, r);
        ASSERT_EQ(expected, rmq_huffman_naive.query(l, r));
        ASSERT_EQ(expected, rmq_huffman.query(l, r));
        ASSERT_EQ(expected, rmq_depth.query(l, r));
        ASSERT_EQ(expected, rmq_breadth.query(l, r));
        ASSERT_EQ(expected, rmq_all_arith.query(l, r));
    }
}

}  // namespace
}  // namespace average_case_optimal_rmq
