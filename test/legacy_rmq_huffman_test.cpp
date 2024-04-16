#include "hyperrmq/legacy_rmq_huffman.hpp"

#include <random>
#include <stack>

#include "gtest/gtest.h"
#include "hyperrmq/rmq_bp.hpp"
#include "hyperrmq/tree_bp.hpp"
#include "hyperrmq/tree_covering.hpp"

namespace hyperrmq {
namespace {

TEST(RMQHuffmanTest, GetMicrotreeInorderTest) {
    const int n = 100000, B = 7;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    LegacyRMQHuffman<128> rmq_huffman(perm, B);
    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree rmm_tree(tree);

    auto [upsilon, microtrees] = legacy_tree_covering(B, tree);
    for (int i = 0; i < microtrees.size(); i++) {
        auto expected = microtrees[i];
        auto actual = rmq_huffman.get_microtree_inorder(i);
        ASSERT_EQ(expected, actual);
    }
}

TEST(RMQHuffmanTest, RestoreBPFromChunk) {
    const int n = 1000, B = 7;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    LegacyRMQHuffman<128> rmq_huffman(perm, B);
    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree rmm_tree(tree);

    auto [upsilon, microtrees] = legacy_tree_covering(B, tree);
    BitArray chunk_concat(2 * n);
    int idx = 0;
    for (int i = 0; i < upsilon.bp.size(); i++) {
        auto [clen, cbp] = rmq_huffman.get_chunk(i);
        for (int i = 0; i < clen; i++) {
            chunk_concat.set(idx + i, (cbp >> i) & 1);
        }
        idx += clen;
    }
    ASSERT_EQ(idx, 2 * n);

    ASSERT_EQ(tree.bp, chunk_concat);
}

TEST(RMQHuffmanTest, Select1Rank1Test) {
    const int n = 1000, B = 7;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    LegacyRMQHuffman<128> rmq_huffman(perm, B);
    TreeBP tree = cartesian_tree_bp(perm);
    RMMTree rmm_tree(tree);

    auto [upsilon, microtrees] = legacy_tree_covering(B, tree);
    std::vector<std::pair<uint32_t, uint32_t>> chunks;
    chunks.reserve(upsilon.bp.size());
    for (int i = 0; i < upsilon.bp.size(); i++) {
        chunks.push_back(rmq_huffman.get_chunk(i));
    }
    int count = 0;
    for (int i = 0; i < chunks.size(); i++) {
        auto [clen, cbp] = chunks[i];
        for (int j = 0; j < clen; j++) {
            auto bit = (cbp >> j) & 1;
            if (bit) {
                std::pair<uint32_t, uint32_t> select_expected = {i, j};
                std::pair<uint32_t, uint32_t> select_actual =
                    rmq_huffman.select1(count);
                ASSERT_EQ(select_expected, select_actual);
            }
            auto rank_actual = rmq_huffman.rank1(i, j);
            auto rank_expected = count;
            ASSERT_EQ(rank_expected, rank_actual);

            count += bit;
        }
    }
}

TEST(RMQHuffmanTest, FewNodes) {
    std::mt19937 engine(0);
    for (int n = 1; n < 20; n++) {
        std::vector<int> perm(n);
        std::generate(perm.begin(), perm.end(), engine);

        LegacyRMQHuffman<128> rmq_huffman(perm);
        RMQBP<> rmq_bp(perm);

        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                ASSERT_EQ(rmq_huffman.query(i, j), rmq_bp.query(i, j));
            }
        }
    }
}

TEST(RMQHuffmanTest, QuerySeventyTest) {
    const int n = 70, B = 6;
    std::vector<int32_t> perm = {
        6,  7,  5,  9,  8,  4,  12, 11, 13, 14, 15, 10, 3,  18, 17, 16, 21, 22,
        20, 24, 25, 23, 19, 27, 26, 30, 29, 28, 33, 32, 31, 2,  35, 34, 36, 1,
        39, 42, 41, 40, 43, 38, 46, 45, 44, 49, 48, 51, 50, 52, 47, 55, 54, 56,
        53, 37, 58, 59, 57, 63, 62, 65, 64, 61, 67, 66, 60, 70, 69, 68};

    LegacyRMQHuffman<128> rmq_huffman(perm, B);
    RMQBP<> rmq_bp(perm);

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            ASSERT_EQ(rmq_huffman.query(i, j), rmq_bp.query(i, j));
        }
    }
}

TEST(RMQHuffmanTest, QueryStressTestSmall) {
    const int n = 200, B = 7;
    std::mt19937 engine(0);
    std::vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    LegacyRMQHuffman<128> rmq_huffman(perm, B);
    RMQBP<> rmq_bp(perm);

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            ASSERT_EQ(rmq_huffman.query(i, j), rmq_bp.query(i, j));
        }
    }
}

}  // namespace
}  // namespace hyperrmq
