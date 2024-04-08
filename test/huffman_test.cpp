#include "hyperrmq/huffman.hpp"

#include <queue>
#include <random>

#include "gtest/gtest.h"
#include "hyperrmq/microtree_array.hpp"
#include "hyperrmq/minimal_cell_array.hpp"
#include "hyperrmq/tree_bp.hpp"

namespace average_case_optimal_rmq {
namespace {

TEST(HuffmanTreeTest, ThreeLeaves) {
    HuffmanTree ht({
        {1, 1},
        {2, 2},
        {3, 4},
    });
    ASSERT_EQ(ht.root->freq, 7);
    ASSERT_EQ(ht.root->left->freq, 3);
    ASSERT_EQ(ht.root->right->freq, 4);
    ASSERT_EQ(ht.root->left->left->freq, 1);
    ASSERT_EQ(ht.root->left->right->freq, 2);

    ASSERT_EQ(ht.root->left->left->key, 1);
    ASSERT_EQ(ht.root->left->right->key, 2);
    ASSERT_EQ(ht.root->right->key, 3);
}

TEST(HuffmanTreeTest, CodeLengthSumStressTest) {
    const int n = 100;
    std::mt19937 engine(0);
    for (int test_case = 0; test_case < 10; test_case++) {
        std::map<uint64_t, uint64_t> value_counts;
        for (int i = 0; i < n; i++) {
            value_counts[i] = engine() % 100 + 1;
        }

        HuffmanTree<> ht(value_counts);
        auto lengths_actual = ht.compute_length();
        uint64_t length_sum_actual = 0;
        for (auto &&[val, length] : lengths_actual) {
            length_sum_actual += length * value_counts[val];
        }

        std::priority_queue<uint64_t, std::vector<uint64_t>,
                            std::greater<uint64_t>>
            pq;
        for (auto &&[val, cnt] : value_counts) {
            pq.push(cnt);
        }
        uint64_t length_sum_expected = 0;
        while (pq.size() > 1) {
            auto c1 = pq.top();
            pq.pop();
            auto c2 = pq.top();
            pq.pop();
            auto c = c1 + c2;
            pq.push(c);
            length_sum_expected += c;
        }

        ASSERT_EQ(length_sum_actual, length_sum_expected);
    }
}

TEST(CanonicalHuffmanCodeTest, OneLeave) {
    std::map<uint64_t, uint64_t> count = {
        {3, 4},
    };
    CanonicalHuffmanCode<> chc(count);
    std::map<uint64_t, std::pair<uint64_t, uint64_t>> expected = {
        {3, {0, 0b0}}};

    ASSERT_EQ(expected, chc.enumerate_alphabet_code_pair());

    BitArray code_seq(0);
    ASSERT_EQ(chc.get_next_length(code_seq, 0), 0);
    ASSERT_EQ(chc.decode(0, 0b0), 3);
}

TEST(CanonicalHuffmanCodeTest, ThreeLeaves) {
    CanonicalHuffmanCode<> chc({
        {1, 1},
        {2, 2},
        {3, 4},
    });
    std::map<uint64_t, std::pair<uint64_t, uint64_t>> expected = {
        {3, {1, 0b0}}, {1, {2, 0b10}}, {2, {2, 0b11}}};

    ASSERT_EQ(expected, chc.enumerate_alphabet_code_pair());

    BitArray code_seq(7);
    code_seq.write_bits(0, 6, 0b010110);
    ASSERT_EQ(chc.get_next_length(code_seq, 0), 1);
    ASSERT_EQ(chc.get_next_length(code_seq, 1), 2);
    ASSERT_EQ(chc.get_next_length(code_seq, 3), 2);
    ASSERT_EQ(chc.get_next_length(code_seq, 5), 1);
    ASSERT_EQ(chc.decode(1, 0b0), 3);
    ASSERT_EQ(chc.decode(2, 0b10), 1);
    ASSERT_EQ(chc.decode(2, 0b11), 2);
}

TEST(CanonicalHuffmanCodeTest, Trees) {
    CanonicalHuffmanCode<TreeBP, EditableMicrotreeArray> chc({
        {TreeBP("()"), 4},
        {TreeBP("()()"), 2},
        {TreeBP("(())"), 1},
    });

    ASSERT_EQ(chc.decode(1, 0b0), TreeBP("()"));
    ASSERT_EQ(chc.decode(2, 0b10), TreeBP("(())"));
    ASSERT_EQ(chc.decode(2, 0b11), TreeBP("()()"));
}

TEST(CanonicalHuffmanCodeTest, ThreeLeavesMinimalCellArray) {
    CanonicalHuffmanCode<uint64_t, MinimalCellArray<uint64_t>> chc({
        {1, 1},
        {2, 2},
        {3, 4},
    });
    std::map<uint64_t, std::pair<uint64_t, uint64_t>> expected = {
        {3, {1, 0b0}}, {1, {2, 0b10}}, {2, {2, 0b11}}};

    ASSERT_EQ(expected, chc.enumerate_alphabet_code_pair());

    BitArray code_seq(7);
    code_seq.write_bits(0, 6, 0b010110);
    ASSERT_EQ(chc.get_next_length(code_seq, 0), 1);
    ASSERT_EQ(chc.get_next_length(code_seq, 1), 2);
    ASSERT_EQ(chc.get_next_length(code_seq, 3), 2);
    ASSERT_EQ(chc.get_next_length(code_seq, 5), 1);
    ASSERT_EQ(chc.decode(1, 0b0), 3);
    ASSERT_EQ(chc.decode(2, 0b10), 1);
    ASSERT_EQ(chc.decode(2, 0b11), 2);
}

template <typename T>
void length_stress_test_helper(std::map<T, uint64_t> value_counts) {
    CanonicalHuffmanCode<T> chc(value_counts);
    auto actual = chc.enumerate_alphabet_code_pair();
    auto expected_pair = HuffmanTree<T>(value_counts).compute_length();
    for (auto &&[k, v] : expected_pair) {
        ASSERT_EQ(actual[k].first, v);
        ASSERT_LT(actual[k].second, 1ull << actual[k].first);
    }
}

TEST(CanonicalHuffmanCodeTest, LengthStressTest) {
    const int n = 100;
    std::mt19937 engine(0);
    for (int test_case = 0; test_case < 10; test_case++) {
        std::map<uint64_t, uint64_t> value_counts;
        for (int i = 0; i < n; i++) {
            value_counts[i] = engine() % 100 + 1;
        }
        length_stress_test_helper(value_counts);
    }
}

TEST(CanonicalHuffmanCodeTest, BiasedLengthStressTest) {
    const int n = 100;
    std::mt19937 engine(0);
    for (int test_case = 0; test_case < 10; test_case++) {
        std::map<uint64_t, uint64_t> value_counts;
        for (int i = 0; i < n; i++) {
            if (engine() % 2) {
                value_counts[i] = engine() % 10000000 + 1;
            } else {
                value_counts[i] = engine() % 10 + 1;
            }
        }
        length_stress_test_helper(value_counts);
    }
}

}  // namespace
}  // namespace average_case_optimal_rmq
