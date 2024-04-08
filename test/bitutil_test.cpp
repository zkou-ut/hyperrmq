#include "hyperrmq/bitutil.hpp"

#include <random>

#include "gtest/gtest.h"

namespace average_case_optimal_rmq {
namespace {

TEST(PopCountTest, Small) {
    ASSERT_EQ(popcount(0), 0);
    ASSERT_EQ(popcount(1), 1);
    ASSERT_EQ(popcount(2), 1);
    ASSERT_EQ(popcount(3), 2);
    ASSERT_EQ(popcount(4), 1);
    ASSERT_EQ(popcount(5), 2);
    ASSERT_EQ(popcount(6), 2);
    ASSERT_EQ(popcount(7), 3);
    ASSERT_EQ(popcount(8), 1);
    ASSERT_EQ(popcount(9), 2);
    ASSERT_EQ(popcount(10), 2);
}

TEST(PopCountTest, Large) { ASSERT_EQ(popcount(0xffffffffffffffff), 64); }

TEST(PopCountTest, StressTest) {
    std::mt19937_64 engine(0);
    std::vector<uint64_t> values(10000);
    std::generate(values.begin(), values.end(), engine);

    for (auto&& v : values) {
        ASSERT_EQ(popcount(v), __builtin_popcountll(v));
    }
}

TEST(CeilLog2Test, Small) {
    ASSERT_EQ(ceil_log2(0), 0);
    ASSERT_EQ(ceil_log2(1), 0);
    ASSERT_EQ(ceil_log2(2), 1);
    ASSERT_EQ(ceil_log2(3), 2);
    ASSERT_EQ(ceil_log2(4), 2);
    ASSERT_EQ(ceil_log2(5), 3);
    ASSERT_EQ(ceil_log2(6), 3);
    ASSERT_EQ(ceil_log2(7), 3);
    ASSERT_EQ(ceil_log2(8), 3);
    ASSERT_EQ(ceil_log2(9), 4);
    ASSERT_EQ(ceil_log2(10), 4);
}

TEST(CeilLog2Test, PowerOfTwoAndPlusOne) {
    for (int shift = 0; shift < 64; shift++) {
        ASSERT_EQ(ceil_log2(1ull << shift), shift);
        ASSERT_EQ(ceil_log2((1ull << shift) + 1), shift + 1);
    }
}

TEST(CeilLog2Test, StressTest) {
    std::mt19937_64 engine(0);
    std::vector<uint64_t> values(10000);
    std::generate(values.begin(), values.end(), engine);

    for (auto&& v : values) {
        while (v) {
            auto c = ceil_log2(v);
            ASSERT_LE(0, c);
            ASSERT_LE(c, 64);
            if (c > 0) {
                ASSERT_LT(1ull << (c - 1), v);
            }
            if (c < 64) {
                ASSERT_LE(v, 1ull << c);
            }
            v >>= 1;
        }
    }
}

}  // namespace
}  // namespace average_case_optimal_rmq
