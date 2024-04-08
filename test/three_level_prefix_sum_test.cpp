#include "three_level_prefix_sum.hpp"

#include <random>

#include "gtest/gtest.h"

namespace average_case_optimal_rmq {
namespace {

TEST(ThreeLevelPrefixSumTest, Small) {
    ThreeLevelPrefixSum<uint32_t, 16, 16> pf({2, 1, 1, 4, 1, 5});

    ASSERT_EQ(pf.size(), 6);

    ASSERT_EQ(pf[0], 2);
    ASSERT_EQ(pf[1], 1);
    ASSERT_EQ(pf[2], 1);
    ASSERT_EQ(pf[3], 4);
    ASSERT_EQ(pf[4], 1);
    ASSERT_EQ(pf[5], 5);

    ASSERT_EQ(pf.sum(0), 0);
    ASSERT_EQ(pf.sum(1), 2);
    ASSERT_EQ(pf.sum(2), 3);
    ASSERT_EQ(pf.sum(3), 4);
    ASSERT_EQ(pf.sum(4), 8);
    ASSERT_EQ(pf.sum(5), 9);
    ASSERT_EQ(pf.sum(6), 14);

    ThreeLevelPrefixSum<uint32_t, 16, 16> pf_eq({2, 1, 1, 4, 1, 5});
    ASSERT_EQ(pf, pf_eq);
    ASSERT_TRUE(pf == pf_eq);
    ASSERT_FALSE(pf != pf_eq);

    ThreeLevelPrefixSum<uint32_t, 16, 16> pf_ne({2, 1, 1, 4, 5, 5});
    ASSERT_NE(pf, pf_ne);
    ASSERT_FALSE(pf == pf_ne);
    ASSERT_TRUE(pf != pf_ne);
}

TEST(ThreeLevelPrefixSumTest, StressTest) {
    std::mt19937 engine(0);
    for (int n = 0; n < 200; n++) {
        std::vector<uint32_t> values(n);
        for (int i = 0; i < n; i++) {
            values[i] = engine() % 1000;
        }

        ThreeLevelPrefixSum<uint32_t, 8, 8> pf(values);

        uint32_t sum = 0;
        ASSERT_EQ(pf.sum(0), sum);
        for (int i = 0; i < n; i++) {
            sum += values[i];
            ASSERT_EQ(pf.sum(i + 1), sum);
        }
    }
}

TEST(ThreeLevelPrefixSumTest, StressTestLarge) {
    const int n = 100000;
    std::mt19937 engine(0);
    std::vector<uint32_t> values(n);
    for (int i = 0; i < n; i++) {
        values[i] = engine() % 1000;
    }

    ThreeLevelPrefixSum<uint32_t, 8, 8> pf(values);

    uint32_t sum = 0;
    ASSERT_EQ(pf.sum(0), sum);
    for (int i = 0; i < n; i++) {
        sum += values[i];
        ASSERT_EQ(pf.sum(i + 1), sum);
    }
}

TEST(ThreeLevelPrefixSumTest, SelectChunkStressTest) {
    for (int n = 10; n <= 100; n += 10) {
        std::mt19937 engine(0);
        std::vector<uint32_t> values(n), csum(n + 1);
        for (int i = 0; i < n; i++) {
            values[i] = engine() % 3;
        }

        ThreeLevelPrefixSum<uint32_t, 8, 8> pf(values);
        int j = 0, s = 0;
        for (int v = 0; v < values.back() + 10; v++) {
            while (j < n && s + values[j] <= v) {
                s += values[j];
                j++;
            }
            ASSERT_EQ(pf.select_chunk(v), j);
        }
    }
}

}  // namespace
}  // namespace average_case_optimal_rmq
