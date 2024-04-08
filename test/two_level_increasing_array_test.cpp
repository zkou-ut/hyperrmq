#include "two_level_increasing_array.hpp"

#include <random>

#include "gtest/gtest.h"

namespace average_case_optimal_rmq {
namespace {

TEST(TwoLevelIncreasingArrayTest, Small) {
    TwoLevelIncreasingArray<> pf({2, 3, 4, 8, 9, 14});

    ASSERT_EQ(pf.size(), 6);

    ASSERT_EQ(pf[0], 2);
    ASSERT_EQ(pf[1], 3);
    ASSERT_EQ(pf[2], 4);
    ASSERT_EQ(pf[3], 8);
    ASSERT_EQ(pf[4], 9);
    ASSERT_EQ(pf[5], 14);

    TwoLevelIncreasingArray<> pf_eq({2, 3, 4, 8, 9, 14});
    ASSERT_EQ(pf, pf_eq);
    ASSERT_TRUE(pf == pf_eq);
    ASSERT_FALSE(pf != pf_eq);

    TwoLevelIncreasingArray<> pf_ne({2, 3, 4, 8, 9});
    ASSERT_NE(pf, pf_ne);
    ASSERT_FALSE(pf == pf_ne);
    ASSERT_TRUE(pf != pf_ne);
}

TEST(TwoLevelIncreasingArrayTest, StressTest) {
    for (int n = 0; n < 100; n++) {
        std::mt19937 engine(0);
        std::vector<uint32_t> values(n);
        for (int i = 0; i < n; i++) {
            values[i] = engine() % 1000;
        }
        for (int i = 0; i < n - 1; i++) {
            values[i + 1] += values[i];
        }

        TwoLevelIncreasingArray<> pf(values);
        for (int i = 0; i < n; i++) {
            ASSERT_EQ(pf[i], values[i]);
        }
    }
}

TEST(TwoLevelIncreasingArrayTest, StressTestLarge) {
    const int n = 100000;
    std::mt19937 engine(0);
    std::vector<uint32_t> values(n);
    for (int i = 0; i < n; i++) {
        values[i] = engine() % 1000;
    }
    for (int i = 0; i < n - 1; i++) {
        values[i + 1] += values[i];
    }

    TwoLevelIncreasingArray<> pf(values);
    for (int i = 0; i < n; i++) {
        ASSERT_EQ(pf[i], values[i]);
    }
}

}  // namespace
}  // namespace average_case_optimal_rmq
