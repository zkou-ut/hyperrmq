#include "hyperrmq/rmq_bp.hpp"

#include <random>

#include "gtest/gtest.h"

namespace hyperrmq {
namespace {

TEST(RMQBPTest, RMQSmall) {
    std::vector<int64_t> perm = {3, 1, 4, 1, 5};

    RMQBP<> rmq_bp(perm);
    ASSERT_EQ(rmq_bp.query(0, 0), 0);
    ASSERT_EQ(rmq_bp.query(0, 1), 1);
    ASSERT_EQ(rmq_bp.query(0, 2), 1);
    ASSERT_EQ(rmq_bp.query(0, 3), 1);
    ASSERT_EQ(rmq_bp.query(0, 4), 1);

    ASSERT_EQ(rmq_bp.query(1, 1), 1);
    ASSERT_EQ(rmq_bp.query(1, 2), 1);
    ASSERT_EQ(rmq_bp.query(1, 3), 1);
    ASSERT_EQ(rmq_bp.query(1, 4), 1);

    ASSERT_EQ(rmq_bp.query(2, 2), 2);
    ASSERT_EQ(rmq_bp.query(2, 3), 3);
    ASSERT_EQ(rmq_bp.query(2, 4), 3);

    ASSERT_EQ(rmq_bp.query(3, 3), 3);
    ASSERT_EQ(rmq_bp.query(3, 4), 3);

    ASSERT_EQ(rmq_bp.query(4, 4), 4);
}

TEST(RMQBPTest, QuerySeventyTest) {
    const int64_t n = 70;
    std::vector<int64_t> perm = {
        6,  7,  5,  9,  8,  4,  12, 11, 13, 14, 15, 10, 3,  18, 17, 16, 21, 22,
        20, 24, 25, 23, 19, 27, 26, 30, 29, 28, 33, 32, 31, 2,  35, 34, 36, 1,
        39, 42, 41, 40, 43, 38, 46, 45, 44, 49, 48, 51, 50, 52, 47, 55, 54, 56,
        53, 37, 58, 59, 57, 63, 62, 65, 64, 61, 67, 66, 60, 70, 69, 68};

    RMQBP<> rmq_bp(perm);
    for (int64_t i = 0; i < n; i++) {
        int64_t min_val = perm[i], min_idx = i;
        for (int64_t j = i; j < n; j++) {
            if (min_val > perm[j]) {
                min_val = perm[j];
                min_idx = j;
            }
            ASSERT_EQ(min_idx, rmq_bp.query(i, j));
        }
    }
}

TEST(RMQBPTest, RMQPermutationStressTest) {
    const int64_t n = 300;
    std::mt19937 engine(0);
    std::vector<int64_t> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    RMQBP<8, 32> rmqbp(perm);
    for (int64_t i = 0; i < n; i++) {
        int64_t min_val = perm[i], min_idx = i;
        for (int64_t j = i; j < n; j++) {
            if (min_val > perm[j]) {
                min_val = perm[j];
                min_idx = j;
            }
            ASSERT_EQ(min_idx, rmqbp.query(i, j));
        }
    }
}

TEST(RMQBPTest, RMQTieStressTest) {
    const int64_t n = 300;
    std::mt19937 engine(0);
    std::vector<int64_t> values(n);

    for (int64_t i = 0; i < n; i++) {
        values[i] = engine() % 10;
    }

    RMQBP<> rmqbp(values);
    for (int64_t i = 0; i < n; i++) {
        int64_t min_val = values[i], min_idx = i;
        for (int64_t j = i; j < n; j++) {
            if (min_val > values[j]) {
                min_val = values[j];
                min_idx = j;
            }
            ASSERT_EQ(min_idx, rmqbp.query(i, j));
        }
    }
}

}  // namespace
}  // namespace hyperrmq
