#include "hyperrmq/runs.hpp"

#include <random>
#include <stack>

#include "gtest/gtest.h"
#include "hyperrmq/arithmetic.hpp"
#include "hyperrmq/bitutil.hpp"

namespace average_case_optimal_rmq {
namespace {

TEST(RunsTest, CountIncreasingRuns) {
    ASSERT_EQ(count_increasing_runs<int>({}), 0);

    ASSERT_EQ(count_increasing_runs<int>({1}), 1);

    ASSERT_EQ(count_increasing_runs<int>({1, 2}), 1);
    ASSERT_EQ(count_increasing_runs<int>({2, 1}), 2);

    ASSERT_EQ(count_increasing_runs<int>({1, 2, 3}), 1);
    ASSERT_EQ(count_increasing_runs<int>({3, 2, 1}), 3);
    ASSERT_EQ(count_increasing_runs<int>({1, 3, 2}), 2);
    ASSERT_EQ(count_increasing_runs<int>({2, 3, 1}), 2);
    ASSERT_EQ(count_increasing_runs<int>({2, 1, 3}), 2);
    ASSERT_EQ(count_increasing_runs<int>({3, 1, 2}), 2);

    ASSERT_EQ(count_increasing_runs<int>({9, 2, 5, 3, 4, 8, 7, 1, 6}), 5);
}

TEST(RunsTest, RandomExactFixedIncreasingRuns) {
    const int n = 100, sqrtn = 10;
    std::vector<int> perm(n);

    std::iota(perm.begin(), perm.end(), 0);

    std::mt19937 engine(0);

    for (int r = 1; r <= sqrtn; r++) {
        random_exact_fixed_incresing_runs(perm, r, false);
        ASSERT_EQ(count_increasing_runs(perm), r);
    }
}

}  // namespace
}  // namespace average_case_optimal_rmq
