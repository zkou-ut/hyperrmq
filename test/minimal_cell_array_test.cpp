#include "minimal_cell_array.hpp"

#include "gtest/gtest.h"

namespace average_case_optimal_rmq {
namespace {

TEST(MinimalCellArrayTest, SmallUnsigned) {
    std::vector<uint64_t> vec_a = {5, 1, 4, 2, 6, 8};
    MinimalCellArray<uint64_t> a(vec_a), b(vec_a, 8, 0);

    ASSERT_EQ(a, b);
    ASSERT_TRUE(a == b);
    ASSERT_FALSE(a != b);

    ASSERT_EQ(vec_a.size(), a.size());
    for (int i = 0; i < vec_a.size(); i++) {
        ASSERT_EQ(vec_a[i], a.get(i));
    }

    ASSERT_EQ(a.size(), b.size());
    for (int i = 0; i < a.size(); i++) {
        ASSERT_EQ(a[i], b[i]);
    }

    auto vec_c = vec_a;
    vec_c[0] = 1;
    MinimalCellArray c(vec_c);

    ASSERT_NE(a, c);
    ASSERT_TRUE(a != c);
    ASSERT_FALSE(a == c);
}

TEST(MinimalCellArrayTest, SmallSigned) {
    std::vector<int64_t> vec_a = {0, 1, -4, 2, -6, 8};
    MinimalCellArray<int64_t> a(vec_a), b(vec_a, 10, -10);

    ASSERT_EQ(a, b);
    ASSERT_TRUE(a == b);
    ASSERT_FALSE(a != b);

    ASSERT_EQ(vec_a.size(), a.size());
    for (int i = 0; i < vec_a.size(); i++) {
        ASSERT_EQ(vec_a[i], a.get(i));
    }

    ASSERT_EQ(a.size(), b.size());
    for (int i = 0; i < a.size(); i++) {
        ASSERT_EQ(a[i], b[i]);
    }

    auto vec_c = vec_a;
    vec_c[0] = 1;
    MinimalCellArray c(vec_c);

    ASSERT_NE(a, c);
    ASSERT_TRUE(a != c);
    ASSERT_FALSE(a == c);
}

}  // namespace
}  // namespace average_case_optimal_rmq
