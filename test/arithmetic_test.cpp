#include "arithmetic.hpp"

#include <queue>
#include <random>

#include "bitutil.hpp"
#include "gtest/gtest.h"

namespace average_case_optimal_rmq {
namespace {

template <bool depth_first>
void arithmetic_test_helper(const TreeBP &bp,
                            const std::vector<uint32_t> &left_seq) {
    ASSERT_EQ(bp_to_left_seq<depth_first>(bp), left_seq);
    ASSERT_EQ(left_seq_to_bp<depth_first>(left_seq).to_string(),
              bp.to_string());

    auto arithmetic = left_seq_to_arithmetic<depth_first>(left_seq);
    ASSERT_EQ(arithmetic_to_left_seq<depth_first>(bp.n, arithmetic), left_seq);
}

TEST(ArithmeticTest, OneNode) {
    TreeBP bp("()");
    std::vector<uint32_t> left_seq = {0};
    // BitArray arithmetic(0);

    arithmetic_test_helper<true>(bp, left_seq);
    arithmetic_test_helper<false>(bp, left_seq);
}

TEST(ArithmeticTest, TwoNodesLeft) {
    TreeBP bp("(())");
    std::vector<uint32_t> left_seq = {1, 0};
    // BitArray arithmetic(1);
    // arithmetic.on(0);

    arithmetic_test_helper<true>(bp, left_seq);
    arithmetic_test_helper<false>(bp, left_seq);
}

TEST(ArithmeticTest, TwoNodesRight) {
    TreeBP bp("()()");
    std::vector<uint32_t> left_seq = {0, 0};
    // BitArray arithmetic(1);

    arithmetic_test_helper<true>(bp, left_seq);
    arithmetic_test_helper<false>(bp, left_seq);
}

TEST(ArithmeticTest, FiveNodes) {
    TreeBP bp("((())())()");
    std::vector<uint32_t> left_seq = {3, 1, 0, 0, 0};

    arithmetic_test_helper<true>(bp, left_seq);
    arithmetic_test_helper<false>(bp, left_seq);
}

TEST(ArithmeticTest, DepthStressTest) {
    std::mt19937 engine(0);

    for (int test_case = 0; test_case < 10000; test_case++) {
        int n = engine() % 1000 + 1;
        std::vector<uint32_t> left_seq;
        left_seq.reserve(n);
        std::stack<uint32_t> subtree_size_st;
        subtree_size_st.push(n);
        for (int i = 0; i < n; i++) {
            uint32_t cur_size = subtree_size_st.top();
            subtree_size_st.pop();
            uint32_t left_size = engine() % cur_size;
            uint32_t right_size = cur_size - left_size - 1;
            left_seq.push_back(left_size);
            if (right_size) {
                subtree_size_st.push(right_size);
            }
            if (left_size) {
                subtree_size_st.push(left_size);
            }
        }

        ASSERT_EQ(left_seq.size(), n);

        ASSERT_EQ(left_seq,
                  bp_to_left_seq<true>(left_seq_to_bp<true>(left_seq)));
        ASSERT_EQ(left_seq, arithmetic_to_left_seq<true>(
                                n, left_seq_to_arithmetic<true>(left_seq)));
    }
}

TEST(ArithmeticTest, BreadthStressTest) {
    std::mt19937 engine(0);

    for (int test_case = 0; test_case < 10000; test_case++) {
        int n = engine() % 1000 + 1;
        std::vector<uint32_t> left_seq;
        left_seq.reserve(n);
        std::queue<uint32_t> subtree_size_st;
        subtree_size_st.push(n);
        for (int i = 0; i < n; i++) {
            uint32_t cur_size = subtree_size_st.front();
            subtree_size_st.pop();
            uint32_t left_size = engine() % cur_size;
            uint32_t right_size = cur_size - left_size - 1;
            left_seq.push_back(left_size);
            if (left_size) {
                subtree_size_st.push(left_size);
            }
            if (right_size) {
                subtree_size_st.push(right_size);
            }
        }

        ASSERT_EQ(left_seq.size(), n);

        ASSERT_EQ(left_seq,
                  bp_to_left_seq<false>(left_seq_to_bp<false>(left_seq)));
        ASSERT_EQ(left_seq, arithmetic_to_left_seq<false>(
                                n, left_seq_to_arithmetic<false>(left_seq)));
    }
}

}  // namespace
}  // namespace average_case_optimal_rmq
