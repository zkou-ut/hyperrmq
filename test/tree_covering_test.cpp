
#include "hyperrmq/tree_covering.hpp"

#include <random>
#include <stack>

#include "gtest/gtest.h"
#include "hyperrmq/bitutil.hpp"
#include "hyperrmq/cartesian_tree.hpp"

namespace hyperrmq {
namespace {

TEST(TreeCoveringTest, LCAOnMicrotreeTest) {
    std::vector<int> perm = {3, 1, 4, 1, 5};
    int n = perm.size();
    auto [upsilon, microtrees] =
        legacy_tree_covering(n, cartesian_tree_bp(perm));
    ASSERT_EQ(microtrees.size(), 1);
    auto microtree = microtrees[0];
    auto [bp, cutpos] = legacy_decode_microtree(microtree);

    auto rmq = [&](int first, int last) -> int {
        if (first > last) {
            std::swap(first, last);
        }
        auto min_val = perm[first];
        auto min_idx = first;
        for (int i = first + 1; i < last + 1; i++) {
            if (min_val > perm[i]) {
                min_val = perm[i];
                min_idx = i;
            }
        }
        return min_idx;
    };

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            ASSERT_EQ(legacy_lca_on_microtree(bp, i, j), rmq(i, j));
        }
    }
}

TEST(TreeCoveringTest, MicrotreeRootsSeventyNodes) {
    const int32_t B = 6;
    std::vector<int32_t> in_to_pre = {
        6,  7,  5,  9,  8,  4,  12, 11, 13, 14, 15, 10, 3,  18, 17, 16, 21, 22,
        20, 24, 25, 23, 19, 27, 26, 30, 29, 28, 33, 32, 31, 2,  35, 34, 36, 1,
        39, 42, 41, 40, 43, 38, 46, 45, 44, 49, 48, 51, 50, 52, 47, 55, 54, 56,
        53, 37, 58, 59, 57, 63, 62, 65, 64, 61, 67, 66, 60, 70, 69, 68};
    std::set<int32_t> roots_in = {5,  11, 12, 15, 18, 22, 24, 27,
                                  31, 35, 41, 50, 55, 58, 63};

    TreeBP ct = cartesian_tree_bp(in_to_pre);
    BitArray roots_pre = microtree_roots_preorder(B, ct);
    for (int i = 0; i < ct.n; i++) {
        ASSERT_EQ(roots_in.count(i), roots_pre.get(in_to_pre[i] - 1));
    }
}

void tree_covering_test_helper(const int32_t B, const std::vector<int> &values,
                               std::set<int32_t> &expected_roots_val) {
    TreeBP ct = cartesian_tree_bp(values);

    std::vector<int> values_inorder = values;
    std::vector<int> values_preorder;
    std::stack<int> st;
    for (int i = ct.bp.size() - 1; i >= 0; i--) {
        if (ct.bp.get(i)) {
            st.push(values_inorder.back());
            values_inorder.pop_back();
        } else {
            values_preorder.push_back(st.top());
            st.pop();
        }
    }
    std::reverse(values_preorder.begin(), values_preorder.end());

    BitArray actual_roots_pre_bitarray = microtree_roots_preorder(B, ct);
    std::set<int32_t> actual_roots_val;
    for (int i = 0; i < actual_roots_pre_bitarray.size(); i++) {
        if (actual_roots_pre_bitarray.get(i)) {
            actual_roots_val.insert(values_preorder[i]);
        }
    }
    ASSERT_EQ(expected_roots_val, actual_roots_val);
}

TEST(TreeCoveringTest, MicrotreeRootsTwoHeavyChildren) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 4, 9, 3, 2, 1, 6, 5, 7, 8};
    std::set<int32_t> roots_val = {0, 1, 2, 5};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(
    TreeCoveringTest,
    MicrotreeRootsPermanentHeavyChildAndTemporaryLightChildAndComponentSizeEqualsBMinusOne) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 3, 2, 1, 6, 5, 7, 8};
    std::set<int32_t> roots_val = {0, 5};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(
    TreeCoveringTest,
    MicrotreeRootsPermanentHeavyChildAndTemporaryLightChildAndComponentSizeEqualsB) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 4, 3, 2, 1, 6, 5, 7, 8};
    std::set<int32_t> roots_val = {0, 1, 5};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(
    TreeCoveringTest,
    MicrotreeRootsTemporaryHeavyChildAndTemporaryLightChildAndComponentSizeEqualsBMinusOne) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 2, 1, 5, 9, 7, 8, 10};
    std::set<int32_t> roots_val = {0, 7};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(
    TreeCoveringTest,
    MicrotreeRootsTemporaryHeavyChildAndTemporaryLightChildAndComponentSizeEqualsB) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 3, 2, 1, 5, 9, 7, 8, 10};
    std::set<int32_t> roots_val = {0, 1, 7};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(TreeCoveringTest, MicrotreeRootsTwoLightChildrenAndComponentSizeEqualsB) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 3, 2, 1, 4};
    std::set<int32_t> roots_val = {0, 1};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(TreeCoveringTest,
     MicrotreeRootsTwoLightChildrenAndComponentSizeEqualsBMinusOne) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 2, 1, 4};
    std::set<int32_t> roots_val = {0};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(TreeCoveringTest, MicrotreeRootsOnePermanentHeavyChild) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 3, 2, 4, 5, 1};
    std::set<int32_t> roots_val = {0, 2};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(TreeCoveringTest,
     MicrotreeRootsOneTemporaryHeavyChildAndComponentSizeEqualsBMinusOne) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 7, 6, 5, 4, 3, 2, 1};
    std::set<int32_t> roots_val = {0, 4};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(TreeCoveringTest,
     MicrotreeRootsOneTemporaryHeavyChildAndComponentSizeEqualsB) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 8, 7, 6, 5, 4, 3, 2, 1};
    std::set<int32_t> roots_val = {0, 1, 5};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(TreeCoveringTest,
     MicrotreeRootsOneLightChildAndComponentSizeEqualsBMinusOne) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 3, 2, 1};
    std::set<int32_t> roots_val = {0};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(TreeCoveringTest, MicrotreeRootsOneLightChildAndComponentSizeEqualsB) {
    const int32_t B = 4;
    std::vector<int32_t> values = {0, 4, 3, 2, 1};
    std::set<int32_t> roots_val = {0, 1};
    tree_covering_test_helper(B, values, roots_val);
}

TEST(TreeCoveringTest, MicrotreeRootsStressTest) {
    const int n = 1000;
    std::mt19937 engine(0);
    for (int test_case = 0; test_case < 10; test_case++) {
        for (int B = 2; B <= 10; B++) {
            std::vector<int> values(n);
            for (int i = 0; i < n; i++) {
                values[i] = engine() % n;
            }
            CartesianTree ct(values);
            TreeBP ctbp = cartesian_tree_bp(values);
            std::set<int> roots_in =
                legacy_microtree_roots_inorder(B, ct.binary_tree);
            BitArray roots_pre = microtree_roots_preorder(B, ctbp);

            std::vector<int> in_to_pre(n);
            auto pre = 0;
            auto dfs = [&](auto self, int v) -> void {
                in_to_pre[v] = pre++;
                int l = ct.binary_tree.left[v];
                if (l != -1) {
                    self(self, l);
                }
                int r = ct.binary_tree.right[v];
                if (r != -1) {
                    self(self, r);
                }
            };
            dfs(dfs, ct.binary_tree.root);
            for (int in = 0; in < n; in++) {
                bool ct_is_root = roots_in.count(in);
                bool bp_is_root = roots_pre.get(in_to_pre[in]);
                ASSERT_EQ(ct_is_root, bp_is_root);
            }
        }
    }
}

TEST(TreeCoveringTest, BuildUpsilonSeventyNodes) {
    const int32_t B = 6;
    std::vector<int32_t> in_to_pre = {
        6,  7,  5,  9,  8,  4,  12, 11, 13, 14, 15, 10, 3,  18, 17, 16, 21, 22,
        20, 24, 25, 23, 19, 27, 26, 30, 29, 28, 33, 32, 31, 2,  35, 34, 36, 1,
        39, 42, 41, 40, 43, 38, 46, 45, 44, 49, 48, 51, 50, 52, 47, 55, 54, 56,
        53, 37, 58, 59, 57, 63, 62, 65, 64, 61, 67, 66, 60, 70, 69, 68};
    TreeBP ct = cartesian_tree_bp(in_to_pre);
    BitArray roots = microtree_roots_preorder(B, ct);
    TreeBP upsilon = legacy_build_upsilon(roots.linear_popcount(), roots, ct);
    ASSERT_EQ(upsilon.to_string(), "((((()))((())(()))))((()))(())");
}

TEST(TreeCoveringTest, BuildUpsilonLeftPackedCheck) {
    auto dfs = [&](auto self, const TreeBP &tree, int &index) -> void {
        auto match = [&](int v) -> void {
            assert(tree.bp.get(index) == v);
            index++;
        };
        auto has_child = [&]() -> bool {
            return index < 2 * tree.n && tree.bp.get(index) == 0;
        };
        match(0);
        bool has_l = has_child();
        if (has_l) {
            self(self, tree, index);
        }

        match(1);
        bool has_r = has_child();
        if (has_r) {
            self(self, tree, index);
        }

        // has_r implies has_l
        ASSERT_TRUE(!has_r || has_l);
    };

    std::mt19937 engine(0);
    for (int shift = 5; shift <= 15; shift++) {
        int n = 1 << shift;
        for (int test_case = 0; test_case < 10; test_case++) {
            for (int B = 2; B <= 10; B += 4) {
                std::vector<int> values(n);
                for (int i = 0; i < n; i++) {
                    values[i] = engine() % n;
                }
                TreeBP ct = cartesian_tree_bp(values);

                BitArray roots = microtree_roots_preorder(B, ct);
                auto microtree_count = roots.linear_popcount();
                TreeBP upsilon =
                    legacy_build_upsilon(microtree_count, roots, ct);

                // It fails with the following code:
                // MicrotreeSplitRankArray microtree_split_rank_array =
                //     enumerate_microtrees(B, microtree_count, roots, ct);
                // TreeBP upsilon = build_upsilon(microtree_count,
                // roots,
                //                                ct,
                //                                microtree_split_rank_array);

                int index = 0;
                dfs(dfs, upsilon, index);
                ASSERT_EQ(index, upsilon.bp.size());
            }
        }
    }
}

TEST(TreeCoveringTest, EnumerateMicrotreesSeventyNodes) {
    const int32_t B = 6;
    std::vector<int32_t> in_to_pre = {
        6,  7,  5,  9,  8,  4,  12, 11, 13, 14, 15, 10, 3,  18, 17, 16, 21, 22,
        20, 24, 25, 23, 19, 27, 26, 30, 29, 28, 33, 32, 31, 2,  35, 34, 36, 1,
        39, 42, 41, 40, 43, 38, 46, 45, 44, 49, 48, 51, 50, 52, 47, 55, 54, 56,
        53, 37, 58, 59, 57, 63, 62, 65, 64, 61, 67, 66, 60, 70, 69, 68};
    TreeBP ct = cartesian_tree_bp(in_to_pre);
    BitArray roots = microtree_roots_preorder(B, ct);
    std::vector<uint64_t> microtrees =
        legacy_enumerate_microtrees(roots.linear_popcount(), roots, ct);
    std::vector<std::pair<std::string, int>> actual;
    for (auto &&m : microtrees) {
        auto [bp_int, cut_pos] = legacy_decode_microtree(m);
        std::string bp_str;
        int length = popcount(bp_int) * 2;
        for (int i = 0; i < length; i++) {
            bp_str += ((((bp_int >> i) & 1) == 0) ? "(" : ")");
        }
        actual.push_back({bp_str, cut_pos});
    }

    std::vector<std::pair<std::string, int>> expected = {
        {"((())()()())", 12},          // 12
        {"((()())(()))", 12},          // 6
        {"()", 1},                     // 13
        {"(()())(()())", 12},          // 19
        {"()", 1},                     // 23
        {"((()))((()))", 12},          // 28
        {"(())", 4},                   // 25
        {"((()))", 6},                 // 16
        {"()(())()", 1},               // 32
        {"()", 1},                     // 36
        {"((())(())())((())())", 20},  // 51
        {"(()((()))())((()))", 18},    // 42
        {"()", 1},                     // 56
        {"((())(()))(())", 14},        // 64
        {"(()())()((()))", 7},         // 59
    };

    ASSERT_EQ(actual, expected);
}

TEST(TreeCoveringTest, EnumerateLargeMicrotreesBPSeventyNodes) {
    const int32_t B = 6;
    std::vector<int32_t> in_to_pre = {
        6,  7,  5,  9,  8,  4,  12, 11, 13, 14, 15, 10, 3,  18, 17, 16, 21, 22,
        20, 24, 25, 23, 19, 27, 26, 30, 29, 28, 33, 32, 31, 2,  35, 34, 36, 1,
        39, 42, 41, 40, 43, 38, 46, 45, 44, 49, 48, 51, 50, 52, 47, 55, 54, 56,
        53, 37, 58, 59, 57, 63, 62, 65, 64, 61, 67, 66, 60, 70, 69, 68};
    TreeBP ct = cartesian_tree_bp(in_to_pre);

    BitArray roots = microtree_roots_preorder(B, ct);
    auto [upsilon, microtree_cut_array] =
        recover_upsilon_and_microtrees(B, roots, ct);
    std::vector<std::pair<std::string, uint64_t>> actual;
    for (int index = 0; index < upsilon.n; index++) {
        const auto &[microtree, cutpos] = microtree_cut_array[index];
        actual.push_back({microtree.to_string(), cutpos});
    }

    std::vector<std::pair<std::string, uint64_t>> expected = {
        {"()", 0},                     // 36
        {"()(())()", 0},               // 32
        {"()", 0},                     // 13
        {"((()())(()))", 6},           // 6
        {"((())()()())", 6},           // 12
        {"((()))", 3},                 // 16
        {"()", 0},                     // 23
        {"(()())(()())", 6},           // 19
        {"(())", 2},                   // 25
        {"((()))((()))", 6},           // 28
        {"()", 0},                     // 56
        {"(()((()))())((()))", 9},     // 42
        {"((())(())())((())())", 10},  // 51
        {"(()())()((()))", 3},         // 59
        {"((())(()))(())", 7},         // 64
    };

    ASSERT_EQ(actual, expected);
}

void tree_covering_restore_bp_test_helper(int B, const TreeBP &ct) {
    const auto [upsilon, microtrees] = legacy_tree_covering(B, ct);

    std::stack<int> st;
    std::vector<int> refs(upsilon.bp.size());
    int rank = 0;
    for (int i = 0; i < upsilon.bp.size(); i++) {
        if (upsilon.bp.get(i) == 0) {
            st.push(i);
        } else {
            refs[st.top()] = refs[i] = rank++;
            st.pop();
        }
    }

    BitArray restored_bp(2 * ct.n);
    int pos = 0;
    for (int i = 0; i < upsilon.bp.size(); i++) {
        auto [mbp, cut_pos] = legacy_decode_microtree(microtrees[refs[i]]);
        if (upsilon.bp.get(i) == 0) {
            for (int j = 0; j < cut_pos; j++) {
                restored_bp.set(pos++, (mbp >> j) & 1);
            }
        } else {
            auto length = popcount(mbp) * 2;
            for (int j = cut_pos; j < length; j++) {
                restored_bp.set(pos++, (mbp >> j) & 1);
            }
        }
    }

    ASSERT_EQ(ct.bp, restored_bp);
}

void tree_covering_large_microtree_restore_bp_test_helper(int B,
                                                          const TreeBP &ct) {
    const auto [upsilon, tree_and_split_rank] = tree_covering(B, ct);

    std::stack<int> st;
    std::vector<int> refs(upsilon.bp.size());
    int rank0 = upsilon.n;
    for (int i = upsilon.bp.size() - 1; i >= 0; i--) {
        if (upsilon.bp.get(i) == 1) {
            st.push(i);
        } else {
            refs[st.top()] = refs[i] = --rank0;
            st.pop();
        }
    }

    BitArray restored_bp(2 * ct.n);
    int pos = 0;
    for (int i = 0; i < upsilon.bp.size(); i++) {
        const auto &[microtree, split_rank] = tree_and_split_rank[refs[i]];
        auto cut = microtree.bp.linear_select1(split_rank);
        if (upsilon.bp.get(i) == 0) {
            restored_bp.write_interval(pos, microtree.bp.read_interval(0, cut));
            pos += cut;
        } else {
            restored_bp.write_interval(
                pos, microtree.bp.read_interval(cut, 2 * microtree.n - cut));
            pos += 2 * microtree.n - cut;
        }
    }

    ASSERT_EQ(ct.bp, restored_bp);
}

TEST(TreeCoveringTest, TreeCoveringRestoreBPSeventyNodes) {
    const int32_t B = 6;
    std::vector<int32_t> in_to_pre = {
        6,  7,  5,  9,  8,  4,  12, 11, 13, 14, 15, 10, 3,  18, 17, 16, 21, 22,
        20, 24, 25, 23, 19, 27, 26, 30, 29, 28, 33, 32, 31, 2,  35, 34, 36, 1,
        39, 42, 41, 40, 43, 38, 46, 45, 44, 49, 48, 51, 50, 52, 47, 55, 54, 56,
        53, 37, 58, 59, 57, 63, 62, 65, 64, 61, 67, 66, 60, 70, 69, 68};
    TreeBP ct = cartesian_tree_bp(in_to_pre);

    tree_covering_restore_bp_test_helper(B, ct);
}

TEST(TreeCoveringTest, TreeCoveringRestoreBPStressTest) {
    std::mt19937 engine(0);
    for (int shift = 5; shift <= 15; shift++) {
        int n = 1 << shift;
        for (int B = 2; B <= 15; B++) {
            std::vector<int> values(n);
            for (int i = 0; i < n; i++) {
                values[i] = engine() % n;
            }
            TreeBP ct = cartesian_tree_bp(values);

            tree_covering_restore_bp_test_helper(B, ct);
        }
    }
}

TEST(TreeCoveringTest, TreeCoveringLargeMicrotreeRestoreBPSeventyNodes) {
    const int32_t B = 6;
    std::vector<int32_t> in_to_pre = {
        6,  7,  5,  9,  8,  4,  12, 11, 13, 14, 15, 10, 3,  18, 17, 16, 21, 22,
        20, 24, 25, 23, 19, 27, 26, 30, 29, 28, 33, 32, 31, 2,  35, 34, 36, 1,
        39, 42, 41, 40, 43, 38, 46, 45, 44, 49, 48, 51, 50, 52, 47, 55, 54, 56,
        53, 37, 58, 59, 57, 63, 62, 65, 64, 61, 67, 66, 60, 70, 69, 68};
    TreeBP ct = cartesian_tree_bp(in_to_pre);

    tree_covering_large_microtree_restore_bp_test_helper(B, ct);
}

TEST(TreeCoveringTest, TreeCoveringLargeMicrotreesRestoreBPStressTest) {
    std::mt19937 engine(0);
    for (int shift = 5; shift <= 15; shift++) {
        int n = 1 << shift;
        for (int i = 1; i <= 10; i++) {
            int B = 1 << i;
            std::vector<int> values(n);
            for (int i = 0; i < n; i++) {
                values[i] = engine() % n;
            }
            TreeBP ct = cartesian_tree_bp(values);

            tree_covering_large_microtree_restore_bp_test_helper(B, ct);
        }
    }
}

TEST(TreeCoveringTest, TreeCoveringRestoreBPLargeBStressTest) {
    std::mt19937 engine(0);
    const int n = 1000000;
    const int B = 100000;

    std::vector<int> values(n);
    for (int i = 0; i < n; i++) {
        values[i] = engine() % n;
    }
    TreeBP ct = cartesian_tree_bp(values);

    tree_covering_large_microtree_restore_bp_test_helper(B, ct);
}

}  // namespace
}  // namespace hyperrmq
