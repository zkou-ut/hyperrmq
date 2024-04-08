#include "hyperrmq/tree_covering.hpp"

#include <cassert>
#include <iostream>
#include <stack>
#include <tuple>

#include "hyperrmq/bitutil.hpp"

namespace average_case_optimal_rmq {

uint64_t legacy_encode_microtree(uint64_t bp, uint64_t cut_pos) {
    assert(bp << 6 >> 6 == bp);
    assert(cut_pos < 64);
    return (bp << 6) | cut_pos;
}

std::pair<uint64_t, uint64_t> legacy_decode_microtree(uint64_t microtree) {
    uint64_t bp = microtree >> 6;
    uint64_t cut_pos = microtree & 0x3f;
    return {bp, cut_pos};
}

uint32_t legacy_lca_on_microtree(uint64_t bp, uint32_t l, uint32_t m) {
    if (l == m) {
        return l;
    }
    if (l > m) {
        std::swap(l, m);
    }
    uint64_t len = popcount(bp) * 2;
    int32_t minval = 1e9, minclose = 0;
    int32_t close = 0;
    l++;
    m++;
    for (int i = 0; i < len; i++) {
        if (l <= close) {
            if (m < close) {
                break;
            }
            int excess = i - 2 * close;
            if (minval > excess) {
                minval = excess;
                minclose = close;
            }
        }
        close += (bp >> i) & 1;
    }
    return minclose - 1;
};

std::set<int32_t> legacy_microtree_roots_inorder(
    int32_t B, const BinaryTree &binary_tree) {
    std::set<int32_t> roots = {binary_tree.root};

    // returns (component size, subtree size)
    auto dfs = [&](auto self, int32_t v) -> std::pair<int32_t, int32_t> {
        int32_t l = binary_tree.left[v], r = binary_tree.right[v];

        int32_t l_sbt_sz = 0, l_cmp_sz = 0;
        if (l != binary_tree.None) {
            std::tie(l_cmp_sz, l_sbt_sz) = self(self, l);
        }

        int32_t r_sbt_sz = 0, r_cmp_sz = 0;
        if (r != binary_tree.None) {
            std::tie(r_cmp_sz, r_sbt_sz) = self(self, r);
        }

        int32_t cmp_sz = l_cmp_sz + r_cmp_sz + 1;
        int32_t sbt_sz = l_sbt_sz + r_sbt_sz + 1;
        if (l_sbt_sz >= B && r_sbt_sz >= B) {
            roots.insert(l);
            roots.insert(r);
            roots.insert(v);
            cmp_sz = 0;
        }
        if (cmp_sz >= B) {
            roots.insert(v);
            cmp_sz = 0;
        }
        return {cmp_sz, sbt_sz};
    };

    dfs(dfs, binary_tree.root);

    return roots;
}

BitArray microtree_roots_preorder(int64_t B, const TreeBP &tree) {
    BitArray res(tree.n);
    res.on(0);

    int64_t pre = 0;
    int64_t index = 0;

    auto match = [&](int v) -> void {
        assert(tree.bp.get(index) == v);
        index++;
    };
    auto has_child = [&]() -> bool {
        return index < 2 * tree.n && tree.bp.get(index) == 0;
    };
    // returns (component size, subtree size)
    auto dfs = [&](auto self) -> std::pair<int64_t, int64_t> {
        int v = pre;
        pre++;

        match(0);
        bool has_l = has_child();
        int l = pre;
        int64_t l_sbt_sz = 0, l_cmp_sz = 0;
        if (has_l) {
            std::tie(l_cmp_sz, l_sbt_sz) = self(self);
        }

        match(1);
        int r = pre;
        int64_t r_sbt_sz = 0, r_cmp_sz = 0;
        if (index < 2 * tree.n && tree.bp.get(index) == 0) {
            std::tie(r_cmp_sz, r_sbt_sz) = self(self);
        }

        int64_t cmp_sz = l_cmp_sz + r_cmp_sz + 1;
        int64_t sbt_sz = l_sbt_sz + r_sbt_sz + 1;
        if (l_sbt_sz >= B && r_sbt_sz >= B) {
            res.on(l);
            res.on(r);
            res.on(v);
            cmp_sz = 0;
        }
        if (cmp_sz >= B) {
            res.on(v);
            cmp_sz = 0;
        }
        return {cmp_sz, sbt_sz};
    };

    dfs(dfs);

    return res;
}

TreeBP legacy_build_upsilon(const uint64_t microtree_count,
                            const BitArray &microtree_roots,
                            const TreeBP &tree) {
    BitArray upsilon_bp(2 * microtree_count);
    uint64_t up_index = 0;

    uint64_t pre = 0, index = 0;
    auto match = [&](uint64_t v) -> void {
        assert(tree.bp.get(index) == v);
        index++;
    };
    auto has_child = [&]() -> bool {
        return index < 2 * tree.n && tree.bp.get(index) == 0;
    };
    auto dfs = [&](auto self) -> void {
        bool is_microtree_root = microtree_roots.get(pre);
        pre++;

        if (is_microtree_root) {
            upsilon_bp.off(up_index++);
        }

        uint64_t last_up_index = up_index;
        match(0);
        bool has_l = has_child();
        if (has_l) {
            self(self);
        }

        auto has_left_microtree = (last_up_index < up_index);

        if (is_microtree_root && has_left_microtree) {
            upsilon_bp.on(up_index++);
        }

        match(1);
        bool has_r = has_child();
        if (has_r) {
            self(self);
        }

        if (is_microtree_root && !has_left_microtree) {
            upsilon_bp.on(up_index++);
        }

        return;
    };

    dfs(dfs);

    assert(up_index == 2 * microtree_count);

    return TreeBP(microtree_count, std::move(upsilon_bp));
}

TreeBP build_upsilon(const uint64_t microtree_count,
                     const BitArray &microtree_roots, const TreeBP &tree,
                     const MicrotreeSplitRankArray &micro_split_rank_array) {
    BitArray upsilon_bp(2 * microtree_count);
    uint64_t up_index = 0, up_pre = 0;

    uint64_t pre = 0, index = 0;
    auto match = [&](uint64_t v) -> void {
        assert(tree.bp.get(index) == v);
        index++;
    };
    auto has_child = [&]() -> bool {
        return index < 2 * tree.n && tree.bp.get(index) == 0;
    };
    auto dfs = [&](auto self) -> void {
        bool is_microtree_root = microtree_roots.get(pre);
        pre++;

        bool early_close = false;
        if (is_microtree_root) {
            upsilon_bp.off(up_index++);
            const auto &[tree, split_rank] = micro_split_rank_array[up_pre];
            early_close = (tree.n == split_rank);
            up_pre++;
        }

        uint64_t last_up_index = up_index;
        match(0);
        bool has_l = has_child();
        if (has_l) {
            self(self);
        }

        auto has_left_microtree = (last_up_index < up_index);

        if (is_microtree_root && (early_close || has_left_microtree)) {
            upsilon_bp.on(up_index++);
        }

        match(1);
        bool has_r = has_child();
        if (has_r) {
            self(self);
        }

        if (is_microtree_root && !early_close && !has_left_microtree) {
            upsilon_bp.on(up_index++);
        }

        return;
    };

    dfs(dfs);

    assert(up_index == 2 * microtree_count);

    return TreeBP(microtree_count, std::move(upsilon_bp));
}

std::vector<uint64_t> legacy_enumerate_microtrees(
    const uint64_t microtree_count, const BitArray &microtree_roots,
    const TreeBP &tree) {
    std::vector<uint64_t> microtrees(microtree_count);

    uint64_t pre = 0, index = 0;
    auto match = [&](uint64_t v) -> void {
        assert(tree.bp.get(index) == v);
        index++;
    };
    auto has_child = [&]() -> bool {
        return index < 2 * tree.n && tree.bp.get(index) == 0;
    };

    uint64_t up_in = 0;
    uint64_t root_count = 0;
    // returns bp and cutpos
    auto dfs = [&](auto self) -> std::pair<uint64_t, uint64_t> {
        bool is_microtree_root = microtree_roots.get(pre);
        root_count += is_microtree_root;
        pre++;

        uint64_t bp = 0;
        uint64_t write_pos = 0, cut_pos = 0;

        uint64_t last_root_count = root_count;
        match(0);
        write_pos++;
        bool has_l = has_child();
        if (has_l) {
            auto [bp_l, cut_pos_l] = self(self);
            if (bp_l == 0 || cut_pos_l != 0) {
                cut_pos = write_pos + cut_pos_l;
            }
            bp |= bp_l << write_pos;
            write_pos += popcount(bp_l) * 2;
        }

        auto has_left_microtree = (last_root_count < root_count);

        uint64_t up_v;
        if (is_microtree_root && has_left_microtree) {
            up_v = up_in++;
        }

        match(1);
        bp |= 1ull << write_pos;
        write_pos++;
        bool has_r = has_child();
        if (has_r) {
            auto [bp_r, cut_pos_r] = self(self);
            if (cut_pos == 0 && (bp_r == 0 || cut_pos_r != 0)) {
                cut_pos = write_pos + cut_pos_r;
            }
            bp |= bp_r << write_pos;
            write_pos += popcount(bp_r) * 2;
        }

        if (is_microtree_root && !has_left_microtree) {
            up_v = up_in++;
        }

        if (is_microtree_root) {
            if (cut_pos == 0) {
                cut_pos = popcount(bp) * 2;
            }
            microtrees[up_v] = legacy_encode_microtree(bp, cut_pos);
            bp = 0;
            cut_pos = 0;
        }
        return {bp, cut_pos};
    };

    dfs(dfs);

    assert(up_in == microtree_count);

    return microtrees;
}

MicrotreeSplitRankArray enumerate_microtrees(const int64_t B,
                                             const uint64_t microtree_count,
                                             const BitArray &microtree_roots,
                                             const TreeBP &tree) {
    EditableMicrotreeArray microtrees(B, microtree_count);
    BitArray split_rank_count_locked(microtree_count);
    std::vector<uint32_t> split_ranks(microtree_count),
        write_pos(microtree_count);

    uint64_t pre = 0, index = 0;
    auto match = [&](uint64_t v) -> void {
        assert(tree.bp.get(index) == v);
        index++;
    };
    auto has_child = [&]() -> bool {
        return index < 2 * tree.n && tree.bp.get(index) == 0;
    };

    uint64_t up_pre = 0;
    auto dfs = [&](auto self, int cur_up_pre = -1) -> void {
        bool is_microtree_root = microtree_roots.get(pre);
        pre++;

        if (is_microtree_root) {
            if (cur_up_pre != -1 && !split_rank_count_locked.get(cur_up_pre)) {
                split_rank_count_locked.on(cur_up_pre);
            }
            cur_up_pre = up_pre;
        }
        up_pre += is_microtree_root;

        match(0);
        microtrees.set(cur_up_pre, write_pos[cur_up_pre]++, false);
        bool has_l = has_child();
        if (has_l) {
            self(self, cur_up_pre);
        }

        match(1);
        microtrees.set(cur_up_pre, write_pos[cur_up_pre]++, true);
        split_ranks[cur_up_pre] += !split_rank_count_locked.get(cur_up_pre);
        bool has_r = has_child();
        if (has_r) {
            self(self, cur_up_pre);
        }

        return;
    };

    dfs(dfs);

    return {microtrees, split_ranks};
}

std::pair<TreeBP, MicrotreeSplitRankArray> recover_upsilon_and_microtrees(
    const int64_t B, const BitArray &microtree_roots, const TreeBP &tree) {
    uint64_t microtree_count = microtree_roots.linear_popcount();

    EditableMicrotreeArray microtrees(B, microtree_count);
    BitArray upsilon_bp(2 * microtree_count);
    BitArray split_rank_count_locked(microtree_count);
    std::vector<uint32_t> split_ranks(microtree_count),
        write_pos(microtree_count);

    uint64_t pre = 0;
    uint64_t up_pre = 0, up_index = 0;

    auto stack_bit = [&](uint64_t up_v, bool b) {
        assert(0 <= up_v && up_v < microtree_count);
        if (b) {
            microtrees.set(up_v, write_pos[up_v], b);
            if (!split_rank_count_locked.get(up_v)) {
                split_ranks[up_v]++;
            }
        }
        write_pos[up_v]++;
    };

    uint64_t excess = 0;
    std::stack<uint64_t> up_pre_stack, excess_stack;
    for (int64_t index = 0; index < tree.n * 2; index++) {
        if (tree.bp.get(index) == 0) {
            if (microtree_roots.get(pre)) {
                if (!up_pre_stack.empty()) {
                    split_rank_count_locked.on(up_pre_stack.top());
                }
                upsilon_bp.off(up_index++);
                up_pre_stack.push(up_pre);
                excess_stack.push(excess);
                up_pre++;
            }
            pre++;

            assert(!up_pre_stack.empty());
            stack_bit(up_pre_stack.top(), 0);
            excess++;
        } else {
            assert(!up_pre_stack.empty());
            stack_bit(up_pre_stack.top(), 1);
            excess--;

            if (!excess_stack.empty() && excess == excess_stack.top()) {
                if (index == 2 * tree.n - 1 || tree.bp.get(index + 1) == 1 ||
                    microtree_roots.get(pre)) {
                    upsilon_bp.on(up_index++);
                    up_pre_stack.pop();
                    excess_stack.pop();
                }
            }
        }
    }

    TreeBP upsilon = TreeBP(microtree_count, std::move(upsilon_bp));
    MicrotreeSplitRankArray microtree_split_rank_array =
        MicrotreeSplitRankArray(microtrees, split_ranks);
    return {upsilon, microtree_split_rank_array};
}

std::pair<TreeBP, std::vector<uint64_t>> legacy_tree_covering(
    const int64_t B, const TreeBP &tree) {
    assert(B <= 15);
    BitArray microtree_roots = microtree_roots_preorder(B, tree);
    uint64_t microtree_count = microtree_roots.linear_popcount();
    TreeBP upsilon =
        legacy_build_upsilon(microtree_count, microtree_roots, tree);
    std::vector<uint64_t> microtrees =
        legacy_enumerate_microtrees(microtree_count, microtree_roots, tree);
    return {upsilon, microtrees};
}

std::pair<TreeBP, MicrotreeSplitRankArray> tree_covering(const int64_t B,
                                                         const TreeBP &tree) {
    return recover_upsilon_and_microtrees(B, microtree_roots_preorder(B, tree),
                                          tree);
}

}  // namespace average_case_optimal_rmq
