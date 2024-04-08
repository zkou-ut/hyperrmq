#pragma once

#include <map>
#include <set>

#include "binary_tree.hpp"
#include "bit_array.hpp"
#include "microtree_array.hpp"
#include "minimal_cell_array.hpp"
#include "tree_bp.hpp"

namespace average_case_optimal_rmq {

uint64_t legacy_encode_microtree(uint64_t bp, uint64_t cut_pos);
std::pair<uint64_t, uint64_t> legacy_decode_microtree(uint64_t microtree);
uint32_t legacy_lca_on_microtree(uint64_t microtree, uint32_t u, uint32_t v);

// Farzan--Munro algorithm
std::set<int32_t> legacy_microtree_roots_inorder(int32_t B,
                                                 const BinaryTree &binary_tree);
BitArray microtree_roots_preorder(int64_t B, const TreeBP &tree);

TreeBP legacy_build_upsilon(const uint64_t microtree_count,
                            const BitArray &microtree_roots,
                            const TreeBP &tree);

TreeBP build_upsilon(const uint64_t microtree_count,
                     const BitArray &microtree_roots, const TreeBP &tree,
                     const MicrotreeSplitRankArray &microtree_split_rank_array);

std::vector<uint64_t> legacy_enumerate_microtrees(
    const uint64_t microtree_count, const BitArray &microtree_roots,
    const TreeBP &tree);

MicrotreeSplitRankArray enumerate_microtrees(const int64_t B,
                                          const uint64_t microtree_count,
                                          const BitArray &microtree_roots,
                                          const TreeBP &tree);

std::pair<TreeBP, MicrotreeSplitRankArray> recover_upsilon_and_microtrees(
    const int64_t B, const BitArray &microtree_roots, const TreeBP &tree);

std::pair<TreeBP, std::vector<uint64_t>> legacy_tree_covering(
    const int64_t B, const TreeBP &tree);

std::pair<TreeBP, MicrotreeSplitRankArray> tree_covering(const int64_t B,
                                                      const TreeBP &tree);

}  // namespace average_case_optimal_rmq
