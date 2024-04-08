#pragma once

#include <cassert>
#include <stack>

#include "hyperrmq/bit_array.hpp"
#include "hyperrmq/bitutil.hpp"
#include "hyperrmq/huffman.hpp"
#include "hyperrmq/memutil.hpp"
#include "hyperrmq/microtree_array.hpp"
#include "hyperrmq/rmm_tree.hpp"
#include "hyperrmq/three_level_prefix_sum.hpp"
#include "hyperrmq/tree_covering.hpp"

namespace average_case_optimal_rmq {

template <uint32_t W = 64, typename CompressedMicrotreeSplitRankArray =
                               CompressedMicrotreeSplitRankArrayHuffman<W>>
struct HyperRMQ {
    static_assert(std::is_base_of_v<CompressedMicrotreeSplitRankArrayInterface,
                                    CompressedMicrotreeSplitRankArray>);

    uint32_t num_of_nodes, num_of_microtrees, num_of_chunks;
    RMMTree<8, 1024> rmm_tree;
    CompressedMicrotreeSplitRankArray compressed_microtree_split_rank_array;
    ThreeLevelPrefixSum<> close_sample;

    HyperRMQ();
    explicit HyperRMQ(const std::vector<int32_t> &values, const int B = 7);

    uint32_t chunk_index_to_microtree_preorder(uint32_t i) const;
    BitArray get_chunk(uint32_t i) const;
    uint64_t get_chunk_popcount(uint32_t i) const;

    // Inorder index to a pair of indices {chunk index, close index}.
    std::pair<uint32_t, uint32_t> inorderselect(uint32_t c) const;

    // A pair of indices {chunk index, close index} to inorder index.
    uint32_t inorder(std::pair<uint32_t, uint32_t> tau) const;

    uint32_t query(uint32_t i, uint32_t j) const;

    uint64_t evaluate_memory_consumption() const;
    std::vector<std::pair<std::string, uint64_t>> memory_table() const;
};

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline HyperRMQ<W, CompressedMicrotreeSplitRankArray>::HyperRMQ() {}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline HyperRMQ<W, CompressedMicrotreeSplitRankArray>::HyperRMQ(
    const std::vector<int32_t> &values, const int B) {
    num_of_nodes = values.size();

    auto [upsilon, raw_microtree_split_rank_array] =
        tree_covering(B, cartesian_tree_bp(values));

    rmm_tree = decltype(rmm_tree)(upsilon);

    compressed_microtree_split_rank_array =
        decltype(compressed_microtree_split_rank_array)(
            raw_microtree_split_rank_array);

    num_of_microtrees = compressed_microtree_split_rank_array.size();
    num_of_chunks = num_of_microtrees * 2;

    std::vector<uint32_t> close_sample_vec;
    close_sample_vec.reserve(num_of_chunks / W);
    uint32_t close = 0;
    int rank0 = 0;
    std::stack<int> opens;
    for (int i = 0; i < num_of_chunks; i++) {
        if (rmm_tree.get_bit(i) == 0) {
            const auto &[tree, split_rank] =
                raw_microtree_split_rank_array[rank0];
            opens.push(rank0);
            rank0++;
            close += split_rank;
        } else {
            const auto &[tree, split_rank] =
                raw_microtree_split_rank_array[opens.top()];
            opens.pop();
            close += tree.n - split_rank;
        }
        if ((i + 1) % W == 0) {
            close_sample_vec.push_back(close);
            close = 0;
        }
    }
    close_sample = decltype(close_sample)(close_sample_vec);
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline uint32_t HyperRMQ<W, CompressedMicrotreeSplitRankArray>::
    chunk_index_to_microtree_preorder(uint32_t i) const {
    uint32_t open_index = (rmm_tree.get_bit(i) == 1) ? rmm_tree.open(i) : i;
    uint32_t open_rank = rmm_tree.rank0(open_index);
    return open_rank;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
BitArray HyperRMQ<W, CompressedMicrotreeSplitRankArray>::get_chunk(
    uint32_t i) const {
    assert(0 <= i && i < num_of_chunks);
    bool isopen = (rmm_tree.get_bit(i) == 0);
    i = chunk_index_to_microtree_preorder(i);
    if (isopen) {
        return compressed_microtree_split_rank_array.get_left_chunk(i);
    } else {
        return compressed_microtree_split_rank_array.get_right_chunk(i);
    }
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline uint64_t
HyperRMQ<W, CompressedMicrotreeSplitRankArray>::get_chunk_popcount(
    uint32_t i) const {
    assert(0 <= i && i < num_of_chunks);
    bool isopen = (rmm_tree.get_bit(i) == 0);
    i = chunk_index_to_microtree_preorder(i);
    if (isopen) {
        return compressed_microtree_split_rank_array.get_left_chunk_popcount(i);
    } else {
        return compressed_microtree_split_rank_array.get_right_chunk_popcount(
            i);
    }
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline std::pair<uint32_t, uint32_t>
HyperRMQ<W, CompressedMicrotreeSplitRankArray>::inorderselect(
    uint32_t c) const {
    assert(0 <= c && c < num_of_nodes);

    auto sample_index = close_sample.select_chunk(c);
    int32_t remain = c - close_sample.sum(sample_index);
    uint32_t chunk_index = sample_index * W;
    while (true) {
        int32_t p = get_chunk_popcount(chunk_index);
        if (remain - p < 0) {
            break;
        }
        remain -= p;
        chunk_index++;
    }

    return {chunk_index, remain};
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline uint32_t HyperRMQ<W, CompressedMicrotreeSplitRankArray>::inorder(
    std::pair<uint32_t, uint32_t> tau) const {
    const auto &[c, k] = tau;
    assert(0 <= c && c < num_of_chunks);

    uint32_t res = close_sample.sum(c / W);
    for (int i = c / W * W; i < c; i++) {
        res += get_chunk_popcount(i);
    }
    res += k;
    return res;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline uint32_t HyperRMQ<W, CompressedMicrotreeSplitRankArray>::query(
    uint32_t i, uint32_t j) const {
    assert(0 <= i && i <= j && j < num_of_nodes);

    auto [i_chunk, i_close] = inorderselect(i);
    auto [j_chunk, j_close] = inorderselect(j);

    uint32_t i_micro = chunk_index_to_microtree_preorder(i_chunk);
    uint32_t j_micro = chunk_index_to_microtree_preorder(j_chunk);

    if (rmm_tree.get_bit(i_chunk) == 1) {
        i_close +=
            compressed_microtree_split_rank_array.get_left_chunk_popcount(
                i_micro);
    }

    if (rmm_tree.get_bit(j_chunk) == 1) {
        j_close +=
            compressed_microtree_split_rank_array.get_left_chunk_popcount(
                j_micro);
    }

    // Recovering lca_chunk is skipped here.
    // Since i_close <= lca_close <= j_close holds, finding lca_chunk can be
    // achieved by choosing appropriate one from i_chunk and j_chunk.
    auto to_chunk_pair = [&](uint32_t lca_close, uint32_t lca_split_rank,
                             uint32_t i_chunk, uint32_t j_chunk) -> uint32_t {
        if (lca_close < lca_split_rank) {
            return inorder({i_chunk, lca_close});
        } else {
            return inorder({j_chunk, lca_close - lca_split_rank});
        }
    };

    if (i_micro == j_micro) {
        auto lca_close = compressed_microtree_split_rank_array.get_lca(
            i_micro, i_close, j_close);
        auto split_rank =
            compressed_microtree_split_rank_array.get_left_chunk_popcount(
                i_micro);
        return to_chunk_pair(lca_close, split_rank, i_chunk, j_chunk);
    }

    // Note that it works even if rmm_tree.get_bit(ic) = 0 or
    // rmm_tree.get_bit(jc) = 0 due to the property of BP.
    uint32_t lca_chunk = rmm_tree.rmq(i_chunk, j_chunk + 1);
    uint32_t lca_micro = chunk_index_to_microtree_preorder(lca_chunk);

    if (lca_micro == i_micro) {
        auto split_rank =
            compressed_microtree_split_rank_array.get_left_chunk_popcount(
                lca_micro);
        j_chunk = lca_chunk;
        j_close = split_rank;
        if (i_close < j_close) {
            j_close--;
        }
        auto lca_close = compressed_microtree_split_rank_array.get_lca(
            lca_micro, i_close, j_close);
        return to_chunk_pair(lca_close, split_rank, i_chunk, j_chunk);
    } else if (lca_micro == j_micro) {
        auto split_rank =
            compressed_microtree_split_rank_array.get_left_chunk_popcount(
                lca_micro);
        i_chunk = lca_chunk;
        i_close = split_rank;
        if (j_close < i_close) {
            i_close--;
        }
        auto lca_close = compressed_microtree_split_rank_array.get_lca(
            lca_micro, i_close, j_close);
        return to_chunk_pair(lca_close, split_rank, i_chunk, j_chunk);
    } else {
        return inorder({lca_chunk, 0});
    }
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline uint64_t HyperRMQ<
    W, CompressedMicrotreeSplitRankArray>::evaluate_memory_consumption() const {
    return rmm_tree.evaluate_memory_consumption() +
           compressed_microtree_split_rank_array.evaluate_memory_consumption() +
           close_sample.evaluate_memory_consumption();
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline std::vector<std::pair<std::string, uint64_t>>
HyperRMQ<W, CompressedMicrotreeSplitRankArray>::memory_table() const {
    return memory_table_combine_children(
        "HyperRMQ", {
                        rmm_tree.memory_table(),
                        compressed_microtree_split_rank_array.memory_table(),
                        close_sample.memory_table(),
                    });
}

}  // namespace average_case_optimal_rmq
