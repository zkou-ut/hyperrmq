#pragma once

#include "bit_array.hpp"
#include "huffman.hpp"
#include "rmm_tree.hpp"
#include "two_level_increasing_array.hpp"

namespace average_case_optimal_rmq {

template <uint32_t W = 64>
struct LegacyRMQHuffman {
    uint32_t num_of_nodes, num_of_microtrees, num_of_chunks;
    RMMTree<8, 1024> rmm_tree;
    CanonicalHuffmanCode<uint64_t> chc;
    BitArray code_seq;
    TwoLevelIncreasingArray<> seq_idx_sample, close_sample;

    LegacyRMQHuffman();
    explicit LegacyRMQHuffman(const std::vector<int32_t> &values,
                              const int B = 7);

    uint32_t to_inorder(uint32_t i);
    uint32_t find_seq_idx(uint32_t i);
    uint64_t get_microtree_inorder(uint32_t i);
    std::pair<uint32_t, uint32_t> get_chunk(uint32_t i);
    std::pair<uint32_t, uint32_t> select1(uint32_t c);
    uint32_t rank1(uint32_t c, uint32_t k);

    uint32_t query(uint32_t i, uint32_t j);

    uint64_t evaluate_memory_consumption() const;
};

}  // namespace average_case_optimal_rmq
