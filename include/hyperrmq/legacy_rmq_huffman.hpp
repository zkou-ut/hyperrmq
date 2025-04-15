#pragma once

#include "hyperrmq/bit_array.hpp"
#include "hyperrmq/huffman.hpp"
#include "hyperrmq/rmm_tree.hpp"
#include "hyperrmq/two_level_increasing_array.hpp"

namespace hyperrmq {

template <uint64_t W = 64>
struct LegacyRMQHuffman {
    uint64_t num_of_nodes, num_of_microtrees, num_of_chunks;
    RMMTree<8, 1024> rmm_tree;
    CanonicalHuffmanCode<uint64_t> chc;
    BitArray code_seq;
    TwoLevelIncreasingArray<> seq_idx_sample, close_sample;

    LegacyRMQHuffman();
    explicit LegacyRMQHuffman(const std::vector<int64_t> &values,
                              const int64_t B = 7);

    uint64_t to_inorder(uint64_t i);
    uint64_t find_seq_idx(uint64_t i);
    uint64_t get_microtree_inorder(uint64_t i);
    std::pair<uint64_t, uint64_t> get_chunk(uint64_t i);
    std::pair<uint64_t, uint64_t> select1(uint64_t c);
    uint64_t rank1(uint64_t c, uint64_t k);

    uint64_t query(uint64_t i, uint64_t j);

    uint64_t evaluate_memory_consumption() const;
};

}  // namespace hyperrmq
