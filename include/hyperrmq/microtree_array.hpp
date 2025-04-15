#pragma once

#include "hyperrmq/array_base.hpp"
#include "hyperrmq/bit_array.hpp"
#include "hyperrmq/huffman.hpp"
#include "hyperrmq/three_level_prefix_sum.hpp"
#include "hyperrmq/tree_bp.hpp"

namespace hyperrmq {

struct EditableMicrotreeArray : FixedLengthCodeArrayBase<TreeBP> {
   public:
    EditableMicrotreeArray();
    EditableMicrotreeArray(uint64_t B, uint64_t microtree_count);
    EditableMicrotreeArray(const std::vector<TreeBP>& tree_vec);

    BitArray encode(const TreeBP& object) const;
    TreeBP decode(const BitArray& code) const;

    using FixedLengthCodeArrayBase::set;
    void set(uint64_t tree_index, uint64_t bit_index, bool bit_value);

    using FixedLengthCodeArrayBase::get;
    bool get(uint64_t tree_index, uint64_t bit_index) const;
};

struct MicrotreeSplitRankArray
    : FixedLengthCodeArrayBase<std::pair<TreeBP, uint64_t>> {
   public:
    MicrotreeSplitRankArray();
    MicrotreeSplitRankArray(
        const std::vector<std::pair<TreeBP, uint64_t>>& vec);
    MicrotreeSplitRankArray(const EditableMicrotreeArray& microtrees,
                            const std::vector<uint64_t>& split_ranks);

    BitArray encode(const std::pair<TreeBP, uint64_t>& object) const;
    std::pair<TreeBP, uint64_t> decode(const BitArray& code) const;

   private:
    uint64_t tree_width, split_rank_width;
};

struct CompressedMicrotreeSplitRankArrayInterface {
    virtual size_t size() const = 0;
    virtual std::pair<TreeBP, uint64_t> operator[](uint64_t index) const = 0;
    virtual uint64_t evaluate_memory_consumption() const = 0;
    virtual std::vector<std::pair<std::string, uint64_t>> memory_table()
        const = 0;

    uint64_t get_node_count(uint64_t index) const;

    BitArray get_left_chunk(uint64_t index) const;
    BitArray get_right_chunk(uint64_t index) const;

    uint64_t get_left_chunk_popcount(uint64_t index) const;
    uint64_t get_right_chunk_popcount(uint64_t index) const;

    std::pair<uint64_t, uint64_t> get_node_count_left_popcount(
        uint64_t index) const;

    uint64_t get_lca(uint64_t index, uint64_t u_inorder,
                     uint64_t v_inorder) const;
};

template <uint64_t W = 64>
struct CompressedMicrotreeSplitRankArrayHuffmanNaive
    : CompressedMicrotreeSplitRankArrayInterface {
   public:
    CompressedMicrotreeSplitRankArrayHuffmanNaive();
    CompressedMicrotreeSplitRankArrayHuffmanNaive(
        const MicrotreeSplitRankArray& microtree_split_rank_array);

    size_t size() const override;
    std::pair<TreeBP, uint64_t> operator[](uint64_t index) const override;
    uint64_t evaluate_memory_consumption() const override;
    std::vector<std::pair<std::string, uint64_t>> memory_table() const override;

    uint64_t length;
    CanonicalHuffmanCode<TreeBP, EditableMicrotreeArray> chc;
    BitArray code_seq;
    ThreeLevelPrefixSum<> code_idx_sample;
    MinimalCellArray<uint64_t> split_ranks;
};

template <uint64_t W = 64>
struct CompressedMicrotreeSplitRankArrayHuffman
    : CompressedMicrotreeSplitRankArrayInterface {
   public:
    CompressedMicrotreeSplitRankArrayHuffman();
    CompressedMicrotreeSplitRankArrayHuffman(
        const MicrotreeSplitRankArray& microtree_split_rank_array);

    size_t size() const override;
    std::pair<TreeBP, uint64_t> operator[](uint64_t index) const override;
    uint64_t evaluate_memory_consumption() const override;
    std::vector<std::pair<std::string, uint64_t>> memory_table() const override;

    uint64_t length;
    CanonicalHuffmanCode<std::pair<TreeBP, uint64_t>, MicrotreeSplitRankArray>
        chc;
    BitArray code_seq;
    ThreeLevelPrefixSum<> code_idx_sample;
};

template <bool depth_first, bool pruning = true>
struct CompressedMicrotreeSplitRankArrayArithmetic
    : CompressedMicrotreeSplitRankArrayInterface {
   public:
    CompressedMicrotreeSplitRankArrayArithmetic();
    CompressedMicrotreeSplitRankArrayArithmetic(
        const MicrotreeSplitRankArray& microtree_split_rank_array);

    size_t size() const override;
    std::pair<TreeBP, uint64_t> operator[](uint64_t index) const override;
    uint64_t evaluate_memory_consumption() const override;
    std::vector<std::pair<std::string, uint64_t>> memory_table() const override;

    uint64_t get_node_count(uint64_t index) const;

    uint64_t get_left_chunk_popcount(uint64_t index) const;
    uint64_t get_right_chunk_popcount(uint64_t index) const;

    std::pair<uint64_t, uint64_t> get_node_count_left_popcount(
        uint64_t index) const;

    uint64_t get_lca(uint64_t index, uint64_t u_inorder,
                     uint64_t v_inorder) const;

    uint64_t length;
    BitArray code_seq;
    MinimalCellArray<uint64_t> node_counts, split_ranks;
    ThreeLevelPrefixSum<uint64_t> code_len;
};

// node_count, split_rank, left_seq
struct CompressedMicrotreeSplitRankArrayAllArithmetic
    : CompressedMicrotreeSplitRankArrayInterface {
   public:
    CompressedMicrotreeSplitRankArrayAllArithmetic();
    CompressedMicrotreeSplitRankArrayAllArithmetic(
        const MicrotreeSplitRankArray& microtree_split_rank_array);

    size_t size() const override;
    std::pair<TreeBP, uint64_t> operator[](uint64_t index) const override;
    uint64_t evaluate_memory_consumption() const override;
    std::vector<std::pair<std::string, uint64_t>> memory_table() const override;

    uint64_t get_node_count(uint64_t index) const;

    uint64_t get_left_chunk_popcount(uint64_t index) const;
    uint64_t get_right_chunk_popcount(uint64_t index) const;

    std::pair<uint64_t, uint64_t> get_node_count_left_popcount(
        uint64_t index) const;

    uint64_t length;
    uint64_t max_tree_n;
    BitArray code_seq;
    ThreeLevelPrefixSum<uint64_t> code_len;
};

}  // namespace hyperrmq
