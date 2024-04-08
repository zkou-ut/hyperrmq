#include "microtree_array.hpp"

#include <algorithm>
#include <cassert>
#include <map>
#include <queue>

#include "arithmetic.hpp"
#include "bitutil.hpp"
#include "memutil.hpp"

namespace average_case_optimal_rmq {

EditableMicrotreeArray::EditableMicrotreeArray() {}

EditableMicrotreeArray::EditableMicrotreeArray(uint32_t B,
                                               uint32_t microtree_count)
    : FixedLengthCodeArrayBase(2 * (2 * B - 1), microtree_count) {}

EditableMicrotreeArray::EditableMicrotreeArray(
    const std::vector<TreeBP> &tree_vec) {
    length = tree_vec.size();
    width = 0;
    for (auto &&tree : tree_vec) {
        if (width < tree.bp.size()) {
            width = tree.bp.size();
        }
    }
    codes = BitArray(width * length);
    for (int tree_index = 0; tree_index < length; tree_index++) {
        codes.write_interval(tree_index * width, tree_vec[tree_index].bp);
    }
}

BitArray EditableMicrotreeArray::encode(const TreeBP &object) const {
    BitArray ret(width);
    ret.write_interval(0, object.bp);
    return ret;
}

TreeBP EditableMicrotreeArray::decode(const BitArray &code) const {
    uint64_t node_count = code.linear_popcount();
    BitArray bp = code.read_interval(0, 2 * node_count);
    return TreeBP(node_count, bp);
}

void EditableMicrotreeArray::set(uint64_t tree_index, uint64_t bit_index,
                                 bool bit_value) {
    assert(0 <= tree_index && tree_index < length);
    assert(0 <= bit_index && bit_index < width);
    codes.set(tree_index * width + bit_index, bit_value);
}

bool EditableMicrotreeArray::get(uint64_t tree_index,
                                 uint64_t bit_index) const {
    assert(0 <= tree_index && tree_index < length);
    assert(0 <= bit_index && bit_index < width);
    return codes.get(tree_index * width + bit_index);
}

MicrotreeSplitRankArray::MicrotreeSplitRankArray() {}

MicrotreeSplitRankArray::MicrotreeSplitRankArray(
    const std::vector<std::pair<TreeBP, uint32_t>> &vec) {
    length = vec.size();
    tree_width = 0;
    uint32_t split_rank_max = 0;
    for (auto &&[tree, split_rank] : vec) {
        if (tree_width < tree.bp.size()) {
            tree_width = tree.bp.size();
        }
        if (split_rank_max < split_rank) {
            split_rank_max = split_rank;
        }
    }
    split_rank_width = ceil_log2(split_rank_max + 1);
    width = tree_width + split_rank_width;
    codes = BitArray(width * length);
    for (int index = 0; index < length; index++) {
        set(index, vec[index]);
    }
}

MicrotreeSplitRankArray::MicrotreeSplitRankArray(
    const EditableMicrotreeArray &microtrees,
    const std::vector<uint32_t> &split_ranks) {
    assert(microtrees.size() == split_ranks.size());
    length = microtrees.size();
    tree_width = microtrees.width;
    split_rank_width = ceil_log2(
        (*std::max_element(split_ranks.begin(), split_ranks.end())) + 1);
    width = tree_width + split_rank_width;
    codes = BitArray(width * length);
    for (int index = 0; index < length; index++) {
        set(index, {microtrees[index], split_ranks[index]});
    }
}

BitArray MicrotreeSplitRankArray::encode(
    const std::pair<TreeBP, uint32_t> &object) const {
    auto &[tree, split_rank] = object;
    BitArray ret(width);
    ret.write_interval(0, tree.bp);
    ret.write_bits(tree_width, split_rank_width, split_rank);
    return ret;
}

std::pair<TreeBP, uint32_t> MicrotreeSplitRankArray::decode(
    const BitArray &code) const {
    uint64_t node_count = code.linear_popcount(0, tree_width);
    BitArray bp = code.read_interval(0, 2 * node_count);
    uint32_t split_rank = code.read_bits(tree_width, split_rank_width);
    return {TreeBP(node_count, std::move(bp)), split_rank};
}

uint32_t CompressedMicrotreeSplitRankArrayInterface::get_node_count(
    uint64_t index) const {
    return (*this)[index].first.n;
}

BitArray CompressedMicrotreeSplitRankArrayInterface::get_left_chunk(
    uint64_t index) const {
    const auto &[tree, split_rank] = (*this)[index];
    auto cut = tree.bp.linear_select1(split_rank);
    return tree.bp.read_interval(0, cut);
}

BitArray CompressedMicrotreeSplitRankArrayInterface::get_right_chunk(
    uint64_t index) const {
    const auto &[tree, split_rank] = (*this)[index];
    auto cut = tree.bp.linear_select1(split_rank);
    return tree.bp.read_interval(cut, tree.bp.size() - cut);
}

uint32_t CompressedMicrotreeSplitRankArrayInterface::get_left_chunk_popcount(
    uint64_t index) const {
    const auto &[tree, split_rank] = (*this)[index];
    return split_rank;
}

uint32_t CompressedMicrotreeSplitRankArrayInterface::get_right_chunk_popcount(
    uint64_t index) const {
    const auto &[tree, split_rank] = (*this)[index];
    return tree.n - split_rank;
}

std::pair<uint32_t, uint32_t>
CompressedMicrotreeSplitRankArrayInterface::get_node_count_left_popcount(
    uint64_t index) const {
    const auto &[tree, split_rank] = (*this)[index];
    return {tree.n, split_rank};
}

uint32_t CompressedMicrotreeSplitRankArrayInterface::get_lca(
    uint64_t index, uint32_t u_inorder, uint32_t v_inorder) const {
    if (u_inorder == v_inorder) {
        return u_inorder;
    }
    const auto &[tree, split_rank] = (*this)[index];
    return tree.naive_lca(u_inorder, v_inorder);
}

template <uint32_t W>
inline CompressedMicrotreeSplitRankArrayHuffmanNaive<
    W>::CompressedMicrotreeSplitRankArrayHuffmanNaive() {}

template <uint32_t W>
CompressedMicrotreeSplitRankArrayHuffmanNaive<W>::
    CompressedMicrotreeSplitRankArrayHuffmanNaive(
        const MicrotreeSplitRankArray &microtree_split_rank_array)
    : length(microtree_split_rank_array.size()) {
    std::map<TreeBP, uint64_t> counter;
    std::vector<uint32_t> split_ranks_vec(length);
    for (int i = 0; i < length; i++) {
        auto [tree, split_rank] = microtree_split_rank_array[i];
        counter[tree]++;
        split_ranks_vec[i] = split_rank;
    }
    chc = decltype(chc)(counter);
    split_ranks = split_ranks_vec;

    auto encode = chc.enumerate_alphabet_code_pair();
    uint32_t code_length_sum = 0;
    for (auto &&[tree, cnt] : counter) {
        code_length_sum += encode[tree].first * cnt;
    }

    std::vector<uint32_t> seq_idx_sample_vec;
    seq_idx_sample_vec.reserve(length / W);
    code_seq = BitArray(code_length_sum);
    uint32_t code_idx = 0, last_code_idx = 0;
    for (int i = 0; i < length; i++) {
        auto [len, code] = encode[microtree_split_rank_array[i].first];
        code_seq.write_bits(code_idx, len, code);
        code_idx += len;
        if ((i + 1) % W == 0) {
            seq_idx_sample_vec.push_back(code_idx - last_code_idx);
            last_code_idx = code_idx;
        }
    }
    assert(code_idx == code_length_sum);
    code_idx_sample = seq_idx_sample_vec;
}

template <uint32_t W>
size_t CompressedMicrotreeSplitRankArrayHuffmanNaive<W>::size() const {
    return length;
}

template <uint32_t W>
std::pair<TreeBP, uint32_t>
CompressedMicrotreeSplitRankArrayHuffmanNaive<W>::operator[](
    uint64_t index) const {
    uint32_t code_idx = code_idx_sample.sum(index / W);
    for (int j = 0; j < index % W; j++) {
        code_idx += chc.get_next_length(code_seq, code_idx);
    }
    uint32_t len = chc.get_next_length(code_seq, code_idx);
    return {chc.decode(len, code_seq.read_bits(code_idx, len)),
            split_ranks[index]};
}

template <uint32_t W>
uint64_t CompressedMicrotreeSplitRankArrayHuffmanNaive<
    W>::evaluate_memory_consumption() const {
    return chc.evaluate_memory_consumption() +
           code_seq.evaluate_memory_consumption() +
           code_idx_sample.evaluate_memory_consumption() +
           split_ranks.evaluate_memory_consumption();
}

template <uint32_t W>
std::vector<std::pair<std::string, uint64_t>>
CompressedMicrotreeSplitRankArrayHuffmanNaive<W>::memory_table() const {
    return memory_table_combine_children(
        "HuffmanArray",
        {
            {{"Sequence", code_seq.evaluate_memory_consumption()}},
            {{"Sampling", code_idx_sample.evaluate_memory_consumption()}},
            {{"Table", chc.evaluate_memory_consumption()}},
            {{"Split Ranks", split_ranks.evaluate_memory_consumption()}},
        });
}

template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<1>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<2>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<4>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<8>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<16>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<32>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<64>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<128>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<256>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<512>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<1024>;
template struct CompressedMicrotreeSplitRankArrayHuffmanNaive<2048>;

template <uint32_t W>
CompressedMicrotreeSplitRankArrayHuffman<
    W>::CompressedMicrotreeSplitRankArrayHuffman() {}

template <uint32_t W>
CompressedMicrotreeSplitRankArrayHuffman<
    W>::CompressedMicrotreeSplitRankArrayHuffman(const MicrotreeSplitRankArray &
                                                     microtree_split_rank_array)
    : length(microtree_split_rank_array.size()) {
    std::map<std::pair<TreeBP, uint32_t>, uint64_t> counter;
    for (int i = 0; i < length; i++) {
        counter[microtree_split_rank_array[i]]++;
    }
    chc = decltype(chc)(counter);

    auto encode = chc.enumerate_alphabet_code_pair();
    uint32_t code_length_sum = 0;
    for (auto &&[tree_and_split_rank, cnt] : counter) {
        code_length_sum += encode[tree_and_split_rank].first * cnt;
    }

    std::vector<uint32_t> seq_idx_sample_vec;
    seq_idx_sample_vec.reserve(length / W);
    code_seq = BitArray(code_length_sum);
    uint32_t code_idx = 0, last_code_idx = 0;
    for (int i = 0; i < length; i++) {
        auto [len, code] = encode[microtree_split_rank_array[i]];
        code_seq.write_bits(code_idx, len, code);
        code_idx += len;
        if ((i + 1) % W == 0) {
            seq_idx_sample_vec.push_back(code_idx - last_code_idx);
            last_code_idx = code_idx;
        }
    }
    assert(code_idx == code_length_sum);
    code_idx_sample = seq_idx_sample_vec;
}

template <uint32_t W>
std::pair<TreeBP, uint32_t>
CompressedMicrotreeSplitRankArrayHuffman<W>::operator[](uint64_t index) const {
    uint32_t code_idx = code_idx_sample.sum(index / W);
    for (int j = 0; j < index % W; j++) {
        code_idx += chc.get_next_length(code_seq, code_idx);
    }
    uint32_t len = chc.get_next_length(code_seq, code_idx);
    return chc.decode(len, code_seq.read_bits(code_idx, len));
}

template <uint32_t W>
uint64_t CompressedMicrotreeSplitRankArrayHuffman<
    W>::evaluate_memory_consumption() const {
    return chc.evaluate_memory_consumption() +
           code_seq.evaluate_memory_consumption() +
           code_idx_sample.evaluate_memory_consumption();
}

template <uint32_t W>
std::vector<std::pair<std::string, uint64_t>>
CompressedMicrotreeSplitRankArrayHuffman<W>::memory_table() const {
    return memory_table_combine_children(
        "HuffmanArray",
        {
            {{"Sequence", code_seq.evaluate_memory_consumption()}},
            {{"Sampling", code_idx_sample.evaluate_memory_consumption()}},
            {{"Table", chc.evaluate_memory_consumption()}},
        });
}

template <uint32_t W>
size_t CompressedMicrotreeSplitRankArrayHuffman<W>::size() const {
    return length;
}

template struct CompressedMicrotreeSplitRankArrayHuffman<1>;
template struct CompressedMicrotreeSplitRankArrayHuffman<2>;
template struct CompressedMicrotreeSplitRankArrayHuffman<4>;
template struct CompressedMicrotreeSplitRankArrayHuffman<8>;
template struct CompressedMicrotreeSplitRankArrayHuffman<16>;
template struct CompressedMicrotreeSplitRankArrayHuffman<32>;
template struct CompressedMicrotreeSplitRankArrayHuffman<64>;
template struct CompressedMicrotreeSplitRankArrayHuffman<128>;
template struct CompressedMicrotreeSplitRankArrayHuffman<256>;
template struct CompressedMicrotreeSplitRankArrayHuffman<512>;
template struct CompressedMicrotreeSplitRankArrayHuffman<1024>;
template struct CompressedMicrotreeSplitRankArrayHuffman<2048>;

template <bool depth_first>
CompressedMicrotreeSplitRankArrayArithmetic<
    depth_first>::CompressedMicrotreeSplitRankArrayArithmetic() {}

template <bool depth_first>
CompressedMicrotreeSplitRankArrayArithmetic<depth_first>::
    CompressedMicrotreeSplitRankArrayArithmetic(
        const MicrotreeSplitRankArray &microtree_split_rank_array)
    : length(microtree_split_rank_array.size()) {
    std::vector<uint32_t> node_counts_vec(length), split_ranks_vec(length);
    std::vector<uint32_t> code_len_vec(length);

    std::map<TreeBP, BitArray> arith_encode;
    for (int i = 0; i < length; i++) {
        const auto &[tree, split_rank] = microtree_split_rank_array[i];
        if (!arith_encode.count(tree)) {
            auto code = left_seq_to_arithmetic<depth_first>(
                bp_to_left_seq<depth_first>(tree));
            arith_encode[tree] = code;
        }

        node_counts_vec[i] = tree.n;
        split_ranks_vec[i] = split_rank;
        code_len_vec[i] = arith_encode[tree].size();
    }

    node_counts = node_counts_vec;
    split_ranks = split_ranks_vec;
    code_len = code_len_vec;

    code_seq = BitArray(code_len.sum(length));
    uint32_t first = 0;
    for (int i = 0; i < length; i++) {
        auto [tree, split_rank] = microtree_split_rank_array[i];
        const auto &code = arith_encode[tree];
        code_seq.write_interval(first, code);
        first += code.size();
    }
    assert(first == code_seq.size());
}

template <bool depth_first>
size_t CompressedMicrotreeSplitRankArrayArithmetic<depth_first>::size() const {
    return length;
}

template <bool depth_first>
std::pair<TreeBP, uint32_t>
CompressedMicrotreeSplitRankArrayArithmetic<depth_first>::operator[](
    uint64_t index) const {
    assert(0 <= index && index < length);
    auto code = code_seq.read_interval(code_len.sum(index), code_len[index]);
    auto tree = left_seq_to_bp<depth_first>(
        arithmetic_to_left_seq<depth_first>(node_counts[index], code));
    return {tree, split_ranks[index]};
}

template <bool depth_first>
uint64_t CompressedMicrotreeSplitRankArrayArithmetic<
    depth_first>::evaluate_memory_consumption() const {
    return code_seq.evaluate_memory_consumption() +
           node_counts.evaluate_memory_consumption() +
           split_ranks.evaluate_memory_consumption() +
           code_len.evaluate_memory_consumption();
}

template <bool depth_first>
std::vector<std::pair<std::string, uint64_t>>
CompressedMicrotreeSplitRankArrayArithmetic<depth_first>::memory_table() const {
    return memory_table_combine_children(
        "Arithmetic",
        {
            {{"Sequence", code_seq.evaluate_memory_consumption()}},
            {{"Length", code_len.evaluate_memory_consumption()}},
            {{"SplitRank", split_ranks.evaluate_memory_consumption()}},
            {{"NodeCount", node_counts.evaluate_memory_consumption()}},
        });
}

template <bool depth_first>
uint32_t
CompressedMicrotreeSplitRankArrayArithmetic<depth_first>::get_node_count(
    uint64_t index) const {
    assert(0 <= index && index < length);
    return node_counts[index];
}

template <bool depth_first>
uint32_t CompressedMicrotreeSplitRankArrayArithmetic<
    depth_first>::get_left_chunk_popcount(uint64_t index) const {
    assert(0 <= index && index < length);
    return split_ranks[index];
}

template <bool depth_first>
uint32_t CompressedMicrotreeSplitRankArrayArithmetic<
    depth_first>::get_right_chunk_popcount(uint64_t index) const {
    assert(0 <= index && index < length);
    return node_counts[index] - split_ranks[index];
}

template <bool depth_first>
std::pair<uint32_t, uint32_t> CompressedMicrotreeSplitRankArrayArithmetic<
    depth_first>::get_node_count_left_popcount(uint64_t index) const {
    assert(0 <= index && index < length);
    return {node_counts[index], split_ranks[index]};
}

template <>
uint32_t CompressedMicrotreeSplitRankArrayArithmetic<true>::get_lca(
    uint64_t index, uint32_t u_inorder, uint32_t v_inorder) const {
    assert(0 <= index && index < length);

    uint32_t node_count = get_node_count(index);
    assert(0 <= u_inorder && u_inorder < node_count);
    assert(0 <= v_inorder && v_inorder < node_count);

    if (u_inorder == v_inorder) {
        return u_inorder;
    }

    if (u_inorder > v_inorder) {
        std::swap(u_inorder, v_inorder);
    }

    auto code = code_seq.read_interval(code_len.sum(index), code_len[index]);
    ArithmeticDecoder decoder(code);

    std::stack<uint32_t> sz_q, offset_q;
    sz_q.push(node_count);
    offset_q.push(0);
    while (!sz_q.empty()) {
        uint32_t sz = sz_q.top();
        sz_q.pop();
        uint32_t offset = offset_q.top();
        offset_q.pop();

        uint32_t lsz = decoder.decode_symbol(sz);
        uint32_t rsz = sz - lsz - 1;

        auto min_inorder = offset + lsz;
        if (u_inorder <= min_inorder && min_inorder <= v_inorder) {
            return min_inorder;
        }

        if (rsz) {
            sz_q.push(rsz);
            offset_q.push(offset + lsz + 1);
        }
        if (lsz) {
            sz_q.push(lsz);
            offset_q.push(offset);
        }
    }

    assert(false);
}

template <>
uint32_t CompressedMicrotreeSplitRankArrayArithmetic<false>::get_lca(
    uint64_t index, uint32_t u_inorder, uint32_t v_inorder) const {
    assert(0 <= index && index < length);

    uint32_t node_count = get_node_count(index);
    assert(0 <= u_inorder && u_inorder < node_count);
    assert(0 <= v_inorder && v_inorder < node_count);

    if (u_inorder == v_inorder) {
        return u_inorder;
    }

    if (u_inorder > v_inorder) {
        std::swap(u_inorder, v_inorder);
    }

    auto code = code_seq.read_interval(code_len.sum(index), code_len[index]);
    ArithmeticDecoder decoder(code);

    std::queue<uint32_t> sz_q, offset_q;
    sz_q.push(node_count);
    offset_q.push(0);
    while (!sz_q.empty()) {
        uint32_t sz = sz_q.front();
        sz_q.pop();
        uint32_t offset = offset_q.front();
        offset_q.pop();

        uint32_t lsz = decoder.decode_symbol(sz);
        uint32_t rsz = sz - lsz - 1;

        auto min_inorder = offset + lsz;
        if (u_inorder <= min_inorder && min_inorder <= v_inorder) {
            return min_inorder;
        }

        if (lsz) {
            sz_q.push(lsz);
            offset_q.push(offset);
        }
        if (rsz) {
            sz_q.push(rsz);
            offset_q.push(offset + lsz + 1);
        }
    }

    assert(false);
}

template struct CompressedMicrotreeSplitRankArrayArithmetic<true>;
template struct CompressedMicrotreeSplitRankArrayArithmetic<false>;

CompressedMicrotreeSplitRankArrayAllArithmetic::
    CompressedMicrotreeSplitRankArrayAllArithmetic()
    : length(0) {}

CompressedMicrotreeSplitRankArrayAllArithmetic::
    CompressedMicrotreeSplitRankArrayAllArithmetic(
        const MicrotreeSplitRankArray &microtree_split_rank_array)
    : length(microtree_split_rank_array.size()) {
    max_tree_n = 1;
    for (int i = 0; i < length; i++) {
        if (max_tree_n < microtree_split_rank_array[i].first.n) {
            max_tree_n = microtree_split_rank_array[i].first.n;
        }
    }

    auto encode = [&](const std::pair<TreeBP, uint32_t> &tree_and_split_rank)
        -> BitArray {
        const auto &[tree, split_rank] = tree_and_split_rank;
        ArithmeticEncoder encoder;
        // 1 <= tree.n <= max_tree_n
        encoder.encode_symbol(tree.n - 1, max_tree_n);
        // 0 <= split_rank <= tree.n
        encoder.encode_symbol(split_rank, tree.n + 1);
        encode_left_seq<false>(bp_to_left_seq<false>(tree), encoder);
        return encoder.terminate();
    };

    std::vector<uint32_t> code_len_vec(length);

    std::map<std::pair<TreeBP, uint32_t>, BitArray> arith_encode;
    for (int i = 0; i < length; i++) {
        auto tree_and_split_rank = microtree_split_rank_array[i];
        if (!arith_encode.count(tree_and_split_rank)) {
            arith_encode[tree_and_split_rank] = encode(tree_and_split_rank);
        }
        code_len_vec[i] = arith_encode[tree_and_split_rank].size();
    }

    code_len = code_len_vec;

    code_seq = BitArray(code_len.sum(length));
    uint32_t first = 0;
    for (int i = 0; i < length; i++) {
        const auto tree_and_split_rank = microtree_split_rank_array[i];
        const auto &code = arith_encode[tree_and_split_rank];
        code_seq.write_interval(first, code);
        first += code.size();
    }
    assert(first == code_seq.size());
}

std::pair<uint32_t, uint32_t>
CompressedMicrotreeSplitRankArrayAllArithmetic::get_node_count_left_popcount(
    uint64_t index) const {
    assert(0 <= index && index < length);
    ArithmeticDecoder decoder(code_seq, code_len.sum(index), code_len[index]);
    uint32_t node_count = decoder.decode_symbol(max_tree_n) + 1;
    uint32_t left_popcount = decoder.decode_symbol(node_count + 1);
    return {node_count, left_popcount};
}

size_t CompressedMicrotreeSplitRankArrayAllArithmetic::size() const {
    return length;
}

std::pair<TreeBP, uint32_t>
CompressedMicrotreeSplitRankArrayAllArithmetic::operator[](
    uint64_t index) const {
    assert(0 <= index && index < length);
    ArithmeticDecoder decoder(code_seq, code_len.sum(index), code_len[index]);
    uint32_t node_count = decoder.decode_symbol(max_tree_n) + 1;
    uint32_t split_rank = decoder.decode_symbol(node_count + 1);
    TreeBP tree =
        left_seq_to_bp<false>(decode_left_seq<false>(node_count, decoder));
    return {tree, split_rank};
}

uint32_t
CompressedMicrotreeSplitRankArrayAllArithmetic::get_left_chunk_popcount(
    uint64_t index) const {
    assert(0 <= index && index < length);
    return get_node_count_left_popcount(index).second;
}

uint32_t
CompressedMicrotreeSplitRankArrayAllArithmetic::get_right_chunk_popcount(
    uint64_t index) const {
    assert(0 <= index && index < length);
    auto [node_count, left_popcount] = get_node_count_left_popcount(index);
    return node_count - left_popcount;
}

uint64_t
CompressedMicrotreeSplitRankArrayAllArithmetic::evaluate_memory_consumption()
    const {
    return code_seq.evaluate_memory_consumption() +
           code_len.evaluate_memory_consumption();
}

std::vector<std::pair<std::string, uint64_t>>
CompressedMicrotreeSplitRankArrayAllArithmetic::memory_table() const {
    return memory_table_combine_children(
        "AllArith", {
                        {{"Sequence", code_seq.evaluate_memory_consumption()}},
                        {{"Length", code_len.evaluate_memory_consumption()}},
                    });
}

uint32_t CompressedMicrotreeSplitRankArrayAllArithmetic::get_node_count(
    uint64_t index) const {
    ArithmeticDecoder decoder(code_seq, code_len.sum(index), code_len[index]);
    uint32_t node_count = decoder.decode_symbol(max_tree_n) + 1;
    return node_count;
}

}  // namespace average_case_optimal_rmq
