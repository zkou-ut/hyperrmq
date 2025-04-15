#pragma once

#include <cassert>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

#include "hyperrmq/bitutil.hpp"
#include "hyperrmq/memutil.hpp"
#include "hyperrmq/minimal_cell_array.hpp"
#include "hyperrmq/table.hpp"
#include "hyperrmq/tree_bp.hpp"

namespace hyperrmq {

template <uint64_t c = 8, uint64_t b = 1024>
struct RMMTree {
    static_assert(c == 8 || c == 16);
    static_assert(b % c == 0);
    using uintc_t = typename std::conditional<c == 8, uint8_t, uint16_t>::type;

    uint64_t num_of_bits, num_of_chunks, num_of_blocks, height, offset;
    std::vector<uintc_t> bp;
    MinimalCellArray<int64_t> excess_sample;
    MinimalCellArray<int64_t> min_tree;

    RMMTree();
    explicit RMMTree(const TreeBP &tree);

    uint64_t get_bit(uint64_t i) const;

    uint64_t idx_to_leaf(uint64_t k) const;
    uint64_t leaf_to_idx(uint64_t v) const;

    bool is_leftmost(uint64_t v) const;
    bool is_rightmost(uint64_t v) const;
    std::pair<uint64_t, uint64_t> vertex_to_interval(uint64_t v) const;
    int64_t node_excess(uint64_t v) const;

    uint64_t excess_prefix_sum(uint64_t n) const;

    uint64_t rank0(uint64_t n) const;
    uint64_t rank1(uint64_t n) const;

    uint64_t select0(uint64_t j) const;
    uint64_t select1(uint64_t j) const;

    std::pair<int64_t, uint64_t> fwdblock(uint64_t i, int64_t d) const;
    uint64_t fwdsearch(uint64_t i, int64_t d) const;
    uint64_t close(uint64_t i) const;

    std::pair<int64_t, uint64_t> bwdblock(uint64_t i, int64_t d) const;
    int64_t bwdsearch(uint64_t i, int64_t d) const;
    uint64_t open(uint64_t i) const;

    std::pair<int64_t, int64_t> minblock(uint64_t i, uint64_t j) const;
    int64_t minexcess(uint64_t i, uint64_t j) const;

    uint64_t rmq(uint64_t i, uint64_t j) const;

    uint64_t evaluate_memory_consumption() const;
    std::vector<std::pair<std::string, uint64_t>> memory_table() const;
};

template <uint64_t c, uint64_t b>
inline RMMTree<c, b>::RMMTree()
    : num_of_bits(0),
      num_of_chunks(0),
      num_of_blocks(0),
      height(0),
      offset(0) {}

template <uint64_t c, uint64_t b>
inline RMMTree<c, b>::RMMTree(const TreeBP &tree)
    : num_of_bits(tree.bp.size()),
      num_of_chunks((num_of_bits + c - 1) / c),
      num_of_blocks((num_of_bits + b - 1) / b),
      height(ceil_log2(num_of_blocks)),
      offset(1ull << height) {
    bp.resize(num_of_chunks);
    for (int64_t i = 0; i < num_of_chunks; i++) {
        bp[i] = tree.bp.read_bits_zero_follow(i * c, c);
    }

    std::vector<int64_t> excess_sample_vec(num_of_blocks + 1),
        min_tree_vec(num_of_blocks * 2 + 1);
    int64_t excess = 0;
    for (int64_t i = 0; i < num_of_blocks; i++) {
        uint64_t l = idx_to_leaf(i);
        for (int64_t j = i * (b / c);
             j < std::min((i + 1) * (b / c), num_of_chunks); j++) {
            min_tree_vec[l] =
                std::min(min_tree_vec[l],
                         min_table<c>[bp[j]] + excess - excess_sample_vec[i]);
            excess += excess_table<c>[bp[j]];
        }
        excess_sample_vec[i + 1] = excess;
    }
    excess_sample = excess_sample_vec;

    for (int64_t i = num_of_blocks - 1; i >= 1; i--) {
        min_tree_vec[i] = std::min(
            min_tree_vec[2 * i], node_excess(2 * i) + min_tree_vec[2 * i + 1]);
    }
    min_tree_vec[2 * num_of_blocks] = num_of_bits + 1;
    min_tree = min_tree_vec;
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::get_bit(uint64_t i) const {
    assert(0 <= i && i < num_of_bits);
    return (bp[i / c] >> (c - 1 - (i % c))) & 1;
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::idx_to_leaf(uint64_t k) const {
    assert(0 <= k && k < num_of_blocks);
    if (k < 2 * num_of_blocks - offset) {
        return k + offset;
    } else {
        return k + offset - num_of_blocks;
    }
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::leaf_to_idx(uint64_t v) const {
    assert(num_of_blocks <= v && v < num_of_blocks * 2);
    if (v >= offset) {
        return v - offset;
    } else {
        return v + num_of_blocks - offset;
    }
}

template <uint64_t c, uint64_t b>
inline bool RMMTree<c, b>::is_leftmost(uint64_t v) const {
    // v is pow of 2
    return (v & (v - 1)) == 0;
}

template <uint64_t c, uint64_t b>
inline bool RMMTree<c, b>::is_rightmost(uint64_t v) const {
    // v + 1 is pow of 2
    return (v & (v + 1)) == 0;
}

template <uint64_t c, uint64_t b>
inline std::pair<uint64_t, uint64_t> RMMTree<c, b>::vertex_to_interval(
    uint64_t v) const {
    assert(1 <= v && v < 2 * num_of_blocks);

    if (v >= num_of_blocks) {
        uint64_t idx = leaf_to_idx(v);
        return {idx, idx + 1};
    }

    auto leftmost_leaf = [&](uint64_t u) -> uint64_t {
        while (u < num_of_blocks) {
            u <<= 1;
        }
        return u;
    };

    if (is_rightmost(v)) {
        return {leaf_to_idx(leftmost_leaf(v)), num_of_blocks};
    } else {
        return {leaf_to_idx(leftmost_leaf(v)),
                leaf_to_idx(leftmost_leaf(v + 1))};
    }
}

template <uint64_t c, uint64_t b>
inline int64_t RMMTree<c, b>::node_excess(uint64_t v) const {
    auto [l, r] = vertex_to_interval(v);
    return excess_sample[r] - excess_sample[l];
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::excess_prefix_sum(uint64_t n) const {
    assert(0 <= n && n <= num_of_bits);
    uint64_t ret = excess_sample[n / b];
    for (int64_t i = n / b * (b / c); i < n / c; i++) {
        ret += excess_table<c>[bp[i]];
    }
    if (n % c) {
        ret += excess_table<c>[bp[n / c] >> (c - (n % c))] - (c - (n % c));
    }
    return ret;
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::rank0(uint64_t n) const {
    return (n + excess_prefix_sum(n)) >> 1;
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::rank1(uint64_t n) const {
    return (n - excess_prefix_sum(n)) >> 1;
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::select0(uint64_t j) const {
    assert(0 <= j && j < num_of_bits / 2);
    uint64_t low = 0, high = num_of_blocks;
    uint64_t mid;
    while (high - low > 1) {
        mid = (low + high) / 2;
        if ((mid * b + excess_sample[mid]) / 2 <= j) {
            low = mid;
        } else {
            high = mid;
        }
    }

    int64_t remain = j - (low * b + excess_sample[low]) / 2;
    uint64_t cidx = low * (b / c);
    while (true) {
        int64_t nxt = (c + excess_table<c>[bp[cidx]]) / 2;
        if (remain - nxt < 0) {
            break;
        }
        remain -= nxt;
        cidx++;
    }

    uint64_t idx = 0;
    while (remain >= 0) {
        remain -= ((~bp[cidx] >> (c - 1 - idx)) & 1);
        idx++;
    }
    return cidx * c + idx - 1;
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::select1(uint64_t j) const {
    assert(0 <= j && j < num_of_bits / 2);
    uint64_t low = 0, high = num_of_blocks;
    uint64_t mid;
    while (high - low > 1) {
        mid = (low + high) / 2;
        if ((mid * b - excess_sample[mid]) / 2 <= j) {
            low = mid;
        } else {
            high = mid;
        }
    }

    int64_t remain = j - (low * b - excess_sample[low]) / 2;
    uint64_t cidx = low * (b / c);
    while (true) {
        int64_t nxt = (c - excess_table<c>[bp[cidx]]) / 2;
        if (remain - nxt < 0) {
            break;
        }
        remain -= nxt;
        cidx++;
    }

    uint64_t idx = 0;
    while (remain >= 0) {
        remain -= ((bp[cidx] >> (c - 1 - idx)) & 1);
        idx++;
    }
    return cidx * c + idx - 1;
}

template <uint64_t c, uint64_t b>
inline std::pair<int64_t, uint64_t> RMMTree<c, b>::fwdblock(uint64_t i,
                                                            int64_t d) const {
    assert(0 <= i && i < num_of_bits);
    assert(d < 0);

    int64_t f = i / c;
    int64_t t = std::min((i / b + 1) * (b / c), num_of_chunks);
    int64_t d_prime = 0;
    for (int64_t j = i % c; j < c; j++) {
        d_prime += 1 - 2 * ((bp[f] >> (c - j - 1)) & 1);
        if (d_prime == d) {
            return {d, f * c + j + 1};
        }
    }

    int64_t p;
    for (p = f + 1; p < t; p++) {
        if (d_prime + min_table<c>[bp[p]] <= d) {
            break;
        }
        d_prime += excess_table<c>[bp[p]];
    }
    if (p >= t) {
        return {d_prime, t * c};
    }

    for (int64_t j = 0; j < c; j++) {
        d_prime += 1 - 2 * ((bp[p] >> (c - j - 1)) & 1);
        if (d_prime == d) {
            return {d, p * c + j + 1};
        }
    }

    assert(false);
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::fwdsearch(uint64_t i, int64_t d) const {
    assert(0 <= i && i < num_of_bits);
    assert(d < 0);

    auto [d_prime, j] = fwdblock(i, d);
    if (d_prime == d) {
        return j;
    }

    int64_t v = idx_to_leaf(i / b);
    while (!is_rightmost(v) && d_prime + min_tree[v + 1] > d) {
        if ((v & 1) == 0) {
            d_prime += node_excess(v + 1);
        }
        v /= 2;
    }
    if (is_rightmost(v)) {
        return num_of_bits + 1;
    }

    v++;
    while (v < num_of_blocks) {
        if (d_prime + min_tree[2 * v] <= d) {
            v = 2 * v;
        } else {
            d_prime += node_excess(2 * v);
            v = 2 * v + 1;
        }
    }

    auto [d_double_prime, j2] = fwdblock(leaf_to_idx(v) * b, d - d_prime);
    assert(d_double_prime == d - d_prime);
    return j2;
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::close(uint64_t i) const {
    assert(get_bit(i) == 0);
    return fwdsearch(i + 1, -1) - 1;
}

template <uint64_t c, uint64_t b>
inline std::pair<int64_t, uint64_t> RMMTree<c, b>::bwdblock(uint64_t i,
                                                            int64_t d) const {
    assert(0 < i && i <= num_of_bits);
    assert(d < 0);

    i--;

    int64_t t = i / c;
    int64_t f = (i / b) * (b / c);
    int64_t d_prime = 0;
    for (int64_t j = i % c; j >= 0; j--) {
        d_prime -= 1 - 2 * ((bp[t] >> (c - j - 1)) & 1);
        if (d_prime == d) {
            return {d, t * c + j};
        }
    }

    int64_t p;
    for (p = t - 1; p >= f; p--) {
        if (d_prime - excess_table<c>[bp[p]] + min_table<c>[bp[p]] <= d) {
            break;
        }
        d_prime -= excess_table<c>[bp[p]];
    }
    if (p < f) {
        return {d_prime, f * c};
    }

    for (int64_t j = c - 1; j >= 0; j--) {
        d_prime -= 1 - 2 * ((bp[p] >> (c - j - 1)) & 1);
        if (d_prime == d) {
            return {d, p * c + j};
        }
    }

    assert(false);
}

template <uint64_t c, uint64_t b>
inline int64_t RMMTree<c, b>::bwdsearch(uint64_t i, int64_t d) const {
    assert(0 < i && i <= num_of_bits);
    assert(d < 0);

    auto [d_prime, j] = bwdblock(i, d);
    if (d_prime == d) {
        return j;
    }

    int64_t v = idx_to_leaf((i - 1) / b);
    while (!is_leftmost(v) &&
           d_prime - node_excess(v - 1) + min_tree[v - 1] > d) {
        if ((v & 1) == 1) {
            d_prime -= node_excess(v - 1);
        }
        v /= 2;
    }
    if (is_leftmost(v)) {
        return -1;
    }

    v--;
    while (v < num_of_blocks) {
        if (d_prime - node_excess(2 * v + 1) + min_tree[2 * v + 1] <= d) {
            v = 2 * v + 1;
        } else {
            d_prime -= node_excess(2 * v + 1);
            v = 2 * v;
        }
    }

    auto [d_double_prime, j2] = bwdblock((leaf_to_idx(v) + 1) * b, d - d_prime);
    assert(d_double_prime == d - d_prime);
    return j2;
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::open(uint64_t i) const {
    assert(get_bit(i) == 1);
    return bwdsearch(i, -1);
}

template <uint64_t c, uint64_t b>
inline std::pair<int64_t, int64_t> RMMTree<c, b>::minblock(uint64_t i,
                                                           uint64_t j) const {
    assert(i <= j);
    assert(i == j || i / b == (j - 1) / b);
    uint64_t k = i / c;
    uint64_t l = j / c;
    int64_t d = 0, m = 0;
    for (int64_t p = i; p < std::min(j, (k + 1) * c); p++) {
        d += 1 - 2 * ((bp[k] >> (c - (p % c) - 1)) & 1);
        m = std::min(m, d);
    }
    if (k == l) {
        return {m, d};
    }

    for (int64_t q = k + 1; q < l; q++) {
        m = std::min(m, d + min_table<c>[bp[q]]);
        d += excess_table<c>[bp[q]];
    }

    for (int64_t p = l * c; p < j; p++) {
        d += 1 - 2 * ((bp[l] >> (c - (p % c) - 1)) & 1);
        m = std::min(m, d);
    }

    return {m, d};
}

template <uint64_t c, uint64_t b>
inline int64_t RMMTree<c, b>::minexcess(uint64_t i, uint64_t j) const {
    assert(0 <= i && i <= j && j <= num_of_bits);
    int64_t k = i / b;
    auto [m, d] = minblock(i, std::min((k + 1) * b, j));
    if (j <= (k + 1) * b) {
        return m;
    }

    int64_t k_prime = j / b;
    int64_t v = idx_to_leaf(k), l = idx_to_leaf(k_prime);

    auto contain = [&](int64_t u) -> bool {
        return (l >> (ceil_log2(l + 1) - ceil_log2(u + 1))) != u;
    };

    while (v + 1 > l || contain(v + 1)) {
        if ((v & 1) == 0) {
            m = std::min(m, d + min_tree[v + 1]);
            d += node_excess(v + 1);
        }
        v >>= 1;
    }

    v++;
    while (v < num_of_blocks) {
        if (d + min_tree[v] >= m) {
            return m;
        }
        if (contain(2 * v)) {
            m = std::min(m, d + min_tree[2 * v]);
            d += node_excess(2 * v);
            v = 2 * v + 1;
        } else {
            v = 2 * v;
        }
    }
    if (d + min_tree[v] >= m) {
        return m;
    }
    auto [m_prime, d_prime] = minblock(k_prime * b, j);
    return std::min(m, d + m_prime);
}

// It searches half-open interval [i, j).
template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::rmq(uint64_t i, uint64_t j) const {
    assert(0 <= i && i < j && j <= num_of_bits);
    int64_t m = minexcess(i, j);
    if (m == 0) {
        return i;
    }
    return fwdsearch(i, m) - 1;
}

template <uint64_t c, uint64_t b>
inline uint64_t RMMTree<c, b>::evaluate_memory_consumption() const {
    return evaluate_vector_memory_consumption(bp) +
           excess_sample.evaluate_memory_consumption() +
           min_tree.evaluate_memory_consumption();
}

template <uint64_t c, uint64_t b>
inline std::vector<std::pair<std::string, uint64_t>>
RMMTree<c, b>::memory_table() const {
    return memory_table_combine_children(
        "RMMTree",
        {
            {{"rmmmain", evaluate_vector_memory_consumption(bp)}},
            {{"excesssample", excess_sample.evaluate_memory_consumption()}},
            {{"mintree", min_tree.evaluate_memory_consumption()}},
        });
}

}  // namespace hyperrmq
