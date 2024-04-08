#pragma once

#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <stack>
#include <vector>

#include "hyperrmq/bit_array.hpp"
#include "hyperrmq/tree_bp.hpp"

namespace average_case_optimal_rmq {

// It encodes interval by using a tag in the interval.
// It saves the last `1` and following zeros.
// The decoder knows it as well.
struct ArithmeticEncoder {
    static constexpr uint32_t m = 32;
    static constexpr uint64_t mask = (1ull << m) - 1;

    ArithmeticEncoder() {}

    void encode_symbol(uint32_t symbol, uint32_t symbol_count) {
        assert(0 <= symbol && symbol < symbol_count);

        uint64_t interval_width = u - l + 1;
        u = l + interval_width * (symbol + 1) / symbol_count - 1;
        l = l + interval_width * symbol / symbol_count;
        while (true) {
            if ((l >> (m - 1)) == (u >> (m - 1))) {
                bool b = l >> (m - 1);
                bs.push(b);
                l = (l << 1) & mask;
                u = ((u << 1) & mask) | 1;
                while (scale3 > 0) {
                    bs.push(!b);
                    scale3--;
                }
            } else if ((l >> (m - 2)) == 0b01 && (u >> (m - 2)) == 0b10) {
                scale3++;
                l = (l << 1) & mask;
                u = ((u << 1) & mask) | 1;
                l ^= (1ull << (m - 1));
                u ^= (1ull << (m - 1));
            } else {
                break;
            }
        }
    }

    BitArray terminate() {
        // It saves the last `1`.
        // The decoder knows it as well.
        // bs.push(1);

        return BitArray(std::move(bs));
    }

    uint64_t l = 0;
    uint64_t u = mask;
    int scale3 = 0;

    BitStack bs;
};

// It decodes the tag in the interval by multiplying specified
// `symbol_count`. It knows that the last `1` and following zeros are
// omitted in the tag.
struct ArithmeticDecoder {
    static constexpr uint32_t m = 32;
    static constexpr uint64_t mask = (1ull << m) - 1;

    ArithmeticDecoder(const BitArray& tag)
        : ArithmeticDecoder(tag, 0, tag.size()) {}
    ArithmeticDecoder(const BitArray& code_seq, uint64_t start, uint64_t width)
        : code_seq(code_seq),
          tag_start(start),
          tag_length(width),
          tag_idx(tag_start + m) {
        t = code_seq.read_bits_zero_follow(start, m);
        if (width < m) {
            t = ((t >> (m - 1 - width)) | 1) << (m - 1 - width);
        }
    }

    uint32_t decode_symbol(uint32_t symbol_count) {
        uint32_t symbol = ((t - l + 1) * symbol_count - 1) / (u - l + 1);

        uint64_t interval_width = u - l + 1;
        u = l + interval_width * (symbol + 1) / symbol_count - 1;
        l = l + interval_width * symbol / symbol_count;

        auto next_bit = [&]() -> uint64_t {
            if (tag_idx < tag_start + tag_length) {
                return code_seq.get(tag_idx++);
            } else if (tag_idx == tag_start + tag_length) {
                tag_idx++;
                return 1;
            } else {
                return 0;
            }
        };

        while (true) {
            if ((l >> (m - 1)) == (u >> (m - 1))) {
                l = (l << 1) & mask;
                u = ((u << 1) & mask) | 1;
                t = ((t << 1) & mask) | next_bit();
            } else if ((l >> (m - 2)) == 0b01 && (u >> (m - 2)) == 0b10) {
                l = (l << 1) & mask;
                u = ((u << 1) & mask) | 1;
                t = ((t << 1) & mask) | next_bit();
                l ^= (1ull << (m - 1));
                u ^= (1ull << (m - 1));
                t ^= (1ull << (m - 1));
            } else {
                break;
            }
        }

        return symbol;
    }

    uint64_t l = 0;
    uint64_t u = mask;
    int scale3 = 0;

    uint64_t t;
    const BitArray& code_seq;
    uint64_t tag_start, tag_length, tag_idx;
};

template <bool depth_first>
std::vector<uint32_t> bp_to_left_seq(const TreeBP& tree) {
    std::vector<uint32_t> depth_first_left_seq(tree.n);

    int64_t pre = 0;
    int64_t index = 0;

    auto match = [&](int v) -> void {
        assert(tree.bp.get(index) == v);
        index++;
    };
    auto has_child = [&]() -> bool {
        return index < 2 * tree.n && tree.bp.get(index) == 0;
    };

    auto dfs = [&](auto self) -> uint32_t {
        int v = pre;
        pre++;

        uint32_t sz = 1;

        match(0);
        if (has_child()) {
            depth_first_left_seq[v] = self(self);
            sz += depth_first_left_seq[v];
        }

        match(1);
        if (has_child()) {
            sz += self(self);
        }

        return sz;
    };

    dfs(dfs);

    if constexpr (depth_first) {
        return depth_first_left_seq;
    } else {
        std::vector<uint32_t> breadth_first_left_seq;
        breadth_first_left_seq.reserve(tree.n);
        std::queue<uint32_t> sz_q, idx_q;
        sz_q.push(tree.n);
        idx_q.push(0);
        for (int i = 0; i < tree.n; i++) {
            uint32_t sz = sz_q.front();
            sz_q.pop();
            uint32_t idx = idx_q.front();
            idx_q.pop();

            uint32_t lsz = depth_first_left_seq[idx];
            uint32_t rsz = sz - lsz - 1;
            breadth_first_left_seq.push_back(lsz);

            if (lsz) {
                sz_q.push(lsz);
                idx_q.push(idx + 1);
            }
            if (rsz) {
                sz_q.push(rsz);
                idx_q.push(idx + lsz + 1);
            }
        }

        return breadth_first_left_seq;
    }
}

template <bool depth_first>
TreeBP left_seq_to_bp(const std::vector<uint32_t>& left_seq);

template <>
inline TreeBP left_seq_to_bp<true>(const std::vector<uint32_t>& left_seq) {
    uint32_t n = left_seq.size();
    BitArray bp(2 * n);
    uint32_t pos = 0;

    int v = -1;
    auto dfs = [&](auto self, uint32_t sz) -> void {
        v++;
        uint32_t lsz = left_seq[v];
        uint32_t rsz = sz - lsz - 1;

        pos++;
        if (lsz) {
            self(self, lsz);
        }
        bp.on(pos++);
        if (rsz) {
            self(self, rsz);
        }
    };

    dfs(dfs, n);
    assert(pos == 2 * n);

    return TreeBP(n, bp);
}

template <>
inline TreeBP left_seq_to_bp<false>(const std::vector<uint32_t>& left_seq) {
    uint32_t n = left_seq.size();
    BitArray bp(2 * n);

    std::queue<int> pre_q, sz_q, write_q;
    pre_q.push(0);
    sz_q.push(n);
    write_q.push(0);
    for (int i = 0; i < n; i++) {
        uint32_t pre = pre_q.front();
        pre_q.pop();
        uint32_t sz = sz_q.front();
        sz_q.pop();
        uint32_t write = write_q.front();
        write_q.pop();

        uint32_t lsz = left_seq[i];
        uint32_t rsz = sz - lsz - 1;

        bp.on(write + lsz * 2 + 1);
        if (lsz) {
            pre_q.push(pre + 1);
            sz_q.push(lsz);
            write_q.push(write + 1);
        }
        if (rsz) {
            pre_q.push(pre + lsz + 1);
            sz_q.push(rsz);
            write_q.push(write + lsz * 2 + 2);
        }
    }

    return TreeBP(n, bp);
}

template <bool depth_first>
void encode_left_seq(const std::vector<uint32_t>& left_seq,
                     ArithmeticEncoder& encoder);

template <>
inline void encode_left_seq<true>(const std::vector<uint32_t>& left_seq,
                                  ArithmeticEncoder& encoder) {
    uint32_t node_count = left_seq.size();

    std::stack<uint32_t> sz_q;
    sz_q.push(node_count);

    int v = 0;
    while (!sz_q.empty()) {
        uint32_t sz = sz_q.top();
        sz_q.pop();

        uint32_t lsz = left_seq[v];
        uint32_t rsz = sz - lsz - 1;
        v++;

        if (rsz) {
            sz_q.push(rsz);
        }
        if (lsz) {
            sz_q.push(lsz);
        }

        encoder.encode_symbol(lsz, sz);
    }
}

template <>
inline void encode_left_seq<false>(const std::vector<uint32_t>& left_seq,
                                   ArithmeticEncoder& encoder) {
    uint32_t node_count = left_seq.size();

    std::queue<uint32_t> sz_q;
    sz_q.push(node_count);

    int v = 0;
    while (!sz_q.empty()) {
        uint32_t sz = sz_q.front();
        sz_q.pop();

        uint32_t lsz = left_seq[v];
        uint32_t rsz = sz - lsz - 1;
        v++;

        if (lsz) {
            sz_q.push(lsz);
        }
        if (rsz) {
            sz_q.push(rsz);
        }

        encoder.encode_symbol(lsz, sz);
    }
}

template <bool depth_first>
BitArray left_seq_to_arithmetic(const std::vector<uint32_t>& left_seq) {
    ArithmeticEncoder encoder;

    encode_left_seq<depth_first>(left_seq, encoder);

    return encoder.terminate();
}

template <bool depth_first>
std::vector<uint32_t> decode_left_seq(const uint32_t node_count,
                                      ArithmeticDecoder& decoder);

template <>
inline std::vector<uint32_t> decode_left_seq<true>(const uint32_t node_count,
                                                   ArithmeticDecoder& decoder) {
    std::vector<uint32_t> left_seq;
    left_seq.reserve(node_count);

    std::stack<uint32_t> sz_q;
    sz_q.push(node_count);

    while (!sz_q.empty()) {
        uint32_t sz = sz_q.top();
        sz_q.pop();

        uint32_t lsz = decoder.decode_symbol(sz);
        uint32_t rsz = sz - lsz - 1;
        left_seq.push_back(lsz);

        if (rsz) {
            sz_q.push(rsz);
        }
        if (lsz) {
            sz_q.push(lsz);
        }
    }

    return left_seq;
}

template <>
inline std::vector<uint32_t> decode_left_seq<false>(
    const uint32_t node_count, ArithmeticDecoder& decoder) {
    std::vector<uint32_t> left_seq;
    left_seq.reserve(node_count);

    std::queue<uint32_t> sz_q;
    sz_q.push(node_count);

    while (!sz_q.empty()) {
        uint32_t sz = sz_q.front();
        sz_q.pop();

        uint32_t lsz = decoder.decode_symbol(sz);
        uint32_t rsz = sz - lsz - 1;
        left_seq.push_back(lsz);

        if (lsz) {
            sz_q.push(lsz);
        }
        if (rsz) {
            sz_q.push(rsz);
        }
    }

    return left_seq;
}

template <bool depth_first>
std::vector<uint32_t> arithmetic_to_left_seq(const uint32_t node_count,
                                             const BitArray& code) {
    ArithmeticDecoder decoder(code);

    return decode_left_seq<depth_first>(node_count, decoder);
}

}  // namespace average_case_optimal_rmq
