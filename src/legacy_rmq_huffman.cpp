#include "legacy_rmq_huffman.hpp"

#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>

#include "bitutil.hpp"
#include "memutil.hpp"
#include "tree_covering.hpp"

namespace average_case_optimal_rmq {

template <uint32_t W>
LegacyRMQHuffman<W>::LegacyRMQHuffman() {}

template <uint32_t W>
LegacyRMQHuffman<W>::LegacyRMQHuffman(const std::vector<int32_t>& values,
                                      const int B) {
    num_of_nodes = values.size();

    auto [upsilon, microtrees] =
        legacy_tree_covering(B, cartesian_tree_bp(values));

    num_of_microtrees = microtrees.size();
    num_of_chunks = num_of_microtrees * 2;

    std::map<uint64_t, uint64_t> counter;
    for (auto&& mt : microtrees) {
        counter[mt]++;
    }

    rmm_tree = RMMTree<8, 1024>(upsilon);
    chc = CanonicalHuffmanCode<>(counter);

    auto encode = chc.enumerate_alphabet_code_pair();
    uint32_t code_length_sum = 0;
    for (auto&& [mt, cnt] : counter) {
        code_length_sum += encode[mt].first * cnt;
    }

    std::vector<uint32_t> seq_idx_sample_vec;
    seq_idx_sample_vec.reserve((num_of_microtrees + W - 1) / W);
    code_seq = BitArray(code_length_sum);
    uint32_t idx = 0;
    for (int i = 0; i < num_of_microtrees; i++) {
        if (i % W == 0) {
            seq_idx_sample_vec.push_back(idx);
        }
        auto [len, code] = encode[microtrees[i]];
        code_seq.write_bits(idx, len, code);
        idx += len;
    }
    assert(idx == code_length_sum);
    seq_idx_sample = seq_idx_sample_vec;

    std::vector<uint32_t> close_sample_vec;
    close_sample_vec.reserve((num_of_chunks + W - 1) / W);
    uint32_t close = 0;
    // reverse the order and then use stack technique to speed up the process
    for (int i = 0; i < num_of_chunks; i++) {
        if (i % W == 0) {
            close_sample_vec.push_back(close);
        }
        if (rmm_tree.get_bit(i) == 0) {
            auto microtree = microtrees[rmm_tree.rank1(rmm_tree.close(i))];
            auto [cbp, cutpos] = legacy_decode_microtree(microtree);
            close += popcount(cbp & ((1ull << cutpos) - 1));
        } else {
            auto microtree = microtrees[rmm_tree.rank1(i)];
            auto [cbp, cutpos] = legacy_decode_microtree(microtree);
            close += popcount(cbp >> cutpos);
        }
    }
    close_sample = close_sample_vec;
}

template <uint32_t W>
uint32_t LegacyRMQHuffman<W>::to_inorder(uint32_t i) {
    uint32_t c = (rmm_tree.get_bit(i) == 0) ? rmm_tree.close(i) : i;
    uint32_t r = rmm_tree.rank1(c);
    return r;
}

template <uint32_t W>
uint32_t LegacyRMQHuffman<W>::find_seq_idx(uint32_t i) {
    assert(0 <= i && i < num_of_microtrees);
    uint32_t idx = seq_idx_sample[i / W];
    for (int j = 0; j < i % W; j++) {
        idx += chc.get_next_length(code_seq, idx);
    }
    return idx;
}

template <uint32_t W>
uint64_t LegacyRMQHuffman<W>::get_microtree_inorder(uint32_t i) {
    assert(0 <= i && i < num_of_microtrees);
    uint32_t idx = find_seq_idx(i);
    uint64_t len = chc.get_next_length(code_seq, idx);
    auto microtree = chc.decode(len, code_seq.read_bits(idx, len));
    return microtree;
}

template <uint32_t W>
std::pair<uint32_t, uint32_t> LegacyRMQHuffman<W>::get_chunk(uint32_t i) {
    assert(0 <= i && i < num_of_chunks);
    auto [bp, cutpos] =
        legacy_decode_microtree(get_microtree_inorder(to_inorder(i)));
    if (rmm_tree.get_bit(i) == 0) {
        return {cutpos, bp & ((1 << cutpos) - 1)};
    } else {
        return {2 * popcount(bp) - cutpos, bp >> cutpos};
    }
}

template <uint32_t W>
std::pair<uint32_t, uint32_t> LegacyRMQHuffman<W>::select1(uint32_t c) {
    assert(0 <= c && c < num_of_nodes);
    uint32_t low = 0, high = close_sample.size();
    uint32_t mid;
    while (high - low > 1) {
        mid = (low + high) / 2;
        if (close_sample[mid] <= c) {
            low = mid;
        } else {
            high = mid;
        }
    }

    int32_t remain = c - close_sample[low];
    uint32_t cidx = low * W;
    while (true) {
        auto [clen, cbp] = get_chunk(cidx);
        int32_t p = popcount(cbp);
        if (remain - p < 0) {
            int idx = 0;
            while (remain >= 0) {
                remain -= (cbp >> idx) & 1;
                idx++;
            }
            return {cidx, idx - 1};
        }
        remain -= p;
        cidx++;
    }
}

template <uint32_t W>
uint32_t LegacyRMQHuffman<W>::rank1(uint32_t c, uint32_t k) {
    assert(0 <= c && c < num_of_chunks);

    auto [clen, cbp] = get_chunk(c);
    assert(0 <= k && k <= clen);

    uint32_t res = close_sample[c / W];
    for (int i = c / W * W; i < c; i++) {
        res += popcount(get_chunk(i).second);
    }
    res += popcount(cbp & ((1ull << k) - 1));
    return res;
}

template <uint32_t W>
uint32_t LegacyRMQHuffman<W>::query(uint32_t i, uint32_t j) {
    assert(0 <= i && i <= j && j < num_of_nodes);

    auto [ic, ik] = select1(i);
    auto [jc, jk] = select1(j);

    uint32_t ii = to_inorder(ic);
    uint32_t ji = to_inorder(jc);

    if (rmm_tree.get_bit(ic) == 1) {
        ik += legacy_decode_microtree(get_microtree_inorder(ii)).second;
    }

    if (rmm_tree.get_bit(jc) == 1) {
        jk += legacy_decode_microtree(get_microtree_inorder(ji)).second;
    }

    auto rmq_on_bp = [&](uint64_t bp, uint32_t l, uint32_t m) -> uint32_t {
        assert(l <= m);
        int32_t minval = 1e9, minidx = 0;
        int32_t excess = 0;
        for (int i = l; i <= m; i++) {
            excess += 1 - 2 * ((bp >> i) & 1);
            if (minval > excess) {
                minval = excess;
                minidx = i;
            }
        }
        return minidx;
    };

    auto to_chunk_pair = [&](uint32_t lk, uint32_t cutpos, uint32_t ic,
                             uint32_t jc) -> uint32_t {
        if (lk < cutpos) {
            return rank1(ic, lk);
        } else {
            return rank1(jc, lk - cutpos);
        }
    };

    if (ii == ji) {
        auto [bp, cutpos] = legacy_decode_microtree(get_microtree_inorder(ii));
        auto lk = rmq_on_bp(bp, ik, jk);
        return to_chunk_pair(lk, cutpos, ic, jc);
    }

    uint32_t lc = rmm_tree.rmq(ic, jc + 1);
    uint32_t li = to_inorder(lc);

    if (li == ii) {
        auto [bp, cutpos] = legacy_decode_microtree(get_microtree_inorder(li));
        jc = lc;
        jk = cutpos - 1;
        if (((bp >> jk) & 1) == 0) {
            jk++;
        }
        auto lk = rmq_on_bp(bp, ik, jk);
        return to_chunk_pair(lk, cutpos, ic, jc);
    } else if (li == ji) {
        auto [bp, cutpos] = legacy_decode_microtree(get_microtree_inorder(li));
        ic = lc;
        ik = cutpos - 1;
        if (((bp >> ik) & 1) == 0) {
            ik++;
        }
        auto lk = rmq_on_bp(bp, ik, jk);
        return to_chunk_pair(lk, cutpos, ic, jc);
    } else {
        return rank1(lc, 0);
    }
}

template <uint32_t W>
uint64_t LegacyRMQHuffman<W>::evaluate_memory_consumption() const {
    return rmm_tree.evaluate_memory_consumption() +
           chc.evaluate_memory_consumption() +
           code_seq.evaluate_memory_consumption() +
           seq_idx_sample.evaluate_memory_consumption() +
           close_sample.evaluate_memory_consumption();
}

template struct LegacyRMQHuffman<32>;
template struct LegacyRMQHuffman<64>;
template struct LegacyRMQHuffman<128>;
template struct LegacyRMQHuffman<256>;
template struct LegacyRMQHuffman<512>;
template struct LegacyRMQHuffman<1024>;
template struct LegacyRMQHuffman<2048>;

}  // namespace average_case_optimal_rmq
