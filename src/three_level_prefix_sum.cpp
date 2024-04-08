#include "hyperrmq/three_level_prefix_sum.hpp"

#include <cassert>

#include "hyperrmq/memutil.hpp"

namespace average_case_optimal_rmq {

template <typename T, uint32_t K, uint32_t L>
ThreeLevelPrefixSum<T, K, L>::ThreeLevelPrefixSum() : length(0) {}

template <typename T, uint32_t K, uint32_t L>
ThreeLevelPrefixSum<T, K, L>::ThreeLevelPrefixSum(const std::vector<T>& arr)
    : length(arr.size()), block(arr) {
    uint32_t superblock_size = length / K;
    std::vector<T> superblock_vec;
    superblock_vec.reserve(superblock_size);
    T s = 0;
    for (int i = 0; i < length; i++) {
        s += arr[i];
        if ((i + 1) % K == 0) {
            superblock_vec.push_back(s);
        }
    }
    superblock = std::move(superblock_vec);
}

template <typename T, uint32_t K, uint32_t L>
size_t ThreeLevelPrefixSum<T, K, L>::size() const {
    return length;
}

template <typename T, uint32_t K, uint32_t L>
T ThreeLevelPrefixSum<T, K, L>::operator[](uint64_t index) const {
    assert(0 <= index && index < length);
    return block[index];
}

template <typename T, uint32_t K, uint32_t L>
T ThreeLevelPrefixSum<T, K, L>::sum(uint64_t index) const {
    assert(0 <= index && index <= length);
    uint64_t q = index / K;
    T ret = q ? superblock[q - 1] : 0;
    for (int i = q * K; i < index; i++) {
        ret += block[i];
    }
    return ret;
}

template <typename T, uint32_t K, uint32_t L>
uint64_t ThreeLevelPrefixSum<T, K, L>::select_chunk(T v) const {
    assert(v >= 0);
    int64_t low = 0, high = superblock.size();
    while (high - low > 1) {
        int64_t mid = (low + high) / 2;
        if (superblock[mid - 1] <= v) {
            low = mid;
        } else {
            high = mid;
        }
    }

    T p = (low == 0 ? 0 : superblock[low - 1]);

    for (int i = low * K; i < block.size(); i++) {
        p += block[i];
        if (v < p) {
            return i;
        }
    }
    return block.size();
}

template <typename T, uint32_t K, uint32_t L>
bool ThreeLevelPrefixSum<T, K, L>::operator==(
    const ThreeLevelPrefixSum& rhs) const {
    if (this->size() != rhs.size()) return false;
    for (int i = 0; i < this->size(); i++) {
        if ((*this)[i] != rhs[i]) return false;
    }
    return true;
}

template <typename T, uint32_t K, uint32_t L>
bool ThreeLevelPrefixSum<T, K, L>::operator!=(
    const ThreeLevelPrefixSum& rhs) const {
    return !(*this == rhs);
}

template <typename T, uint32_t K, uint32_t L>
uint64_t ThreeLevelPrefixSum<T, K, L>::evaluate_memory_consumption() const {
    return block.evaluate_memory_consumption() +
           superblock.evaluate_memory_consumption();
}

template <typename T, uint32_t K, uint32_t L>
std::vector<std::pair<std::string, uint64_t>>
ThreeLevelPrefixSum<T, K, L>::memory_table() const {
    return memory_table_combine_children(
        "ThreeLevel", {
                          {{"Block", block.evaluate_memory_consumption()}},
                          {{"Super", superblock.evaluate_memory_consumption()}},
                      });
}

template struct ThreeLevelPrefixSum<uint32_t, 8, 8>;
template struct ThreeLevelPrefixSum<uint32_t, 16, 16>;
template struct ThreeLevelPrefixSum<uint32_t, 32, 32>;
template struct ThreeLevelPrefixSum<uint32_t, 64, 64>;

}  // namespace average_case_optimal_rmq
