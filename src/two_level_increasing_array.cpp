#include "hyperrmq/two_level_increasing_array.hpp"

#include <cassert>

namespace hyperrmq {

template <typename T, uint64_t K>
TwoLevelIncreasingArray<T, K>::TwoLevelIncreasingArray() : length(0) {}

template <typename T, uint64_t K>
TwoLevelIncreasingArray<T, K>::TwoLevelIncreasingArray(
    const std::vector<T>& arr)
    : length(arr.size()) {
    uint64_t superblock_size = (arr.size() + K - 1) / K,
             block_size = arr.size() - superblock_size;
    std::vector<T> block_vec(block_size), superblock_vec(superblock_size);
    for (int64_t i = 0; i < arr.size(); i++) {
        int64_t q = i / K;
        int64_t r = i % K;
        if (r == 0) {
            superblock_vec[q] = arr[i];
        } else {
            block_vec[q * (K - 1) + r - 1] = arr[i] - superblock_vec[q];
        }
    }
    block = std::move(block_vec);
    superblock = std::move(superblock_vec);
}

template <typename T, uint64_t K>
size_t TwoLevelIncreasingArray<T, K>::size() const {
    return length;
}

template <typename T, uint64_t K>
T TwoLevelIncreasingArray<T, K>::get(uint64_t index) const {
    assert(0 <= index && index < size());
    int64_t q = index / K, r = index % K;
    return superblock[q] + (r ? block[q * (K - 1) + r - 1] : 0);
}

template <typename T, uint64_t K>
T TwoLevelIncreasingArray<T, K>::operator[](uint64_t index) const {
    return get(index);
}

template <typename T, uint64_t K>
bool TwoLevelIncreasingArray<T, K>::operator==(
    const TwoLevelIncreasingArray& rhs) const {
    if (this->size() != rhs.size()) return false;
    for (int64_t i = 0; i < this->size(); i++) {
        if (this->get(i) != rhs.get(i)) return false;
    }
    return true;
}

template <typename T, uint64_t K>
bool TwoLevelIncreasingArray<T, K>::operator!=(
    const TwoLevelIncreasingArray& rhs) const {
    return !(*this == rhs);
}

template <typename T, uint64_t K>
uint64_t TwoLevelIncreasingArray<T, K>::evaluate_memory_consumption() const {
    return block.evaluate_memory_consumption() +
           superblock.evaluate_memory_consumption();
}

template struct TwoLevelIncreasingArray<uint64_t, 8>;
template struct TwoLevelIncreasingArray<uint64_t, 16>;
template struct TwoLevelIncreasingArray<uint64_t, 32>;
template struct TwoLevelIncreasingArray<uint64_t, 64>;

}  // namespace hyperrmq
