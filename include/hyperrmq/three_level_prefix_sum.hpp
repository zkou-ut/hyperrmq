#pragma once

#include "hyperrmq/two_level_increasing_array.hpp"

namespace hyperrmq {

template <typename T = uint32_t, uint32_t K = 32, uint32_t L = 32>
struct ThreeLevelPrefixSum {
   public:
    ThreeLevelPrefixSum();
    ThreeLevelPrefixSum(const std::vector<T>& arr);

    size_t size() const;

    T operator[](uint64_t index) const;
    T sum(uint64_t index) const;

    // It returns first j s.t. v < sum(j) (same as std::upper_bound).
    uint64_t select_chunk(T v) const;

    bool operator==(const ThreeLevelPrefixSum& rhs) const;
    bool operator!=(const ThreeLevelPrefixSum& rhs) const;

    uint64_t evaluate_memory_consumption() const;
    std::vector<std::pair<std::string, uint64_t>> memory_table() const;

    //    private:
    uint64_t length;
    MinimalCellArray<T> block;
    TwoLevelIncreasingArray<T, L> superblock;
};

}  // namespace hyperrmq
