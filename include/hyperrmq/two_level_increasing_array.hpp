#pragma once

#include "hyperrmq/minimal_cell_array.hpp"

namespace average_case_optimal_rmq {

template <typename T = uint32_t, uint32_t K = 32>
struct TwoLevelIncreasingArray {
   public:
    TwoLevelIncreasingArray();
    TwoLevelIncreasingArray(const std::vector<T>& arr);

    size_t size() const;

    T get(uint64_t index) const;
    T operator[](uint64_t index) const;

    bool operator==(const TwoLevelIncreasingArray& rhs) const;
    bool operator!=(const TwoLevelIncreasingArray& rhs) const;

    uint64_t evaluate_memory_consumption() const;

   private:
    uint64_t length;
    MinimalCellArray<T> block, superblock;
};

}  // namespace average_case_optimal_rmq
