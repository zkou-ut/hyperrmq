#pragma once

#include "hyperrmq/bit_array.hpp"

namespace average_case_optimal_rmq {

template <typename T>
struct MinimalCellArray {
   public:
    MinimalCellArray();
    MinimalCellArray(const std::vector<T>& arr, uint32_t cell_width = 0,
                     T cell_offset = 0);

    size_t size() const;

    T get(uint64_t index) const;
    T operator[](uint64_t index) const;

    bool operator==(const MinimalCellArray& rhs) const;
    bool operator!=(const MinimalCellArray& rhs) const;

    uint64_t evaluate_memory_consumption() const;

   private:
    uint32_t width;
    uint64_t length;
    T offset;
    BitArray bit_array;
};

}  // namespace average_case_optimal_rmq
