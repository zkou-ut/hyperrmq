#include "hyperrmq/minimal_cell_array.hpp"

#include <algorithm>
#include <cassert>

#include "hyperrmq/bitutil.hpp"

namespace hyperrmq {

template <typename T>
MinimalCellArray<T>::MinimalCellArray() : width(0), length(0), offset(0) {}

template <typename T>
MinimalCellArray<T>::MinimalCellArray(const std::vector<T>& arr,
                                      uint32_t cell_width, T cell_offset)
    : width(cell_width), length(arr.size()), offset(cell_offset) {
    if (length == 0) {
        return;
    }
    if (width == 0 && offset == 0) {
        auto [min_it, max_it] = std::minmax_element(arr.begin(), arr.end());
        offset = *min_it;
        width = ceil_log2((*max_it) - offset + 1);
    }
    bit_array = BitArray(width * length);
    if (width) {
        for (int i = 0; i < length; i++) {
            bit_array.write_bits(i * width, width, arr[i] - offset);
        }
    }
}

template <typename T>
size_t MinimalCellArray<T>::size() const {
    return length;
}

template <typename T>
T MinimalCellArray<T>::get(uint64_t index) const {
    assert(0 <= index && index < length);
    return (width ? bit_array.read_bits(index * width, width) : 0) + offset;
}

template <typename T>
T MinimalCellArray<T>::operator[](uint64_t index) const {
    return get(index);
}

template <typename T>
bool MinimalCellArray<T>::operator==(const MinimalCellArray& rhs) const {
    if (this->size() != rhs.size()) return false;
    for (int i = 0; i < this->size(); i++) {
        if (this->get(i) != rhs.get(i)) return false;
    }
    return true;
}

template <typename T>
bool MinimalCellArray<T>::operator!=(const MinimalCellArray& rhs) const {
    return !(*this == rhs);
}

template <typename T>
uint64_t MinimalCellArray<T>::evaluate_memory_consumption() const {
    return bit_array.evaluate_memory_consumption();
}

template struct MinimalCellArray<uint8_t>;
template struct MinimalCellArray<uint16_t>;
template struct MinimalCellArray<uint32_t>;
template struct MinimalCellArray<uint64_t>;

template struct MinimalCellArray<int8_t>;
template struct MinimalCellArray<int16_t>;
template struct MinimalCellArray<int32_t>;
template struct MinimalCellArray<int64_t>;

}  // namespace hyperrmq
