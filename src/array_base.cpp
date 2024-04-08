#include "hyperrmq/array_base.hpp"

#include <cassert>

#include "hyperrmq/tree_bp.hpp"

namespace average_case_optimal_rmq {

template <typename T>
FixedLengthCodeArrayBase<T>::FixedLengthCodeArrayBase()
    : width(0), length(0), codes() {}

template <typename T>
FixedLengthCodeArrayBase<T>::FixedLengthCodeArrayBase(uint32_t width,
                                                      uint64_t length)
    : width(width), length(length), codes(width * length) {}

template <typename T>
FixedLengthCodeArrayBase<T>::FixedLengthCodeArrayBase(uint32_t width,
                                                      std::vector<T> vec)
    : width(width), length(vec.size()), codes(width * length) {
    for (int index = 0; index < length; index++) {
        set(index, vec[index]);
    }
}

template <typename T>
size_t FixedLengthCodeArrayBase<T>::size() const {
    return length;
}

template <typename T>
void FixedLengthCodeArrayBase<T>::set(uint64_t index, T object) {
    assert(0 <= index && index < length);
    auto code = encode(object);
    assert(code.size() == width);
    codes.write_interval(width * index, code);
}

template <typename T>
T FixedLengthCodeArrayBase<T>::get(uint64_t index) const {
    assert(0 <= index && index < length);
    return decode(codes.read_interval(width * index, width));
}

template <typename T>
T FixedLengthCodeArrayBase<T>::operator[](uint64_t index) const {
    return get(index);
}

template <typename T>
bool FixedLengthCodeArrayBase<T>::operator==(
    const FixedLengthCodeArrayBase &rhs) const {
    if (this->size() != rhs.size()) return false;
    for (int i = 0; i < this->size(); i++) {
        if (this->get(i) != rhs.get(i)) return false;
    }
    return true;
}

template <typename T>
bool FixedLengthCodeArrayBase<T>::operator!=(
    const FixedLengthCodeArrayBase &rhs) const {
    return !(*this == rhs);
}

template <typename T>
uint64_t FixedLengthCodeArrayBase<T>::evaluate_memory_consumption() const {
    return codes.evaluate_memory_consumption();
}

template struct FixedLengthCodeArrayBase<TreeBP>;
template struct FixedLengthCodeArrayBase<std::pair<TreeBP, uint32_t>>;

}  // namespace average_case_optimal_rmq
