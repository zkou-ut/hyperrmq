#pragma once

#include "bit_array.hpp"

namespace average_case_optimal_rmq {

// T: type of encoded object
template <typename T>
struct FixedLengthCodeArrayBase {
   public:
    FixedLengthCodeArrayBase();
    FixedLengthCodeArrayBase(uint32_t width, uint64_t length);
    FixedLengthCodeArrayBase(uint32_t width, std::vector<T> vec);

    // It returns a BitArray of width W.
    virtual BitArray encode(const T& object) const = 0;

    // It takes a BitArray of width W.
    virtual T decode(const BitArray& code) const = 0;

    size_t size() const;

    void set(uint64_t index, T object);

    T get(uint64_t index) const;
    T operator[](uint64_t index) const;

    bool operator==(const FixedLengthCodeArrayBase& rhs) const;
    bool operator!=(const FixedLengthCodeArrayBase& rhs) const;

    uint64_t evaluate_memory_consumption() const;

    uint32_t width;
    uint64_t length;
    BitArray codes;
};

}  // namespace average_case_optimal_rmq
