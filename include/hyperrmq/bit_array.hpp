#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace hyperrmq {

struct BitStack {
    uint64_t num_of_bits;
    std::vector<uint64_t> cells;

    BitStack();

    void push(uint64_t b);
};

// It packs bits from the top.
struct BitArray {
   public:
    BitArray();
    explicit BitArray(uint64_t length_in_bits);
    explicit BitArray(uint64_t num_of_bits, const std::vector<uint64_t>& cells);
    explicit BitArray(uint64_t num_of_bits, std::vector<uint64_t>&& cells);
    explicit BitArray(const std::string& bp);
    BitArray(BitStack&& bs);

    size_t size() const;

    void on(uint64_t i);
    void off(uint64_t i);

    void set(uint64_t i, bool v);
    uint64_t get(uint64_t i) const;

    uint64_t get_raw_cell(uint64_t cell_index) const;

    void write_bits(uint64_t first, uint64_t width, uint64_t bits);
    uint64_t read_bits(uint64_t first, uint64_t width) const;
    uint64_t read_bits_zero_follow(uint64_t first, uint64_t width) const;

    void write_interval(uint64_t first, const BitArray& interval);
    BitArray read_interval(uint64_t first, uint64_t width) const;

    uint64_t linear_popcount() const;
    uint64_t linear_popcount(uint64_t first, uint64_t width) const;

    uint64_t linear_select1(uint64_t j) const;

    std::string to_string(std::string zero = "0", std::string one = "1") const;

    bool operator==(const BitArray& rhs) const;
    bool operator!=(const BitArray& rhs) const;

    bool operator<(const BitArray& rhs) const;

    uint64_t evaluate_memory_consumption() const;

   private:
    uint64_t num_of_bits;
    uint64_t cell_count;
    std::vector<uint64_t> cells;
};

std::ostream& operator<<(std::ostream& stream, const BitArray& bit_array);

}  // namespace hyperrmq
