#include "hyperrmq/bit_array.hpp"

#include <cassert>
#include <climits>

#include "hyperrmq/bitutil.hpp"
#include "hyperrmq/memutil.hpp"

namespace hyperrmq {

BitStack::BitStack() : num_of_bits(0), cells() {}

void BitStack::push(uint64_t b) {
    assert(0 <= b && b <= 1);
    uint64_t write_pos = num_of_bits;
    num_of_bits++;
    if (((num_of_bits + 63) >> 6) > cells.size()) {
        cells.push_back(0);
    }
    cells[write_pos >> 6] |= b << (~write_pos & 0x3f);
}

BitArray::BitArray() : num_of_bits(0), cell_count(0), cells() {}

BitArray::BitArray(uint64_t length_in_bits)
    : num_of_bits(length_in_bits),
      cell_count((num_of_bits + 63) >> 6),
      cells(cell_count) {}

BitArray::BitArray(uint64_t num_of_bits, const std::vector<uint64_t>& cells)
    : num_of_bits(num_of_bits),
      cell_count((num_of_bits + 63) >> 6),
      cells(cells) {
    assert(cells.size() == cell_count);
    assert(cells.back() << (num_of_bits & 0x3f) == 0);
}

BitArray::BitArray(uint64_t num_of_bits, std::vector<uint64_t>&& cells)
    : num_of_bits(num_of_bits),
      cell_count((num_of_bits + 63) >> 6),
      cells(cells) {
    assert(cells.size() == cell_count);
    assert(cells.back() << (num_of_bits & 0x3f) == 0);
}

BitArray::BitArray(const std::string& bp)
    : num_of_bits(bp.size()),
      cell_count((num_of_bits + 63) >> 6),
      cells(cell_count) {
    for (int i = 0; i < num_of_bits; i++) {
        assert(bp[i] == '(' || bp[i] == ')');
        set(i, (bp[i] == ')'));
    }
}

BitArray::BitArray(BitStack&& bs)
    : num_of_bits(bs.num_of_bits), cell_count((num_of_bits + 63) >> 6) {
    bs.cells.shrink_to_fit();
    cells = std::move(bs.cells);
}

size_t BitArray::size() const { return num_of_bits; }

void BitArray::on(uint64_t i) {
    assert(0 <= i && i < num_of_bits);
    cells[i >> 6] |= (1ull << (~i & 0x3f));
}

void BitArray::off(uint64_t i) {
    assert(0 <= i && i < num_of_bits);
    cells[i >> 6] &= ~(1ull << (~i & 0x3f));
}

void BitArray::set(uint64_t i, bool v) {
    assert(0 <= i && i < num_of_bits);
    cells[i >> 6] &= ~(1ull << (~i & 0x3f));
    cells[i >> 6] |= (uint64_t(v) << (~i & 0x3f));
}

uint64_t BitArray::get(uint64_t i) const {
    assert(0 <= i && i < num_of_bits);
    return (cells[i >> 6] >> (~i & 0x3f)) & 1;
}

uint64_t BitArray::get_raw_cell(uint64_t cell_index) const {
    assert(0 <= cell_index && cell_index < cell_count);
    return cells[cell_index];
}

void BitArray::write_bits(uint64_t first, uint64_t width, uint64_t bits) {
    assert(0 <= width && width <= 64);
    assert(0 <= first && first + width <= num_of_bits);
    assert(width == 64 || (bits < (1ull << width)));
    if (width == 0) {
        return;
    }

    uint64_t last = first + width - 1;

    uint64_t fb = first >> 6, fi = (~first & 0x3f), top_mask = ((-2ull) << fi);
    uint64_t lb = last >> 6, li = (~last & 0x3f),
             bottom_mask = ((1ull << li) - 1ull);
    if (fb == lb) {
        cells[fb] &= top_mask | bottom_mask;
        cells[fb] |= bits << li;
    } else {
        cells[fb] &= top_mask;
        cells[fb] |= bits >> (64 - li);
        cells[lb] &= bottom_mask;
        cells[lb] |= bits << li;
    }
}

uint64_t BitArray::read_bits(uint64_t first, uint64_t width) const {
    assert(0 <= width && width <= 64);
    assert(0 <= first && first + width <= num_of_bits);
    if (width == 0) {
        return 0;
    }

    uint64_t last = first + width - 1;

    uint64_t fb = first >> 6, fi = (~first & 0x3f),
             top_mask = (2ull << fi) - 1ull;
    uint64_t lb = last >> 6, li = (~last & 0x3f), bottom_mask = (-1ull) << li;
    if (fb == lb) {
        return (cells[fb] & top_mask & bottom_mask) >> li;
    } else {
        return ((cells[fb] & top_mask) << (64 - li)) |
               ((cells[lb] & bottom_mask) >> li);
    }
}

uint64_t BitArray::read_bits_zero_follow(uint64_t first, uint64_t width) const {
    assert(0 <= width && width <= 64);
    assert(0 <= first && first <= num_of_bits);
    if (width == 0 || first == num_of_bits) {
        return 0;
    }

    uint64_t last = first + width - 1;

    uint64_t fb = first >> 6, fi = (~first & 0x3f),
             top_mask = (2ull << fi) - 1ull;
    uint64_t lb = last >> 6, li = (~last & 0x3f), bottom_mask = (-1ull) << li;
    if (fb == lb) {
        return (cells[fb] & top_mask & bottom_mask) >> li;
    } else {
        uint64_t bottom = (lb < cell_count) ? (cells[lb] & bottom_mask) : 0;
        return ((cells[fb] & top_mask) << (64 - li)) | ((bottom) >> li);
    }
}

void BitArray::write_interval(uint64_t first, const BitArray& interval) {
    assert(0 <= first && first + interval.num_of_bits <= num_of_bits);
    uint64_t interval_seek = 0;
    while (interval_seek < interval.num_of_bits) {
        uint64_t copy_width =
            std::min(64ul, interval.num_of_bits - interval_seek);
        this->write_bits(first + interval_seek, copy_width,
                         interval.read_bits(interval_seek, copy_width));
        interval_seek += copy_width;
    }
}

BitArray BitArray::read_interval(uint64_t first, uint64_t width) const {
    assert(0 <= first && first + width <= num_of_bits);
    BitArray interval(width);
    uint64_t interval_seek = 0;
    while (interval_seek < interval.num_of_bits) {
        uint64_t copy_width =
            std::min(64ul, interval.num_of_bits - interval_seek);
        interval.write_bits(interval_seek, copy_width,
                            this->read_bits(first + interval_seek, copy_width));
        interval_seek += copy_width;
    }
    return interval;
}

uint64_t BitArray::linear_popcount() const {
    uint64_t res = 0;
    for (int i = 0; i < cell_count; i++) {
        res += popcount(cells[i]);
    }
    return res;
}

uint64_t BitArray::linear_popcount(uint64_t first, uint64_t width) const {
    assert(0 <= first && first + width <= num_of_bits);
    if (width == 0) {
        return 0;
    }

    uint64_t last = first + width - 1;

    uint64_t fb = first >> 6, fi = (~first & 0x3f),
             top_mask = (2ull << fi) - 1ull;
    uint64_t lb = last >> 6, li = (~last & 0x3f), bottom_mask = (-1ull) << li;
    if (fb == lb) {
        return popcount((cells[fb] & top_mask & bottom_mask) >> li);
    } else {
        uint64_t res =
            popcount(cells[fb] & top_mask) + popcount(cells[lb] & bottom_mask);
        for (int b = fb + 1; b < lb; b++) {
            res += popcount(cells[b]);
        }
        return res;
    }
}

uint64_t BitArray::linear_select1(uint64_t j) const {
    uint64_t r = 0;
    uint64_t cell_idx = 0;
    while (cell_idx < cell_count) {
        uint64_t p = popcount(cells[cell_idx]);
        if (r + p > j) {
            break;
        }
        r += p;
        cell_idx++;
    }

    if (cell_idx == cell_count) {
        return num_of_bits;
    }

    uint64_t idx = 0;
    while (r <= j) {
        r += (cells[cell_idx] >> (63 - idx)) & 1;
        idx++;
    }
    return (cell_idx << 6) + idx - 1;
}

std::string BitArray::to_string(std::string zero, std::string one) const {
    std::string res;
    for (int i = 0; i < num_of_bits; i++) {
        res += get(i) ? one : zero;
    }
    return res;
}

bool BitArray::operator==(const BitArray& rhs) const {
    return this->num_of_bits == rhs.num_of_bits && this->cells == rhs.cells;
}

bool BitArray::operator!=(const BitArray& rhs) const {
    return this->num_of_bits != rhs.num_of_bits || this->cells != rhs.cells;
}

bool BitArray::operator<(const BitArray& rhs) const {
    return (this->num_of_bits < rhs.num_of_bits) ||
           (this->num_of_bits == rhs.num_of_bits && this->cells < rhs.cells);
}

uint64_t BitArray::evaluate_memory_consumption() const {
    return evaluate_vector_memory_consumption(cells);
}

std::ostream& operator<<(std::ostream& stream, const BitArray& bit_array) {
    return stream << bit_array.to_string();
}

}  // namespace hyperrmq
