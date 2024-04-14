#include "hyperrmq/tree_bp.hpp"

#include <cassert>
#include <stack>

#include "hyperrmq/table.hpp"

namespace hyperrmq {

TreeBP::TreeBP() : n(0), bp() {}

TreeBP::TreeBP(uint64_t n, const BitArray &bp) : n(n), bp(bp) {
    assert(bp.size() == 2 * n);
    assert(bp.linear_popcount() == n);
}

TreeBP::TreeBP(uint64_t n, BitArray &&bp) : n(n), bp(bp) {
    assert(bp.size() == 2 * n);
    assert(bp.linear_popcount() == n);
}

TreeBP::TreeBP(const std::string &bp_str) : n(bp_str.size() / 2), bp(bp_str) {
    assert(bp.size() == 2 * n);
    assert(bp.linear_popcount() == n);
}

bool TreeBP::operator==(const TreeBP &rhs) const {
    return this->n == rhs.n && this->bp == rhs.bp;
}
bool TreeBP::operator!=(const TreeBP &rhs) const {
    return this->n != rhs.n || this->bp != rhs.bp;
}

bool TreeBP::operator<(const TreeBP &rhs) const {
    return (this->n < rhs.n) || (this->n == rhs.n && this->bp < rhs.bp);
}

std::string TreeBP::to_string() const { return bp.to_string("(", ")"); }

uint64_t TreeBP::evaluate_memory_consumption() const {
    return bp.evaluate_memory_consumption();
}

uint32_t TreeBP::naive_fwdsearch(uint32_t j, int32_t d) const {
    constexpr int c = 16;
    constexpr int mask = (1 << c) - 1;
    constexpr int part_count = 64 / c;

    int32_t d_prime = 0;

    int s = j / c + 1;
    int t = (bp.size() + c - 1) / c;

    if (s * c > bp.size()) {
        for (int i = j; i < bp.size(); i++) {
            d_prime += 1 - 2 * bp.get(i);
            if (d_prime == d) {
                return i + 1;
            }
        }
        return bp.size() + 1;
    }

    for (int i = j; i < s * c; i++) {
        d_prime += 1 - 2 * bp.get(i);
        if (d_prime == d) {
            return i + 1;
        }
    }

    for (int i = s; i < t; i++) {
        uint64_t cell = bp.get_raw_cell(i / part_count);
        uint64_t cell_part = (cell >> (64 - c * (i % part_count + 1))) & mask;
        if (d_prime + min_table<c>[cell_part] <= d) {
            for (int k = c - 1; k >= 0; k--) {
                d_prime += 1 - 2 * ((cell_part >> k) & 1);
                if (d_prime == d) {
                    return i * c + (c - k);
                }
            }
        }
        d_prime += excess_table<c>[cell_part];
    }

    return bp.size() + 1;
}

int32_t TreeBP::naive_bwdsearch(uint32_t j, int32_t d) const {
    constexpr int c = 16;
    constexpr int mask = (1 << c) - 1;
    constexpr int part_count = 64 / c;

    int32_t d_prime = 0;

    int s = (j - 1) / c;

    for (int i = j - 1; i >= s * c; i--) {
        d_prime -= 1 - 2 * bp.get(i);
        if (d_prime == d) {
            return i;
        }
    }

    for (int i = s - 1; i >= 0; i--) {
        uint64_t cell = bp.get_raw_cell(i / part_count);
        uint64_t cell_part = (cell >> (64 - c * (i % part_count + 1))) & mask;
        if (d_prime - excess_table<c>[cell_part] + min_table<c>[cell_part] <=
            d) {
            for (int k = 0; k < c; k++) {
                d_prime -= 1 - 2 * ((cell_part >> k) & 1);
                if (d_prime == d) {
                    return i * c + (c - k) - 1;
                }
            }
        }
        d_prime -= excess_table<c>[cell_part];
    }

    return -1;
}

int32_t TreeBP::naive_minexcess(uint32_t j, int32_t k) const {
    constexpr int c = 16;
    constexpr int mask = (1 << c) - 1;
    constexpr int part_count = 64 / c;

    int32_t d = 0, m = 0;

    int s = j / c + 1;
    int t = k / c;

    if (s > t) {
        for (int i = j; i < k; i++) {
            d += 1 - 2 * bp.get(i);
            m = std::min(m, d);
        }
        return m;
    }

    for (int i = j; i < s * c; i++) {
        d += 1 - 2 * bp.get(i);
        m = std::min(m, d);
    }

    for (int i = s; i < t; i++) {
        uint64_t cell = bp.get_raw_cell(i / part_count);
        uint64_t cell_part = (cell >> (64 - c * (i % part_count + 1))) & mask;
        m = std::min(m, d + min_table<c>[cell_part]);
        d += excess_table<c>[cell_part];
    }

    for (int i = t * c; i < k; i++) {
        d += 1 - 2 * bp.get(i);
        m = std::min(m, d);
    }

    return m;
}

uint32_t TreeBP::naive_open(uint32_t index) const {
    assert(0 <= index && index < bp.size());
    assert(bp.get(index) == 1);
    return naive_bwdsearch(index, -1);
}

uint32_t TreeBP::naive_close(uint32_t index) const {
    assert(0 <= index && index < bp.size());
    assert(bp.get(index) == 0);
    return naive_fwdsearch(index + 1, -1) - 1;
}

uint32_t TreeBP::naive_lca(uint32_t u_inorder, uint32_t v_inorder) const {
    assert(0 <= u_inorder && u_inorder < n);
    assert(0 <= v_inorder && v_inorder < n);
    if (u_inorder > v_inorder) {
        std::swap(u_inorder, v_inorder);
    }
    auto i = bp.linear_select1(u_inorder);
    auto j = bp.linear_select1(v_inorder);
    int32_t m = naive_minexcess(i, j + 1);
    if (m == 0) {
        return u_inorder;
    }
    return bp.linear_popcount(0, naive_fwdsearch(i, m) - 1);
}

TreeBP cartesian_tree_bp(const std::vector<int32_t> &values) {
    uint64_t n = values.size();
    BitArray bp(2 * n);

    std::stack<uint64_t> st;
    uint64_t pos = 2 * n;
    for (int64_t i = n - 1; i >= 0; i--) {
        while (!st.empty() && values[st.top()] >= values[i]) {
            st.pop();
            bp.off(--pos);
        }
        st.push(i);
        bp.on(--pos);
    }
    while (!st.empty()) {
        st.pop();
        bp.off(--pos);
    }
    assert(pos == 0);

    return TreeBP(n, std::move(bp));
}

}  // namespace hyperrmq
