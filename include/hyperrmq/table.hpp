#pragma once

#include <array>
#include <cstdint>

namespace hyperrmq {

template <uint32_t c>
static std::pair<std::array<int32_t, (1 << c)>, std::array<int32_t, (1 << c)>>
construct_tables() {
    std::array<int32_t, 1 << c> excess_table{}, min_table{};

    for (int i = 0; i < (1 << c); i++) {
        excess_table[i] = 0;
        min_table[i] = 0;
        for (int p = c - 1; p >= 0; p--) {
            excess_table[i] += 1 - 2 * ((i >> p) & 1);
            if (min_table[i] > excess_table[i]) {
                min_table[i] = excess_table[i];
            }
        }
    }

    return {excess_table, min_table};
}

template <uint32_t c>
static const std::pair<std::array<int32_t, (1 << c)>,
                       std::array<int32_t, (1 << c)>>
    tables = construct_tables<c>();

template <uint32_t c>
const std::array<int32_t, (1 << c)> excess_table = tables<c>.first;

template <uint32_t c>
const std::array<int32_t, (1 << c)> min_table = tables<c>.second;

}  // namespace hyperrmq
