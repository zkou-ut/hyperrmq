#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace hyperrmq {

const std::string memory_table_tab("   ", 3);

std::vector<std::pair<std::string, uint64_t>> memory_table_combine_children(
    std::string parent_name,
    const std::vector<std::vector<std::pair<std::string, uint64_t>>> &children);

template <typename T>
uint64_t evaluate_vector_memory_consumption(const std::vector<T> &vec);

}  // namespace hyperrmq
