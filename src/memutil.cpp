#include "hyperrmq/memutil.hpp"

#include <climits>
#include <string>

namespace average_case_optimal_rmq {

template <typename T>
uint64_t evaluate_vector_memory_consumption(const std::vector<T>& vec) {
    return CHAR_BIT * sizeof(T) * vec.capacity();
}

template uint64_t evaluate_vector_memory_consumption<uint8_t>(
    const std::vector<uint8_t>&);
template uint64_t evaluate_vector_memory_consumption<uint16_t>(
    const std::vector<uint16_t>&);
template uint64_t evaluate_vector_memory_consumption<uint32_t>(
    const std::vector<uint32_t>&);
template uint64_t evaluate_vector_memory_consumption<uint64_t>(
    const std::vector<uint64_t>&);

template uint64_t evaluate_vector_memory_consumption<int8_t>(
    const std::vector<int8_t>&);
template uint64_t evaluate_vector_memory_consumption<int16_t>(
    const std::vector<int16_t>&);
template uint64_t evaluate_vector_memory_consumption<int32_t>(
    const std::vector<int32_t>&);
template uint64_t evaluate_vector_memory_consumption<int64_t>(
    const std::vector<int64_t>&);

std::vector<std::pair<std::string, uint64_t>> memory_table_combine_children(
    std::string parent_name,
    const std::vector<std::vector<std::pair<std::string, uint64_t>>>&
        children) {
    std::vector<std::pair<std::string, uint64_t>> result = {{parent_name, 0}};

    for (auto&& child : children) {
        for (auto&& [key, val] : child) {
            result.push_back({memory_table_tab + key, val});
        }
        result.front().second += child.front().second;
    }

    return result;
}

}  // namespace average_case_optimal_rmq
