#include "hyperrmq/bitutil.hpp"

namespace average_case_optimal_rmq {

uint64_t popcount(uint64_t b) {
    b -= (b >> 1) & 0x5555555555555555;
    b = (b & 0x3333333333333333) + ((b >> 2) & 0x3333333333333333);
    return (((b + (b >> 4)) & 0x0F0F0F0F0F0F0F0F) * 0x0101010101010101) >> 56;
}

uint64_t ceil_log2(uint64_t b) {
    uint64_t lg = 0;
    while (lg < 64 && (1ull << lg) < b) {
        lg++;
    }
    return lg;
}

}  // namespace average_case_optimal_rmq
