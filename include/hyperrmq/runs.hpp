#pragma once

#include <cstdint>
#include <vector>

namespace hyperrmq {

template <typename T>
uint32_t count_increasing_runs(const std::vector<T> &perm);

// It works well if `run_count` is small compared to `perm` size.
// For example, if `perm` size is 1e7, it works for run_count < 1e4.
template <typename T>
void random_roughly_fixed_incresing_runs(std::vector<T> &perm,
                                         uint32_t run_count, uint32_t seed = 0);

// It repeats trying random_roughly_fixed_increasing_runs
// until it gives the demanded run_count.
template <typename T>
void random_exact_fixed_incresing_runs(std::vector<T> &perm, uint32_t run_count,
                                       bool verbose = true);

}  // namespace hyperrmq
