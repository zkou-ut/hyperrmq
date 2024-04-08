#include "runs.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <set>

namespace average_case_optimal_rmq {

template <typename T>
uint32_t count_increasing_runs(const std::vector<T>& perm) {
    if (perm.size() == 0) return 0;
    T answer = 1;
    for (int i = 0; i < perm.size() - 1; i++) {
        if (perm[i] > perm[i + 1]) {
            answer++;
        }
    }
    return answer;
}

template <typename T>
void random_roughly_fixed_incresing_runs(std::vector<T>& perm,
                                         uint32_t run_count, uint32_t seed) {
    assert(1 <= run_count && run_count <= perm.size());

    std::mt19937_64 engine(seed);

    std::shuffle(perm.begin(), perm.end(), engine);

    std::set<uint64_t> split_indices = {0, perm.size()};
    while (split_indices.size() < run_count + 1) {
        uint64_t index = engine() % (perm.size() - 1) + 1;
        split_indices.insert(index);
    }
    auto it = split_indices.begin();
    while (true) {
        auto start = *it;
        it++;
        if (it == split_indices.end()) break;
        auto end = *it;
        std::sort(perm.begin() + start, perm.begin() + end);
    }
}

template <typename T>
void random_exact_fixed_incresing_runs(std::vector<T>& perm, uint32_t run_count,
                                       bool verbose) {
    if (verbose && uint64_t(run_count) * run_count > perm.size()) {
        std::cerr << "Warning: Large run_count " << run_count
                  << " compared to permutation size " << perm.size() << ".\n";
        std::cerr << "It may be time consuming." << std::endl;
    }

    for (uint32_t seed = 0;; seed++) {
        random_roughly_fixed_incresing_runs(perm, run_count, seed);
        if (count_increasing_runs(perm) == run_count) {
            if (verbose) {
                std::cerr << "Found, seed: " << seed << std::endl;
            }
            return;
        }
    }
}

template uint32_t count_increasing_runs(const std::vector<int32_t>&);
template uint32_t count_increasing_runs(const std::vector<uint32_t>&);

template void random_roughly_fixed_incresing_runs(std::vector<int32_t>&,
                                                  uint32_t, uint32_t);
template void random_roughly_fixed_incresing_runs(std::vector<uint32_t>&,
                                                  uint32_t, uint32_t);

template void random_exact_fixed_incresing_runs(std::vector<int32_t>&, uint32_t,
                                                bool);
template void random_exact_fixed_incresing_runs(std::vector<uint32_t>&,
                                                uint32_t, bool);

}  // namespace average_case_optimal_rmq
