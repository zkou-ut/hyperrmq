#include <bits/stdc++.h>

#include "hyperrmq/hyper_rmq.hpp"
#include "hyperrmq/runs.hpp"

using namespace std;
using namespace hyperrmq;

void output_benchmark(
    pair<vector<pair<string, uint64_t>>, vector<pair<string, int>>> bm_result) {
    auto [memory_result, time_result] = bm_result;
    cout << "memory:" << endl;
    for (auto&& [label, val] : memory_result) {
        cout << left << setw(25) << (label + ":") << right << setw(20) << val
             << endl;
    }

    cout << "time:" << endl;
    for (auto&& [label, val] : time_result) {
        cout << setw(20) << label << ":" << setw(20) << val << " ms" << endl;
    }
}

template <typename RMQ>
void benchmark(const vector<int>& perm, int B,
               const vector<pair<int, int>>& queries) {
    cout << typeid(RMQ).name() << endl;
    cout << "B = " << B << endl;

    vector<pair<string, int>> time_result;
    auto start = chrono::system_clock::now();
    auto timer = [&]() -> int {
        auto now = chrono::system_clock::now();
        return static_cast<double>(
            chrono::duration_cast<chrono::microseconds>(now - start).count() /
            1000.0);
    };
    auto time_keeper = [&](string message, bool end) -> void {
        if (end) {
            time_result.push_back({message, timer()});
        } else {
            start = chrono::system_clock::now();
        }
    };

    time_keeper("construction", 0);
    RMQ rmq(perm, B);
    time_keeper("construction", 1);

    time_keeper("query", 0);
    for (auto&& [i, j] : queries) {
        rmq.query(i, j);
    }
    time_keeper("query", 1);

    vector<pair<string, uint64_t>> memory_result = rmq.memory_table();

    output_benchmark({memory_result, time_result});
};

int main(int argc, char const* argv[]) {
    constexpr int Q = 1e6;

    if (argc < 3) {
        cout << "Specify the number of elements N and the expected number of "
                "increasing runs r as arguments"
             << endl;
        cout << "usage : ./build/experiment/increasing_runs 10000 100" << endl;
        return 0;
    }

    int N = atoi(argv[1]);
    int r = atoi(argv[2]);

    auto timeinfo =
        chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << std::ctime(&timeinfo) << endl;
    cout << "N = " << N << endl;
    cout << "Q = " << Q << endl;
    cout << "r(expected) = " << r << endl;

    mt19937 engine(0);
    vector<int32_t> perm(N);
    iota(perm.begin(), perm.end(), 0);

    random_roughly_fixed_incresing_runs(perm, r, 0);
    cout << "r = " << count_increasing_runs(perm) << endl;

    vector<pair<int, int>> queries(Q);
    for (int i = 0; i < Q; i++) {
        queries[i] = {engine() % N, engine() % N};
        if (queries[i].first > queries[i].second) {
            swap(queries[i].first, queries[i].second);
        }
    }

    for (int B = 5; B <= 100; B += 5) {
        benchmark<HyperRMQ<16>>(perm, B, queries);
    }

    return 0;
}
