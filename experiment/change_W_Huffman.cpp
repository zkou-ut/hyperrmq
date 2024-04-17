#include <bits/stdc++.h>

#include "hyperrmq/hyper_rmq.hpp"

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

template <int i, int j>
struct bm {
    bm(const vector<int>& perm, int B, const vector<pair<int, int>>& queries) {
        constexpr int si = 1 << i;
        constexpr int sj = 1 << j;
        cout << si << " " << sj << endl;
        benchmark<HyperRMQ<si, CompressedMicrotreeSplitRankArrayHuffman<sj>>>(
            perm, B, queries);
    }
};

int main() {
    cout << fixed << setprecision(5);

    constexpr int N = 1e9, Q = 1e6;

    auto timeinfo =
        chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << std::ctime(&timeinfo) << endl;
    cout << "N = " << N << endl;
    cout << "Q = " << Q << endl;

    mt19937 engine(0);
    vector<int32_t> perm(N);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    vector<pair<int, int>> queries(Q);
    for (int i = 0; i < Q; i++) {
        queries[i] = {engine() % N, engine() % N};
        if (queries[i].first > queries[i].second) {
            swap(queries[i].first, queries[i].second);
        }
    }

    const int B = 8;

    bm<0, 0>(perm, B, queries);
    bm<1, 1>(perm, B, queries);
    bm<2, 2>(perm, B, queries);
    bm<3, 3>(perm, B, queries);
    bm<4, 4>(perm, B, queries);
    bm<5, 5>(perm, B, queries);
    bm<6, 6>(perm, B, queries);
    bm<7, 7>(perm, B, queries);
    bm<8, 8>(perm, B, queries);

    return 0;
}
