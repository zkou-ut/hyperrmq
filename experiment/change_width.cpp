#include <bits/stdc++.h>

#include "RMQRMM64.h"
#include "hyperrmq/hyper_rmq.hpp"
#include "sdsl/memory_management.hpp"
#include "sdsl/rmq_support.hpp"

using namespace std;
using namespace average_case_optimal_rmq;
using namespace sdsl;

using HyperRMQHuffman =
    HyperRMQ<16, CompressedMicrotreeSplitRankArrayHuffman<16>>;
using HyperRMQBreadthFirstArithmetic =
    HyperRMQ<16, CompressedMicrotreeSplitRankArrayArithmetic<false>>;

template <typename HyperRMQ>
vector<double> benchmarkHyperRMQ(const vector<int>& perm, int B,
                                 const vector<vector<pair<int, int>>>& queries,
                                 string RMQname) {
    const int N = perm.size();

    cout << RMQname << endl;

    mt19937 engine(0);

    HyperRMQ rmq(perm, B);

    cout << rmq.evaluate_memory_consumption() << endl;

    vector<double> time_result;

    for (auto&& queries_row : queries) {
        auto start = chrono::system_clock::now();
        for (auto&& [l, r] : queries_row) {
            rmq.query(l, r);
        }
        auto end = chrono::system_clock::now();
        double time =
            static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(end - start)
                    .count()) /
            queries_row.size();
        time_result.push_back(time);
    }

    return time_result;
}

vector<double> benchmarkFerrada(const vector<int>& perm,
                                const vector<vector<pair<int, int>>>& queries) {
    cout << "FerradaRMQ" << endl;
    const int N = perm.size();

    mt19937 engine(0);

    auto perm_reverse = new int[N];
    copy(perm.rbegin(), perm.rend(), perm_reverse);
    RMQRMM64 rmq(perm_reverse, N);

    cout << 8 * rmq.getSize() << endl;

    vector<double> time_result;

    for (auto&& queries_row : queries) {
        auto start = chrono::system_clock::now();
        for (auto&& [l, r] : queries_row) {
            rmq.queryRMQ(N - 1 - r, N - 1 - l);
        }
        auto end = chrono::system_clock::now();
        double time =
            static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(end - start)
                    .count()) /
            queries_row.size();
        time_result.push_back(time);
    }

    delete[] perm_reverse;

    return time_result;
}

vector<double> benchmarkSDSLNEW(const vector<int>& perm,
                                const vector<vector<pair<int, int>>>& queries) {
    cout << "SDSLNEWRMQ" << endl;

    const int N = perm.size();

    mt19937 engine(0);

    rmq_succinct_rec_new<true, 2048, 1024, 128, 0> rmq(&perm);

    cout << 8 * size_in_bytes(rmq) << endl;

    vector<double> time_result;

    for (auto&& queries_row : queries) {
        auto start = chrono::system_clock::now();
        for (auto&& [l, r] : queries_row) {
            rmq(l, r);
        }
        auto end = chrono::system_clock::now();
        double time =
            static_cast<double>(
                chrono::duration_cast<chrono::microseconds>(end - start)
                    .count()) /
            queries_row.size();
        time_result.push_back(time);
    }

    return time_result;
}

int main() {
    cout << fixed << setprecision(6);

    constexpr int N = 1e9;
    constexpr int Q = 1e5;

    auto timeinfo =
        chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << std::ctime(&timeinfo) << endl;
    cout << "N = " << N << endl;
    cout << "Q = " << Q << endl;

    mt19937 engine(0);
    vector<int32_t> perm(N);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    vector<vector<pair<int, int>>> queries;
    vector<int> ws;
    for (int w = 2; w <= N; w <<= 1) {
        ws.push_back(w);
        queries.push_back(vector<pair<int, int>>(Q));
        for (int i = 0; i < Q; i++) {
            uint32_t l = engine() % (N - w + 1);
            uint32_t r = l + w - 1;
            queries.back()[i] = {l, r};
        }
    }

    auto output = [&](const vector<double>& time_result) -> void {
        assert(time_result.size() == ws.size());
        for (int i = 0; i < ws.size(); i++) {
            cout << ws[i] << " " << time_result[i] << endl;
        }
    };

    output(benchmarkFerrada(perm, queries));

    output(benchmarkSDSLNEW(perm, queries));

    int B_huf = 8;
    output(benchmarkHyperRMQ<
           HyperRMQ<16, CompressedMicrotreeSplitRankArrayHuffman<16>>>(
        perm, B_huf, queries, "HyperRMQHuffman " + to_string(B_huf)));

    for (int B = 64; B <= 1024; B <<= 1) {
        output(benchmarkHyperRMQ<HyperRMQ<
                   16, CompressedMicrotreeSplitRankArrayArithmetic<false>>>(
            perm, B, queries, "HyperRMQBreadth " + to_string(B)));
    }

    int B_depth = 512;
    output(benchmarkHyperRMQ<
           HyperRMQ<16, CompressedMicrotreeSplitRankArrayArithmetic<true>>>(
        perm, B_depth, queries, "HyperRMQDepth " + to_string(B_depth)));

    return 0;
}
