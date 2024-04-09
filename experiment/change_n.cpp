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
void benchmarkHyperRMQ(const vector<int>& perm, int B,
                       const vector<pair<int, int>>& queries, string RMQname) {
    cout << RMQname << endl;

    vector<pair<string, double>> benchmark_result;
    auto start = chrono::system_clock::now();
    auto timer = [&]() -> double {
        auto now = chrono::system_clock::now();
        return static_cast<double>(
            chrono::duration_cast<chrono::microseconds>(now - start).count() /
            1000.0);
    };
    auto time_keeper = [&](string message, bool end) -> void {
        if (end) {
            benchmark_result.push_back({message, timer()});
        } else {
            start = chrono::system_clock::now();
        }
    };

    time_keeper("construction", 0);
    HyperRMQ rmq(perm, B);
    time_keeper("construction", 1);

    time_keeper("query", 0);
    for (auto&& [i, j] : queries) {
        rmq.query(i, j);
    }
    time_keeper("query", 1);

    benchmark_result.push_back(
        {"memory", double(rmq.evaluate_memory_consumption()) / perm.size()});

    for (auto&& [label, value] : benchmark_result) {
        cout << left << setw(20) << (label + ":") << right << setw(20) << value
             << endl;
    }
}

void benchmarkFerrada(const vector<int>& perm,
                      const vector<pair<int, int>>& queries) {
    cout << "FerradaRMQ" << endl;

    const int N = perm.size();

    vector<pair<string, double>> benchmark_result;
    auto start = chrono::system_clock::now();
    auto timer = [&]() -> double {
        auto now = chrono::system_clock::now();
        return static_cast<double>(
            chrono::duration_cast<chrono::microseconds>(now - start).count() /
            1000.0);
    };
    auto time_keeper = [&](string message, bool end) -> void {
        if (end) {
            benchmark_result.push_back({message, timer()});
        } else {
            start = chrono::system_clock::now();
        }
    };

    auto perm_reverse = new int[N];
    copy(perm.rbegin(), perm.rend(), perm_reverse);

    time_keeper("construction", 0);
    RMQRMM64 rmq(perm_reverse, N);
    time_keeper("construction", 1);

    time_keeper("query", 0);
    for (auto&& [i, j] : queries) {
        rmq.queryRMQ(N - 1 - j, N - 1 - i);
    }
    time_keeper("query", 1);

    benchmark_result.push_back({"memory", 8.0 * rmq.getSize() / perm.size()});

    for (auto&& [label, value] : benchmark_result) {
        cout << left << setw(20) << (label + ":") << right << setw(20) << value
             << endl;
    }

    delete[] perm_reverse;
}

void benchmarkSDSLNEW(const vector<int>& perm,
                      const vector<pair<int, int>>& queries) {
    cout << "SDSLNEWRMQ" << endl;

    vector<pair<string, double>> benchmark_result;
    auto start = chrono::system_clock::now();
    auto timer = [&]() -> double {
        auto now = chrono::system_clock::now();
        return static_cast<double>(
            chrono::duration_cast<chrono::microseconds>(now - start).count() /
            1000.0);
    };
    auto time_keeper = [&](string message, bool end) -> void {
        if (end) {
            benchmark_result.push_back({message, timer()});
        } else {
            start = chrono::system_clock::now();
        }
    };

    time_keeper("construction", 0);
    rmq_succinct_rec_new<true, 2048, 1024, 128, 0> rmq(&perm);
    time_keeper("construction", 1);

    time_keeper("query", 0);
    for (auto&& [i, j] : queries) {
        rmq(i, j);
    }
    time_keeper("query", 1);

    benchmark_result.push_back(
        {"memory", 8.0 * size_in_bytes(rmq) / perm.size()});

    for (auto&& [label, value] : benchmark_result) {
        cout << left << setw(20) << (label + ":") << right << setw(20) << value
             << endl;
    }
}

int main() {
    cout << fixed << setprecision(6);
    int Q = 1e6;

    auto timeinfo =
        chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << std::ctime(&timeinfo) << endl;
    cout << "Q = " << Q << endl;

    int B_huf = 3;
    for (int64_t N = 1e4; N <= int64_t(1e9); N *= 10) {
        cout << "N = " << N << endl;
        B_huf++;

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

        benchmarkFerrada(perm, queries);

        benchmarkSDSLNEW(perm, queries);

        benchmarkHyperRMQ<HyperRMQHuffman>(
            perm, B_huf, queries, "HyperRMQHuffman B = " + to_string(B_huf));

        for (int B = 64; B <= 1024; B <<= 1) {
            benchmarkHyperRMQ<HyperRMQBreadthFirstArithmetic>(
                perm, B, queries, "HyperRMQHuffman B = " + to_string(B));
        }
    }

    return 0;
}
