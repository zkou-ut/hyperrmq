#include <bits/stdc++.h>

#include "RMQRMM64.h"
#include "hyperrmq/hyper_rmq.hpp"
#include "hyperrmq/memutil.hpp"
#include "sdsl/rmq_succinct_rec_new.hpp"

using namespace std;
using namespace hyperrmq;
using namespace sdsl;

using HyperRMQHuffman =
    HyperRMQ<16, CompressedMicrotreeSplitRankArrayHuffman<16>>;
using HyperRMQBreadthFirstArithmetic =
    HyperRMQ<16, CompressedMicrotreeSplitRankArrayArithmetic<false>>;

void test_dfs(const vector<int64_t>& values, int64_t Q = 1'000'000,
              int64_t B_huffman = 8, int64_t B_arithmetic = 512) {
    const int64_t N = values.size();

    cerr << " - Start DFS Test" << endl;

    HyperRMQHuffman hyper_rmq_huffman(values, B_huffman);
    cerr << " - Finished constructing Huffman" << endl;

    HyperRMQBreadthFirstArithmetic hyper_rmq_arith(values, B_arithmetic);
    cerr << " - Finished constructing Arithmetic" << endl;

    auto values_array = new int64_t[N];
    // reverse input array and queries to find the leftmost minimum.
    copy(values.rbegin(), values.rend(), values_array);
    RMQRMM64 rmq_rmm(values_array, N);
    delete[] values_array;
    cerr << " - Finished constructing RMM" << endl;

    auto rmq_rmm_query = [&](int64_t l, int64_t r) -> int64_t {
        return N - 1 - rmq_rmm.queryRMQ(N - 1 - r, N - 1 - l);
    };

    rmq_succinct_rec_new<true, 2048, 1024, 128, 0> rmq_new(&values);
    cerr << " - Finished constructing REC" << endl;

    auto check_query = [&](int64_t l, int64_t r) -> int64_t {
        vector<int64_t> answers = {hyper_rmq_huffman.query(l, r),
                                   hyper_rmq_arith.query(l, r),
                                   rmq_rmm_query(l, r), rmq_new(l, r)};

        for (int64_t k = 0; k < answers.size(); k++) {
            if (answers[k] != answers[0]) {
                cout << "Test for rmq failed." << endl;
                if (N <= 20) {
                    cout << "Input array:\n[";
                    for (int64_t i = 0; i < N; i++) {
                        cout << values[i] << (i == N - 1 ? "]" : ", ");
                    }
                    cout << endl;
                }
                cout << "query: " << l << ", " << r << endl;
                cout << "answers:\n[";
                for (int64_t i = 0; i < answers.size(); i++) {
                    cout << answers[i]
                         << (i == answers.size() - 1 ? "]" : ", ");
                }
                cout << endl;
                exit(-1);
            }
        }

        return answers[0];
    };

    cerr << " - Query info:" << endl;
    stack<pair<int64_t, int64_t>> queries;
    queries.push({0, N - 1});
    for (int64_t i = 0; i < Q && !queries.empty(); i++) {
        auto [l, r] = queries.top();
        queries.pop();

        if ((i & -i) == i) {
            cerr << right << setw(12) << i << ": [" << l << ", " << r << "]"
                 << endl;
        }

        auto m = check_query(l, r);
        if (l < m - 1) queries.push({l, m - 1});
        if (m + 1 < r) queries.push({m + 1, r});
    }
}

int main(int argc, char const* argv[]) {
    if (argc < 2) {
        cout << "Specify the number of elements N as an argument" << endl;
        cout << "usage : ./build/overflow-test 2200000000" << endl;
        return 0;
    }
    const int64_t N = atol(argv[1]);

    const int64_t Q = 1'000'000;
    mt19937 engine(0);

    cout << "Overflow test for rmq data structures..." << endl;

    cout << "[Permutation start]" << endl;
    {
        vector<int64_t> perm(N);
        iota(perm.begin(), perm.end(), 0L);
        shuffle(perm.begin(), perm.end(), engine);

        test_dfs(perm, Q);
    }
    cout << "[Permutation OK]" << endl;

    cout << "[Uniformly random start]" << endl;
    {
        vector<int64_t> values(N);
        for (int64_t i = 0; i < N; i++) {
            values[i] = engine() % 10000;
        }

        test_dfs(values, Q);
    }
    cout << "[Uniformly random OK]" << endl;

    cout << "All the tests have successfully finished." << endl;

    return 0;
}
