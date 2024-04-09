#include <bits/stdc++.h>

#include "RMQRMM64.h"
#include "hyperrmq/hyper_rmq.hpp"
#include "hyperrmq/memutil.hpp"
#include "sdsl/rmq_succinct_rec_new.hpp"

using namespace std;
using namespace average_case_optimal_rmq;
using namespace sdsl;

using HyperRMQHuffman =
    HyperRMQ<16, CompressedMicrotreeSplitRankArrayHuffman<16>>;
using HyperRMQBreadthFirstArithmetic =
    HyperRMQ<16, CompressedMicrotreeSplitRankArrayArithmetic<false>>;

void test_exhaustive(const vector<int>& values, int B_huffman = 4,
                     int B_arithmetic = 64) {
    const int N = values.size();

    HyperRMQHuffman hyper_rmq_huffman(values, B_huffman);
    HyperRMQBreadthFirstArithmetic hyper_rmq_arith(values, B_arithmetic);

    auto values_array = new int[N];
    // reverse input array and queries to find the leftmost minimum.
    copy(values.rbegin(), values.rend(), values_array);
    RMQRMM64 rmq_rmm(values_array, N);
    auto rmq_rmm_query = [&](int l, int r) -> int {
        return N - 1 - rmq_rmm.queryRMQ(N - 1 - r, N - 1 - l);
    };

    rmq_succinct_rec_new<true, 2048, 1024, 128, 0> rmq_new(&values);

    auto check_query = [&](int l, int r) {
        vector<int> answers = {hyper_rmq_huffman.query(l, r),
                               hyper_rmq_arith.query(l, r), rmq_rmm_query(l, r),
                               rmq_new(l, r)};

        for (int k = 0; k < answers.size(); k++) {
            if (answers[k] != answers[0]) {
                cout << "Test for rmq failed." << endl;
                if (N <= 20) {
                    cout << "Input array:\n[";
                    for (int i = 0; i < N; i++) {
                        cout << values[i] << (i == N - 1 ? "]" : ", ");
                    }
                    cout << endl;
                }
                cout << "query: " << l << ", " << r << endl;
                cout << "answers:\n[";
                for (int i = 0; i < answers.size(); i++) {
                    cout << answers[i]
                         << (i == answers.size() - 1 ? "]" : ", ");
                }
                cout << endl;
                exit(-1);
            }
        }
    };

    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            check_query(i, j);
        }
    }

    delete[] values_array;
}

int main() {
    cout << "Simple test for rmq data structures..." << endl;

    mt19937 engine(0);

    for (int N = 1; N <= 100; N++) {
        vector<int32_t> perm(N);
        iota(perm.begin(), perm.end(), 0);
        shuffle(perm.begin(), perm.end(), engine);

        test_exhaustive(perm);
    }

    for (int N = 1; N <= 100; N++) {
        vector<int32_t> values(N);
        for (int i = 0; i < N; i++) {
            values[i] = engine() % 10;
        }

        test_exhaustive(values);
    }

    cout << "All the tests have successfully finished." << endl;

    return 0;
}
