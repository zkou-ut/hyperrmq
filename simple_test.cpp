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

template <typename T>
void check_equal(T left, T right) {
    if (left != right) {
        cout << "check_equal failed" << endl;
        cout << "left:  " << left << endl;
        cout << "right: " << right << endl;
        exit(-1);
    }
}

void test_exhaustive(const vector<int>& values, int B_huffman = 4,
                     int B_arithmetic = 64) {
    const int N = values.size();

    HyperRMQHuffman hyper_rmq_huffman(values, B_huffman);
    HyperRMQBreadthFirstArithmetic hyper_rmq_arith(values, B_arithmetic);

    RMQRMM64 rmq_rmm(const_cast<int*>(values.data()), N);

    rmq_succinct_rec_new<true, 2048, 1024, 128, 0> rmq_new(&values);

    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            vector<int> answers = {hyper_rmq_huffman.query(i, j),
                                   hyper_rmq_arith.query(i, j),
                                   rmq_rmm.queryRMQ(i, j), rmq_new(i, j)};

            for (int k = 0; k < answers.size(); k++) {
                check_equal(answers[k], answers[0]);
            }
        }
    }
}

int main() {
    cout << "Simple test for rmq data structures..." << endl;

    for (int N = 1; N <= 100; N++) {
        mt19937 engine(0);
        vector<int32_t> perm(N);
        iota(perm.begin(), perm.end(), 0);
        shuffle(perm.begin(), perm.end(), engine);

        test_exhaustive(perm);
    }

    cout << "All the tests have successfully finished." << endl;

    return 0;
}
