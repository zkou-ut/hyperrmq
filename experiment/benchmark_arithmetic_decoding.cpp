#include <bits/stdc++.h>

#include "hyperrmq/arithmetic.hpp"

using namespace std;
using namespace hyperrmq;

int main() {
    cout << fixed << setprecision(6);

    auto timeinfo =
        chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << std::ctime(&timeinfo) << endl;

    std::mt19937_64 engine(0);
    auto random_arith = [&](int64_t tree_n) -> BitArray {
        std::vector<uint64_t> left_seq;
        left_seq.reserve(tree_n);
        std::queue<uint64_t> subtree_size_st;
        subtree_size_st.push(tree_n);
        for (int64_t i = 0; i < tree_n; i++) {
            uint64_t cur_size = subtree_size_st.front();
            subtree_size_st.pop();
            uint64_t left_size = engine() % cur_size;
            uint64_t right_size = cur_size - left_size - 1;
            left_seq.push_back(left_size);
            if (left_size) {
                subtree_size_st.push(left_size);
            }
            if (right_size) {
                subtree_size_st.push(right_size);
            }
        }
        return left_seq_to_arithmetic<false>(left_seq);
    };

    for (int64_t s = 1; s <= 20; s++) {
        int64_t n = 1 << s;
        int64_t Q = 1 << (25 - s);
        vector<BitArray> arith_list;
        for (int64_t i = 0; i < Q; i++) {
            arith_list.push_back(random_arith(n));
        }

        auto start = chrono::system_clock::now();
        for (auto&& arith : arith_list) {
            left_seq_to_bp<false>(arithmetic_to_left_seq<false>(n, arith));
        }
        auto end = chrono::system_clock::now();
        double time = static_cast<double>(
            chrono::duration_cast<chrono::microseconds>(end - start).count());

        cout << n << " " << Q << " " << time / Q << endl;
    }

    return 0;
}
