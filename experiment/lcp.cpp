
#include <bits/stdc++.h>

#include "RMQRMM64.h"
#include "hyperrmq/hyper_rmq.hpp"
#include "hyperrmq/hypersuccinct_binary_tree.hpp"
#include "sdsl/rmq_support.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/suffix_trees.hpp"
#include "sdsl/construct_lcp.hpp"
#include "sdsl/construct_bwt.hpp"

// This file is based on the file from the following repository:
// https://github.com/kittobi1992/rmq-experiments/blob/master/executer/lcp_experiment.cpp

// usage:
// ./build/experiment/lcp ./pizza/&chilli/dna.200MB ./pizza\&chilli/build

using namespace std;
using namespace sdsl;
using namespace hyperrmq;

using ll = long long;

struct state {
    ll i, j, l;
    state(ll i, ll j, ll l) : i(i), j(j), l(l) {}
};

using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

HighResClockTimepoint s, e;

inline HighResClockTimepoint time() {
    return std::chrono::high_resolution_clock::now();
}

double seconds() {
    std::chrono::duration<double> elapsed_seconds = e - s;
    return elapsed_seconds.count();
}

using HyperRMQHuffman =
    HyperRMQ<16, CompressedMicrotreeSplitRankArrayHuffman<16>>;
using HyperRMQBreadth =
    HyperRMQ<16, CompressedMicrotreeSplitRankArrayArithmetic<false>>;

using HypersuccinctBinaryTreeHuffman =
    HypersuccinctBinaryTree<16, CompressedMicrotreeSplitRankArrayHuffman<16>>;
using HypersuccinctBinaryTreeBreadth =
    HypersuccinctBinaryTree<16,
                            CompressedMicrotreeSplitRankArrayArithmetic<false>>;

void construct_lcp(cache_config& test_config, string& test_file) {
    {
        cout << "Load text..." << endl;
        int_vector<8> text;
        load_vector_from_file(text, test_file, 1);
        append_zero_symbol(text);
        store_to_cache(text, conf::KEY_TEXT, test_config);

        cout << "Construct Suffix Array..." << endl;
        int_vector<> sa(text.size(), 0, bits::hi(text.size()) + 1);
        algorithm::calculate_sa((const unsigned char*)text.data(), text.size(),
                                sa);
        store_to_cache(sa, conf::KEY_SA, test_config);
    }

    {
        cout << "Construct LCP Array..." << endl;
        construct_lcp_PHI<8>(test_config);
    }
}

template <class RMQ>
size_t traverseSuffixTree(RMQ& rmq, int_vector<> lcp) {
    size_t num_queries = 0;
    size_t N = rmq.size();
    stack<state> s;
    s.emplace(0, N - 1, 0);

    while (!s.empty()) {
        state cur = s.top();
        s.pop();
        //                    if(cur.i != cur.j) cout << "Internal Node: (" <<
        //                    cur.i << "," << cur.j << ") - " << cur.l << endl;
        if (cur.i == cur.j) {
            //                           cout << "Leaf: " << cur.i << endl;
            continue;
        }

        ll cur_i = cur.i;
        while (cur_i < cur.j) {
            num_queries++;
            size_t min_i = rmq(cur_i + 1, cur.j);
            //                           cout << cur_i << " " << cur.j << " " <<
            //                           min_i << endl;
            ll ii = cur_i;
            ll jj = cur.j - 1;
            if (lcp[min_i] == cur.l && min_i < cur.j)
                jj = min_i - 1;
            else if (lcp[min_i] != cur.l)
                jj = cur.j;
            if (ii + 1 <= jj) {
                num_queries++;
                size_t l_idx = rmq(ii + 1, jj);
                s.emplace(ii, jj, lcp[l_idx]);
            } else if (ii == jj)
                s.emplace(ii, ii, lcp[ii]);
            cur_i = jj + 1;
        }
    }

    return num_queries;
}

template <class RMQ>
size_t traverseSuffixTreeHyper(RMQ & rmq, int_vector<> lcp) {
    size_t num_queries = 0;
    size_t N = lcp.size();
    stack<state> s;
    s.emplace(0, N - 1, 0);

    while (!s.empty()) {
        state cur = s.top();
        s.pop();
        //                    if(cur.i != cur.j) cout << "Internal Node: (" <<
        //                    cur.i << "," << cur.j << ") - " << cur.l << endl;
        if (cur.i == cur.j) {
            //                           cout << "Leaf: " << cur.i << endl;
            continue;
        }

        ll cur_i = cur.i;
        while (cur_i < cur.j) {
            num_queries++;
            size_t min_i = rmq.query(cur_i + 1, cur.j);
            //                           cout << cur_i << " " << cur.j << " " <<
            //                           min_i << endl;
            ll ii = cur_i;
            ll jj = cur.j - 1;
            if (lcp[min_i] == cur.l && min_i < cur.j)
                jj = min_i - 1;
            else if (lcp[min_i] != cur.l)
                jj = cur.j;
            if (ii + 1 <= jj) {
                num_queries++;
                size_t l_idx = rmq.query(ii + 1, jj);
                s.emplace(ii, jj, lcp[l_idx]);
            } else if (ii == jj)
                s.emplace(ii, ii, lcp[ii]);
            cur_i = jj + 1;
        }
    }

    return num_queries;
}

// template <class HypersuccinctT>
// size_t traverseSuffixTreeHyper(RMQ& rmq, int_vector<> lcp) {
//     size_t num_queries = 0;
//     size_t N = lcp.size();
//     stack<state> s;
//     s.emplace(0, N - 1, 0);

//     while (!s.empty()) {
//         state cur = s.top();
//         s.pop();
//         //                    if(cur.i != cur.j) cout << "Internal Node: (" <<
//         //                    cur.i << "," << cur.j << ") - " << cur.l << endl;
//         if (cur.i == cur.j) {
//             //                           cout << "Leaf: " << cur.i << endl;
//             continue;
//         }

//         ll cur_i = cur.i;
//         while (cur_i < cur.j) {
//             num_queries++;
//             size_t min_i = rmq.query(cur_i + 1, cur.j);
//             //                           cout << cur_i << " " << cur.j << " " <<
//             //                           min_i << endl;
//             ll ii = cur_i;
//             ll jj = cur.j - 1;
//             if (lcp[min_i] == cur.l && min_i < cur.j)
//                 jj = min_i - 1;
//             else if (lcp[min_i] != cur.l)
//                 jj = cur.j;
//             if (ii + 1 <= jj) {
//                 num_queries++;
//                 size_t l_idx = rmq.query(ii + 1, jj);
//                 s.emplace(ii, jj, lcp[l_idx]);
//             } else if (ii == jj)
//                 s.emplace(ii, ii, lcp[ii]);
//             cur_i = jj + 1;
//         }
//     }

//     return num_queries;
// }

size_t traverseSuffixTreeFerrada(RMQRMM64& rmq, int_vector<> lcp) {
    size_t num_queries = 0;
    size_t N = lcp.size();
    stack<state> s;
    s.emplace(0, N - 1, 0);

    while (!s.empty()) {
        state cur = s.top();
        s.pop();
        //         if(cur.i != cur.j) cout << "Internal Node: (" << cur.i << ","
        //         << cur.j << ") - " << cur.l << endl;
        if (cur.i == cur.j) {
            //             cout << "Leaf: " << cur.i << endl;
            continue;
        }

        ll cur_i = cur.i;
        while (cur_i < cur.j) {
            num_queries++;
            size_t min_i =
                N - rmq.queryRMQ(N - cur.j - 1, N - (cur_i + 1) - 1) - 1;
            //                           cout << cur_i << " " << cur.j << " " <<
            //                           min_i << endl;
            ll ii = cur_i;
            ll jj = cur.j - 1;
            if (lcp[min_i] == cur.l && min_i < cur.j)
                jj = min_i - 1;
            else if (lcp[min_i] != cur.l)
                jj = cur.j;
            if (ii + 1 <= jj) {
                num_queries++;
                size_t l_idx =
                    N - rmq.queryRMQ(N - jj - 1, N - (ii + 1) - 1) - 1;
                s.emplace(ii, jj, lcp[l_idx]);
            } else if (ii == jj)
                s.emplace(ii, ii, lcp[ii]);
            cur_i = jj + 1;
        }
    }
    return num_queries;
}

int main(int argc, char* argv[]) {
    string test_file, temp_dir, test_id;

    test_file = argv[1];
    temp_dir = argv[2];
    test_id = test_file.substr(test_file.find_last_of("/\\") + 1);

    cache_config test_config = cache_config(false, temp_dir, test_id);

    string lcp_file = cache_file_name(conf::KEY_LCP, test_config);
    int_vector<> lcp;
    if (!load_from_file(lcp, lcp_file)) {
        construct_lcp(test_config, test_file);
        load_from_file(lcp, lcp_file);
    }

    string algo1 = "RMQ_SDSL_REC_ST";
    string algo2 = "RMQ_SUCCINCT_SCT";
    string algo3 = "RMQ_FERRADA";
    string algo4 = "RMQ_HUFFMAN";
    string algo5 = "RMQ_BREADTH_B_";
    string algo6 = "HBT_HUFFMAN";
    string algo7 = "HBT_BREADTH_B_";

    {
        rmq_succinct_rec_new<true, 0, 1024, 128, 0> rmq(&lcp);
        cout << "Start Suffix-Tree Traversion for RMQ " << algo1 << "..."
             << endl;
        s = time();
        size_t num_queries =
            traverseSuffixTree<rmq_succinct_rec_new<true, 0, 1024, 128, 0>>(
                rmq, lcp);
        e = time();
        // double percentage_avoided =
        // (static_cast<double>(rmq.num_avoided_selects)/static_cast<double>(rmq.num_queries));
        // std::cout << rmq.num_avoided_selects << " out of " << rmq.num_queries
        // << " queries (" << percentage_avoided << "%) avoids second select" <<
        // std::endl;
        double t = seconds();
        std::cout << "LCP_RESULT Benchmark=" << test_id << " Algo=" << algo1
                  << " Time=" << t << " NumQueries=" << num_queries
                  << std::endl;
    }

    {
        rmq_succinct_rec_new<true, 2048, 1024, 128, 0> rmq(&lcp);
        cout << "Start Suffix-Tree Traversion for RMQ " << algo2 << "..."
             << endl;
        s = time();
        size_t num_queries =
            traverseSuffixTree<rmq_succinct_rec_new<true, 2048, 1024, 128, 0>>(
                rmq, lcp);
        e = time();
        // double percentage_avoided =
        // (static_cast<double>(rmq.num_avoided_selects)/static_cast<double>(rmq.num_queries));
        // std::cout << rmq.num_avoided_selects << " out of " << rmq.num_queries
        // << " queries (" << percentage_avoided << "%) avoids second select" <<
        // std::endl;
        double t = seconds();
        std::cout << "LCP_RESULT Benchmark=" << test_id << " Algo=" << algo2
                  << " Time=" << t << " NumQueries=" << num_queries
                  << std::endl;
    }

    {
        size_t N = lcp.size();
        long int* B = new long int[N];
        for (size_t i = 0; i < N; ++i) {
            B[N - i - 1] = lcp[i];
        }
        RMQRMM64 rmq(B, N);
        delete[] B;
        cout << "Start Suffix-Tree Traversion for RMQ " << algo3 << "..."
             << endl;
        s = time();
        size_t num_queries = traverseSuffixTreeFerrada(rmq, lcp);
        e = time();
        double t = seconds();
        std::cout << "LCP_RESULT Benchmark=" << test_id << " Algo=" << algo3
                  << " Time=" << t << " NumQueries=" << num_queries
                  << std::endl;
    }

    {
        size_t N = lcp.size();
        std::vector<int> C(N);
        for (size_t i = 0; i < N; ++i) {
            C[i] = lcp[i];
        }

        {
            HyperRMQHuffman rmq(C, int(round(log2(N))));
            cout << "Start Suffix-Tree Traversion for RMQ " << algo4 << "..."
                 << endl;
            s = time();
            size_t num_queries = traverseSuffixTreeHyper(rmq, lcp);
            e = time();
            double t = seconds();
            std::cout << "LCP_RESULT Benchmark=" << test_id << " Algo=" << algo4
                      << " Time=" << t << " NumQueries=" << num_queries
                      << std::endl;
        }
        for (int B = 64; B <= 1024; B <<= 1) {
            HyperRMQBreadth rmq(C, B);
            cout << "Start Suffix-Tree Traversion for RMQ " << algo5 << B
                 << "..." << endl;
            s = time();
            size_t num_queries = traverseSuffixTreeHyper(rmq, lcp);
            e = time();
            double t = seconds();
            std::cout << "LCP_RESULT Benchmark=" << test_id << " Algo=" << algo5
                      << " Time=" << t << " NumQueries=" << num_queries
                      << std::endl;
        }
    }

    // {
    //     size_t N = lcp.size();
    //     std::vector<int> C(N);
    //     for (size_t i = 0; i < N; ++i) {
    //         C[i] = lcp[i];
    //     }

    //     {
    //         HypersuccinctBinaryTreeHuffman rmq(N, round(log2(N)));
    //         cout << "Start Suffix-Tree Traversion for RMQ " << algo6 << "..."
    //              << endl;
    //         s = time();
    //         size_t num_queries = traverseSuffixTreeHyper(rmq, lcp);
    //         e = time();
    //         double t = seconds();
    //         std::cout << "LCP_RESULT Benchmark=" << test_id << " Algo=" << algo6
    //                   << " Time=" << t << " NumQueries=" << num_queries
    //                   << std::endl;
    //     }
    //     for (int B = 64; B <= 1024; B <<= 1) {
    //         HypersuccinctBinaryTreeBreadth rmq(C, B);
    //         cout << "Start Suffix-Tree Traversion for RMQ " << algo7 << B
    //              << "..." << endl;
    //         s = time();
    //         size_t num_queries = traverseSuffixTreeHyper(rmq, lcp);
    //         e = time();
    //         double t = seconds();
    //         std::cout << "LCP_RESULT Benchmark=" << test_id << " Algo=" << algo7
    //                   << " Time=" << t << " NumQueries=" << num_queries
    //                   << std::endl;
    //     }
    // }
}
