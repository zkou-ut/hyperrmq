
#include <bits/stdc++.h>

#include "RMQRMM64.h"
#include "hyperrmq/hyper_rmq.hpp"
#include "hyperrmq/hypersuccinct_binary_tree.hpp"
#include "sdsl/construct_bwt.hpp"
#include "sdsl/construct_lcp.hpp"
#include "sdsl/rmq_support.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/suffix_trees.hpp"

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

using HighResClockTimepoint =
    std::chrono::time_point<std::chrono::high_resolution_clock>;

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
size_t traverseSuffixTreeHyper(RMQ& rmq, int_vector<> lcp) {
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

// template <class HBT>
// size_t traverseSuffixTreeHBT(HBT& hbt, int_vector<> lcp) {
//     using Node = typename HBT::Node;

//     size_t num_queries = 0;
//     size_t N = lcp.size();
//     stack<state> s;
//     s.emplace(0, N - 1, 0);
//     stack<Node> nodes;
//     nodes.push(hbt.root());

//     while (!s.empty()) {
//         state cur = s.top();
//         s.pop();
//         Node v = nodes.top();
//         nodes.pop();

//         if (cur.i == cur.j) {
//             cout << "Leaf: " << cur.i << endl;
//             continue;
//         }

//         ll cur_i = cur.i;
//         while (cur_i < cur.j) {
//             num_queries++;
//             v = hbt.right_child(v);
//             auto l = hbt.left_child(v);
//             size_t min_i = hbt.node_to_inorder(v);
//             cout << cur_i << " " << cur.j << " " << min_i << endl;
//             ll ii = cur_i;
//             ll jj = cur.j - 1;
//             if (lcp[min_i] == cur.l && min_i < cur.j) {
//                 jj = min_i - 1;
//             } else if (lcp[min_i] != cur.l) {
//                 jj = cur.j;
//                 l = v;
//             }
//             if (ii + 1 <= jj) {
//                 num_queries++;
//                 size_t l_idx = hbt.node_to_inorder(l);
//                 s.emplace(ii, jj, lcp[l_idx]);
//                 nodes.push(l);
//             }
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
    if (argc < 3) {
        cout << "usage : ./build/experiment/lcp ./pizza/&chilli/dna.200MB "
                "./pizza\\&chilli/build"
             << endl;
        return 0;
    }
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

    auto timeinfo =
        chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << std::ctime(&timeinfo) << endl;
    cout << test_id << endl;

    {
        rmq_succinct_rec_new<true, 2048, 1024, 128, 0> rmq(&lcp);
        cout << "RMQ_SDSL_REC_ST" << endl;
        s = time();
        size_t num_queries =
            traverseSuffixTree<rmq_succinct_rec_new<true, 2048, 1024, 128, 0>>(
                rmq, lcp);
        e = time();
        double t = seconds();
        cout << "Time = " << t << endl;
        cout << "Memory = " << size_in_bytes(rmq) * 8 << endl;
        cout << "# of queries=" << num_queries << std::endl;
    }

    {
        size_t N = lcp.size();
        long int* B = new long int[N];
        for (size_t i = 0; i < N; ++i) {
            B[N - i - 1] = lcp[i];
        }
        RMQRMM64 rmq(B, N);
        delete[] B;

        cout << "RMQ_FERRADA" << endl;
        s = time();
        size_t num_queries = traverseSuffixTreeFerrada(rmq, lcp);
        e = time();
        double t = seconds();
        cout << "Time = " << t << endl;
        cout << "Memory = " << rmq.getSize() * 8 << endl;
        cout << "# of queries=" << num_queries << std::endl;
    }

    {
        size_t N = lcp.size();
        std::vector<int> C(N);
        for (size_t i = 0; i < N; ++i) {
            C[i] = lcp[i];
        }

        {
            HyperRMQHuffman rmq(C, int(round(log2(N) / 4)));
            cout << "RMQ_Huffman" << endl;
            s = time();
            size_t num_queries = traverseSuffixTreeHyper(rmq, lcp);
            e = time();
            double t = seconds();
            cout << "Time = " << t << endl;
            cout << "Memory = " << rmq.evaluate_memory_consumption() << endl;
            cout << "# of queries=" << num_queries << std::endl;
        }
        for (int B = 64; B <= 1024; B <<= 1) {
            HyperRMQBreadth rmq(C, B);
            cout << "RMQ_Breadth_" << B << endl;
            s = time();
            size_t num_queries = traverseSuffixTreeHyper(rmq, lcp);
            e = time();
            double t = seconds();
            cout << "Time = " << t << endl;
            cout << "Memory = " << rmq.evaluate_memory_consumption() << endl;
            cout << "# of queries=" << num_queries << std::endl;
        }

        // {
        //     HypersuccinctBinaryTreeHuffman hbt(cartesian_tree_bp(C),
        //                                        int(round(log2(N) / 4)));
        //     cout << "HBT_Huffman" << endl;
        //     s = time();
        //     size_t num_queries = traverseSuffixTreeHBT(hbt, lcp);
        //     e = time();
        //     double t = seconds();
        //     cout << "Time = " << t << endl;
        //     cout << "Memory = " << hbt.evaluate_memory_consumption() << endl;
        //     cout << "# of queries=" << num_queries << std::endl;
        // }
        // for (int B = 64; B <= 1024; B <<= 1) {
        //     HypersuccinctBinaryTreeBreadth hbt(cartesian_tree_bp(C), B);
        //     cout << "HBT_Breadth_" << B << endl;
        //     s = time();
        //     size_t num_queries = traverseSuffixTreeHBT(hbt, lcp);
        //     e = time();
        //     double t = seconds();
        //     cout << "Time = " << t << endl;
        //     cout << "Memory = " << hbt.evaluate_memory_consumption() << endl;
        //     cout << "# of queries=" << num_queries << std::endl;
        // }
    }
}
