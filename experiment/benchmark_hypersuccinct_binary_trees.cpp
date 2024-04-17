#include <bits/stdc++.h>

#include "hyperrmq/hypersuccinct_binary_tree.hpp"

using namespace std;
using namespace hyperrmq;

void output_benchmark(
    pair<vector<pair<string, uint64_t>>, vector<pair<string, int>>> bm_result) {
    auto [memory_result, time_result] = bm_result;
    cout << "memory:" << endl;
    for (auto&& [label, val] : memory_result) {
        cout << setw(20) << label << ":" << setw(20) << val << endl;
    }

    cout << "time:" << endl;
    for (auto&& [label, val] : time_result) {
        cout << setw(20) << label << ":" << setw(20) << val << " ms" << endl;
    }
}

template <int N, int Q, int B, int W,
          typename CompressedMicrotreeSplitRankArray>
void benchmark() {
    using Node = typename HypersuccinctBinaryTree<
        W, CompressedMicrotreeSplitRankArray>::Node;

    auto timeinfo =
        chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << std::ctime(&timeinfo) << endl;
    cout << "N = " << N << endl;
    cout << "Q = " << Q << endl;
    cout << "B = " << B << endl;
    cout << "W = " << W << endl;
    cout << typeid(CompressedMicrotreeSplitRankArray).name() << endl;

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

    mt19937 engine(0);
    vector<int32_t> perm(N);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), engine);

    time_keeper("cnst.", 0);
    HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray> hbt(
        cartesian_tree_bp(perm), B);
    time_keeper("cnst.", 1);

    std::vector<int> inorders(Q + 1);
    for (int i = 0; i < Q + 1; i++) {
        inorders[i] = engine() % N;
    }

    std::vector<Node> nodes;
    time_keeper("to node", 0);
    for (int i = 0; i < Q; i++) {
        nodes.push_back(hbt.inorder_to_node(inorders[i]));
    }
    time_keeper("to node", 1);
    nodes.push_back(hbt.inorder_to_node(inorders.back()));

    time_keeper("to in.", 0);
    for (int i = 0; i < Q; i++) {
        hbt.node_to_inorder(nodes[i]);
    }
    time_keeper("to in.", 1);

    time_keeper("root", 0);
    for (int i = 0; i < Q; i++) {
        hbt.root();
    }
    time_keeper("root", 1);

    time_keeper("parent", 0);
    for (int i = 0; i < Q; i++) {
        hbt.parent(nodes[i]);
    }
    time_keeper("parent", 1);

    time_keeper("left", 0);
    for (int i = 0; i < Q; i++) {
        hbt.left_child(nodes[i]);
    }
    time_keeper("left", 1);

    time_keeper("right", 0);
    for (int i = 0; i < Q; i++) {
        hbt.right_child(nodes[i]);
    }
    time_keeper("right", 1);

    time_keeper("is leaf", 0);
    for (int i = 0; i < Q; i++) {
        hbt.is_leaf(nodes[i]);
    }
    time_keeper("is leaf", 1);

    time_keeper("child_label", 0);
    for (int i = 0; i < Q; i++) {
        hbt.child_label(nodes[i]);
    }
    time_keeper("child_label", 1);

    time_keeper("subtree", 0);
    for (int i = 0; i < Q; i++) {
        hbt.subtree_size(nodes[i]);
    }
    time_keeper("subtree", 1);

    time_keeper("is anc.", 0);
    for (int i = 0; i < Q; i++) {
        hbt.is_ancestor(nodes[i], nodes[i + 1]);
    }
    time_keeper("is anc.", 1);

    time_keeper("ldesc.", 0);
    for (int i = 0; i < Q; i++) {
        hbt.leftmost_desc(nodes[i]);
    }
    time_keeper("ldesc.", 1);

    time_keeper("rdesc.", 0);
    for (int i = 0; i < Q; i++) {
        hbt.rightmost_desc(nodes[i]);
    }
    time_keeper("rdesc.", 1);

    time_keeper("lca", 0);
    for (int i = 0; i < Q; i++) {
        hbt.lca(nodes[i], nodes[i + 1]);
    }
    time_keeper("lca", 1);

    vector<pair<string, uint64_t>> memory_result;

    memory_result.push_back({"total", hbt.evaluate_memory_consumption()});
    memory_result.push_back(
        {"rmm tree", hbt.rmm_tree.evaluate_memory_consumption()});
    memory_result.push_back(
        {"code seq", hbt.compressed_microtree_split_rank_array
                         .evaluate_memory_consumption()});
    memory_result.push_back(
        {"close_sample", hbt.close_sample.evaluate_memory_consumption()});

    output_benchmark({memory_result, time_result});
};

int main() {
    cout << fixed << setprecision(5);

    constexpr int N = 1e9, Q = 1e6;

    benchmark<N, Q, 8, 16, CompressedMicrotreeSplitRankArrayHuffman<16>>();
    benchmark<N, Q, 64, 16,
              CompressedMicrotreeSplitRankArrayArithmetic<true>>();
    benchmark<N, Q, 128, 16,
              CompressedMicrotreeSplitRankArrayArithmetic<true>>();
    benchmark<N, Q, 256, 16,
              CompressedMicrotreeSplitRankArrayArithmetic<true>>();
    benchmark<N, Q, 512, 16,
              CompressedMicrotreeSplitRankArrayArithmetic<true>>();
    benchmark<N, Q, 1024, 16,
              CompressedMicrotreeSplitRankArrayArithmetic<true>>();

    return 0;
}
