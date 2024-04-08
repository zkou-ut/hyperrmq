#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <stack>

#include "hyperrmq/bit_array.hpp"
#include "hyperrmq/bitutil.hpp"
#include "hyperrmq/huffman.hpp"
#include "hyperrmq/memutil.hpp"
#include "hyperrmq/microtree_array.hpp"
#include "hyperrmq/rmm_tree.hpp"
#include "hyperrmq/three_level_prefix_sum.hpp"
#include "hyperrmq/tree_covering.hpp"

namespace average_case_optimal_rmq {

namespace {
// I gave up writing Node class with lazy evaluation.

// struct Node {
//     // `chunk_index == -1` represents that the node does not exist.

//     Node(const HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>
//              &tree_ref,
//          int32_t chunk_index, int32_t local_inorder_or_index,
//          bool is_local_inorder = true)
//         : tree_ref(tree_ref), chunk_index(chunk_index) {
//         if (!is_valid()) return;

//         if (tree_ref.rmm_tree.get_bit(chunk_index) == 0) {
//             left_chunk_index = chunk_index;
//         } else {
//             right_chunk_index = chunk_index;
//         }
//         if (is_local_inorder) {
//             local_inorder = local_inorder_or_index;
//         } else {
//             local_index = local_inorder_or_index;
//         }
//     }

//     // inexpensive operations
//     bool is_valid() { return chunk_index != -1; }

//     int32_t get_chunk_index() { return chunk_index; }

//     int32_t get_left_chunk_index() {
//         if (!is_valid()) return -1;

//         if (left_chunk_index == -1) {
//             left_chunk_index = tree_ref.rmm_tree.open(right_chunk_index);
//         }
//         return left_chunk_index;
//     }

//     int32_t get_right_chunk_index() {
//         if (!is_valid()) return -1;

//         if (right_chunk_index == -1) {
//             right_chunk_index = tree_ref.rmm_tree.close(left_chunk_index);
//         }
//         return right_chunk_index;
//     }

//     int32_t get_microtree_preorder() {
//         if (!is_valid()) return -1;

//         if (microtree_preorder == -1) {
//             microtree_preorder =
//                 tree_ref.rmm_tree.rank0(get_left_chunk_index());
//         }
//         return microtree_preorder;
//     }

//     int32_t get_microtree_node_count() {
//         if (!is_valid()) return -1;

//         if (node_count == -1) {
//             node_count =
//                 tree_ref.compressed_microtree_split_rank_array.get_node_count(
//                     get_microtree_preorder());
//         }
//         return node_count;
//     }

//     int32_t get_microtree_split_rank() {
//         if (!is_valid()) return -1;

//         if (split_rank == -1) {
//             split_rank = tree_ref.compressed_microtree_split_rank_array
//                          .get_left_chunk_popcount(get_microtree_preorder());
//         }
//         return split_rank;
//     }

//     // possibly expensive operations
//     std::shared_ptr<TreeBP> get_microtree() {
//         if (!is_valid()) return nullptr;

//         if (microtree_ptr == nullptr) {
//             auto [microtree, split_rank_] =
//                 tree_ref.compressed_microtree_split_rank_array
//                     [get_microtree_preorder()];
//             microtree_ptr = std::make_shared<TreeBP>(microtree);
//             split_rank = split_rank_;
//         }
//         return microtree_ptr;
//     }

//     int32_t get_cut() {
//         if (!is_valid()) return -1;

//         if (cut == -1) {
//             auto ptr = get_microtree();
//             cut = ptr->bp.linear_select1(ptr->second);
//         }
//         return cut;
//     }

//     int32_t get_local_inorder() {
//         if (!is_valid()) return -1;

//         if (local_inorder == -1) {
//             local_inorder = get_microtree()->bp.linear_popcount(0,
//             local_index);
//         }
//         return local_inorder;
//     }

//     int32_t get_local_index() {
//         if (!is_valid()) return -1;

//         if (local_index == -1) {
//             local_index = get_microtree()->bp.linear_select1(local_inorder);
//         }
//         return local_index;
//     }

//     void inc() {
//         if (!is_valid()) return;

//         // Note that either of these is defined and the equality does not
//         // hold for the undefined one.
//         bool reach_end = (local_index == get_microtree_node_count() * 2 - 1)
//         ||
//                          (local_inorder == get_microtree_node_count() - 1);

//         if (reach_end) {
//             this = Node(right_chunk_index + 1, 0, false);
//         } else {
//             if (microtree_ptr == nullptr) {
//                 local_inorder = -1;
//             } else {
//                 local_inorder = microtree_ptr->bp.get(local_index);
//             }
//             local_index++;
//         }
//     }

//     void dec() {
//         if (!is_valid()) return;

//         // Note that either of these is defined and the equality does not
//         // hold for the undefined one.
//         if (get_local_index() == 0) {
//             if (chunk_index == 0) {
//                 this = Node(tree_ref, -1, -1);
//                 return;
//             }
//             while (true) {
//                 left_chunk_index--;
//                 Node tmp = Node(left_chunk_index, 0);
//                 if (tmp.get_microtree_node_count()) {
//                     tmp.local_index = tmp.get_microtree_node_count() * 2 - 1;
//                     tmp.local_inorder = tmp.get_microtree_node_count() - 1;
//                 }
//             }
//         } else {
//             local_index--;
//             if (microtree_ptr == nullptr) {
//                 local_inorder = -1;
//             } else {
//                 local_inorder = microtree_ptr->bp.get(local_index);
//             }
//         }
//     }

//     void open() {
//         local_index = get_microtree()->naive_open(get_local_index());
//         local_inorder = -1;
//     }

//     void close() {
//         local_index = get_microtree()->naive_close(get_local_index());
//         local_inorder = -1;
//     }

//     void enclose() {}

//    private:
//     const HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>
//     &tree_ref;

//     int32_t chunk_index = -1, left_chunk_index = -1, right_chunk_index = -1,
//             microtree_preorder = -1;
//     int32_t local_inorder = -1, local_index = -1;
//     int32_t node_count = -1, split_rank = -1, cut = -1;
//     std::shared_ptr<TreeBP> microtree_ptr = nullptr;
// };
}  // namespace

template <uint32_t W = 64, typename CompressedMicrotreeSplitRankArray =
                               CompressedMicrotreeSplitRankArrayHuffman<W>>
struct HypersuccinctBinaryTree {
    static_assert(std::is_base_of_v<CompressedMicrotreeSplitRankArrayInterface,
                                    CompressedMicrotreeSplitRankArray>);

    struct Node {
        // `chunk_index == -1` represents that the node does not exist.

        friend HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>;

       private:
        const HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>
            *tree_ptr = nullptr;

        Node(const HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>
                 *tree_ptr)
            : tree_ptr(tree_ptr) {}

        void set_chunk_index(int32_t chunk_index_) {
            chunk_index = chunk_index_;
            if (!is_valid()) return;

            if (tree_ptr->rmm_tree.get_bit(chunk_index) == 0) {
                left_chunk_index = chunk_index;
                right_chunk_index = tree_ptr->rmm_tree.close(chunk_index);
            } else {
                left_chunk_index = tree_ptr->rmm_tree.open(chunk_index);
                right_chunk_index = chunk_index;
            }
            microtree_preorder = tree_ptr->rmm_tree.rank0(left_chunk_index);

            node_count =
                tree_ptr->compressed_microtree_split_rank_array.get_node_count(
                    microtree_preorder);
            split_rank = tree_ptr->compressed_microtree_split_rank_array
                             .get_left_chunk_popcount(microtree_preorder);
        }

        void set_tree() {
            auto [microtree, split_rank_] =
                tree_ptr
                    ->compressed_microtree_split_rank_array[microtree_preorder];
            microtree_ptr = std::make_shared<TreeBP>(microtree);

            cut = microtree.bp.linear_select1(split_rank);
        }

        void set_local(int32_t local_inorder_or_index,
                       bool is_local_inorder = true, bool reset_chunk = false) {
            if (!is_valid()) return;

            if (is_local_inorder) {
                local_inorder = local_inorder_or_index;
                local_index = microtree_ptr->bp.linear_select1(local_inorder);
            } else {
                local_index = local_inorder_or_index;
                local_inorder =
                    microtree_ptr->bp.linear_popcount(0, local_index);
            }

            if (reset_chunk) {
                if (local_index < cut) {
                    chunk_index = left_chunk_index;
                } else {
                    chunk_index = right_chunk_index;
                }
            }
        }

        // Returns node with set_chunk already called.
        Node search_nonzero_chunk(int32_t chunk_index_start,
                                  int32_t step) const {
            auto next_index = chunk_index_start;
            while (0 <= next_index && next_index < tree_ptr->num_of_chunks) {
                auto tmp = Node(tree_ptr);
                tmp.set_chunk_index(next_index);
                if (tmp.chunk_index == tmp.left_chunk_index ||
                    tmp.split_rank < tmp.node_count) {
                    return tmp;
                }
                next_index += step;
            }
            return Node(tree_ptr);
        }

       public:
        // The following variables are assigned by set_chunk_index.
        int32_t chunk_index = -1, left_chunk_index = -1, right_chunk_index = -1,
                microtree_preorder = -1;
        int32_t node_count = -1, split_rank = -1;
        // The following variables are assigned by set_tree.
        std::shared_ptr<TreeBP> microtree_ptr = nullptr;
        int32_t cut = -1;
        // The following variables are assigned by set_local.
        int32_t local_inorder = -1, local_index = -1;

        Node(const HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>
                 *tree_ptr,
             int32_t chunk_index_, int32_t local_inorder_or_index,
             bool is_local_inorder = true)
            : tree_ptr(tree_ptr) {
            if (chunk_index_ != -1) {
                set_chunk_index(chunk_index_);
                set_tree();
            }
            if (local_inorder_or_index != -1) {
                set_local(local_inorder_or_index, is_local_inorder);
            }
        }

        bool is_valid() const { return chunk_index != -1; }

        int access() const {
            if (!is_valid()) return -1;
            return microtree_ptr->bp.get(local_index);
        }

        void inc() {
            if (!is_valid()) return;

            if (local_index == microtree_ptr->n * 2 - 1 ||
                local_index == cut - 1) {
                auto tmp = search_nonzero_chunk(chunk_index + 1, 1);
                if (tmp.is_valid()) {
                    tmp.set_tree();
                    tmp.set_local(
                        tmp.chunk_index == tmp.left_chunk_index ? 0 : tmp.cut,
                        false);
                }
                *this = tmp;
            } else {
                local_inorder += microtree_ptr->bp.get(local_index);
                local_index++;
            }
        }

        void dec() {
            if (!is_valid()) return;

            if (local_index == 0 || local_index == cut) {
                auto tmp = search_nonzero_chunk(chunk_index - 1, -1);
                if (tmp.is_valid()) {
                    tmp.set_tree();
                    tmp.set_local(tmp.chunk_index == tmp.left_chunk_index
                                      ? tmp.cut - 1
                                      : tmp.node_count * 2 - 1,
                                  false);
                }
                *this = tmp;
            } else {
                local_index--;
                local_inorder -= microtree_ptr->bp.get(local_index);
            }
        }

        void open() {
            if (!is_valid()) return;

            set_local(microtree_ptr->naive_open(local_index), false, true);
        }

        void close() {
            if (!is_valid()) return;

            set_local(microtree_ptr->naive_close(local_index), false, true);
        }

        void rightmost_desc() {
            if (!is_valid()) return;
            assert(access() == 1);

            auto pos = microtree_ptr->naive_fwdsearch(local_index, -2) - 1;
            if (pos == node_count * 2) {
                auto next_index =
                    tree_ptr->rmm_tree.fwdsearch(right_chunk_index, -2) - 2;
                Node tmp(tree_ptr);
                if (next_index != tree_ptr->num_of_chunks) {
                    tmp.set_chunk_index(next_index);
                    tmp.set_tree();
                    tmp.set_local(tmp.node_count - 1);
                }
                *this = tmp;
            } else {
                set_local(pos, false, true);
                dec();
            }
        }

        bool operator==(const Node &rhs) const {
            return this->chunk_index == rhs.chunk_index &&
                   (this->chunk_index == -1 ||
                    this->local_index == rhs.local_index);
        }
        bool operator!=(const Node &rhs) const { return !(*this == rhs); }
        bool operator<(const Node &rhs) const {
            return this->chunk_index < rhs.chunk_index ||
                   (this->chunk_index == rhs.chunk_index &&
                    this->local_index < rhs.local_index);
        }
        bool operator<=(const Node &rhs) const { return !(rhs < *this); }

        friend std::ostream &operator<<(std::ostream &stream,
                                        const Node &node) {
            return stream << "Node { " << node.chunk_index << ", "
                          << node.local_index << " }";
        }
    };

    uint32_t num_of_nodes, num_of_microtrees, num_of_chunks;
    RMMTree<8, 1024> rmm_tree;
    CompressedMicrotreeSplitRankArray compressed_microtree_split_rank_array;
    ThreeLevelPrefixSum<> close_sample;

    HypersuccinctBinaryTree();
    explicit HypersuccinctBinaryTree(const TreeBP &tree, const int B);
    explicit HypersuccinctBinaryTree(
        uint32_t n, const std::pair<TreeBP, MicrotreeSplitRankArray>
                        &upsilon_and_microtree_split_rank_array);

    uint32_t chunk_index_to_microtree_preorder(uint32_t i) const;
    BitArray get_chunk(uint32_t i) const;
    uint64_t get_chunk_popcount(uint32_t i) const;

    Node inorder_to_node(uint32_t inorder) const;

    int32_t node_to_inorder(Node node) const;

    Node root() const;
    Node parent(Node v) const;
    Node left_child(Node v) const;
    Node right_child(Node v) const;

    bool is_leaf(Node v) const;
    int32_t child_label(Node v) const;

    int32_t subtree_size(Node v) const;
    bool is_ancestor(Node u, Node v) const;

    Node leftmost_desc(Node v) const;
    Node rightmost_desc(Node v) const;

    Node lca(Node u, Node v) const;

    uint64_t evaluate_memory_consumption() const;
};

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline HypersuccinctBinaryTree<
    W, CompressedMicrotreeSplitRankArray>::HypersuccinctBinaryTree() {}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::
    HypersuccinctBinaryTree(const TreeBP &tree, const int B)
    : HypersuccinctBinaryTree(tree.n, tree_covering(B, tree)) {}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::
    HypersuccinctBinaryTree(uint32_t n,
                            const std::pair<TreeBP, MicrotreeSplitRankArray>
                                &upsilon_and_microtree_split_rank_array)
    : num_of_nodes(n) {
    const auto &[upsilon, microtree_split_rank_array] =
        upsilon_and_microtree_split_rank_array;

    rmm_tree = decltype(rmm_tree)(upsilon);
    compressed_microtree_split_rank_array =
        decltype(compressed_microtree_split_rank_array)(
            microtree_split_rank_array);

    num_of_microtrees = compressed_microtree_split_rank_array.size();
    num_of_chunks = num_of_microtrees * 2;

    std::vector<uint32_t> close_sample_vec;
    close_sample_vec.reserve(num_of_chunks / W);
    uint32_t close = 0;
    int rank0 = 0;
    std::stack<int> opens;
    for (int i = 0; i < num_of_chunks; i++) {
        if (rmm_tree.get_bit(i) == 0) {
            const auto &[tree, split_rank] = microtree_split_rank_array[rank0];
            opens.push(rank0);
            rank0++;
            close += split_rank;
        } else {
            const auto &[tree, split_rank] =
                microtree_split_rank_array[opens.top()];
            opens.pop();
            close += tree.n - split_rank;
        }
        if ((i + 1) % W == 0) {
            close_sample_vec.push_back(close);
            close = 0;
        }
    }
    close_sample = decltype(close_sample)(close_sample_vec);
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline uint32_t HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::
    chunk_index_to_microtree_preorder(uint32_t i) const {
    uint32_t open_index = (rmm_tree.get_bit(i) == 1) ? rmm_tree.open(i) : i;
    uint32_t open_rank = rmm_tree.rank0(open_index);
    return open_rank;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
BitArray
HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::get_chunk(
    uint32_t i) const {
    assert(0 <= i && i < num_of_chunks);
    bool isopen = (rmm_tree.get_bit(i) == 0);
    i = chunk_index_to_microtree_preorder(i);
    if (isopen) {
        return compressed_microtree_split_rank_array.get_left_chunk(i);
    } else {
        return compressed_microtree_split_rank_array.get_right_chunk(i);
    }
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline uint64_t HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::
    get_chunk_popcount(uint32_t i) const {
    assert(0 <= i && i < num_of_chunks);
    bool isopen = (rmm_tree.get_bit(i) == 0);
    i = chunk_index_to_microtree_preorder(i);
    if (isopen) {
        return compressed_microtree_split_rank_array.get_left_chunk_popcount(i);
    } else {
        return compressed_microtree_split_rank_array.get_right_chunk_popcount(
            i);
    }
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline typename HypersuccinctBinaryTree<W,
                                        CompressedMicrotreeSplitRankArray>::Node
HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::inorder_to_node(
    uint32_t inorder) const {
    assert(0 <= inorder && inorder < num_of_nodes);

    auto sample_index = close_sample.select_chunk(inorder);
    int32_t remain = inorder - close_sample.sum(sample_index);
    uint32_t chunk_index = sample_index * W;

    while (true) {
        int32_t p = get_chunk_popcount(chunk_index);
        if (remain - p < 0) {
            break;
        }
        remain -= p;
        chunk_index++;
    }

    Node node(this);
    node.set_chunk_index(chunk_index);
    node.set_tree();
    node.set_local(remain + (node.chunk_index == node.right_chunk_index
                                 ? node.split_rank
                                 : 0));

    return node;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline int32_t
HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::node_to_inorder(
    Node node) const {
    if (!node.is_valid()) return -1;
    auto chunk_index = node.chunk_index;
    assert(0 <= chunk_index && chunk_index < num_of_chunks);

    uint32_t res = close_sample.sum(chunk_index / W);
    for (int i = chunk_index / W * W; i < chunk_index; i++) {
        res += get_chunk_popcount(i);
    }
    res += node.local_inorder;
    if (node.chunk_index == node.right_chunk_index) {
        res -= node.split_rank;
    }

    return res;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline
    typename HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::Node
    HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::root()
        const {
    Node ret(this, 0, 0, false);
    ret.close();
    return ret;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline
    typename HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::Node
    HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::parent(
        Node v) const {
    v.open();
    v.dec();
    if (v.access() == 0) {
        v.close();
    }
    return v;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline
    typename HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::Node
    HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::left_child(
        Node v) const {
    v.open();
    v.inc();
    if (v.access() == 1) {
        return Node(this);
    }
    v.close();
    return v;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline
    typename HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::Node
    HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::right_child(
        Node v) const {
    v.inc();
    if (v.access() == 1) {
        return Node(this);
    }
    v.close();
    return v;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline bool
HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::is_leaf(
    Node v) const {
    auto u = v;
    u.dec();
    v.inc();
    return u.access() != 1 && v.access() != 0;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline int32_t
HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::child_label(
    Node v) const {
    v.open();
    v.dec();
    return v.access();
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline int32_t
HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::subtree_size(
    Node v) const {
    if (!v.is_valid()) return -1;
    auto u = v;
    u.open();
    v.rightmost_desc();
    // node_to_inorder works as a rank if u.access() == 0, i.e., '('.
    return (node_to_inorder(v) - node_to_inorder(u) + 1);
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline bool
HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::is_ancestor(
    Node u, Node v) const {
    return leftmost_desc(u) <= v && v <= rightmost_desc(u);
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline
    typename HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::Node
    HypersuccinctBinaryTree<
        W, CompressedMicrotreeSplitRankArray>::leftmost_desc(Node v) const {
    // Naive approach:
    v.open();
    while (v.access() == 0) {
        v.inc();
    }
    return v;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline
    typename HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::Node
    HypersuccinctBinaryTree<
        W, CompressedMicrotreeSplitRankArray>::rightmost_desc(Node v) const {
    v.rightmost_desc();
    return v;
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline
    typename HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::Node
    HypersuccinctBinaryTree<W, CompressedMicrotreeSplitRankArray>::lca(
        Node u, Node v) const {
    if (u.right_chunk_index == v.right_chunk_index) {
        int32_t lca_local_inorder =
            u.microtree_ptr->naive_lca(u.local_inorder, v.local_inorder);
        u.set_local(lca_local_inorder, true, true);
        return u;
    }

    auto [rmq_first, rmq_last] =
        std::minmax(u.right_chunk_index, v.right_chunk_index);
    auto lca_right_chunk_index = rmm_tree.rmq(rmq_first, rmq_last + 1);
    if (lca_right_chunk_index == u.right_chunk_index) {
        int32_t v_anc_inorder =
            u.split_rank - u.microtree_ptr->bp.get(u.cut - 1);
        int32_t lca_local_inorder =
            u.microtree_ptr->naive_lca(u.local_inorder, v_anc_inorder);
        u.set_local(lca_local_inorder, true, true);
        return u;
    } else if (lca_right_chunk_index == v.right_chunk_index) {
        int32_t u_anc_inorder =
            v.split_rank - v.microtree_ptr->bp.get(v.cut - 1);
        int32_t lca_local_inorder =
            v.microtree_ptr->naive_lca(u_anc_inorder, v.local_inorder);
        v.set_local(lca_local_inorder, true, true);
        return v;
    } else {
        return Node(this, lca_right_chunk_index, 0);
    }
}

template <uint32_t W, typename CompressedMicrotreeSplitRankArray>
inline uint64_t HypersuccinctBinaryTree<
    W, CompressedMicrotreeSplitRankArray>::evaluate_memory_consumption() const {
    return rmm_tree.evaluate_memory_consumption() +
           compressed_microtree_split_rank_array.evaluate_memory_consumption() +
           close_sample.evaluate_memory_consumption();
}

}  // namespace average_case_optimal_rmq
