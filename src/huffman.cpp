#include "huffman.hpp"

#include <algorithm>
#include <cassert>

#include "memutil.hpp"
#include "microtree_array.hpp"
#include "minimal_cell_array.hpp"

namespace average_case_optimal_rmq {

template <typename T>
HuffmanTree<T>::HuffmanTree(const std::map<T, uint64_t> &alphabet_count)
    : number_of_alphabet(alphabet_count.size()) {
    std::vector<NodePtr> leaves;
    leaves.reserve(number_of_alphabet);
    for (auto &&[key, freq] : alphabet_count) {
        leaves.push_back(std::make_shared<Node>(key, freq));
    }
    std::sort(leaves.begin(), leaves.end(),
              [](const NodePtr &a, const NodePtr &b) -> bool {
                  return a->freq < b->freq;
              });
    for (int i = 0; i < number_of_alphabet - 1; i++) {
        leaves[i]->next = leaves[i + 1];
    }
    NodePtr insert_pos = leaves[0];
    NodePtr look_pos = leaves[0];
    while (look_pos->next != nullptr) {
        NodePtr n =
            std::make_shared<Node>(look_pos->freq + look_pos->next->freq);
        n->left = look_pos;
        n->right = look_pos->next;
        while (insert_pos->next != nullptr &&
               insert_pos->next->freq <= n->freq) {
            insert_pos = insert_pos->next;
        }
        n->next = insert_pos->next;
        insert_pos->next = n;
        look_pos = look_pos->next->next;
    }
    root = look_pos;

    auto delete_next = [](auto self, NodePtr n) -> void {
        n->next = nullptr;
        if (n->left != nullptr) {
            self(self, n->left);
        }
        if (n->right != nullptr) {
            self(self, n->right);
        }
    };
    delete_next(delete_next, root);
}

template <typename T>
std::vector<std::pair<T, uint64_t>> HuffmanTree<T>::compute_length() {
    std::vector<std::pair<T, uint64_t>> res;
    res.reserve(number_of_alphabet);

    auto dfs = [&](auto self, NodePtr n, uint64_t d = 0) -> void {
        if (n->left == nullptr) {
            res.push_back({n->key, d});
        } else {
            self(self, n->left, d + 1);
            self(self, n->right, d + 1);
        }
    };
    dfs(dfs, root);

    return res;
}

template <typename T, typename AlphabetArray>
CanonicalHuffmanCode<T, AlphabetArray>::CanonicalHuffmanCode()
    : number_of_alphabet(0) {}

template <typename T, typename AlphabetArray>
CanonicalHuffmanCode<T, AlphabetArray>::CanonicalHuffmanCode(
    const std::map<T, uint64_t> &alphabet_count)
    : number_of_alphabet(alphabet_count.size()) {
    HuffmanTree<T> ht(alphabet_count);
    auto alphabet_and_code_lengths = ht.compute_length();
    sort(alphabet_and_code_lengths.begin(), alphabet_and_code_lengths.end(),
         [](const std::pair<T, uint64_t> &a, const std::pair<T, uint64_t> &b)
             -> bool { return a.second < b.second; });

    std::vector<T> alphabets_vec;
    alphabets_vec.reserve(number_of_alphabet);
    for (auto &&[alphabet, code_length] : alphabet_and_code_lengths) {
        alphabets_vec.push_back(alphabet);
    }
    alphabets = std::move(alphabets_vec);

    maximum_code_length = alphabet_and_code_lengths.back().second;
    assert(maximum_code_length <= 64);

    std::vector<uint64_t> codes;
    codes.reserve(number_of_alphabet);
    codes.push_back(0);
    for (int i = 1; i < number_of_alphabet; i++) {
        codes.push_back((codes.back() + 1)
                        << (alphabet_and_code_lengths[i].second -
                            alphabet_and_code_lengths[i - 1].second));
    }

    first_index.resize(maximum_code_length + 1, number_of_alphabet);
    first_code.resize(maximum_code_length + 1);
    for (int i = number_of_alphabet - 1; i >= 0; i--) {
        first_index[alphabet_and_code_lengths[i].second] = i;
        first_code[alphabet_and_code_lengths[i].second] = codes[i];
    }
    for (int l = maximum_code_length - 1; l >= 1; l--) {
        if (first_index[l] == number_of_alphabet) {
            first_index[l] = first_index[l + 1];
            first_code[l] = first_code[l + 1] >> 1;
        }
    }
}

template <typename T, typename AlphabetArray>
std::map<T, std::pair<uint64_t, uint64_t>>
CanonicalHuffmanCode<T, AlphabetArray>::enumerate_alphabet_code_pair() const {
    std::map<T, std::pair<uint64_t, uint64_t>> result;
    uint64_t code_length = 0, code = 0;
    for (int i = 0; i < number_of_alphabet; i++) {
        while (code_length < maximum_code_length &&
               first_code[code_length + 1] <= (code << 1)) {
            code_length++;
            code <<= 1;
        }
        result[alphabets[i]] = {code_length, code};
        code++;
    }
    return result;
}

template <typename T, typename AlphabetArray>
uint64_t CanonicalHuffmanCode<T, AlphabetArray>::get_next_length(
    const BitArray &code_seq, uint64_t index) const {
    uint64_t next = code_seq.read_bits_zero_follow(index, maximum_code_length);
    int l = 0, r = maximum_code_length + 1;
    int m;
    while (r - l > 1) {
        m = (l + r) / 2;
        if ((first_code[m] << (maximum_code_length - m)) <= next) {
            l = m;
        } else {
            r = m;
        }
    }
    return l;
}

template <typename T, typename AlphabetArray>
T CanonicalHuffmanCode<T, AlphabetArray>::decode(uint64_t code_length,
                                                 uint64_t code) const {
    return alphabets[first_index[code_length] + code - first_code[code_length]];
}

template <typename T, typename AlphabetArray>
uint64_t CanonicalHuffmanCode<T, AlphabetArray>::evaluate_memory_consumption()
    const {
    uint64_t alphabet_size = 0;
    if constexpr (std::is_same_v<AlphabetArray, std::vector<T>>) {
        alphabet_size = evaluate_vector_memory_consumption(alphabets);
    } else {
        alphabet_size = alphabets.evaluate_memory_consumption();
    }
    return alphabet_size + evaluate_vector_memory_consumption(first_code) +
           evaluate_vector_memory_consumption(first_index);
}

template struct HuffmanTree<uint64_t>;
template struct HuffmanTree<TreeBP>;
template struct HuffmanTree<std::pair<TreeBP, uint32_t>>;

template struct CanonicalHuffmanCode<uint64_t, std::vector<uint64_t>>;
template struct CanonicalHuffmanCode<uint64_t, MinimalCellArray<uint64_t>>;
template struct CanonicalHuffmanCode<TreeBP, EditableMicrotreeArray>;
template struct CanonicalHuffmanCode<std::pair<TreeBP, uint32_t>,
                                     MicrotreeSplitRankArray>;

}  // namespace average_case_optimal_rmq
