#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <vector>

#include "hyperrmq/bit_array.hpp"

namespace hyperrmq {

template <typename T = uint64_t>
struct HuffmanTree {
    struct Node;
    using NodePtr = std::shared_ptr<Node>;
    struct Node {
        T key;
        uint64_t freq;
        NodePtr next, left, right;
        Node(uint64_t freq) : key(), freq(freq), next(), left(), right() {}
        Node(T key, uint64_t freq)
            : key(key), freq(freq), next(), left(), right() {}
    };

    uint64_t number_of_alphabet;
    NodePtr root;

    explicit HuffmanTree(const std::map<T, uint64_t> &alphabet_count);

    std::vector<std::pair<T, uint64_t>> compute_length();
};

// `AlphabetArray` requires:
//     - constructor from std::vector<T>
//     - access by []
template <typename T = uint64_t, typename AlphabetArray = std::vector<T>>
struct CanonicalHuffmanCode {
    uint64_t number_of_alphabet, maximum_code_length;

    CanonicalHuffmanCode();
    explicit CanonicalHuffmanCode(const std::map<T, uint64_t> &alphabet_count);

    std::map<T, std::pair<uint64_t, uint64_t>> enumerate_alphabet_code_pair()
        const;

    uint64_t get_next_length(const BitArray &code_seq, uint64_t index) const;
    T decode(uint64_t code_length, uint64_t code) const;

    uint64_t evaluate_memory_consumption() const;

   private:
    AlphabetArray alphabets;
    std::vector<uint64_t> first_index, first_code;
};

}  // namespace hyperrmq
