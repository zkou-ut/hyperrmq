#include "hyperrmq/cartesian_tree.hpp"

#include <stack>

namespace hyperrmq {

CartesianTree::CartesianTree(const std::vector<int32_t> &values) {
    n = values.size();
    std::stack<int32_t> st;
    std::vector<int32_t> parent(n);
    int c;

    for (int i = 0; i < n; i++) {
        c = i;
        while (!st.empty() && values[st.top()] > values[i]) {
            c = st.top();
            st.pop();
        }
        parent[c] = i;
        parent[i] = st.empty() ? -1 : st.top();
        st.push(i);
    }

    binary_tree = BinaryTree(parent);
}

}  // namespace hyperrmq
