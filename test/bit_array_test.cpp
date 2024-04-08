#include "hyperrmq/bit_array.hpp"

#include <random>
#include <set>

#include "gtest/gtest.h"

namespace average_case_optimal_rmq {
namespace {

TEST(BitArrayTest, Small) {
    BitArray a(3);
    a.set(0, 1);
    a.set(1, 0);
    a.set(2, 1);
    ASSERT_EQ(a.get(0), 1);
    ASSERT_EQ(a.get(1), 0);
    ASSERT_EQ(a.get(2), 1);

    ASSERT_EQ(a.read_bits(0, 0), 0b0);
    ASSERT_EQ(a.read_bits(0, 1), 0b1);
    ASSERT_EQ(a.read_bits(0, 2), 0b10);
    ASSERT_EQ(a.read_bits(0, 3), 0b101);

    ASSERT_EQ(a.read_bits(1, 0), 0b0);
    ASSERT_EQ(a.read_bits(1, 1), 0b0);
    ASSERT_EQ(a.read_bits(1, 2), 0b01);

    ASSERT_EQ(a.read_bits(2, 0), 0b0);
    ASSERT_EQ(a.read_bits(2, 1), 0b1);

    ASSERT_EQ(a.read_bits(3, 0), 0b0);

    ASSERT_EQ(a.size(), 3);
    ASSERT_EQ(a.linear_popcount(), 2);

    ASSERT_EQ(a.linear_select1(0), 0);
    ASSERT_EQ(a.linear_select1(1), 2);
    ASSERT_EQ(a.linear_select1(2), 3);

    ASSERT_EQ(a.read_interval(0, 3), a);
}

TEST(BitArrayTest, ZeroLength) {
    BitArray a;
    ASSERT_EQ(a.read_bits(0, 0), 0b0);
    for (int i = 0; i <= 64; i++) {
        ASSERT_EQ(a.read_bits_zero_follow(0, i), 0b0);
    }

    ASSERT_EQ(a.size(), 0);
    ASSERT_EQ(a.linear_popcount(), 0);
    ASSERT_EQ(a.linear_select1(0), 0);
    ASSERT_EQ(a.read_interval(0, 0), a);
}

TEST(BitArrayTest, BuildFromVector) {
    BitArray a(3, {0b101ull << 61});
    ASSERT_EQ(a.get(0), 1);
    ASSERT_EQ(a.get(1), 0);
    ASSERT_EQ(a.get(2), 1);

    ASSERT_EQ(a.read_bits(0, 1), 0b1);
    ASSERT_EQ(a.read_bits(0, 2), 0b10);
    ASSERT_EQ(a.read_bits(0, 3), 0b101);

    ASSERT_EQ(a.read_bits(1, 1), 0b0);
    ASSERT_EQ(a.read_bits(1, 2), 0b01);

    ASSERT_EQ(a.read_bits(2, 1), 0b1);

    ASSERT_EQ(a.size(), 3);
    ASSERT_EQ(a.linear_popcount(), 2);
}

#ifdef DEBUG
TEST(BitArrayDeathTest, BuildFromVectorOutOfRange) {
    ASSERT_DEATH_IF_SUPPORTED({ BitArray(3, {-1ull}); },
                              "Assertion .* failed\\.");
    ASSERT_DEATH_IF_SUPPORTED({ BitArray(3, {1ull << 60}); },
                              "Assertion .* failed\\.");
}
#endif

TEST(BitArrayTest, BuildFromBP) {
    BitArray a(")()");
    ASSERT_EQ(a.get(0), 1);
    ASSERT_EQ(a.get(1), 0);
    ASSERT_EQ(a.get(2), 1);

    ASSERT_EQ(a.read_bits(0, 1), 0b1);
    ASSERT_EQ(a.read_bits(0, 2), 0b10);
    ASSERT_EQ(a.read_bits(0, 3), 0b101);

    ASSERT_EQ(a.read_bits(1, 1), 0b0);
    ASSERT_EQ(a.read_bits(1, 2), 0b01);

    ASSERT_EQ(a.read_bits(2, 1), 0b1);

    ASSERT_EQ(a.size(), 3);
    ASSERT_EQ(a.linear_popcount(), 2);
}

TEST(BitArrayTest, StressTest) {
    const int n = 100;
    const int q = 100000;
    BitArray a(n);
    std::vector<int> expected(n);
    std::mt19937 mt(0);
    for (int i = 0; i < n; i++) {
        int val = mt() % 2;
        a.set(i, val);
        expected[i] = val;
    }
    for (int i = 0; i < q; i++) {
        uint64_t idx = mt() % n;
        int query_type = mt() % 3;
        if (query_type == 0) {
            int val = mt() % 2;
            a.set(idx, val);
            expected[idx] = val;
        } else if (query_type == 1) {
            int val = mt() % 2;
            if (val) {
                a.on(idx);
            } else {
                a.off(idx);
            }
            expected[idx] = val;
        } else {
            ASSERT_EQ(a.get(idx), expected[idx]);
        }
    }
}

TEST(BitArrayTest, LinearPopCount) {
    const int n = 1000;
    BitArray ba(n);
    std::mt19937 mt(0);
    for (int i = 0; i < n; i++) {
        ba.set(i, mt() % 2);
    }
    for (int l = 0; l <= n; l++) {
        ASSERT_EQ(ba.linear_popcount(l, 0), 0);

        uint32_t expected = 0;
        for (int r = l; r < n; r++) {
            expected += ba.get(r);
            ASSERT_EQ(expected, ba.linear_popcount(l, r - l + 1));
        }
    }

    ASSERT_EQ(ba.linear_popcount(), ba.linear_popcount(0, ba.size()));
}

TEST(BitArrayTest, SelectRankIdentity) {
    const int n = 1000;
    BitArray ba(n);
    std::mt19937 mt(0);
    for (int i = 0; i < n; i++) {
        ba.set(i, mt() % 2);
    }

    for (int i = 0; i < ba.size(); i++) {
        if (ba.get(i)) {
            ASSERT_EQ(i, ba.linear_select1(ba.linear_popcount(0, i)));
        }
    }
}

TEST(BitArrayTest, ToString) {
    BitArray a(4);
    a.set(0, 0);
    a.set(1, 1);
    a.set(2, 0);
    a.set(3, 1);
    ASSERT_EQ(a.to_string(), "0101");
}

TEST(BitArrayTest, EqualityAndInequality) {
    BitArray a(2), b(2), c(2), d(1);
    a.set(0, 1);
    a.set(1, 0);

    b.set(0, 1);
    b.set(1, 0);
    ASSERT_TRUE(a == b);
    ASSERT_FALSE(a != b);
    ASSERT_FALSE(a < b);
    ASSERT_FALSE(b < a);

    c.set(0, 0);
    c.set(1, 0);
    ASSERT_FALSE(a == c);
    ASSERT_TRUE(a != c);
    ASSERT_FALSE(a < c);
    ASSERT_TRUE(c < a);

    d.set(0, 1);
    ASSERT_FALSE(a == d);
    ASSERT_TRUE(a != d);
    ASSERT_FALSE(a < d);
    ASSERT_TRUE(d < a);

    std::set<BitArray> s = {a, b, c, d};
    ASSERT_EQ(s.size(), 3);
}

TEST(BitArrayTest, ReadWriteBitsStressTest) {
    const int n = 1000;
    const int q = 10000;
    BitArray actual(n), expected(n);
    std::mt19937_64 mt(0);
    for (int i = 0; i < q; i++) {
        int query_type = mt() % 2;
        uint64_t width = mt() % 65;
        uint64_t start = mt() % (n - width + 1);
        if (query_type == 0) {
            uint64_t bits = mt();
            if (width < 64) {
                bits &= (1ull << width) - 1;
            }
            for (int i = 0; i < width; i++) {
                expected.set(start + i, (bits >> (width - i - 1)) & 1);
            }
            actual.write_bits(start, width, bits);
        } else {
            uint64_t expected_res = 0;
            for (int i = 0; i < width; i++) {
                expected_res |= expected.get(start + i) << (width - i - 1);
            }
            ASSERT_EQ(actual.read_bits(start, width), expected_res);
            ASSERT_EQ(actual.read_bits_zero_follow(start, width), expected_res);
        }
        ASSERT_EQ(expected, actual);
    }
}

TEST(BitArrayTest, ReadWriteBitsStressTest64) {
    const int n = 1000;
    const int q = 10000;
    BitArray actual(n), expected(n);
    std::mt19937_64 mt(0);
    for (int i = 0; i < q; i++) {
        int query_type = mt() % 2;
        uint64_t width = 64;
        uint64_t start = mt() % (n - width + 1);
        if (query_type == 0) {
            uint64_t bits = mt();
            for (int i = 0; i < width; i++) {
                expected.set(start + i, (bits >> (width - i - 1)) & 1);
            }
            actual.write_bits(start, width, bits);
        } else {
            uint64_t expected_res = 0;
            for (int i = 0; i < width; i++) {
                expected_res |= expected.get(start + i) << (width - i - 1);
            }
            ASSERT_EQ(actual.read_bits(start, width), expected_res);
            ASSERT_EQ(actual.read_bits_zero_follow(start, width), expected_res);
        }
        ASSERT_EQ(expected, actual);
    }
}

TEST(BitArrayTest, ReadBitsZeroFollowStressTest) {
    BitArray ba(32), padded(96);
    const int t = 10;
    std::mt19937 mt(0);
    for (int i = 0; i < t; i++) {
        uint64_t bits = mt();
        ba.write_bits(0, 32, bits);
        padded.write_bits(0, 32, bits);

        for (int first = 0; first < 32; first++) {
            for (int width = 1; width <= 64; width++) {
                ASSERT_EQ(ba.read_bits_zero_follow(first, width),
                          padded.read_bits(first, width));
            }
        }
    }
}

TEST(BitArrayTest, ReadWriteIntervalStressTest) {
    const int n = 1000;
    const int q = 10000;
    BitArray actual(n), expected(n);
    std::mt19937_64 mt(0);
    for (int i = 0; i < q; i++) {
        int query_type = mt() % 2;
        uint64_t width = mt() % (n + 1);
        uint64_t start = mt() % (n - width + 1);
        if (query_type == 0) {
            BitArray interval(width);
            for (int i = 0; i < width; i++) {
                interval.set(i, mt() % 2);
            }

            for (int i = 0; i < width; i++) {
                expected.set(start + i, interval.get(i));
            }
            actual.write_interval(start, interval);
        } else {
            BitArray expected_res(width);
            for (int i = 0; i < width; i++) {
                expected_res.set(i, expected.get(start + i));
            }
            ASSERT_EQ(actual.read_interval(start, width), expected_res);
        }
        ASSERT_EQ(expected, actual);
    }
}

}  // namespace
}  // namespace average_case_optimal_rmq
