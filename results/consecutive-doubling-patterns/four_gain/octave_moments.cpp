// Exact integer moments of prime-ideal norms in one dyadic octave.
#include <primesieve/iterator.hpp>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;
constexpr int bits = 15;
constexpr unsigned cells = 1u << bits;

struct Moment {
    i64 sum = 0;
    i64 displacement = 0;
    u64 absolute = 0;
    void add(int weight, u64 offset) {
        sum += weight;
        displacement += weight * static_cast<i64>(offset);
        absolute += static_cast<unsigned>(std::abs(weight));
    }
};

template <unsigned q> std::array<int, q> character() {
    std::array<int, q> result{};
    result.fill(-1);
    result[0] = 0;
    for (unsigned i = 1; i < q; ++i)
        result[i * i % q] = 1;
    return result;
}

const auto au = character<547>();
const auto av = character<727>();
const auto bu = character<571>();
const auto bv = character<599>();

bool bad_a(u64 p) { return p == 547 || p == 727 || p == 343057911569 || p == 343070636977; }
bool bad_b(u64 p) { return p == 571 || p == 599 || p == 413718278401 || p == 603339156001; }

int main(int argc, char **argv) {
    if (argc != 3)
        return 2;
    unsigned k = std::stoul(argv[1]), threads = std::stoul(argv[2]);
    // Even if every norm is present, the total absolute displacement in one cell
    // is bounded by 4 * width^2. At k=44 this fits signed 64-bit arithmetic.
    if (k < bits || k > 44 || threads == 0 || threads > 64)
        return 2;
    u64 lo = u64(1) << k, hi = 2 * lo, width = lo >> bits;
    std::vector<std::array<Moment, 2>> data(cells);
    std::vector<u64> prime_totals(cells);
    constexpr unsigned batch = 128;
    std::atomic<unsigned> next{0};
    auto worker = [&]() {
        for (;;) {
            unsigned first = next.fetch_add(batch);
            if (first >= cells)
                break;
            u64 start = lo + first * width, end = start + batch * width;
            primesieve::iterator it(start, end - 1);
            for (u64 p = it.next_prime(); p < end; p = it.next_prime()) {
                unsigned cell = static_cast<unsigned>((p - lo) / width);
                ++prime_totals[cell];
                if (((p + 1) & 15) > 2)
                    continue; // odd p is 1 or 15 modulo 16
                u64 offset = (p - lo) & (width - 1);
                int a = au[p % 547], b = bu[p % 571];
                if (a == av[p % 727] && !bad_a(p))
                    data[cell][0].add(-4 * a, offset);
                if (b == bv[p % 599] && !bad_b(p))
                    data[cell][1].add(-4 * b, offset);
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned i = 0; i < threads; ++i)
        pool.emplace_back(worker);
    for (auto &th : pool)
        th.join();

    primesieve::iterator small(3);
    for (u64 p = small.next_prime(); p <= (hi - 1) / p; p = small.next_prime()) {
        u64 square = p * p;
        bool class1 = (p % 16 == 1 || p % 16 == 15);
        bool class2 = (square % 16 == 1 || square % 16 == 15);
        auto insert = [&](unsigned field, bool omitted, bool equal) {
            if (omitted || (class1 && equal))
                return;
            if (!class2 && square > (hi - 1) / square)
                return;
            unsigned degree = class2 ? 2 : 4;
            u64 norm = class2 ? square : square * square;
            if (norm < lo || norm >= hi)
                return;
            unsigned cell = static_cast<unsigned>((norm - lo) / width);
            data[cell][field].add(-4 / static_cast<int>(degree), (norm - lo) & (width - 1));
        };
        insert(0, bad_a(p), au[p % 547] == av[p % 727]);
        insert(1, bad_b(p), bu[p % 571] == bv[p % 599]);
    }

    std::cout << "# octave " << k << " cells " << cells << '\n';
    for (unsigned j = 0; j < cells; ++j) {
        std::cout << j << ' ' << prime_totals[j];
        for (const auto &m : data[j])
            std::cout << ' ' << m.sum << ' ' << m.displacement << ' ' << m.absolute;
        std::cout << '\n';
    }
}
