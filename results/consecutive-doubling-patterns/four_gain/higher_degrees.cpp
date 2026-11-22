// Reuse the two published fixed-point logarithm enclosures at rational primes p <= 2^29.
#define THETA_REUSE_LOG_BOUNDS
#if defined(INDEPENDENT)
#include "../word_sieve_recheck.cpp"
#else
#include "../word_sieve.cpp"
#include <primesieve/iterator.hpp>
#endif

constexpr std::uint64_t norm_limit = std::uint64_t(1) << 58;
constexpr std::uint64_t root_limit = std::uint64_t(1) << 29;

int main(int argc, char **argv) {
    using Num = std::uint64_t;
    if (argc != 2)
        return 2;
    std::string host = argv[1];
    if (host != "547" && host != "571")
        return 2;
    unsigned u = host == "547" ? 547 : 571, v = host == "547" ? 727 : 599;
#if defined(INDEPENDENT)
    auto left = character(u), right = character(v);
    std::vector<U> mesh(subdivisions + 1);
    for (auto &val : mesh)
        if (!(std::cin >> val))
            return 2;
    std::vector<Bin> output(59);
#else
    auto left = legendre(u), right = legendre(v);
    Sums output;
    output.logTable.resize(bins + 1);
    output.sums.resize(59);
    for (auto &val : output.logTable)
        if (!(std::cin >> val))
            return 2;
#endif
    for (unsigned i = 0; i <= 58; ++i) {
        Num val;
        if (!(std::cin >> val) || val != (Num(1) << i))
            return 2;
    }
    auto emit = [&](Num p) {
        if (p == 2 || p == u || p == v)
            return;
        const Num square = p * p;
        Num n;
#if defined(INDEPENDENT)
        const unsigned ord16 =
            (p % 16 == 1 || p % 16 == 15) ? 1 : ((p % 16 == 7 || p % 16 == 9) ? 2 : 4);
        const unsigned ord2 = left[p % u] * right[p % v] == 1 ? 1 : 2;
        unsigned f = std::max(ord16, ord2);
        if (f == 1 || (f == 4 && square > norm_limit / square))
            return;
        n = f == 2 ? square : square * square;
#else
        if ((p % 16 == 1 || p % 16 == 15) && left[p % u] == right[p % v])
            return;
        if (square % 16 == 1 || square % 16 == 15)
            n = square;
        else {
            if (square > norm_limit / square)
                return;
            n = square * square;
        }
#endif
        unsigned bucket = 64u - static_cast<unsigned>(__builtin_clzll(n - 1));
#if defined(INDEPENDENT)
        auto bounds = logBounds(p, mesh);
        Bin &out = output[bucket];
        out.low -= Signed(bounds[1]);
        out.high -= Signed(bounds[0]);
        out.reciprocalLow -= Signed(((bounds[1] << 40) + n - 1) / n);
        out.reciprocalHigh -= Signed((bounds[0] << 40) / n);
#else
        output.add(bucket, p, n, -1);
#endif
    };
#if defined(INDEPENDENT)
    std::vector<unsigned> base;
    for (unsigned p = 3; Num(p) * p <= root_limit; p += 2) {
        bool ok = true;
        for (Num q : base) {
            if (q * q > p)
                break;
            if (p % q == 0) {
                ok = false;
                break;
            }
        }
        if (ok)
            base.push_back(p);
    }
    constexpr Num step = 1u << 20;
    std::vector<char> flags(step / 2);
    for (Num l = 3; l <= root_limit; l += step) {
        Num last = std::min(l + step - 1, root_limit);
        std::size_t size = (last - l) / 2 + 1;
        std::fill_n(flags.begin(), size, 1);
        for (Num p : base) {
            if (p * p > last)
                break;
            Num start = std::max(p * p, ((l + p - 1) / p) * p);
            if (start % 2 == 0)
                start += p;
            for (Num j = (start - l) / 2; j < size; j += p)
                flags[j] = 0;
        }
        for (std::size_t j = 0; j < size; ++j)
            if (flags[j])
                emit(l + 2 * j);
    }
#else
    primesieve::iterator it(3, root_limit);
    for (Num p = it.next_prime(); p <= root_limit; p = it.next_prime())
        emit(p);
#endif
    for (unsigned i = 0; i <= 58; ++i) {
#if defined(INDEPENDENT)
        const auto &v = output[i];
        std::cout << (Num(1) << i) << ' ' << number(v.low) << ' ' << number(v.high) << ' '
                  << number(v.reciprocalLow) << ' ' << number(v.reciprocalHigh) << '\n';
#else
        const auto &v = output.sums[i];
        std::cout << (Num(1) << i) << ' ' << decimal(v.a) << ' ' << decimal(v.b) << ' '
                  << decimal(v.c) << ' ' << decimal(v.d) << '\n';
#endif
    }
}
