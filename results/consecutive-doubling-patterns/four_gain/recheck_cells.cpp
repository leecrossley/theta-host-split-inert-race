// Independently sieve selected cells using a segmented sieve over every odd integer.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>
using U = std::uint64_t;
using I = std::int64_t;

U power(U base, U exponent, U modulus) {
    U res = 1;
    for (; exponent; exponent >>= 1, base = base * base % modulus)
        if (exponent & 1)
            res = res * base % modulus;
    return res;
}

std::vector<int> signs(unsigned q) {
    std::vector<int> out(q);
    for (unsigned v = 0; v < q; ++v) {
        U e = power(v, (q - 1) / 2, q);
        out[v] = e == q - 1 ? -1 : static_cast<int>(e);
    }
    return out;
}

int main() {
    std::vector<std::pair<unsigned, unsigned>> cases;
    unsigned k, j;
    U max_high = 0;
    while (std::cin >> k >> j) {
        if (k < 15 || k > 44 || j >= (1u << 15))
            return 2;
        cases.emplace_back(k, j);
        max_high = std::max(max_high, (U(1) << k) + (U(j + 1) << (k - 15)));
    }
    if (cases.empty())
        return 2;
    unsigned root = static_cast<unsigned>(std::sqrt(static_cast<double>(max_high)));
    while (U(root + 1) * (root + 1) <= max_high)
        ++root;
    while (U(root) * root > max_high)
        --root;
    std::vector<bool> prime(root + 1, true);
    prime[0] = prime[1] = false;
    std::vector<unsigned> primes;
    for (unsigned p = 2; p <= root; ++p)
        if (prime[p]) {
            primes.push_back(p);
            for (U q = U(p) * p; q <= root; q += p)
                prime[q] = false;
        }
    std::array<unsigned, 2> u{547, 571}, v{727, 599};
    std::array<std::vector<int>, 2> su{signs(u[0]), signs(u[1])}, sv{signs(v[0]), signs(v[1])};
    const std::array<std::array<U, 2>, 2> barred{
        {{343057911569, 343070636977}, {413718278401, 603339156001}}};
    for (auto [exponent, cell] : cases) {
        U width = U(1) << (exponent - 15);
        U left = (U(1) << exponent) + cell * width, end = left + width;
        std::array<std::array<I, 3>, 2> moments{};
        U number = 0;
        auto add = [&](U p, bool first) {
            if (p == 2)
                return;
            for (unsigned f = 0; f < 2; ++f) {
                if (p == u[f] || p == v[f] || p == barred[f][0] || p == barred[f][1])
                    continue;
                int a = su[f][p % u[f]], b = sv[f][p % v[f]];
                unsigned degree = 1;
                while (degree <= 4) {
                    U residue = power(p % 16, degree, 16);
                    if ((residue == 1 || residue == 15) && (degree % 2 == 0 || a * b == 1))
                        break;
                    degree *= 2;
                }
                if (degree > 4 || (degree == 1) != first)
                    continue;
                U norm = 1;
                bool outside = false;
                for (unsigned e = 0; e < degree; ++e) {
                    if (norm > (end - 1) / p) {
                        outside = true;
                        break;
                    }
                    norm *= p;
                }
                if (outside || norm < left)
                    continue;
                I w = degree == 1 ? -4 * a : -4 / static_cast<I>(degree);
                moments[f][0] += w;
                moments[f][1] += w * static_cast<I>(norm - left);
                moments[f][2] += std::abs(w);
            }
        };
        for (U p : primes) {
            if (p * p >= end)
                break;
            add(p, false);
        }
        constexpr U chunk = U(1) << 20;
        std::vector<unsigned char> mask(chunk / 2);
        for (U start = left; start < end; start += chunk) {
            U limit = std::min(end, start + chunk), odd = start | 1;
            auto length = static_cast<std::size_t>((limit - odd + 1) / 2);
            std::fill_n(mask.begin(), length, 1);
            for (U p : primes) {
                if (p == 2)
                    continue;
                if (p * p >= limit)
                    break;
                U first = std::max(p * p, ((odd + p - 1) / p) * p);
                if (!(first % 2))
                    first += p;
                for (U n = first; n < limit; n += 2 * p)
                    mask[(n - odd) / 2] = 0;
            }
            for (std::size_t index = 0; index < length; ++index)
                if (mask[index]) {
                    U p = odd + 2 * index;
                    ++number;
                    if (p % 16 == 1 || p % 16 == 15)
                        add(p, true);
                }
        }
        std::cout << exponent << ' ' << cell << ' ' << number;
        for (auto &row : moments)
            for (I val : row)
                std::cout << ' ' << val;
        std::cout << '\n';
    }
}
