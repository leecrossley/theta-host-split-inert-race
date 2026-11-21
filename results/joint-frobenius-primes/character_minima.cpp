#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <vector>

// Write the least partial sum on every integer block, for each primitive
// nonprincipal real character inducing a character of the joint modulus.
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    constexpr int limit = 1 << 26;
    constexpr std::array<int, 4> primes = {547, 727, 571, 599};
    constexpr std::array<int, 8> minus_four = {0, 1, 0, -1, 0, 1, 0, -1};
    constexpr std::array<int, 8> plus_eight = {0, 1, 0, -1, 0, -1, 0, 1};
    std::ofstream output(argv[1]);
    if (!output) return 3;
    std::vector<int> cuts = {1};
    while (cuts.back() < limit) {
        int left = cuts.back();
        cuts.push_back(std::min(limit, left + std::max(1, left / 192)));
    }
    std::array<std::vector<std::int8_t>, 4> symbols;
    for (int i = 0; i < 4; ++i) {
        int p = primes[i];
        symbols[i].assign(p, -1);
        symbols[i][0] = 0;
        for (int a = 1; a < p; ++a) symbols[i][a * a % p] = 1;
    }
    for (int mask = 0; mask < 16; ++mask) {
        long long odd_discriminant = 1;
        std::vector<std::int8_t> values(limit, 1);
        for (int i = 0; i < 4; ++i) {
            if ((mask & (1 << i)) == 0) continue;
            int p = primes[i];
            odd_discriminant *= -p;
            for (int start = 0; start < limit; start += p) {
                int count = std::min(p, limit - start);
                for (int k = 0; k < count; ++k) values[start + k] *= symbols[i][k];
            }
        }
        for (int d : {1, -4, 8, -8}) {
            if (mask == 0 && d == 1) continue;
            std::array<int, 8> two_part;
            two_part.fill(1);
            for (int j = 0; j < 8; ++j) {
                if (d == -4 || d == -8) two_part[j] *= minus_four[j];
                if (d == 8 || d == -8) two_part[j] *= plus_eight[j];
            }
            int sum = 0;
            output << odd_discriminant * d;
            for (std::size_t i = 1; i < cuts.size(); ++i) {
                int minimum = limit;
                for (int n = cuts[i - 1]; n < cuts[i]; ++n) {
                    sum += values[n] * two_part[n % 8];
                    minimum = std::min(minimum, sum);
                }
                output << '\t' << minimum;
            }
            output << '\n';
        }
    }
    return output.good() ? 0 : 4;
}
