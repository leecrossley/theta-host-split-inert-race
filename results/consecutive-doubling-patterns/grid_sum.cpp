// Evaluate the finite Fourier sums with signed fixed-point weights and modular phases.
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

struct Mode {
    std::uint64_t half_step;
    std::int64_t weight;
};

int main() {
    constexpr unsigned bits = 15;
    constexpr std::int64_t scale = std::int64_t{1} << 24;
    std::vector<std::int64_t> sine(1U << bits);
    for (auto& s : sine) {
        if (!(std::cin >> s) || std::abs(s) > scale) return 1;
    }
    unsigned count = 0;
    if (!(std::cin >> count) || count > 20000) return 2;
    std::vector<Mode> modes(count);
    std::int64_t absolute_sum = 0;
    for (auto& mode : modes) {
        if (!(std::cin >> mode.half_step >> mode.weight) || std::abs(mode.weight) > 100*scale) return 3;
        absolute_sum += std::abs(mode.weight);
    }
    if (absolute_sum > std::numeric_limits<std::int64_t>::max()/scale) return 4;
    unsigned length;
    std::uint64_t code, exponent;
    while (std::cin >> length >> code >> exponent) {
        if (length == 0 || length > 12 || exponent > (std::uint64_t{1} << 40)) return 5;
        std::vector<std::int64_t> total(length, 0);
        for (const auto& mode : modes) {
            auto phase = mode.half_step*(2*exponent+1);
            for (auto& x : total) {
                const auto k = (phase + (std::uint64_t{1} << (63-bits))) >> (64-bits);
                x += mode.weight*sine[k];
                phase += 2*mode.half_step;
            }
        }
        std::cout << length << ' ' << code << ' ' << exponent;
        for (const auto x : total) std::cout << ' ' << x;
        std::cout << '\n';
    }
    return std::cin.eof() ? 0 : 6;
}
