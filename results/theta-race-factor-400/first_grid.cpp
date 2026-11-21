// Integer Fourier evaluation for the full grid of candidate prime-race starts.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

using Integer = std::int64_t;
using Wide = __int128;
#ifndef GRID_SIZE
#define GRID_SIZE 32768
#endif
constexpr int size = GRID_SIZE;
static_assert(size == 32768 || size == 16384);

constexpr int degree = 16;
constexpr int steps_per_unit = 8;
constexpr int last_block = 1620550*steps_per_unit/size;
constexpr Integer scale = Integer{1} << 48;

struct Complex {
    Integer real = 0;
    Integer imag = 0;
};

Complex multiply(Complex a, Complex b) {
    Wide real = Wide{a.real} * b.real - Wide{a.imag} * b.imag;
    Wide imag = Wide{a.real} * b.imag + Wide{a.imag} * b.real;
    return {static_cast<Integer>(real / scale), static_cast<Integer>(imag / scale)};
}

struct Frequency {
    int bin;
    Complex advance;
    Complex phase{scale, 0};
    std::array<Complex, degree + 1> coefficients;
};

// Unnormalized transform with a positive exponential. Each division in the
// butterfly multiplication truncates toward zero and is charged in the proof.
void transform(std::vector<Complex>& values, const std::vector<Complex>& roots) {
    if constexpr (size == 16384) {
        // Decimation in frequency. The unnormalized positive transform
        // emerges in bit-reversed order, which the final permutation restores.
        for (int length = size; length >= 2; length /= 2) {
            const int half = length/2;
            const int stride = size/length;
            for (int start = 0; start < size; start += length) {
                for (int j = 0; j < half; ++j) {
                    const Complex u = values[start+j], v = values[start+j+half];
                    values[start+j] = {u.real+v.real, u.imag+v.imag};
                    values[start+j+half] = multiply(roots[j*stride],
                                                     {u.real-v.real, u.imag-v.imag});
                }
            }
        }
        for (int i = 1, j = 0; i < size; ++i) {
            int bit = size >> 1;
            for (; (j & bit) != 0; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) std::swap(values[i], values[j]);
        }
        return;
    }
    for (int i = 1, j = 0; i < size; ++i) {
        int bit = size >> 1;
        for (; (j & bit) != 0; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(values[i], values[j]);
    }
    for (int length = 2; length <= size; length *= 2) {
        const int half = length / 2;
        const int stride = size / length;
        for (int start = 0; start < size; start += length) {
            for (int j = 0; j < half; ++j) {
                const Complex u = values[start+j];
                const Complex v = multiply(roots[j*stride], values[start+j+half]);
                values[start+j] = {u.real+v.real, u.imag+v.imag};
                values[start+j+half] = {u.real-v.real, u.imag-v.imag};
            }
        }
    }
}

struct Run {
    int first = -1, last = -1;
    std::size_t count = 0;

    void finish(int sign) const {
        if (first >= 0) std::cout << sign << " " << first << " " << last << "\n";
    }
    void record(int tick, int sign) {
        if (first < 0 || tick-last > 42) {
            finish(sign);
            first = tick;
        }
        last = tick;
        ++count;
    }
};

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int total;
    if (!(std::cin >> total) || total != 13946) return EXIT_FAILURE;
    std::vector<Complex> roots(size/2);
    for (auto& root : roots) std::cin >> root.real >> root.imag;
    std::vector<Frequency> data(total);
    for (auto& term : data) {
        std::cin >> term.bin >> term.advance.real >> term.advance.imag;
        if (term.bin < 0 || term.bin >= size) return EXIT_FAILURE;
        for (auto& c : term.coefficients) std::cin >> c.real >> c.imag;
    }
    if (!std::cin) return EXIT_FAILURE;

    Run gain, loss;
    std::vector<Complex> bins(size);
    std::vector<Integer> values(size);
    for (int block = 0; block <= last_block; ++block) {
        if ((block+1)*size > 12494*steps_per_unit) {
            std::fill(values.begin(), values.end(), 0);
            for (int order = degree; order >= 0; --order) {
                std::fill(bins.begin(), bins.end(), Complex{});
                for (const auto& term : data) {
                    const Complex z = multiply(term.coefficients[order], term.phase);
                    bins[term.bin].real += z.real;
                    bins[term.bin].imag += z.imag;
                }
                transform(bins, roots);
                for (int i = 0; i < size; ++i) {
                    const int k = i-size/2;
                    const Integer old = static_cast<Integer>(Wide{values[i]}*k/(size/2));
                    values[i] = old+bins[(i+size/2)%size].imag;
                }
            }
            for (int i = 0; i < size; ++i) {
                const int tick = block*size+i;
                if (tick < 12494*steps_per_unit || tick > 1620550*steps_per_unit) continue;
                if (Wide{100}*values[i] > Wide{276}*scale) gain.record(tick, 1);
                if (-Wide{100}*values[i] > Wide{276}*scale) loss.record(tick, -1);
            }
        }
        for (auto& term : data) term.phase = multiply(term.phase, term.advance);
    }
    gain.finish(1);
    loss.finish(-1);
    std::cout << "0 " << gain.count << " " << loss.count << "\n";
    return EXIT_SUCCESS;
}
