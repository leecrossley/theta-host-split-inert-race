#include <gmpxx.h>
#include <fstream>
#include <iostream>

// Check every supplied modular equality exactly. A successful test forces each
// prime divisor to be 1 modulo the full power of two dividing p-1.
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::ifstream input(argv[1]);
    if (!input) return 3;
    mpz_class p, a;
    unsigned long count = 0;
    while (input >> p) {
        if (!(input >> a) || p < 3 || mpz_even_p(p.get_mpz_t()) || a < 2 || a >= p) return 4;
        mpz_class m = p - 1;
        auto n = mpz_scan1(m.get_mpz_t(), 0);
        mpz_class k = m >> n;
        mpz_class step = mpz_class(1) << n;
        if (k >= step) return 5;
        mpz_class exponent = m / 2, value;
        mpz_powm(value.get_mpz_t(), a.get_mpz_t(), exponent.get_mpz_t(), p.get_mpz_t());
        if (value != m) return 6;
        ++count;
    }
    if (!input.eof()) return 7;
    std::cout << count << '\n';
    return 0;
}
