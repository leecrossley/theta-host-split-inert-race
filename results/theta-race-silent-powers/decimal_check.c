/* Independent survivor sieve using decimal fixed-point phases and 128-bit products. */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

#define BLOCK 8192
#define LIMIT UINT32_C(1000000000)

int main(void) {
    const uint64_t modulus = UINT64_C(10000000000000000);
    const uint64_t edge = 17 * (modulus / 100);
    uint64_t phases[25];
    for (size_t j = 0; j < 25; j++) {
        if (scanf("%" SCNu64, &phases[j]) != 1 || phases[j] >= modulus) return 2;
    }
    uint32_t candidates[BLOCK];
    unsigned char hits[BLOCK];
    for (uint32_t start = 1; start <= LIMIT; start += BLOCK) {
        size_t count = LIMIT - start + 1 < BLOCK ? LIMIT - start + 1 : BLOCK;
        for (size_t k = 0; k < count; k++) {
            candidates[k] = start + k;
            hits[k] = 0;
        }
        for (size_t j = 25; j > 0 && count > 0; j--) {
            size_t remaining = 0;
            for (size_t k = 0; k < count; k++) {
                uint32_t n = candidates[k];
                uint64_t phase = ((__uint128_t)n * phases[j-1]) % modulus;
                if (edge <= phase && phase <= modulus-edge) hits[n-start]++;
                if (hits[n-start] < 2) candidates[remaining++] = n;
            }
            count = remaining;
        }
        if (count != 0) {
            fprintf(stderr, "Uncovered exponent %" PRIu32 "\n", candidates[0]);
            return 1;
        }
    }
    puts("Independent decimal sieve finds two witnesses for all 1000000000 exponents in blocks of 8192.");
    return 0;
}
