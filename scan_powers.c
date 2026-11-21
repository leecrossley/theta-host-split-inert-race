/* Check a phase witness for every exponent from 1 through one billion. */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    FILE *file = fopen(argv[1], "r");
    if (file == NULL) return 2;
    uint64_t phase[25];
    for (size_t i = 0; i < 25; i++) {
        unsigned row;
        if (fscanf(file, "%u %" SCNu64, &row, &phase[i]) != 2) return 2;
    }
    char extra;
    if (fscanf(file, " %c", &extra) != EOF) return 2;
    fclose(file);
    const uint64_t edge = UINT64_C(11) << 58;
    const uint64_t opposite = UINT64_C(0) - edge;
    uint64_t deepest = 0;
    for (uint64_t n = 1; n <= UINT64_C(1000000000); n++) {
        size_t j = 0;
        unsigned hits = 0;
        for (; j < 25; j++) {
            /* Unsigned multiplication gives the phase numerator modulo 2^64. */
            uint64_t residue = n * phase[j];
            if (edge <= residue && residue <= opposite && ++hits == 2) break;
        }
        if (j == 25) {
            fprintf(stderr, "Uncovered exponent: %" PRIu64 "\n", n);
            return 1;
        }
        if (j + 1 > deepest) deepest = j + 1;
    }
    printf("All 1000000000 positive exponents have two certified phase witnesses; at most %" PRIu64
           " of the 25 phases were needed per exponent.\n", deepest);
    return 0;
}
