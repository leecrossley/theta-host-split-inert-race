/* Certify closed dyadic cells by adding nonnegative, downward-rounded deficits. */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

enum { MODES = 100, DEPTH = 6, BINS = 256 };
static uint64_t step[MODES], count[DEPTH+1];
static uint16_t price[DEPTH+1][MODES][BINS];

static int cover(uint64_t middle, unsigned level) {
    ++count[level];
    unsigned lower = 0;
    for (unsigned j = 0; j < MODES; ++j) {
        /* The unsigned multiplication takes the residue modulo 2^64. */
        unsigned bin = (unsigned)((middle*step[j]) >> 56);
        lower += price[level][j][bin];
        if (lower > 5600) return 1;
    }
    if (level == DEPTH) {
        fprintf(stderr, "Uncovered dyadic cell at tick %" PRIu64 "\n", middle);
        return 0;
    }
    uint64_t offset = UINT64_C(32) >> level;
    return cover(middle-offset, level+1) && cover(middle+offset, level+1);
}

int main(void) {
    for (unsigned j = 0; j < MODES; ++j) {
        if (scanf("%" SCNu64, &step[j]) != 1) return 2;
        for (unsigned level = 0; level <= DEPTH; ++level) {
            for (unsigned k = 0; k < BINS; ++k) {
                unsigned amount;
                if (scanf("%u", &amount) != 1 || amount > UINT16_MAX) return 2;
                price[level][j][k] = (uint16_t)amount;
            }
        }
    }
    char extra;
    if (scanf(" %c", &extra) != EOF) return 2;
    /* Start with [1/2,1], then cover [n,n+1] through the right endpoint. */
    if (!cover(96, 1)) return 1;
    for (uint64_t n = 1; n < UINT64_C(1000000000); ++n)
        if (!cover(128*n+64, 0)) return 1;
    printf("Every real log displacement from 1/2 to 1000000000 has deficit > 56.\n");
    printf("Dyadic cells tested by depth:");
    for (unsigned k = 0; k <= DEPTH; ++k) printf(" %" PRIu64, count[k]);
    puts("");
    return 0;
}
