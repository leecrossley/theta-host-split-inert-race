/* Verify a second cover on a translated three-quarter grid using a Mersenne modulus. */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

enum { MODES = 128, DEPTH = 5, DISTANCES = 129 };
static const uint64_t modulus = (UINT64_C(1)<<61)-1;
static uint64_t numerator[MODES], visited[DEPTH+1];
static uint16_t reward[DEPTH+1][MODES][DISTANCES];

static uint64_t phase(uint64_t ticks, uint64_t multiplier) {
    __uint128_t product = (__uint128_t)ticks*multiplier;
    uint64_t reduced = ((uint64_t)product & modulus) + (uint64_t)(product >> 61);
    return reduced >= modulus ? reduced-modulus : reduced;
}

static int inspect(uint64_t left, uint64_t width, unsigned depth) {
    uint64_t centre = left+width/2;
    unsigned total = 0;
    ++visited[depth];
    for (unsigned j = 0; j < MODES; ++j) {
        uint64_t residue = phase(centre, numerator[j]);
        uint64_t distance = residue < modulus-residue ? residue : modulus-residue;
        unsigned bin = (unsigned)(distance >> 53);
        total += reward[depth][j][bin];
        if (total > 5600) return 1;
    }
    if (depth == DEPTH) {
        fprintf(stderr, "Uncovered three-quarter cell at tick %" PRIu64 "\n", left);
        return 0;
    }
    return inspect(left, width/2, depth+1) && inspect(left+width/2, width/2, depth+1);
}

int main(void) {
    for (unsigned j = 0; j < MODES; ++j) {
        if (scanf("%" SCNu64, &numerator[j]) != 1 || numerator[j] >= modulus) return 2;
        for (unsigned l = 0; l <= DEPTH; ++l) {
            for (unsigned k = 0; k < DISTANCES; ++k) {
                unsigned value;
                if (scanf("%u", &value) != 1 || value > UINT16_MAX) return 2;
                reward[l][j][k] = (uint16_t)value;
            }
        }
    }
    char extra;
    if (scanf(" %c", &extra) != EOF) return 2;
    for (uint64_t left = 128; left < UINT64_C(256000000000); left += 192)
        if (!inspect(left, 192, 0)) return 1;
    printf("The independent three-quarter grid certifies the full displacement range.\n");
    printf("Three-quarter cells tested by depth:");
    for (unsigned l = 0; l <= DEPTH; ++l) printf(" %" PRIu64, visited[l]);
    puts("");
    return 0;
}
