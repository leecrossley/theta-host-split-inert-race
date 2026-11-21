"""Directed Abel summation and rational bounds for the joint progression modulus."""
from fractions import Fraction as F
from itertools import combinations
from math import isqrt, prod
from pathlib import Path
import sys

PRIMES = (547, 727, 571, 599)
LIMIT = 1 << 26
GRID = 1 << 40


def cuts():
    points = [1]
    while points[-1] < LIMIT:
        left = points[-1]
        points.append(min(LIMIT, left + max(1, left // 192)))
    return points


def root_enclosures(n):
    low = high = n * GRID
    for _ in range(8):
        low = isqrt(low * GRID)
        high = isqrt(high * GRID - 1) + 1
    return (GRID, low), (GRID, high, (high * high + GRID - 1) // GRID)


def kronecker_two(d, n):
    if d == 1:
        return 1
    if n % 2 == 0:
        return 0
    a = 1 if n % 4 == 1 else -1
    b = 1 if n % 8 in (1, 7) else -1
    return {-4: a, 8: b, -8: a * b}[d]


def independent_prefix(primes, d, points):
    # Cross-check the compiled stream at the first thousand terms by Euler's criterion.
    total = 0
    first = []
    for l, r in zip(points, points[1:]):
        if r > 1025:
            break
        minimum = LIMIT
        for n in range(l, r):
            value = kronecker_two(d, n)
            for p in primes:
                residue = pow(n, (p - 1) // 2, p)
                value *= -1 if residue == p - 1 else residue
            total += value
            minimum = min(minimum, total)
        first.append(minimum)
    return first


def log_free_check(path):
    points = cuts()
    bounds = {x: root_enclosures(x) for x in points}
    labels = {d * prod(-p for p in pp): (pp, d)
              for n in range(5) for pp in combinations(PRIMES, n)
              for d in (1, -4, 8, -8) if pp or d != 1}
    assert len(labels) == 63
    lows = {}
    for line in path.read_text().splitlines():
        disc, *minima = map(int, line.split())
        assert disc not in lows and disc in labels and len(minima) == len(points) - 1
        pp, d = labels[disc]
        prefix = independent_prefix(pp, d, points)
        assert minima[:len(prefix)] == prefix
        q = abs(disc)
        pv10 = (isqrt(q) + 1) * (10 + 7 * q.bit_length())
        cell_lows = []
        for k in range(2):
            value = 0
            for left, right, m in zip(points, points[1:], minima):
                lo = (255-k) * bounds[left][0][k]
                hi = (256-k) * bounds[right][1][k+1]
                scale = lo if m >= 0 else hi
                value += m * scale * (right-left) // (256 * left * right)
            tail = pv10 * bounds[LIMIT][1][k+1]
            value -= (tail + 10*LIMIT - 1) // (10*LIMIT)
            cell_lows.append(value)
        lows[disc] = min(cell_lows)
    assert set(lows) == set(labels)
    disc = min(lows, key=lows.get)
    assert disc == -1_088_114_643_208 and lows[disc] * 1000 > 288 * GRID
    print(f"All 63 real primitive nonprincipal factors exceed 0.288 on [127/128,1]. "
          f"Two rational cells, {len(points)-1:,} Abel blocks and a uniform character-sum tail per factor.")


def verify_constants():
    modulus = 16 * prod(PRIMES)
    assert modulus == 2_176_229_286_416
    assert modulus * 3**30 < 8**30
    assert F(366, 24) + 11 < 366
    assert 10**7 * modulus * 366**12 < 2**173 < 2**183
    error = F(1, 10**7) * (2 + F(1, 366**10))
    assert F(1, 65) < F(1, 64) - error
    assert F(1, 64) + error < F(1, 63)
    assert 366**10 > 3 and F(1, 198) < F(1, 3*65) < F(1, 3*63) < F(1, 186)
    print("Checked the rational margins in Lee's bound under Dirichlet GRH for log L >= 366. For widths from L/(log L)^10 the prime weight of each sign pattern lies between width/65 and width/63.")


if __name__ == "__main__":
    log_free_check(Path(sys.argv[1]))
    verify_constants()
