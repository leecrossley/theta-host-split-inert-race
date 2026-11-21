"""Directed Abel summation and rational bounds for the joint progression modulus."""
from fractions import Fraction as F
from itertools import combinations
from math import factorial, isqrt, prod
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


def exp_upper(x):
    terms = [F(1)]
    for i in range(1, 81):
        terms.append(terms[-1] * x / i)
    return sum(terms) + terms[-1] * x / 81 / (1 - x / 82)


def verify_constants():
    modulus = 16 * prod(PRIMES)
    assert modulus == 2_176_229_286_416 < 1097**5
    assert exp_upper(F(7)) < 1097
    assert exp_upper(F(128000, 6397)) < modulus
    assert sum(F(7, 10)**j / factorial(j) for j in range(5)) > 2
    # log 2 > 693/1000 by the positive atanh series at 1/3.
    assert 2 * sum(F(1, (2*j+1)*3**(2*j+1)) for j in range(4)) > F(693, 1000)
    assert 1585 * F(693, 1000) > 1097
    assert F(1097, 2)**10 / factorial(10) > 4000 * 1097**3
    error = F(283, 6400) + F(1, 2000)
    assert error < F(45, 1000)
    assert F(1, 64) - 7*error/1097 > F(1, 66)
    assert F(1, 64) + 7*error/1097 < F(1, 62)
    assert F(1, 64) - error*(F(2, 6) + F(1, 1097)) > F(1, 1800)
    # Bordignon, Table 6, p. 1435. The triples below are (width factor, lower denominator, upper denominator).
    assert exp_upper(F(87, 10)) < 6004 and modulus <= 6004**5
    assert F(6004, 2)**28 / factorial(28) > (2 * 10**7) * 6004**12
    assert modulus <= 12959**3 and exp_upper(F(83, 10)) < 12959
    for log_start, table_constant, bounds in (
        (6004, F(19, 1000), [(1, 67, 61), (3, 65, 63)]),
        (12959, F(32, 10000), [(1, 65, 63)]),
    ):
        endpoint = table_constant / 64 + F(1, 10**7)
        for width_factor, lower, upper in bounds:
            interval = endpoint * (F(2, width_factor) + F(1, log_start**10))
            assert F(1, lower) < F(1, 64) - interval
            assert F(1, 64) + interval < F(1, upper)
    print("Bordignon's published Table 6 gives each of 16 patterns log mass between L/198 and L/186 in (L,4L/3] for log L >= 1,097.")
    print("With t = log L >= 1,097 the mass exceeds width/1,800 for widths at least 6L/t. At t >= 6,004, widths from L/t^10 give mass between width/67 and width/61.")
    print("The bounds width/65 and width/63 hold from width 3L/t^10 at t >= 6,004 and from width L/t^10 at t >= 12,959.")


if __name__ == "__main__":
    log_free_check(Path(sys.argv[1]))
    verify_constants()
