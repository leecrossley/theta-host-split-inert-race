"""Enclose all real starts for the first full two-doubling sign catalogue."""

import argparse
from fractions import Fraction
from functools import lru_cache
from math import isqrt
from pathlib import Path
import subprocess

from flint import arb, ctx

F = Fraction
A = arb(2).log
WITNESSES = (1469, 9, 89, 33, 17, 15, 20, 28, 45, 8, 39, 150, 56, 248, 629)
PEAK = (F("46695.21"), F("46695.22"))
JUMP = (F("20288514.03"), F("20288514.09"))
FIRST = (F("20306533.58"), F("20306533.61"))
LIMIT = F(3947, 2000)
EARLY = F(45700)
EARLY_END = F(46701)
LATE_BEGIN = F(20288400)
LATE_END = F(20309000)


def mesh(left, right, step):
    points = [left + i * step for i in range(1 + (right - left) // step)]
    if points[-1] != right:
        points.append(right)
    return points


def real(value):
    value = F(value)
    return arb(value.numerator) / value.denominator


def subdivision(n):
    upper = FIRST[0]
    points = {F(1), F(6), EARLY, EARLY_END, LATE_BEGIN, upper, *PEAK, *JUMP, *FIRST}
    for power in range(25):
        start = 2**power
        if start > upper:
            break
        points.update(F(start * (n + j), n) for j in range(n + 1) if start * (n + j) <= upper * n)
    for centre in (46695, 20288514, 20306533):
        points.update(F(centre + i) for i in range(-1400, 1401, 3))
        points.update(F(100 * centre + i, 100) for i in range(-250, 251))
        points.update(F(1000 * centre + i, 1000) for i in range(20, 251))
    points.update(F(y) for y in range(20000000, 20307000, 128))
    return sorted(t for t in points if 1 <= t <= FIRST[1])


def prefix_tables(binary, mode, edges):
    count, bits, recip = (2**15, 58, 96) if mode == "primary" else (2**12, 62, 102)
    loggrid = []
    for j in range(count + 1):
        lower = ((1 + arb(j) / count).log() * 2**bits).floor().unique_fmpz()
        assert lower is not None
        loggrid.append(int(lower))
    data = "".join(str(i) + "\n" for i in [*loggrid, *edges])
    fields = []
    for label in ("547", "571"):
        result = subprocess.run(
            [str(binary), label], input=data, text=True, capture_output=True, check=True
        )
        totals = [0, 0, 0, 0]
        sums = {}
        for e, line in zip(edges, result.stdout.splitlines(), strict=True):
            k, *v = map(int, line.split())
            assert k == e and len(v) == 4 and v[0] <= v[1] and v[2] <= v[3]
            totals = [x + y for x, y in zip(totals, v)]
            sums[e] = tuple(totals)
        fields.append(sums)
    return fields, bits, recip


def small_prime_power(n):
    """Whether an integer can be a rational prime power; only used on tiny local bands."""
    d = 2
    while d * d <= n and n % d:
        d += 1 if d == 2 else 2
    if d * d > n:
        return True
    while n % d == 0:
        n //= d
    return n == 1


def check_band(lo, hi):
    for k in (1, 2, 4, 8, 16):
        for n in range(int(k * lo) + 1, int(k * hi) + 1):
            assert not small_prime_power(n), (lo, hi, n)


class Race:
    def __init__(self, sums, bits, recip):
        self.table = sums
        self.bits = bits
        self.recip = recip
        self.a = A()
        self.r = arb(2).sqrt()
        self.norm = 2 * self.r * (self.r - 1) ** 2
        self.factor = 8 / self.norm

    def band(self, i, j):
        a, b, c, d = (v - u for u, v in zip(self.table[int(i)], self.table[int(j)]))
        # The cumulative endpoint differences are sums of the enclosed disjoint bins.
        return arb(a + b, b - a) / 2 ** (self.bits + 1), arb(c + d, d - c) / 2 ** (self.recip + 1)

    @lru_cache(maxsize=None)
    def level(self, y):
        p, q = self.band(y, 2 * y)
        s, t = self.band(2 * y, 4 * y)
        v = real(y)
        return self.factor * (2 * p - 2 * v * q + 4 * v * t - s) / v.sqrt()

    def coefficients(self, y, k):
        z = k * y
        p, q = self.band(z, 2 * z)
        s, t = self.band(2 * z, 4 * z)
        return self.factor * arb(k).sqrt() * (4 * t - 2 * q), self.factor * (2 * p - s) / arb(
            k
        ).sqrt()

    def move_band(self, y, k):
        a, b = self.coefficients(y, k)
        c, d = self.coefficients(y, 2 * k)
        return (c - a) / self.a, (d - b) / self.a

    def move(self, y, k):
        return (self.level(2 * k * y) - self.level(k * y)) / self.a


def peak_data(races):
    check_band(*PEAK)
    a, b = races
    left, right = map(real, PEAK)
    u, v = b.move_band(PEAK[0], 1)
    s, t = b.move_band(PEAK[0], 2)
    y = (t - v) / (u - s)
    value = u * y.sqrt() + v / y.sqrt()
    assert left < y < right
    whole = left.union(right)
    assert u - v / whole > 0 and s - t / whole < 0
    for k in (1, 2):
        U, V = a.move_band(PEAK[0], k)
        centre = real(sum(PEAK) / 2)
        base = U * centre.sqrt() + V / centre.sqrt()
        derivative = (U - V / whole) / (2 * whole.sqrt())
        assert base + derivative * (whole - centre) > 3 > value
    assert arb("1.68244928") < value < arb("1.68244930")
    return y, value


def onset(races, interval, threshold, *, ceiling=None, other_bounds=None):
    check_band(*interval)
    x, z = map(real, interval)
    whole = x.union(z)
    a, b = races
    coeff = [race.move_band(interval[0], k) for k in (1, 2) for race in (a, b)]
    u, v = coeff[0]
    assert (u - v / whole) > 0

    def value(y):
        return u * y.sqrt() + v / y.sqrt() - threshold

    assert value(x) < 0 < value(z)
    if ceiling is not None:
        assert u * z.sqrt() + v / z.sqrt() < ceiling
    # Stop when the directed log enclosure straddles the threshold at a test point.
    for _ in range(70):
        m = (x + z) / 2
        res = value(m)
        if res < 0:
            x = m
        elif res > 0:
            z = m
        else:
            break
    for (U, V), bound in zip(coeff[1:], other_bounds or [threshold] * 3, strict=True):
        centre = real(sum(interval) / 2)
        deriv = (U - V / whole) / (2 * whole.sqrt())
        assert U * centre.sqrt() + V / centre.sqrt() + deriv * (whole - centre) > bound
        if other_bounds is not None:
            assert bound > threshold
    return x.union(z)


def log_errors(lo, hi):
    """Bound raw level secant error in log height, including all derivative jumps."""
    root = arb(2).sqrt()
    a = A()
    B = 2 * root * (root - 1) ** 2
    l = real(lo)
    h = real(hi - lo) / l

    def level(k):
        L = 2 * k * l
        U = real(2 * k * hi)
        weight = 8 * (2 * U).log()
        smooth = (2 * U + 1) * weight * 4 * (a / 2).sinh() / (B * (L / 2).sqrt())
        jumps = arb(0)
        for s, c in ((arb(1) / 2, 2 / B), (1, 4 * (a / 2).cosh() / B), (2, 2 / B)):
            jumps += c * weight * (1 + s * U * h) / (s * L).sqrt()
        return smooth * h * h / 32 + jumps * h / 4

    e = [level(k) for k in (1, 2, 4)]
    return ((e[0] + e[1]) / a, (e[1] + e[2]) / a)


def linear_errors(lo, hi):
    """Bound the concave secant excess of sqrt(y) times the threshold-subtracted move."""
    r = arb(2).sqrt()
    B = 2 * r * (r - 1) ** 2
    l = real(lo)
    h = real(hi - lo)

    def jumps(k):
        total = arb(0)
        for s, c in ((1, 2), (2, 6), (4, 4)):
            L = s * k * l
            weight = (real(s * k * hi)).log()
            total += c * (1 + s * k * h) * weight / L
        return 8 * arb(k).sqrt() * total / B

    e = [jumps(k) for k in (1, 2, 4)]
    return ((e[0] + e[1]) * h / (4 * A()), (e[1] + e[2]) * h / (4 * A()))


def verify_branches(races, early, late, integer_margin):
    # On each cell H(y)=sqrt(y)D(ky) is piecewise affine. Its right slope
    # may change by at most the total variation of its interior slope jumps.
    for points, rules in (
        (early, ((1, 1, "rising"), (1, 2, "falling"), (0, 1, "above"), (0, 2, "above"))),
        (late, ((0, 1, "rising"), (1, 1, "above"), (0, 2, "above"), (1, 2, "above"))),
    ):
        for left, right in zip(points, points[1:]):
            width = real(right - left)
            cell = real(left).union(real(right))
            offset = arb(0).union(width)
            errors = linear_errors(left, right)
            for i, k, rule in rules:
                u, v = races[i].move_band(left, k)
                jumps = 4 * errors[k - 1] / width
                slope = u + arb(0, jumps.upper())
                value = u * real(left) + v + slope * offset
                if rule == "above":
                    assert value > 2 * cell.sqrt(), (left, right, i, k)
                else:
                    derivative_sign = 2 * cell * slope - value
                    assert derivative_sign > 0 if rule == "rising" else derivative_sign < 0, (
                        left,
                        right,
                        i,
                        k,
                    )
    assert races[1].move(early[0], 1) < 0
    assert races[1].move(early[-1], 2) < integer_margin
    assert races[0].move(late[0], 1) < integer_margin
    assert races[0].move(late[-1], 1) > 2
    assert races[1].move(F(46696), 2) < integer_margin
    # Below six the first 571/599 move vanishes since both scores use norms at most 48.
    possible = []
    for p in range(3, 50, 2):
        if any(p % d == 0 for d in range(2, 1 + isqrt(p))):
            continue
        same = pow(p, 285, 571) == pow(p, 299, 599) == 1 or (
            pow(p, 285, 571) == 570 and pow(p, 299, 599) == 598
        )
        f = next(f for f in (1, 2, 4) if pow(p, f, 16) in (1, 15) and (same or f % 2 == 0))
        possible.append(p**f)
    assert min(possible) == 49


def threshold_bracket(race, interval, target):
    lo, hi = interval
    assert race.move(lo, 1) < target < race.move(hi, 1)
    for _ in range(50):
        mid = (lo + hi) / 2
        value = race.move(mid, 1)
        if value < target:
            lo = mid
        elif value > target:
            hi = mid
        else:
            break
    return real(lo).union(real(hi))


def word(races, y):
    return tuple(race.move(y, k) for k in (1, 2) for race in races)


def coverage(races, nodes, mode, integer_margin, threshold):
    count = [0, 0, 0]
    previous = None
    for lo, hi in zip(nodes, nodes[1:]):
        if (
            hi <= 6
            or (EARLY <= lo and hi <= EARLY_END)
            or (LATE_BEGIN <= lo and hi <= JUMP[1])
            or lo >= FIRST[0]
        ):
            previous = None
            continue
        # After the second crossing of the early peak, only the larger benchmark is tested.
        target = arb(0) if hi <= EARLY else integer_margin if hi <= LATE_BEGIN else threshold
        previous = word(races, lo) if previous is None else previous
        top = word(races, hi)
        if mode == "primary":
            err = log_errors(lo, hi)
            upper = [
                max(x.upper(), y.upper()) + err[j // 2] - target
                for j, (x, y) in enumerate(zip(previous, top))
            ]
        else:
            err = linear_errors(lo, hi)
            l = real(lo).sqrt()
            h = real(hi).sqrt()
            upper = [
                max((l * (x - target)).upper(), (h * (y - target)).upper()) + err[j // 2]
                for j, (x, y) in enumerate(zip(previous, top))
            ]
        if not any(v < 0 for v in upper):
            raise ArithmeticError((lo, hi, min(upper), mode))
        count[(0 if hi <= EARLY else 1 if hi <= JUMP[0] else 2)] += 1
        previous = top
    return count


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("method", choices=("primary", "independent"))
    parser.add_argument("sieve", type=Path)
    args = parser.parse_args()
    ctx.prec = 152 if args.method == "primary" else 188
    nodes = subdivision(4096 if args.method == "primary" else 3500)
    early = mesh(EARLY, EARLY_END, F(1, 4) if args.method == "primary" else F(1, 5))
    late = mesh(LATE_BEGIN, LATE_END, F(6) if args.method == "primary" else F(5))
    extra = ((F(45700), F(45800)), (F(46600), F(46630)), (F(20307500), F(20309000)))
    observations = {*nodes, *early, *late, *map(F, WITNESSES), F(20306534), F(46695), F(46696)}
    edges = sorted({int(k * y) for k in (1, 2, 4, 8, 16) for y in observations})
    edges = sorted(
        set(edges).union(
            *(
                range(int(k * l), int(k * r) + 1)
                for l, r in (*extra, (LATE_BEGIN, JUMP[1]))
                for k in (1, 2, 4, 8)
            )
        )
    )
    assert edges[-1] == 324944000 and edges[-1] < 343057911569
    fields, bits, recip = prefix_tables(args.sieve.resolve(), args.method, edges)
    races = [Race(table, bits, recip) for table in fields]
    prime = 183263
    assert all(prime % d for d in range(2, 1 + isqrt(prime)))
    assert prime % 16 in (1, 15) and all(pow(prime, (q - 1) // 2, q) == 1 for q in (571, 599))
    assert 4 * extra[0][1] < prime < 8 * extra[0][0]
    prime = 380207
    assert all(prime % d for d in range(2, 1 + isqrt(prime)))
    assert prime % 16 == 15 and all(pow(prime, (q - 1) // 2, q) == 1 for q in (571, 599))
    assert 8 * PEAK[1] < prime < 16 * PEAK[0]
    peak_height, peak = peak_data(races)
    benchmark = real(LIMIT)
    jump = onset(races, JUMP, peak, ceiling=benchmark)
    first = onset(races, FIRST, benchmark, other_bounds=(20, 48, 38))
    for code, y in enumerate(WITNESSES):
        v = word(races, F(y))
        assert all((1 if code & (1 << (3 - j)) else -1) * x > 2 for j, x in enumerate(v))
    assert all(x > benchmark for x in word(races, F(20306534)))
    assert first / peak_height > 434 and jump / peak_height > 434
    integer_margin = races[1].move(F(46695), 1)
    assert real(F(22639, 10**8)) < peak - integer_margin < real(F(22640, 10**8))
    verify_branches(races, early, late, integer_margin)
    integer_jump = threshold_bracket(races[0], (LATE_BEGIN, JUMP[1]), integer_margin)
    assert integer_jump / peak_height > 434
    crossings = [
        threshold_bracket(race, interval, target)
        for race, interval, target in zip(
            (races[1], races[1], races[0]), extra, (real(0), real(F(8, 5)), real(2)), strict=True
        )
    ]
    counts = coverage(races, nodes, args.method, integer_margin, benchmark)
    print(
        args.method,
        "continuous cover:",
        *counts,
        "closed arcs at zero, the early integer maximum and the benchmark.",
    )
    print(
        "Two strictly increasing first-passage branches:",
        len(early) - 1,
        len(late) - 1,
        "slope cells. Every margin from 0 through 2 is covered.",
    )
    print("The 15 other words occur at integer starts through 1469, all four margins > 2.")
    for label, value, digits in (
        ("early peak height", peak_height, 8),
        ("early peak margin", peak, 10),
        ("first return above peak", jump, 5),
        ("full catalogue onset", first, 5),
        ("integer peak", integer_margin, 10),
        ("late root at the integer peak", integer_jump, 5),
        ("zero margin", crossings[0], 5),
        ("margin 1.6", crossings[1], 5),
        ("margin 2", crossings[2], 5),
    ):
        scale = 10**digits
        q = (value * scale).floor().unique_fmpz()
        assert q is not None, (label, value)
        print(f"{label}: ({q}, {q + 1}) / 10^{digits}")
    print("The first successful integer start for four gains at 1.9735 is 20306534.")


if __name__ == "__main__":
    main()
