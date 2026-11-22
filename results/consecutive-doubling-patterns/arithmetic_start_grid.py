"""Locate the first four-gain word on prime and square starting heights."""

import argparse
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

from flint import arb, ctx

from word_cover import (
    LIMIT,
    Race,
    linear_errors,
    log_errors,
    mesh,
    prefix_tables,
    real,
    subdivision,
    threshold_bracket,
    word,
)

EARLY = (F(45700), F(46740))
LATE = (F(20284000), F(20314000))
PRIME_STARTS = (1471, 11, 89, 127, 17, 4073, 23, 29, 47, 157, 41, 25057, 59, 251, 631)
SQUARE_STARTS = (1521, 9, 100, 144, 9604, 4096, 529, 324, 1089, 169, 81, 25281, 64, 256, 676)


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def next_prime(n):
    p = n + 1
    while not is_prime(p):
        p += 1
    return p


def next_square(n):
    return (isqrt(n) + 1) ** 2


def check_monotone(races, points, early):
    rules = (
        ((1, 1, 1), (1, 2, -1), (0, 1, 0), (0, 2, 0))
        if early
        else ((0, 1, 1), (1, 1, 0), (0, 2, 0), (1, 2, 0))
    )
    for left, right in zip(points, points[1:]):
        width = real(right - left)
        heights = real(left).union(real(right))
        offsets = arb(0).union(width)
        jumps = [e * 4 / width for e in linear_errors(left, right)]
        for field, k, orientation in rules:
            slope, intercept = races[field].move_band(left, k)
            slopes = slope + arb(0, jumps[k - 1].upper())
            affine = slope * real(left) + intercept + slopes * offsets
            if orientation:
                assert orientation * (2 * heights * slopes - affine) > 0, (left, field, k)
            else:
                assert affine > 2 * heights.sqrt(), (left, field, k)


def check_exclusion(races, nodes, threshold, method):
    previous = word(races, nodes[0])
    for left, right in zip(nodes, nodes[1:]):
        top = word(races, right)
        if method == "primary":
            error = log_errors(left, right)
            bounds = [
                max(a.upper(), b.upper()) + error[j // 2] - threshold
                for j, (a, b) in enumerate(zip(previous, top))
            ]
        else:
            error = linear_errors(left, right)
            lo, hi = real(left).sqrt(), real(right).sqrt()
            bounds = [
                max((lo * (a - threshold)).upper(), (hi * (b - threshold)).upper()) + error[j // 2]
                for j, (a, b) in enumerate(zip(previous, top))
            ]
        assert any(x < 0 for x in bounds), (left, right, method)
        previous = top


def show_bracket(label, x, digits):
    v = (x * 10**digits).floor().unique_fmpz()
    assert v is not None, (label, x)
    print(f"{label}: ({v}, {v + 1}) / 10^{digits}.")


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("method", choices=("primary", "independent"))
    parser.add_argument("sieve", type=Path)
    args = parser.parse_args()
    primary = args.method == "primary"
    ctx.prec = 160 if primary else 196
    early = mesh(*EARLY, F(1, 4) if primary else F(1, 5))
    late = mesh(*LATE, F(6) if primary else F(5))
    nodes = sorted(
        {
            EARLY[1],
            LATE[0],
            *(y for y in subdivision(4096 if primary else 3500) if EARLY[1] <= y <= LATE[0]),
        }
    )
    next_p = next_prime(46695)
    assert next_prime(46690) == 46691 and next_prime(46691) == next_p
    assert next_square(46655) == 46656 and next_square(46656) == 47089
    extra = (
        set(PRIME_STARTS) | set(SQUARE_STARTS) | {46691, 46695, 46656, next_p, 4506**2, 4507**2}
    )
    roots = ((F(20284000), F(20288000)), (F(20287900), F(20288500)))
    points = set(early) | set(late) | set(nodes) | set(map(F, extra))
    for y in (
        next_prime(20306533),
        next_square(20306533),
        next_prime(20308190),
        next_square(20308190),
    ):
        points.add(F(y))
    edges = {int(k * y) for y in points for k in (1, 2, 4, 8, 16)}
    for l, r in roots:
        for k in (1, 2, 4, 8):
            edges.update(range(int(k * l), int(k * r) + 1))
    assert max(edges) == 325024000 < 343057911569
    fields, bits, recip = prefix_tables(args.sieve.resolve(), args.method, sorted(edges))
    races = [Race(d, bits, recip) for d in fields]
    for classification, starts in (
        (is_prime, PRIME_STARTS),
        (lambda n: isqrt(n) ** 2 == n, SQUARE_STARTS),
    ):
        for code, y in enumerate(starts):
            assert classification(y)
            assert all(
                (1 if code & (1 << (3 - j)) else -1) * v > 2
                for j, v in enumerate(word(races, F(y)))
            )
    square_margin = races[1].move(F(46656), 1)
    prime_margin = races[1].move(F(46691), 1)
    integer_margin = races[1].move(F(46695), 1)
    assert 1 < square_margin < prime_margin < integer_margin < 2
    assert races[1].move(F(next_p), 2) < prime_margin
    check_monotone(races, early, True)
    check_monotone(races, late, False)
    assert races[1].move(EARLY[0], 1) < 0 and races[1].move(EARLY[1], 2) < square_margin
    assert races[0].move(LATE[0], 1) < square_margin and races[0].move(LATE[1], 1) > 2
    check_exclusion(races, nodes, square_margin, args.method)
    positive_roots = [
        threshold_bracket(races[0], interval, q)
        for interval, q in zip(roots, (square_margin, prime_margin), strict=True)
    ]
    assert positive_roots[0] > 434 * 46691 and positive_roots[1] > 434 * 46695
    grid_floor = [int(x.floor().unique_fmpz()) for x in positive_roots]
    for q, start in ((real(LIMIT), 20306533), (real(2), 20308190)):
        for succ in (next_prime, next_square):
            y = succ(start)
            assert LATE[0] < y < LATE[1] and all(v > q for v in word(races, F(y)))
    print(
        f"{args.method}: {len(nodes) - 1} inter-window exclusions and {len(early) - 1} early / {len(late) - 1} late slope cells."
    )
    print(
        "All fifteen other words have prime starts through 25057 and square starts through 25281, with four strict margins above 2."
    )
    for label, x, d in (
        ("Square threshold", square_margin, 10),
        ("Prime threshold", prime_margin, 10),
        ("Late square-threshold root", positive_roots[0], 5),
        ("Late prime-threshold root", positive_roots[1], 5),
    ):
        show_bracket(label, x, d)
    print(
        "At their exact thresholds the first square and prime starters are",
        next_square(grid_floor[0]),
        "and",
        next_prime(grid_floor[1]),
        "respectively.",
    )
    print(
        f"At 1.9735 the first prime is {next_prime(20306533)} and the first square is {next_square(20306533)}."
    )
    print(
        f"At 2 the first prime is {next_prime(20308190)} and the first square is {next_square(20308190)}."
    )


if __name__ == "__main__":
    main()
