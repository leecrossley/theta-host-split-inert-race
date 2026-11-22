"""Verify the disappearing coefficients by factoring norms and testing the exact three-point formula."""

from fractions import Fraction
from math import isqrt

from flint import arb, ctx

from clock_phases_recheck import LocalWeights

Y = (216**2, 46691, 46695)
LIMIT = 46993


def prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))


def weight(norm, height):
    return max(0, min(2 - 2 * height / norm, 4 * height / norm - 1))


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    # Each good norm is a first, second or fourth power. Work backwards from that norm.
    square_residues = {u: {x * x % u for x in range(1, u)} for u in (547, 727, 571, 599)}

    def degree(p, u, v):
        if p % u == 0 or p % v == 0:
            return None
        equal = (p % u in square_residues[u]) == (p % v in square_residues[v])
        if p % 16 in (1, 15) and equal:
            return 1
        return 2 if p * p % 16 in (1, 15) else 4

    rows = []
    for n in range(min(Y) + 1, LIMIT + 1):
        if n % 16 not in (1, 15):
            continue
        for d, p in ((1, n), (2, isqrt(n)), (4, isqrt(isqrt(n)))):
            if p**d != n or not prime(p) or d != degree(p, 571, 599):
                continue
            a = degree(p, 547, 727)
            if a is not None and p**a > LIMIT:
                continue
            entries = tuple(n * weight(n, Fraction(y)) / 2 for y in Y)
            assert all(weight(n, Fraction(2 * y)) == 0 for y in Y)
            if any(entries):
                rows.append((p, entries))
    assert rows == [(46687, (31, 0, 0)), (46703, (47, 12, 8)), (46993, (337, 302, 298))]
    print("Independent factor-by-norm enumeration found exactly the same three nonzero rows.")

    # Check every possible residue-class knot inside the half-unit stencil. Then test
    # the left and right affine pieces exactly, separately for both terms of the move.
    for p, _ in rows:
        step = Fraction(1, 2)
        t = (p - step, Fraction(p), p + step)
        for k in (1, 2, 4, 8):
            assert [
                j
                for j in range(int(k * t[0]), int(k * t[2]) + 1)
                if k * t[0] < j < k * t[2] and j % 16 in (1, 15)
            ] == ([p] if k == 1 else [])
        for k in (1, 2):
            vals = [weight(p, k * s) for s in t]
            delta = (vals[0] - 2 * vals[1] + vals[2]) / step
            assert delta == (Fraction(2, p) if k == 1 else 0)
    print(
        "Rational second differences at all three last knots give exactly the predicted coefficient."
    )

    count = 0
    a_only = []
    b_only = []
    for n in range(min(Y) + 1, 8 * max(Y)):
        if n % 16 not in (1, 15):
            continue
        for d, p in ((1, n), (2, isqrt(n)), (4, isqrt(isqrt(n)))):
            if p**d != n or not prime(p) or degree(p, 571, 599) != d:
                continue
            relevant = any(y < n < 8 * y for y in Y)
            if not relevant:
                continue
            # A rational nonzero weight cannot cancel sqrt(2) times another rational weight.
            a = degree(p, 547, 727)
            assert a is not None
            count += 1
            if p**a <= min(Y) + 31:
                a_only.append((p**a, n))
            if n <= LIMIT:
                b_only.append((n, p**a))
    assert count == 3327
    assert sorted(a_only) == [
        (223, 49729),
        (271, 73441),
        (353, 124609),
        (367, 134689),
        (463, 214369),
        (46687, 46687),
    ]
    assert sorted(b_only) == [(46687, 46687), (46703, 46703), (46769, 46769**2), (46993, 46993)]

    # Directly check both sides of the entrance knot as exact rational numbers.
    n = 46769**2
    centre = Fraction(n, 8)
    step = Fraction(1, 16)
    for k in (1, 2):
        vals = [weight(n, k * t) for t in (centre - step, centre, centre + step)]
        delta = (vals[0] - 2 * vals[1] + vals[2]) / step
        assert delta == (0 if k == 1 else Fraction(8, n))
    assert centre == Fraction(2187339361, 8)
    print(
        "Independent enumeration recovered all 3327 active primes; both single-field tables and the exact entrance slope at 46769 squared agree."
    )

    for q in (46447, 46559):
        assert prime(q) and degree(q, 547, 727) == 2 and degree(q, 571, 599) == 1
        assert q * q < 46769**2 < 2 * q * q
    assert Fraction(46447**2, 4) == 539330952 + Fraction(1, 4)
    ctx.prec = 256
    field = LocalWeights()
    for name, y in zip(("square", "prime", "integer"), Y, strict=True):
        q = field.move(Fraction(y), 1)
        lower = (q * 10**18).floor().unique_fmpz()
        assert lower is not None
        print(f"Direct 256-bit {name} threshold in ({lower}, {lower + 1}) / 10^18.")


if __name__ == "__main__":
    main()
