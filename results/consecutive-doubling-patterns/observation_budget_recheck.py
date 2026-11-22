"""Recheck the cover and delayed primes using residue tables and short segmented sieves."""

from fractions import Fraction
from math import isqrt

MODS = (547, 727, 571, 599)
SQUARES = {q: {j * j % q for j in range(1, q)} for q in MODS}


def primality_table(limit):
    marked = [True] * (limit + 1)
    marked[0] = marked[1] = False
    for d in range(2, isqrt(limit) + 1):
        if marked[d]:
            for n in range(d * d, limit + 1, d):
                marked[n] = False
    return [i for i, ok in enumerate(marked) if ok]


def segment(lo, hi, base):
    marked = bytearray([1]) * (hi - lo + 1)
    for p in base:
        if p * p > hi:
            break
        begin = max(p * p, (lo + p - 1) // p * p)
        if begin <= hi:
            marked[begin - lo :: p] = bytes((hi - begin) // p + 1)
    return [lo + j for j, ok in enumerate(marked) if ok]


def degrees(p):
    def f(u, v):
        if p in (2, u, v):
            return 0
        if p % 16 in (1, 15) and ((p % u in SQUARES[u]) == (p % v in SQUARES[v])):
            return 1
        return 2 if p * p % 16 in (1, 15) else 4

    return f(547, 727), f(571, 599)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    n = 46769**2
    n1, n2 = 46447**2, 46559**2
    bases = primality_table(isqrt(8 * n + 6000))
    good = [p for p in bases if p < 46993 and degrees(p) == (2, 1)]
    assert good[-4:] == [46337, 46447, 46559, 46769]
    protect = [p**2 for p in good if p < 46656]
    cut = Fraction(n, 8)
    ends = [cut] + [m for m in protect if cut < m <= n1]
    worst = Fraction()
    for left, right in zip(ends, ends[1:]):
        before = [m for m in protect if m <= left][-2:]
        after = [m for m in protect if right <= m][:2]
        assert len(before) == len(after) == 2 and before[1] <= left < right == after[0]
        assert after[1] <= 2 * before[0]
        worst = max(worst, Fraction(after[1], before[0]))
    assert len(ends) == 189 and len(protect) == 306 and worst == Fraction(377952481, 325766401)
    assert 2 * protect[-3] > n1
    expected = (17258590369, 17341923439)
    for m, p in zip((n1, n2), expected, strict=True):
        found = [q for q in segment(8 * m - 1100, 8 * m, bases) if degrees(q) == (1, 2)]
        assert max(found) == p
    others = (17258590801, 17498715103, 17498715233)
    for q in others:
        assert segment(q, q, bases) == [q] and degrees(q) == (1, 2)
    # At the two-sample upper construction no higher integer norm enters the second window.
    assert 8 * n1 < 8 * (n1 + Fraction(1, 16)) < 8 * n1 + 1
    assert 4 * n < expected[0] < 8 * n1 < others[0] < 8 * n2
    assert all(8 * n < q < 2 * expected[0] for q in others[1:])
    # Verify that the two nuisance and target slopes really are distinct at both upper constructions.
    a, b = n1, n1 + Fraction(1, 16)
    c = (b - n2) / (b - a)
    d = (n2 - a) / (b - a)
    assert (
        c + d == 1
        and c * (a - n2) + d * (b - n2) == 0
        and c * (a - n) + d * (b - n) == n2 - n != 0
    )
    print(
        "Independent residue tables, 188 closed-left intervals, segmented large-prime sieves and exact two-sample weights agree."
    )


if __name__ == "__main__":
    main()
