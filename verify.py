"""Verify the four Hardy sign changes, the Riesz-product frequencies and the margin."""
from fractions import Fraction
from itertools import product

from flint import arb, ctx, dirichlet_char

SCALE = 100_000_000
# Conrey label and two endpoints, in units of 10^-8.
BRACKETS = (
    (9584, 5989, 12294899, 12294900),
    (4568, 1141, 27941486, 27941487),
    (4792, 1197, 31952267, 31952268),
    (9136, 1141, 37800727, 37800728),
)


def characters():
    chars = {}
    for p in (571, 599):
        for m, n in ((1, 1), (8, 5), (16, 5), (16, 13)):
            q = p * m
            a = next(a for a in range(p - 1, q, p) if a % m == n % m)
            chi = dirichlet_char(q, a)
            assert chi.is_primitive() and int(chi.conductor()) == q and chi.parity() == 1
            chars[q, a] = chi
    assert len(chars) == 8
    return chars


def hardy(chi, x):
    z = chi.hardy_z(x)
    assert z.imag.contains(0)
    return z.real


def sign(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    raise ArithmeticError(f"unresolved sign {x}")


def main():
    ctx.prec = 100
    chars = characters()
    for (q, a), chi in chars.items():
        print(f"central value of Hardy Z for {q}.{a} has sign {sign(hardy(chi, 0)):+d}")
    zeros = []
    for q, a, lo, hi in BRACKETS:
        chi = chars[q, a]
        left = sign(hardy(chi, arb(lo) / SCALE))
        right = sign(hardy(chi, arb(hi) / SCALE))
        assert left == -right
        zeros.append(arb(lo + hi, hi - lo) / (2 * SCALE))
        print(f"zero in ({lo}, {hi}) / {SCALE} for {q}.{a}; endpoint signs {left:+d}, {right:+d}")
    for i, a in enumerate(zeros):
        assert a > 0 and all(not a.overlaps(b) for b in zeros[:i])
    tested = 0
    for coeff in product((-1, 0, 1), repeat=len(zeros)):
        terms = [c for c in coeff if c]
        if len(terms) < 2 or terms[0] < 0:
            continue
        freq = abs(sum(c * g for c, g in zip(coeff, zeros)))
        assert freq > 0
        for chi in chars.values():
            sign(hardy(chi, freq))
        tested += 1
    assert tested == 36
    print("36 mixed frequencies up to sign exclude zero and all eight zero spectra (288 interval tests)")
    margin = sum(Fraction(1, 2) / (Fraction(1, 4) + Fraction(hi, SCALE) ** 2)
                 for _, _, _, hi in BRACKETS)
    assert margin > Fraction(61027, 10000)
    print("sum of four coefficient lower bounds strictly exceeds 61027/10000")


if __name__ == "__main__":
    main()
