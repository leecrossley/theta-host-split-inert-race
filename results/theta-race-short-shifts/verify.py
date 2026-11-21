"""Verify five Hardy sign changes, the signed frequencies and two uniform race bounds."""
from fractions import Fraction
from itertools import product
from math import factorial

from flint import arb, ctx, dirichlet_char

if not __debug__:
    raise RuntimeError("The interval checks require Python assertions.")

SCALE = 1_000_000_000
# Conrey label and two endpoints, in units of 10^-9.
BRACKETS = (
    (11632, 7269, 48536362, 48536363),
    (8752, 1093, 248572181, 248572182),
    (4376, 1093, 434216664, 434216665),
    (547, 546, 446011529, 446011530),
    (11632, 7269, 579615600, 579615601),
)


def characters():
    chars = {}
    for p in (547, 727):
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
    assert tested == 116
    print("116 mixed frequencies up to sign exclude zero and all eight zero spectra (928 interval tests)")
    upper_log2 = Fraction(69315, 100000)
    assert sum(upper_log2**k / factorial(k) for k in range(10)) > 2
    level = Fraction(0)
    swing = Fraction(0)
    for _, _, lower, upper in BRACKETS:
        lo, hi = Fraction(lower, SCALE), Fraction(upper, SCALE)
        coefficient = Fraction(1, 2) / (Fraction(1, 4) + hi**2)
        level += coefficient
        swing += coefficient * lo * (1 - (hi*upper_log2)**2 / 24)
    assert level > Fraction(66921, 10000)
    assert swing > Fraction(19735, 10000)
    print("Rational lower bounds: level amplitude > 6.6921; logarithmic shift slope > 1.9735")


if __name__ == "__main__":
    main()
