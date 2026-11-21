"""Compare the standalone first-loss certificate with the earlier family certificate."""

from collections import Counter
from decimal import Decimal, localcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from itertools import product
from pathlib import Path
import lzma
import sys

from flint import arb, ctx

RESULT = Path(__file__).resolve().parent


def check_prime_coefficients():
    moduli = (547, 727)
    residues = {q: {j*j % q for j in range(1, q)} for q in moduli}
    count = Counter()
    for residue in (1, 15):
        for n in range(residue, 16*547*727, 16):
            if n % 547 and n % 727:
                count[n % 547 in residues[547], n % 727 in residues[727]] += 1
    assert len(count) == 4 and set(count.values()) == {198198}

    for residue, u, v in product(range(1, 16, 2), (-1, 1), (-1, 1)):
        degree = next(k for k in range(1, 5)
                      if (u*v)**k == 1 and pow(residue, k, 16) in (1, 15))
        for k in range(1, 65):
            odd_trace = -4*(u**k + v**k)*int(pow(residue, k, 16) in (1, 15))
            prime_trace = -8*u**degree*int(k == degree)
            difference = prime_trace - odd_trace
            if k == 1:
                assert difference == 0
            if k == 2:
                assert difference == 8*(degree == 1)
            if k == 3:
                assert abs(difference) == 8*(degree == 1)
            assert abs(difference) <= 8

    slopes = (0, 4, -2, 0)
    knots = (4, 2, 1)
    jumps = [Fraction(abs(b-a), k) for a, b, k in zip(slopes, slopes[1:], knots)]
    assert jumps == [1, 3, 2]
    assert 4 + 4*sum(k*d for k, d in zip(knots, jumps)) == 52
    print("All four degree-one residue types have 198,198 classes. All 64 exponents and the three one-sided tent jumps give the stated correction and smoothing constants.")


def check_moments(data):
    from cell_check import even_moments
    from early_400 import higher_even_moments

    ctx.prec = 148
    a, root = arb(2).log(), arb(2).sqrt()
    v0, v1 = even_moments(a, root)
    v2, v3 = higher_even_moments(a, root)
    total = (71*v0 - 42*v1 - 6*v2 - v3)/(54*a)
    for line in data.splitlines():
        _, _, lo, hi = map(int, line.split())
        g = (arb(lo)/10**14).union(arb(hi)/10**14)
        u = (a*g/2).sin()**2
        weight = ((root-1)**2 + 4*root*u)/(2*(root-1)**2*(arb(1)/4 + g*g))
        polynomial = (16*u**3 - 48*u**2 + 75*u + 11)/54
        total -= 4*weight*polynomial/a
    assert total < arb(751)/1000
    print("Independent real gamma quadratures at four shifts and interval sine weights put the signed unlisted majorant below 0.751. This looser bound still leaves the full score above -1.96.")


def check_integer_endpoints():
    with localcontext() as context:
        context.prec = 95
        root = Decimal(8).sqrt()
        lower = (Decimal('40.0936').exp()/root).to_integral_value(rounding=ROUND_CEILING)
        upper = (Decimal('42.2909').exp()/(9*root)).to_integral_value(rounding=ROUND_FLOOR)
        assert (lower, upper) == (91386963248971966, 91393856157112503)
        assert upper-lower+1 == 6892908140538
        assert Decimal('19.95')*root/Decimal('40.0507').exp() < Decimal('2.28e-16')
    print("An independent 95-place decimal evaluation reproduces the integer endpoints and the real-versus-integer factor separation.")


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    base = (Path(sys.argv[1]).resolve() if len(sys.argv) > 1
            else RESULT.parent/'theta-race-factor-400')
    if not (base/'check_zeros.py').is_file():
        raise SystemExit("Supply the path to the theta-race-factor-400 result in the family checkout.")
    sys.path.insert(0, str(base))
    from check_zeros import check_family

    check_family()
    data = lzma.decompress((RESULT/'evidence/onset-zeros.tsv.xz').read_bytes())
    assert data == (base/'evidence/thousand-brackets.tsv').read_bytes()
    print("The 13,946 brackets match the earlier complete certificate byte for byte. Every Conrey character matches its field twist on all unit classes.")
    check_prime_coefficients()
    check_moments(data)
    check_integer_endpoints()


if __name__ == '__main__':
    main()
