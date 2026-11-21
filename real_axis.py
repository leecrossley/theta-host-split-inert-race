"""Cover the real interval [1/2,1] for each primitive factor with certified balls."""
from fractions import Fraction
from math import isqrt
from flint import arb, acb, ctx, dirichlet_char

from independent import ball, phase

FACTORS = (
    (571, (570, 1141, 1141, 5709)),
    (599, (598, 5989, 1197, 1197)),
)


def partial_sum_bound(p, j, chi):
    a = b = largest = 0
    roots = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    for n in range(1, int(chi.modulus()) + 1):
        k = phase(p, j, n)
        x, y = (0, 0) if k is None else roots[k]
        # The possible values are 0 and the four unit roots; this disc
        # distinguishes each of them and also checks the Conrey labels.
        assert abs(chi(n) - acb(x, y)) < arb(1)/4
        a, b = a + x, b + y
        largest = max(largest, a*a + b*b)
    assert (a, b) == (0, 0)
    return isqrt(largest) + 1


def cover(chi, bound):
    m = arb(bound)
    derivative = 4*m.sqrt() + 4*m.sqrt()*m.log() + 4
    pending = [(Fraction(1, 2), Fraction(1))]
    accepted = 0
    while pending:
        lo, hi = pending.pop()
        mid = (lo + hi) / 2
        value = abs(chi.l_function(acb(ball(mid))))
        if value > derivative * ball((hi - lo)/2):
            accepted += 1
        else:
            if hi - lo < Fraction(1, 2**20):
                raise ArithmeticError(f"unresolved real interval {lo} to {hi}")
            pending.extend(((lo, mid), (mid, hi)))
    return accepted


def main():
    ctx.prec = 96
    for p, labels in FACTORS:
        for j, a in enumerate(labels):
            q = p * (1, 16, 8, 16)[j]
            chi = dirichlet_char(q, a)
            assert int(chi.order()) <= 4 and chi.is_primitive()
            m = partial_sum_bound(p, j, chi)
            print(f"{q}.{a}: real interval [1/2,1] covered by {cover(chi, m)} cells; partial sums at most {m}")


if __name__ == "__main__":
    main()
