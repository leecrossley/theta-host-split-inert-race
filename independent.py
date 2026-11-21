"""Check the spectrum via residue classes and Hurwitz zeta, without a Dirichlet L evaluator."""
from fractions import Fraction
from itertools import product
from math import isqrt

from flint import acb, arb, ctx

# A character is a Legendre symbol at p times the j-th power of the even
# quartic character modulo 16 taking 5 to i. The primitive modulus is used.
TARGETS = (
    (599, 1, 12294899, 12294900),
    (571, 2, 27941486, 27941487),
    (599, 2, 31952267, 31952268),
    (571, 1, 37800727, 37800728),
)
SCALE = 10**8
DEGREE = 54


def ball(q):
    return arb(q.numerator) / q.denominator


def symbol(p, n):
    if n % p == 0:
        return 0
    return 1 if pow(n, (p - 1) // 2, p) == 1 else -1


def phase(p, j, n):
    legendre = symbol(p, n)
    if not legendre or (j and n % 2 == 0):
        return None
    quartic = {1: 0, 15: 0, 5: 1, 11: 1, 9: 2, 7: 2, 13: 3, 3: 3}
    return ((2 if legendre < 0 else 0) + (quartic[n % 16] * j if j else 0)) % 4


def check_group():
    # Characters odd on complex conjugation vs. Frobenius classes over F.
    fourth_roots = ((1, 0), (0, 1), (-1, 0), (0, -1))
    for a, d, k in product((1, -1), (1, -1), range(4)):
        f = next(r for r in (1, 2, 4) if (a*d)**r == 1 and (k*r) % 4 == 0)
        for n in range(1, 17):
            powers = [fourth_roots[(n*(j*k + (2 if v == -1 else 0))) % 4]
                      for v in (a, d) for j in range(4)]
            factors = tuple(map(sum, zip(*powers)))
            local = 8 * a**n if n % f == 0 else 0
            assert factors == (local, 0)
    print("Euler factors agree on all 16 Frobenius classes through a full period")


def check_windows():
    limit = 20_000
    sieve = bytearray([1]) * (limit + 1)
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = bytes(len(range(p*p, limit + 1, p)))
    b = arb(2)
    drift = 2*b.sqrt()*(b.sqrt() - 1)**2
    for y, lo, hi in ((50, 7512, 7513), (5000, -38475, -38474)):
        total = arb(0)
        for p in range(3, 4*y + 1, 2):
            if not sieve[p] or p in (571, 599):
                continue
            a, d = symbol(571, p), symbol(599, p)
            f = next(r for r in (1, 2, 4)
                     if pow(p, r, 16) in (1, 15) and (a*d)**r == 1)
            n = p**f
            if y < n < 4*y:
                u = arb(y)/n
                weight = 2*(1-u) if n <= 2*y else 4*u-1
                # (8/f) ideals, each with log norm f*log(p).
                total -= 8*a**f*arb(p).log()*weight
        normalized = total / (drift * arb(y).sqrt())
        assert arb(lo)/1000 < normalized < arb(hi)/1000
        print(f"Direct b=2 window at y={y}: normalized race in ({lo}/1000, {hi}/1000)")


class Evaluator:
    def __init__(self, p, j):
        self.q = q = p * (1, 16, 8, 16)[j]
        self.first = []
        moments = [[0, 0] for _ in range(DEGREE + 1)]
        for a in range(1, q + 1):
            k = phase(p, j, a)
            if k is None:
                continue
            sg = 1 if k < 2 else -1
            self.first.append((arb(a).log(), sg / arb(a).sqrt(), k % 2))
            offset = 2 * a - q
            power = 1
            for row in moments:
                row[k % 2] += sg * power
                power *= offset
        assert moments[0] == [0, 0]
        self.moments = [acb(re, im) / (2*q)**k for k, (re, im) in enumerate(moments)]
        # For |s| < 2, Re(s)=1/2, the discarded terms are bounded by
        # 4 sqrt(q) sum_{k>K} (k+1)/3^k = sqrt(q) (2K+5)/3^K.
        self.error = arb(q).sqrt() * (2*DEGREE + 5) / 3**DEGREE

    def value(self, t):
        # t is an exact rational ordinate evaluated as a ball.
        s = acb(arb(1)/2, t)
        assert abs(s) < 2
        sums = [acb(0), acb(0)]
        for log_n, weight, idx in self.first:
            angle = -t * log_n
            sn, cs = angle.sin_cos()
            sums[idx] += acb(cs, sn) * weight
        tail = acb(0)
        coef = acb(1)
        for k in range(1, DEGREE + 1):
            coef *= -(s + k - 1) / k
            tail += coef * self.moments[k] * (s+k).zeta(arb(3)/2)
        err = arb(0, self.error)
        return sums[0] + acb(0, 1) * sums[1] + tail * acb(self.q)**(-s) + acb(err, err)

    def completed(self, t):
        u = acb(arb(3)/4, t/2)
        return self.value(t) * (self.q / arb.pi())**u * u.gamma()

    def nonzero_on(self, lo, hi):
        mid = ball((lo + hi) / 2)
        radius = ball((hi - lo) / 2)
        q = arb(self.q)
        # Abel summation, |sum_{n<=x} chi(n)| <= min(x,q).
        derivative = 4*q.sqrt() + 2*(4*q.sqrt()*q.log() + 4)
        return abs(self.value(mid)) > radius * derivative


def main():
    ctx.prec = 100
    check_group()
    check_windows()
    evaluators = {(p, j): Evaluator(p, j) for p in (571, 599) for j in range(4)}
    for ev in evaluators.values():
        assert abs(ev.value(arb(0))) > 0
    intervals = []
    for p, j, a, b in TARGETS:
        ev = evaluators[p, j]
        ref = ev.completed(arb(0)).conjugate()
        lo, hi = (Fraction(v, SCALE) for v in (a, b))
        values = [ev.completed(ball(t)) * ref for t in (lo, hi)]
        assert all(z.imag.contains(0) for z in values)
        assert values[0].real * values[1].real < 0
        intervals.append((lo, hi))
    count = 0
    for cs in product((-1, 0, 1), repeat=4):
        nonzero = [c for c in cs if c]
        if len(nonzero) < 2 or nonzero[0] < 0:
            continue
        los, his = zip(*[sorted((c*l, c*h)) for c, (l, h) in zip(cs, intervals)])
        lo, hi = sum(los), sum(his)
        assert lo*hi > 0
        lo, hi = sorted((abs(lo), abs(hi)))
        for ev in evaluators.values():
            assert ev.nonzero_on(lo, hi), (ev.q, lo, hi)
            count += 1
    assert count == 288
    print("Residue-class computation confirms 8 centres, 4 sign changes and 288 exclusions")


if __name__ == "__main__":
    main()
