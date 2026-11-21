"""Check the spectrum via residue classes and Hurwitz zeta, without a Dirichlet L evaluator."""
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, lcm

from flint import acb, arb, ctx, nmod_poly

if not __debug__:
    raise RuntimeError("The interval checks require Python assertions.")

# A character is a Legendre symbol at p times the j-th power of the even
# quartic character modulo 16 taking 5 to i. The primitive modulus is used.
TARGETS = (
    (727, 1, 48536362, 48536363),
    (547, 1, 248572181, 248572182),
    (547, 2, 434216664, 434216665),
    (547, 0, 446011529, 446011530),
    (727, 1, 579615600, 579615601),
)
SCALE = 10**9
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
    # Directly count base prime ideals for b=r=2, omitting rational ramification.
    limit = 400
    sieve = bytearray([1]) * (limit + 1)
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = bytes(len(range(p*p, limit + 1, p)))
    ideals = []
    checked = 0
    for p in range(3, limit+1, 2):
        if not sieve[p] or p in (547, 727):
            continue
        a, d = symbol(547, p), symbol(727, p)
        f = next(k for k in (1, 2, 4) if pow(p, k, 16) in (1, 15) and (a*d)**k == 1)

        def degrees(coefficients):
            factors = nmod_poly(coefficients, p).factor()[1]
            assert all(e == 1 for _, e in factors)
            return [poly.degree() for poly, _ in factors]

        # Tensor the residue fields of the real quartic and the real quadratic.
        polynomial_profile = []
        extension_degree = degrees([547, 0, 1])[0]
        for j, k in product(degrees([2, 0, -4, 0, 1]), degrees([-397669, 0, 1])):
            degree = lcm(j, k)
            splitting_sign = 1 if degree % extension_degree == 0 else -1
            polynomial_profile.extend([(degree, splitting_sign)] * gcd(j, k))
        assert polynomial_profile == [(f, a**f)] * (8//f)
        checked += 1
        if p**f <= limit:
            ideals.append((p**f, -8*a**f*arb(p).log()))
    assert checked == 77
    print("Polynomial factorizations confirm all 77 unramified rational-prime profiles through 400")
    drift = 2*arb(2).sqrt()*(arb(2).sqrt()-1)**2

    def score(y):
        total = arb(0)
        for n, coefficient in ideals:
            if y < n < 4*y:
                u = arb(y)/n
                total += coefficient * (2*(1-u) if n <= 2*y else 4*u-1)
        return total/(drift * arb(y).sqrt())

    for y, lower, upper in ((25, -12294, -12293), (50, 34499, 34500)):
        change = (score(2*y)-score(y))/arb(2).log()
        assert arb(lower)/1000 < change < arb(upper)/1000
        print(f"Direct b=r=2 slope from {y} to {2*y} in ({lower}/1000, {upper}/1000)")


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
    evaluators = {(p, j): Evaluator(p, j) for p in (547, 727) for j in range(4)}
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
    for cs in product((-1, 0, 1), repeat=len(TARGETS)):
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
    assert count == 928
    print("Residue-class computation confirms 8 centres, 5 sign changes and 928 exclusions")
    delta = arb(2).log()
    level = arb(0)
    swing = arb(0)
    for lower, upper in intervals:
        a, b = ball(lower), ball(upper)
        coefficient = 1 / (2*(arb(1)/4 + b*b))
        level += coefficient
        swing += coefficient * 2*(a*delta/2).sin()/delta
    assert level > arb(66921)/10000 and swing > arb(19735)/10000
    print("Independent sine and level bounds exceed 1.9735 and 6.6921")


if __name__ == "__main__":
    main()
