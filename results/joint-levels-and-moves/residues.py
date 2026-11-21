"""Periodic character sum evaluated via a convergent Hurwitz expansion."""
from flint import acb, arb

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
        # For |s| < 4 use (k+3)^3 as an upper bound for |(s)_k/k!|.
        self.error = 8*arb(q).sqrt()*(DEGREE+4)**3 / 3**DEGREE

    def value(self, t):
        # t is an exact rational ordinate evaluated as a ball.
        s = acb(arb(1)/2, t)
        assert abs(s) < 4
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
