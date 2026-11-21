"""Quadratic and quartic character phases for the two fields."""
from flint import arb


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
