"""Exact harmonic counts and point counts for first-crossing fields on both hosts."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import isqrt

INFINITY = None
F = Fraction


def legendre(points):
    finite = sorted(x for x in points if x is not INFINITY)
    if len(finite) == 3:
        a, b, c = finite
        return (b - a) / (c - a)
    a, b, c, d = finite
    return (b - a) * (d - c) / ((c - a) * (d - b))


def invariant(value):
    return 256 * (1 - value + value**2)**3 / (value * (1 - value))**2


def harmonic_subsets(squares):
    return [group for group in combinations([INFINITY, *squares], 4)
            if invariant(legendre(group)) == 1728]


def completions(triple):
    finite = sorted(x for x in triple if x is not INFINITY)
    if len(finite) == 2:
        a, b = finite
        return {2*a - b, 2*b - a, (a + b)/2}
    a, b, c = finite
    # Send a, b, c to 0, 1, infinity, then pull the three harmonic values back.
    result = set()
    for z in (F(-1), F(1, 2), F(2)):
        den = z * (b - a) - (b - c)
        result.add(INFINITY if den == 0 else (z*c*(b-a) - a*(b-c)) / den)
    return result


def count_over_eleven(elements, roots):
    # i^2 = -1 gives the quadratic field since -1 has no square root modulo 11.
    def multiply(a, b):
        return ((a[0]*b[0] - a[1]*b[1]) % 11, (a[0]*b[1] + a[1]*b[0]) % 11)

    square_counts = Counter(multiply(y, y) for y in elements)
    total = 1  # Odd degree gives exactly one point at infinity.
    for x, y in elements:
        value = (1, 0)
        for root in roots:
            value = multiply(value, (x - root, y))
        total += square_counts[value]
    return total


def check():
    points = [INFINITY, F(0), F(1), F(2)]
    alternatives = []
    for triple in combinations(points, 3):
        missing = next(x for x in points if x not in triple)
        values = completions(triple)
        assert missing in values
        alternatives.extend(values - {missing})
    assert len(alternatives) == len(set(alternatives)) == 8
    assert set(alternatives) == {F(-2), F(-1), F(1,2), F(2,3), F(4,3), F(3,2), F(3), F(4)}
    branch_set = {1 + 24*x for x in alternatives}
    exceptions = (13, 17, 33, 37, 73, 97)
    assert branch_set == set(exceptions) | {-47, -23}
    base = [F(1,4), F(25,4), F(49,4)]
    for value in exceptions:
        assert 0 < value < 4 * 6**2 and isqrt(value)**2 != value
        assert len(harmonic_subsets(base + [F(value,4)])) == 2
    examples = (([1,2,3,4],0), ([1,5,7,8],1), ([7,13,17,23],2))
    for values, count in examples:
        squares = [F(x,4)**2 for x in values]
        assert len(harmonic_subsets(squares)) == count
    print("Rational four-target examples attain ranks 0, 2 and 4. The eight harmonic completions outside a fixed harmonic four-set are distinct.")
    print("For the fixed triple (1/2,5/2,7/2), the six positive exceptional values of 4t^2 with t^2 rational are 13, 17, 33, 37, 73 and 97.")

    roots = [1, 25, 49, 625, 1225]
    assert len({x % 11 for x in roots}) == 5
    groups = harmonic_subsets(list(map(F, roots)))
    assert groups == [(None,1,25,49), (None,25,625,1225)]
    fields = ([(a,0) for a in range(11)], list(product(range(11), repeat=2)))
    counts = [count_over_eleven(field, roots) for field in fields]
    fixed = [count_over_eleven(field, [-1,0,1]) for field in fields]
    assert counts == [16,118] and fixed == [12,144]
    traces = [len(field) + 1 - count for field, count in zip(fields, counts)]
    assert traces == [-4,4] and traces[0]**2 - traces[1] == 2*6
    print("For targets (1,5,7,25,35)/10 only two of the fifteen elliptic subcovers are harmonic. The genus-two curve has 16 and 118 points over F11 and F121; E has 12 and 144. The full rank is four.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Use the verifier with assertions enabled.")
    check()
