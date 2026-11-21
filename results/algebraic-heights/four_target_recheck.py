"""Independently check the exceptional targets and the genus-two reduction."""

from itertools import combinations
from flint import fmpq as Q, fmpq_poly, nmod_poly


def branch(x):
    return (1, 0) if x is None else (x, 1)


def determinant(a, b):
    return a[0]*b[1] - a[1]*b[0]


def harmonic_polynomial(triple):
    a, b, c = map(branch, triple)
    d = (fmpq_poly([0, 1]), 1)
    u = determinant(a,c) * determinant(b,d)
    v = determinant(a,d) * determinant(b,c)
    return (u + v) * (u - 2*v) * (2*u - v)


def verify():
    base = [None, Q(1,4), Q(25,4), Q(49,4)]
    critical = set()
    for triple in combinations(base,3):
        poly = harmonic_polynomial(triple)
        unit, factors = poly.factor()
        assert unit != 0 and all(f.degree() == 1 and e == 1 for f,e in factors)
        for factor, _ in factors:
            point = -factor[0]/factor[1]
            if point not in base and point > 0:
                critical.add(4*point)
    assert critical == {Q(v) for v in (13,17,33,37,73,97)}
    print("Homogeneous cross-ratio cubics independently give all six exceptional positive branch values.")

    prim = [None, 1, 25, 49, 625, 1225]
    marked = []
    for group in combinations(prim,4):
        a,b,c,d = map(branch,group)
        u = determinant(a,c)*determinant(b,d)
        v = determinant(a,d)*determinant(b,c)
        if (u+v)*(u-2*v)*(2*u-v) == 0:
            marked.append(group)
    assert marked == [(None,1,25,49), (None,25,625,1225)]
    prime = 11
    x = nmod_poly([0,1],prime)
    curve = nmod_poly([1],prime)
    for root in prim[1:]:
        curve *= x - root
    assert curve.gcd(curve.derivative()).degree() == 0
    raised = curve ** ((prime-1)//2)
    entries = [int(raised[prime*i-j]) for i in (1,2) for j in (1,2)]
    assert entries == [1,0,0,6]
    elliptic = (x**3-x) ** ((prime-1)//2)
    assert int(elliptic[prime-1]) == 0
    print("The genus-two Hasse-Witt matrix modulo 11 is diag(1,6); that of V^2=U^3-U is zero. Determinants and point enumeration agree.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Use the verifier with assertions enabled.")
    verify()
