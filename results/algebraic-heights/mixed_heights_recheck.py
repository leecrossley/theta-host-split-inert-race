"""Reconstruct the prime bands from two original tent evaluations over Q(sqrt(2))."""
from itertools import combinations

from flint import fmpq, fmpq_poly


X = fmpq_poly([0, 1])
MODULUS = X**2 - 2
ZERO = fmpq_poly()
ONE = fmpq_poly([1])


def quotient(a, b):
    gcd, inverse, _ = b.xgcd(MODULUS)
    assert gcd == ONE
    return (a * inverse) % MODULUS


def tent(y):
    if y < 1 < 2*y:
        return 2-2*y
    if 2*y < 1 < 4*y:
        return 4*y-1
    return fmpq(0)


def value(y):
    return X * (tent(2*y)/2) - fmpq_poly([tent(y)])


def recover(a, b):
    slope = (value(b)-value(a))/(b-a)
    return slope, value(a)-slope*a


def abs_key(z):
    a,b = z[0],z[1]
    if a*b >= 0:
        positive = a+b > 0
    else:
        positive = (a*a > 2*b*b) == (a > 0)
    sign = 1 if positive else -1
    return sign*a, sign*b


def step_rows(left,right):
    a,b = 110*left+right, left+110*right
    return recover(a/111,b/111)


def main():
    bounds = [fmpq(k,8) for k in (1,2,4,8)]
    bands = [step_rows(a,b) for a,b in zip(bounds,bounds[1:])]
    expected = (
        {(fmpq(1,8),0),(fmpq(1,4),0),(0,fmpq(1,4)),(fmpq(3,4),fmpq(1,2))},
        {(fmpq(3,2),-1),(0,fmpq(1,4)),(fmpq(1,2),0),(1,0)},
    )
    for (early,late),table in zip(zip(bands,bands[1:]),expected):
        types = [(ZERO,ZERO),early,late,tuple(-v for v in early),tuple(-v for v in late)]
        actual = {abs_key(quotient(b-d,a-c)) for (a,b),(c,d) in combinations(types,2)}
        assert actual == table
        assert quotient(quotient(late[1],late[0]),quotient(early[1],early[0])) == 2*X
    print("Exact FLINT polynomial arithmetic recovered both four-value secant tables by interpolation at rational points of the original tents.")


if __name__ == "__main__":
    main()
