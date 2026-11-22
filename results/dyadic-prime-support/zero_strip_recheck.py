"""An independent derivative estimate for exceptional zeros and a second spectral mass bound."""

from fractions import Fraction as F
from itertools import combinations
from math import factorial

M = 65536
BASE = 1 << 42


def euler(q):
    return [0 if n == 0 else (1 if pow(n, (q-1)//2, q) == 1 else -1) for n in range(q)]


def derivative_check(u, v, twist):
    a, b = euler(u), euler(v)
    modulus = u*v*(8 if twist else 1)
    partial = maximum = lower = derivative = 0
    for n in range(1, modulus+1):
        value = a[n%u]*b[n%v]
        if twist:
            value *= (0 if n%2 == 0 else 1 if n%8 in (1,7) else -1)
        partial += value
        maximum = max(maximum, abs(partial))
        if n <= M:
            lower += value*BASE//n
        if n < M:
            # log(n+1) < 0.7 bit_length(n).
            numerator = abs(partial)*BASE*(10+7*n.bit_length())
            denominator = 10*n*(n+1)
            derivative += (numerator+denominator-1)//denominator
    assert partial == 0
    s = F(99,100)
    tail = F(9*maximum,8*M)*((1+16*F(7,10))/s + 1/s**2)
    slope = F(9,8)*F(derivative,BASE) + tail
    return F(lower,BASE) - F(2*maximum,M) - slope/100


def exp_below(x, n):
    return sum(x**j / factorial(j) for j in range(n+1))


def alternate_mass():
    # Use 640 as a different cut and 11.6 as a different Fourier envelope.
    q,t,b = 7000000,640,F(58,5)
    assert F(8,1)/F(69,100) < b
    assert exp_below(F(183,10),90) > q*b
    assert exp_below(F(223,10),100) > q*(t+2)
    assert exp_below(F(2000,223),70) > 7700
    # N(U) <= 0.42 U log(qmax U) for U>=11.6, from the same zero-count theorem.
    assert F(1,3) + (F(247,1000) + F(8949,16000))/b < F(42,100)
    mass = F(4,3)*F(42,100)*b*(F(183,10)+F(1,3))
    tail = F(4,3)*F(42,100)*b**4/t**3*(F(223,10)+F(1,3))
    assert mass/7700 + tail + F(2,10**6) < F(1,50)
    print("The independent zero-count integral with cutoff 640 gives absolute error < 1/50.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    vals = [derivative_check(u, v, t) for u,v in combinations((547,727,571,599),2) for t in (False,True)]
    assert min(vals)>F(1,4)
    print("Direct values at 1 and uniform derivative bounds give L(s) > 1/4 on all 12 intervals.")
    alternate_mass()
