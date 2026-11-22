"""Exclude real exceptional zeros on the twelve real pair characters and price the strip tail."""

from fractions import Fraction as F
from itertools import combinations
from math import factorial

CONDUCTORS = (547, 727, 571, 599)
MOD8 = (0, 1, 0, -1, 0, -1, 0, 1)
CUT = 2**16
UNIT = 2**36


def quadratic(q):
    values = [-1] * q
    for j in range(1, q):
        values[j*j % q] = 1
    values[0] = 0
    return values


def positivity(u, v, twisted):
    a, b = quadratic(u), quadratic(v)
    twist = MOD8 if twisted else (1,)*8
    q = u*v*(8 if twisted else 1)
    partial = bound = positive = negative = 0
    for n in range(1, q+1):
        partial += a[n % u]*b[n % v]*twist[n % 8]
        bound = max(bound, abs(partial))
        if n < CUT:
            denominator = n*(n+1)
            if partial > 0:
                positive += partial*UNIT//denominator
            elif partial < 0:
                negative += (-partial*UNIT + denominator-1)//denominator
    assert partial == 0
    # s integral A(t) t^(-s-1) dt, separating the signs of A on each integer interval.
    lower = F(99,100)*F(positive,UNIT) - F(9,8)*F(negative,UNIT) - F(9*bound,8*CUT)
    return q, lower


def exp_lower(x, last=100):
    return sum(x**j/factorial(j) for j in range(last+1))


def strip_estimate():
    q, t, b = 7000000, 512, 12
    assert max(16*u*v for u,v in combinations(CONDUCTORS, 2)) < q
    assert F(9,8)**100 > CUT
    assert exp_lower(F(19)) > b*q
    assert exp_lower(F(45,2)) > q*(t+2)
    assert exp_lower(F(80,9)) > 7000
    # BMOR Corollary 1.2 gives 0.247 log(qU) + 6.894, hence the looser 8.949.
    # log(qmax U)>16 and pi>3, so the discarded main-term denominator only improves the bound.
    assert F(1,3) + (F(247,1000) + F(8949,16000))/12 < F(1,2)
    total_mass = F(2*b,3)*(19+F(1,3))
    high_mass = F(2*b**4,3*t**3)*(F(45,2)+F(1,3))
    assert total_mass < 155 and high_mass < F(1,400)
    # Include both the trivial-zero series and removal of proper prime powers.
    assert total_mass/7000 + high_mass + F(2,10**6) < F(1,40)
    assert F(249,2000)-F(3,80) > F(1,12)
    print("At log X >= 2000 every pair-character prime sum has absolute value < 1/40.")
    print("Each unit projection has weighted square sum > 1/12 and more than X/(144 log(2X)) visible primes.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    results = [positivity(u, v, t) for u,v in combinations(CONDUCTORS,2) for t in (False,True)]
    assert all(value > F(46,100) for _,value in results)
    print("All 12 real pair characters have L(s) > 0.46 throughout [0.99, 1].")
    print("Floors of 1000 times the uniform lower bounds: " + " ".join(str(int(1000*z)) for _,z in results))
    strip_estimate()
