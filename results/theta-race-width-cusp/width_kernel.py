"""Complete prime-power and gamma evaluations in multiplicative coordinates."""
from bisect import bisect_right
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from pathlib import Path
from flint import arb

EVIDENCE = Path(__file__).resolve().parent / "evidence"
CONDUCTORS = (547, 4376, 8752, 8752, 727, 5816, 11632, 11632)


def check_inputs():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    assert sha256((EVIDENCE/"curvature-check.txt").read_bytes()).hexdigest() == "fdb7b01a06e080bcf003817b9ac9537ed5257343f51bf921bd305044ff7eb2ab"


def ball(x):
    x = Fraction(x)
    return arb(x.numerator)/x.denominator


def trace(p, j):
    symbols = [pow(p, (q-1)//2, q) for q in (547,727)]
    if p > 2 and pow(p, j, 16) not in (1,15):
        return 0
    return (1 if p == 2 else 4)*sum((s if s <= 1 else -1)**j for s in symbols)


class PrimeTable:
    def __init__(self, bound):
        sieve = bytearray(b"\1")*(bound+1)
        sieve[:2] = b"\0\0"
        for p in range(2, isqrt(bound)+1):
            if sieve[p]:
                sieve[p*p::p] = bytes((bound-p*p)//p+1)
        rows = []
        for p in range(2, bound+1):
            if sieve[p]:
                n, j = p, 1
                while n <= bound:
                    t = trace(p,j)
                    if t:
                        rows.append((n, p, t))
                    n, j = n*p, j+1
        self.bound = bound
        self.powers = []
        self.weights = [arb(0)]
        self.divided = [arb(0)]
        for n,p,t in sorted(rows):
            self.powers.append(n)
            w = arb(p).log()*t
            self.weights.append(self.weights[-1]+w)
            self.divided.append(self.divided[-1]+w/n)

    def at(self, b, z, c):
        assert 0 < z < 2 and b*z <= self.bound
        k = bisect_right(self.powers, b*z)
        sq = ball(z).sqrt()
        return (ball(b)*sq*self.divided[k]-self.weights[k]/sq)/c


@lru_cache(maxsize=None)
def gamma_series(z):
    z = min(z, 1/z)
    if z == 1:
        return arb(2).log()
    x = ball(z)
    return ((1+x)*(1+x).log()+(1-x)*(1-x).log())/(2*x.sqrt())


class Window:
    def __init__(self, b, table):
        self.b, self.table = Fraction(b), table
        self.r = ball(b).sqrt()
        self.c = 2*(self.r-1)**2
        self.constant = sum((q/arb.pi()).log() for q in CONDUCTORS)-8*(arb.const_euler()+arb(2).log())

    def gamma(self, z):
        b, r, c = self.b, self.r, self.c
        u = b*min(z,1/z)
        k = ball(u).sqrt()
        g = r/c*(k-1/k) if u > 1 else arb(0)
        return (self.constant*g-8/c*(r*(gamma_series(z*b)+gamma_series(z/b))-
                                      (ball(b)+1)*gamma_series(z)))

    def score_profile(self, z):
        return self.gamma(z)-self.table.at(self.b,z,self.c)-self.table.at(self.b,1/z,self.c)

    def measurements(self, ratio):
        z = Fraction(ratio)
        assert 1 < z < 2
        d = (ball(z).sqrt()+1/ball(z).sqrt())/2
        origin, shifted = self.score_profile(Fraction(1)), self.score_profile(z)
        curvature = 2*(origin-shifted)
        weighted = 2*(d*(1+origin)-(1+shifted))
        background = 2*(d-1)+2*(d*self.gamma(Fraction(1))-self.gamma(z))
        return curvature, weighted, background
