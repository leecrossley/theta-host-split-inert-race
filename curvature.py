"""Archimedean correction and exact rational sampling locations."""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from flint import arb

EVIDENCE = Path(__file__).resolve().parent/"evidence"
DEN = 10**36
CASES = (2, 3, 5, 7, 11, 17, 31, 47, 547, 727)
CONDUCTORS = (547, 4376, 8752, 8752, 727, 5816, 11632, 11632)


def check_inputs():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    assert sha256((EVIDENCE/"peak-check.txt").read_bytes()).hexdigest() == "a8892af30c7063d10f575e9a3ca1f671ef5695f614e6b3cf03559f33a7f1c285"


def ball(x):
    q = Fraction(x)
    return arb(q.numerator)/q.denominator


def locations(n):
    step = Fraction(4*n+1, 4*n)
    return Fraction(n)/step, Fraction(n), n*step


class GammaTerms:
    def __init__(self):
        self.a = arb(2).log()
        self.root = arb(2).sqrt()
        self.c = 2*(self.root-1)**2
        self.b = self.root*self.c
        self.constant = sum((q/arb.pi()).log() for q in CONDUCTORS)-8*(arb.const_euler()+self.a)
        self.origin = self.at(Fraction(1))

    def tent(self, x):
        z = min(x, 1/x)
        if z <= Fraction(1, 2):
            return arb(0)
        a = (2*ball(z)).sqrt()
        return 2*(a-1/a)/self.b

    def series(self, x):
        ratio = min(x, 1/x)
        if ratio == 1:
            return self.a
        z = ball(ratio)
        return ((1+z)*(1+z).log()+(1-z)*(1-z).log())/(2*z.sqrt())

    def at(self, x):
        return (self.constant*self.tent(x)-8/self.c *
                (self.root*(self.series(2*x)+self.series(x/2))-3*self.series(x)))


def read_samples():
    rows = {}
    for line in (EVIDENCE/"curvatures.tsv").read_text().splitlines():
        p, m, l, c, u = map(int, line.split())
        assert (p, m) not in rows
        rows[p,m] = [arb(2*v+1, 1)/(2*DEN) for v in (l, c, u)]
    expected = {(p, m) for p in CASES for m in ((1,2) if p in (2,547,727) else (1,2,4))}
    assert set(rows) == expected
    return rows
