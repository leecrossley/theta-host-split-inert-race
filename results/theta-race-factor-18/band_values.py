"""Normalized move intervals at rational abscissae from the integer segmented sieve."""
from flint import arb
from finite import sieve, box


def values(binary, points, denominator, *, fine=False):
    edges = sorted({(m*k)//denominator for k in points for m in (1, 2, 4, 8)})
    totals = sieve(binary, edges, fine=fine)
    root = arb(2).sqrt()
    normal = 8/(2*root*(root-1)**2*arb(2).log())
    def weights(lo, hi):
        a, b, c, d = [r-l for l, r in zip(totals[lo], totals[hi])]
        return box(a, b, 58 if fine else 40), box(c, d, 96 if fine else 72)
    def score(k):
        edges = [(m*k)//denominator for m in (1, 2, 4, 8)]
        y = arb(k)/denominator
        first, middle, last = [weights(a, b) for a, b in zip(edges, edges[1:])]
        def tent(y, inside, outside):
            return 2*inside[0]-2*y*inside[1]+4*y*outside[1]-outside[0]
        return normal*(tent(2*y,middle,last)/root-tent(y,first,middle))/y.sqrt()
    return {k: score(k) for k in points}
