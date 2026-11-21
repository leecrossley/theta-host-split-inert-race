"""Compare forward and backward schedules on every centilog cell, including both endpoints."""
from array import array
from fractions import Fraction
from math import ceil, floor
from flint import arb, ctx

from ordered_averages import supports

FIRST, LAST, EDGE = 104, 1250104, 1254900
BAD = EDGE+1
EXPECTED = (483,738,953)
COVERED = (1249999,1241977,977328)


def forward(tables):
    ends = [array("i", [BAD])*(EDGE+2) for _ in range(2)]
    for top, sign, m, k in tables:
        radius = Fraction(top+82,16984)+Fraction(2,k)
        x = Fraction(m,64)
        left, right = floor(100*(x-radius)),ceil(100*(x+radius))
        assert FIRST <= left < right <= EDGE
        target = ends[int(sign == 1)]
        target[left] = min(target[left],right)
    for target in ends:
        for j in range(EDGE-1,-1,-1):
            target[j] = min(target[j],target[j+1])
        for j in range(FIRST,2177):
            target[j] = min(target[j],j+220)
    neg,pos = ends
    both = array("i",map(max,neg,pos))
    def earliest(j):
        twice = both[both[j]]
        pnp,npn = pos[neg[pos[j]]],neg[pos[neg[j]]]
        return both[j],min(twice,pnp,npn),min(both[twice],neg[pnp],pos[npn])
    return earliest


def backward(tables):
    starts = [array("i",[0])*(EDGE+2) for _ in range(2)]
    for top,sign,m,k in tables:
        denominator = 64*16984*k
        midpoint = 100*m*16984*k
        offset = 6400*((top+82)*k+2*16984)
        left, right = (midpoint-offset)//denominator, -(-(midpoint+offset)//denominator)
        target = starts[(sign+1)//2]
        target[right] = max(target[right],left)
    for t in starts:
        for j in range(FIRST+220,2177+220):
            t[j] = max(t[j],j-220)
        for j in range(1,EDGE+1):
            t[j] = max(t[j],t[j-1])
    neg,pos = starts
    both = array("i",map(min,neg,pos))
    latest = [array("i") for _ in range(3)]
    for j in range(EDGE+1):
        n,p = neg[j],pos[j]
        npn,pnp = neg[pos[n]],pos[neg[p]]
        twice = both[both[j]]
        latest[0].append(both[j])
        latest[1].append(max(twice,npn,pnp))
        latest[2].append(max(both[twice],pos[npn],neg[pnp]))
    return latest


def check():
    if not __debug__:
        raise RuntimeError("Use the verifier with assertions enabled.")
    ctx.prec = 144
    offset = 3*arb(2).log()/2
    assert 103 < 100*offset < FIRST and FIRST/arb(100) > offset
    assert LAST-1 < 100*(12500+offset) < LAST
    assert arb(2176)/100-offset < arb(10**9).log()
    assert arb(9).log() < arb(220)/100
    assert arb(529)/100 < arb(200).log()
    tables = supports()
    earliest = forward(tables)
    latest = backward(tables)
    worst, covered = [0]*3,[0]*3
    for i in range(FIRST,LAST+1):
        for mode,j in enumerate(earliest(i)):
            assert i < j <= EDGE
            # The first deadline in the backward schedule is exactly this one.
            assert latest[mode][j-1] < i <= latest[mode][j]
            span = j-i+1
            worst[mode] = max(worst[mode],span)
            if FIRST < i < LAST and span <= 529:
                covered[mode] += 1
    assert tuple(worst) == EXPECTED and tuple(covered) == COVERED
    assert all(2*(arb(d)/100).exp() < c for d,c in zip(worst,(251,3208,27534)))
    print("Forward earliest finishes and backward latest starts agree on all 1,250,001 centilog cutoffs.")
    print("The three maximal certified log starting-widths are 4.83, 7.38 and 9.53, giving completion factors 251, 3208 and 27534.")
    print("At factor 400 the stronger two conclusions cover at least 1,241,977 and 977,328 whole centilog cells in a log range of length 12,500.")


if __name__ == "__main__":
    check()
