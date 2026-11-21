"""Bound the prime-power correction using measured degree-one Chebyshev sums."""

from functools import lru_cache
from math import isqrt
from pathlib import Path
import subprocess
import sys

from flint import arb, ctx

CAP = 2_700_000_000
OMITTED = (2, 547, 727, 343057911569, 343070636977)
# Thresholds include slack for a bin that meets a boundary.
STRATA = ((10001, 70), (32768, 62), (262200, 32), (1048800, 8), (4200000, 7),
          (33560000, 3), (67200000, 2), (135000000, 1))
SCALE = 1000


def edges():
    ends = [2]
    while ends[-1] < CAP:
        left = ends[-1]
        ends.append(min(CAP, left + max(1, left//8192)))
    return ends


def certify_counts(iterator):
    ctx.prec = 128
    lower, upper = arb(0), arb(0)
    start = 1
    bins = edges()
    for end in bins:
        got, count = next(iterator)
        assert got == end and 0 <= count <= end-start
        low_before = lower
        lower += 8*count*arb(start+1).log()
        upper += 8*count*arb(end).log()
        # Between two endpoints the Chebyshev function cannot decrease.
        if end >= 1000:
            assert upper < 2*start
        for cutoff, num in STRATA:
            if end >= cutoff:
                eps = arb(num)/SCALE
                assert low_before > (1-eps)*end and upper < (1+eps)*start
        start = end
    assert next(iterator, None) is None
    assert len(bins) == 117010
    return len(bins)



def ordinary_theta():
    # The fourth and higher powers need only these small rational primes.
    top = 100000
    assert (arb(42) + 3*arb(2).log()).exp() < arb(top)**4
    mask = bytearray([1])*(top+1)
    mask[:2] = bytes(2)
    for p in range(2, isqrt(top)+1):
        if mask[p]:
            mask[p*p::p] = bytes(1 + (top-p*p)//p)
    value = arb(0)
    for p, yes in enumerate(mask):
        if yes:
            value += arb(p).log()
            assert value < arb(102)*p/100


def verify(binary, counts_file=None):
    ordinary_theta()
    payload = ''.join(str(n)+'\n' for n in edges())
    output = subprocess.check_output([binary, '--degree-one-counts'], input=payload, text=True)
    rows = (tuple(map(int,row.split())) for row in output.splitlines())
    count = certify_counts(rows)
    if counts_file is not None:
        Path(counts_file).write_text(output)
    print(f"The degree-one splitting counts cover every real endpoint in {count:,} bins through 2,700,000,000. The relative error for eight times the Chebyshev sum is below 0.001 above 135,000,000.")


@lru_cache(maxsize=None)
def correction(ell):
    """Uniform error when exp(ell) <= y and the square window stays below CAP."""
    x = arb(ell)
    a, r = arb(2).log(), arb(2).sqrt()
    B = 2*r*(r-1)**2
    Z = (x/2).exp()
    applicable = [arb(n)/SCALE for cut,n in STRATA if Z >= cut]
    assert applicable and ((x+3*a)/3).exp() > 1000
    square = 2*(2+r)*applicable[-1]
    cube = 2*((3*a+x)/3-x/2).exp()
    # For y <= exp(42), no rational exponent exceeds 64.
    assert (arb(42)+3*a)/a < 64
    higher = 8*arb(102)/100*sum((((3*a+x)/k-x/2).exp() for k in range(4,65)), arb(0))
    omissions = 8*sum((arb(p).log()+3*a for p in OMITTED),arb(0))/Z
    return (square + cube + higher + omissions)/(B*a)


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    verify(*sys.argv[1:])
