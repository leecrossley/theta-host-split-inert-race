"""Independently enumerate all odd primes and check the degree-one square correction."""

from bisect import bisect_left
from itertools import compress
from math import isqrt
import sys

from splitting_square_bound import CAP, certify_counts, edges


def counts():
    end = edges()
    result = [0]*len(end)
    limit = isqrt(CAP)
    base = bytearray([1])*(limit+1)
    base[:2] = bytes(2)
    for p in range(2,isqrt(limit)+1):
        if base[p]:
            base[p*p::p] = bytes(1 + (limit-p*p)//p)
    divisors = [p for p in range(3,limit+1,2) if base[p]]
    euler = {q: [pow(a,(q-1)//2,q) == 1 for a in range(q)]
             for q in (547,727)}
    size = 1 << 21
    for low in range(3,CAP+1,2*size):
        high = min(low+2*size, CAP+1)
        flags = bytearray([1])*((high-low+1)//2)
        for p in divisors:
            if p*p >= high:
                break
            first = max(p*p, ((low+p-1)//p)*p)
            if first % 2 == 0:
                first += p
            index = (first-low)//2
            if index < len(flags):
                flags[index::p] = bytes(1+(len(flags)-1-index)//p)
        for residue in (1,15):
            first = low+(residue-low)%16
            offset = (first-low)//2
            for p in compress(range(first,high,16),flags[offset::8]):
                if p not in (547,727) and euler[547][p%547] == euler[727][p%727]:
                    result[bisect_left(end,p)] += 1
    return list(zip(end,result))


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError("Use the verifier with assertions enabled.")
    data = counts()
    certify_counts(iter(data))
    if len(sys.argv) == 2:
        with open(sys.argv[1]) as inp:
            other = [tuple(map(int,row.split())) for row in inp]
        assert data == other
    print("The separate odd-integer sieve with Euler’s criterion reproduced every degree-one prime-count bin and certified the Chebyshev bounds.")
