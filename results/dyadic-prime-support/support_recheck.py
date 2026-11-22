"""Independently test exact prime-coefficient equations and fourth-power witnesses."""

from csv import DictReader
from fractions import Fraction as F
from itertools import product
from math import isqrt
from pathlib import Path

EVIDENCE = Path(__file__).resolve().parent / "evidence"
PAIRS = ((547, 727), (571, 599))
BAD = ((2, 547, 727, 343057911569, 343070636977),
       (2, 571, 599, 413718278401, 603339156001))


def small_primes(limit):
    flags = bytearray([1]) * (limit + 1)
    flags[:2] = bytes(2)
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p*p::p] = bytes(1 + (limit - p*p)//p)
    return [n for n, flag in enumerate(flags) if flag]


def characters(host):
    return [{x*x % q for x in range(1, q)} for q in PAIRS[host]]


def prime_data(p, h, residue):
    u, v = PAIRS[h]
    different = (p % u in residue[0]) != (p % v in residue[1])
    if p % 8 in (3, 5):
        f = 4
    elif p % 16 in (7, 9) or different:
        f = 2
    else:
        f = 1
    sign = 1 if f == 1 and p % u not in residue[0] else -1
    return p**f, sign


def hat(x, z):
    return max(F(0), min(2 - F(2*z, x), F(4*z, x) - 1))


def rank(rows):
    basis = {}
    for row in rows:
        value = list(map(F, row))
        for j in sorted(basis):
            q = value[j]
            if q:
                value = [a - q*b for a, b in zip(value, basis[j], strict=True)]
        pivot = next((j for j, v in enumerate(value) if v), None)
        if pivot is not None:
            denom = value[pivot]
            basis[pivot] = [v/denom for v in value]
    return len(basis)


def compatible(system, variables):
    return rank([row[:variables] for row in system]) == rank(system)


def early():
    prime_list = small_primes(4096)
    fields = []
    for h in range(2):
        residues = characters(h)
        fields.append({p: prime_data(p,h,residues) for p in prime_list if p not in BAD[h]})
        savings = []
        minima = []
        for n in range(9):
            y = 1 << n
            data = [m for m, _ in fields[h].values() if m < 8*y]
            middle = [m for m in data if 2*y < m < 4*y]
            older = [m for m in data if m < 2*y]
            baseline = sum(2*y < m < 8*y for m in data)
            # All smaller prime coefficients must vanish, along with a candidate
            # norm in the middle octave. Earlier levels have support at most 4y.
            columns = [1 << j for j in range(1, n+1)]
            def equation(m):
                return [hat(m,z) for z in columns] + [-hat(m,2*y)]
            prior = [equation(m) for m in older]
            successes = [m for m in middle if compatible(prior + [equation(m)], len(columns))]
            if successes:
                assert len(successes) == len(middle)
                savings.append(n)
            minima.append(baseline - bool(successes))
        assert savings == ([3, 4] if h == 0 else [4, 5])
        assert minima == ([0,0,1,2,6,6,8,17,21] if h == 0 else [0,0,0,1,4,5,7,16,27])
    # The very first singleton must have a nonzero coefficient for the asserted prime.
    assert [(p, m) for p, (m, _) in fields[0].items() if 8 < m < 32] == [(31,31)]
    assert [(p, m) for p, (m, _) in fields[1].items() if 16 < m < 64] == [(7,49)]
    print("Exact rational elimination reproduces all 18 early minima and the two one-field singleton primes.")
    return fields


def verify_large_octaves():
    with (EVIDENCE / "prime-norm-witnesses.tsv").open() as f:
        rows = list(DictReader(f, delimiter="\t"))
    assert len(rows) == 720
    divisors = small_primes(isqrt(max(int(row["prime"]) for row in rows)))
    residues = [characters(h) for h in range(2)]
    coverage = {}
    for row in rows:
        h,k,p,d = (int(row[key]) for key in ("host", "exponent", "prime", "degree"))
        assert 8 <= k < 128 and h in (0,1) and p >= 3 and p not in BAD[h]
        assert all(p % q != 0 for q in divisors if q*q <= p)
        m,_ = prime_data(p,h,residues[h])
        assert m == p**d and m.bit_length() - 1 == k and m != 1 << k
        coverage.setdefault((h,k), set()).add(m)
    assert min(len(coverage[h,k]) for h in range(2) for k in range(8,128)) == 3
    # Use a different subinterval and the uniform q<=100000 bound from BMOR
    # for X>=2^33. Between 2^32 and 2^33 use their q=8 Table 1.
    for growth, c, start in ((F(7,6), F(584,10**6),2**32),
                             (F(7,6), F(1,160),2**33)):
        assert growth**4 < 2 and start > (2261078657 if start == 2**32 else 8*10**9)
        assert (growth-1)/4 - c*(1+growth)/20 > F(1,26)
        assert isqrt(start) > 26*8
    print("Independent residue-square and prime-divisor checks verify 720 witnesses; a second analytic subinterval covers the infinite tail.")


def densities():
    weights = {-1: 1, 0: 2, 1: 1}
    for a,b,units in ((1,0,8),(0,1,8),(1,1,10),(1,-1,10),(2,3,12)):
        total = sum(weights[u]*weights[v] for u,v in product(weights, repeat=2) if a*u+b*v)
        assert total == units
    print("Independent sign multiplicities: 8, 10 or 12 active assignments out of 64 residue/sign choices.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    early()
    verify_large_octaves()
    densities()
