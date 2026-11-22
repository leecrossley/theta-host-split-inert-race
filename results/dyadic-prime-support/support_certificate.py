"""Verify dyadic prime-norm witnesses, early minima and the two-field density constants."""

from csv import DictReader
from fractions import Fraction
from itertools import product
from math import comb, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parent
HOSTS = ((547, 727, 343057911569, 343070636977),
         (571, 599, 413718278401, 603339156001))


def is_prime(n):
    return n > 2 and n % 2 == 1 and all(n % d for d in range(3, isqrt(n) + 1, 2))


def degree(host, p):
    u, v, *extra = HOSTS[host]
    assert p not in (2, u, v, *extra)
    cu = 1 if pow(p, (u - 1)//2, u) == 1 else -1
    cv = 1 if pow(p, (v - 1)//2, v) == 1 else -1
    return next(d for d in (1, 2, 4) if pow(p, d, 16) in (1, 15) and (cu*cv)**d == 1)


def early_counts(host, largest):
    counts = [0] * (largest + 1)
    for p in range(3, 1 << (largest + 1), 2):
        if p in HOSTS[host] or not is_prime(p):
            continue
        norm = p ** degree(host, p)
        k = norm.bit_length() - 1
        if k <= largest:
            counts[k] += 1
    return counts


def occupancy_test():
    witnesses = {}
    with (ROOT / "evidence" / "prime-norm-witnesses.tsv").open() as f:
        for row in DictReader(f, delimiter="\t"):
            h, k, p, d = (int(row[key]) for key in ("host", "exponent", "prime", "degree"))
            assert h in (0, 1) and 8 <= k < 128
            assert p not in HOSTS[h] and is_prime(p) and degree(h, p) == d
            assert (1 << k) < p**d < (1 << (k + 1))
            group = witnesses.setdefault((h, k), set())
            assert p not in group
            group.add(p)
    assert all(len(witnesses[h, k]) == 3 for h in (0, 1) for k in range(8, 128))
    print("Three distinct good norms in all 240 field/octave pairs from exponent 8 through 127.")

    # Analytic continuation. If X >= 2^32, restrict to primes 3 mod 8 in (X, 19X/16].
    # Theorem 1.2 and Table 1 of Bennett et al. bound theta(x;8,3) by x/4
    # with error 0.0005840*x/log(x) above 2261078657.
    c = Fraction(584, 10**6)
    t = Fraction(19, 16)
    assert 2**32 > 2261078657 and t**4 < 2
    assert Fraction(t-1, 4) - c*(1+t)/20 > Fraction(1, 24)
    # log X > 20 and log(tX) <= sqrt(X) for X >= 2^32.
    # sqrt(X)/24 > 8, leaving three primes after even five omissions.
    assert 2**16 > 24*8
    print("For every exponent at least 128, primes 3 mod 8 supply at least three degree-four norms.")


def exact_exceptions():
    for h, label in enumerate(("547/727", "571/599")):
        counts = early_counts(h, 10)
        exceptions = []
        minima = []
        for n in range(9):
            k = n
            while k >= 0 and counts[k] == 1:
                k -= 1
            erase = counts[n+1] > 0 and (k < 0 or counts[k] == 0)
            baseline = counts[n+1] + counts[n+2]
            minima.append(baseline - int(erase))
            if erase:
                exceptions.append(n)
        assert exceptions == ([3, 4] if h == 0 else [4, 5])
        print(f"{label} octave counts through 10: " + " ".join(map(str, counts)))
        print(f"{label} minimum support through exponent 8: " + " ".join(map(str, minima)))
        print(f"{label} saves one prime exactly at exponents " + ", ".join(map(str, exceptions)) + ".")


def two_field_densities():
    signs = tuple(product((-1, 1), repeat=4))
    for a, b, expected in ((1, 0, Fraction(1, 8)), (0, 1, Fraction(1, 8)),
                           (1, 1, Fraction(5, 32)), (1, -1, Fraction(5, 32)),
                           (1, 2, Fraction(3, 16))):
        fraction = Fraction(sum(a*(u+v) + b*(s+t) != 0 for u,v,s,t in signs), 4*len(signs))
        assert fraction == expected
    assert Fraction(comb(4, 2), 4**2) == Fraction(3, 8)
    print("Six times the three prime densities gives the sharp constants 3/4, 15/16 and 9/8.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    occupancy_test()
    exact_exceptions()
    two_field_densities()
