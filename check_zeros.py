"""Verify the eight characters and the sign-changing rational ordinate intervals."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from math import gcd
import os
from pathlib import Path

from flint import arb, ctx, dirichlet_char

DEN = 10**11
TABLE = Path(__file__).parent / "evidence" / "zero-brackets.tsv"
# Modulus, Conrey label, odd prime, exponent of the even quartic character.
FAMILY = ((547, 546, 547, 0), (4376, 1093, 547, 2),
          (8752, 1093, 547, 1), (8752, 5469, 547, 3),
          (727, 726, 727, 0), (5816, 1453, 727, 2),
          (11632, 7269, 727, 1), (11632, 1453, 727, 3))


def read_table():
    rows = [tuple(map(int, line.split())) for line in TABLE.read_text().splitlines()]
    labels = {(q, n) for q, n, _, _ in FAMILY}
    if len(rows) != 5200 or rows != sorted(rows) or {(q,n) for q,n,_,_ in rows} != labels:
        raise ValueError("Unexpected bracket table or family")
    previous = {}
    for q, n, lo, hi in rows:
        if not previous.get((q,n), 0) < lo < hi:
            raise ValueError(f"Overlapping or nonpositive bracket: {(q,n,lo,hi)}")
        previous[q,n] = hi
    return rows


def check_family():
    quartic = {sgn * pow(5, j, 16) % 16: j for j in range(4) for sgn in (1,-1)}
    for q, label, p, j in FAMILY:
        chi = dirichlet_char(q, label)
        if not (chi.is_primitive() and chi.conductor() == q and chi.parity() == 1):
            raise ValueError(f"Conductor or parity mismatch for {q}, {label}")
        exponent = int(chi.group().exponent())
        for n in range(1, 16*p, 2):
            if gcd(n, p) != 1:
                continue
            legendre = 0 if pow(n,(p-1)//2,p) == 1 else 2
            expected = legendre + j*quartic[n % 16]
            if (4*int(chi.chi_exponent(n)) - exponent*expected) % (4*exponent):
                raise ValueError(f"Character mismatch at {(q,label,n)}")


def check_batch(rows):
    ctx.prec = 96
    chars = {}
    for q, n, lo, hi in rows:
        if (q,n) not in chars:
            chars[q,n] = dirichlet_char(q,n)
        chi = chars[q,n]
        left = chi.hardy_z(arb(lo)/DEN)
        right = chi.hardy_z(arb(hi)/DEN)
        if not (left.imag.contains(0) and right.imag.contains(0) and left.real*right.real < 0):
            raise ArithmeticError(f"Uncertified Hardy sign change at {(q,n,lo,hi)}")
    return len(rows)


def main():
    check_family()
    rows = read_table()
    work = [rows[i:i+24] for i in range(0, len(rows), 24)]
    workers = int(os.environ.get("THETA_WORKERS", min(12, os.cpu_count() or 1)))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        if sum(pool.map(check_batch, work)) != len(rows):
            raise ArithmeticError("Missing sign check")
    counts = Counter((q,n) for q,n,_,_ in rows)
    print("All eight primitive odd characters match the field characters on every unit class.")
    print("5200 disjoint-within-factor Hardy sign changes at 96 bits (10400 endpoints).")
    print("Positive brackets by conductor/label: " + ", ".join(
        f"{q}/{n}={counts[q,n]}" for q,n,_,_ in FAMILY))


if __name__ == "__main__":
    main()
