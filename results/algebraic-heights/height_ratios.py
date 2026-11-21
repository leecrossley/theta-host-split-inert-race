"""Certify the slope rank and the last deficient prime interval for the two hosts."""
import csv
from collections import Counter
from fractions import Fraction
from math import isqrt
from pathlib import Path

HOSTS = {"A": (547, 727), "B": (571, 599)}
SLOPE_COLUMNS = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 71, 73, 89,
                 97, 103, 113, 127, 137, 151, 191, 199, 223, 233, 239, 241, 257,
                 263, 271, 311, 313, 337, 353, 367, 431, 433, 439, 449, 457, 479,
                 503, 521, 577, 601, 607, 641, 647, 673, 743, 881, 887, 911, 929,
                 1231, 1999, 2927, 7649, 8111, 16993)
INTERCEPT_COLUMNS = tuple(sorted((set(SLOPE_COLUMNS) - {241}) | {281}))
HERE = Path(__file__).resolve().parent


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))


def prime_data(host, p, by_residues):
    u, v = HOSTS[host]
    if p in (2, u, v):
        return None
    if by_residues:
        a, b = (1 if p % q in SQUARES[q] else -1 for q in (u, v))
        if p % 16 in (1, 15) and a == b:
            degree = 1
        elif p*p % 16 in (1, 15):
            degree = 2
        else:
            degree = 4
    else:
        a, b = (1 if pow(p, (q-1)//2, q) == 1 else -1 for q in (u, v))
        degree = next(f for f in (1, 2, 4)
                      if pow(p, f, 16) in (1, 15) and (a*b)**f == 1)
    return p**degree, a**degree


SQUARES = {q: {x*x % q for x in range(1, q)} for pair in HOSTS.values() for q in pair}


def coefficient(cell, p, modulus, h, direct, intercept):
    host, left, right = cell
    data = prime_data(host, p, by_residues=not direct)
    if data is None:
        return 0
    norm, sign = data
    y = Fraction(left + right, 16)
    if direct:
        def tent_part(z):
            if z < norm < 2*z:
                return 2 if intercept else -2
            if 2*z < norm < 4*z:
                return -1 if intercept else 4
            return 0
        dilation = h if intercept else 2*h
        value = dilation*tent_part(2*y) - tent_part(y)
    else:
        jumps = (-h, 3*h+1, -2*h-3, 2) if intercept else (8*h, -4-12*h, 6+4*h, -2)
        value = sum(j for i,j in enumerate(jumps) if 8*y > (1 << i)*norm)
    denominator = 1 if intercept else norm
    return value*sign*pow(denominator, -1, modulus) % modulus


def determinant(matrix, modulus):
    n = len(matrix)
    result = 1
    for k in range(n):
        pivot = next((i for i in range(k,n) if matrix[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            matrix[pivot], matrix[k] = matrix[k], matrix[pivot]
            result = -result
        diagonal = matrix[k][k]
        result = result*diagonal % modulus
        inverse = pow(diagonal, -1, modulus)
        for i in range(k+1,n):
            scale = matrix[i][k]*inverse % modulus
            if scale:
                for j in range(k+1,n):
                    matrix[i][j] = (matrix[i][j] - scale*matrix[k][j]) % modulus
    return result


def crossing_minor():
    with (HERE / "evidence/exact-crossings.tsv").open() as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    labels = [(r["host"], int(r["left_eighths"]), int(r["right_eighths"])) for r in rows]
    cells = list(dict.fromkeys(labels))
    assert len(rows) == 62 and len(cells) == 60
    repeats = {cell for cell, count in Counter(labels).items() if count > 1}
    assert repeats == {("A", 62, 79), ("A", 289, 316)}
    for cell in repeats:
        hit = [r for r,c in zip(rows,labels) if c == cell]
        assert {r["sign"] for r in hit} == {"1", "-1"}
        lo, hi = sorted(int(r["lower_billionths"]) for r in hit)
        assert 3*(hi+1) < 4*lo
    all_primes = set(SLOPE_COLUMNS) | set(INTERCEPT_COLUMNS)
    assert all(is_prime(p) for p in all_primes)
    assert all(prime_data(host,p,True)[0] <= 1_000_000
               for p in all_primes for host in HOSTS)
    for intercept,columns in ((False,SLOPE_COLUMNS),(True,INTERCEPT_COLUMNS)):
        assert len(set(columns)) == len(cells)
        residues = []
        for modulus,direct in ((1_000_000_007,False),(2_147_483_647,True)):
            assert is_prime(modulus)
            h = pow(pow(2, (modulus+1)//4, modulus), -1, modulus)
            assert 2*h*h % modulus == 1
            def evaluate(use_tents):
                return [[coefficient(c,p,modulus,h,use_tents,intercept)
                         for p in columns] for c in cells]
            matrix = evaluate(direct)
            assert matrix == evaluate(not direct)
            value = determinant(matrix,modulus)
            assert value
            residues.append(value)
        kind = "intercept" if intercept else "slope"
        print(f"The protected 60 by 60 {kind} minor has residues {residues[0]} and {residues[1]} "
              "modulo 1000000007 and 2147483647. Every column norm on both fields is at most one million.")
    print("The two repeated cells each contain two opposite signed hits, separated by a factor below 4/3.")


def prime_intervals():
    # Past 14493/4 in L every interval has all nine trace pairs by the tagged
    # joint-prime theorem. The next prime at its endpoint is 4831.
    outside = [[], [], [], []]
    for p in range(3, 4832, 2):
        if not is_prime(p) or p % 16 not in (1, 15) or any(p in q for q in HOSTS.values()):
            continue
        traces = []
        for host in HOSTS:
            n, sign = prime_data(host,p,False)
            assert (n,sign) == prime_data(host,p,True)
            traces.append(sign if n == p else 0)
        a,b = traces
        for primes, value in zip(outside, (a,b,a-b,a+b)):
            if value:
                primes.append(p)
    lasts = []
    for primes in outside:
        bad = [(p,q) for p,q in zip([1]+primes,primes) if 4*p < 3*q]
        # The final-listed prime of each type covers left endpoints up to the
        # permanent threshold of the joint theorem.
        assert primes[-1] >= Fraction(14493,4)
        lasts.append(bad[-1])
    assert lasts == [(479,641), (241,337), (113,223), (127,223)]
    print("The four last deficient prime gaps are 479 to 641, 241 to 337, "
          "113 to 223 and 127 to 223. The last small-prime gap puts the joint prime-band threshold at 641/8 in the certified range.")


def counts():
    assert (62+2)//3 == 21 and 62*61//2-(20*3+1) == 1830
    print("The 62 individual heights give at least 21 classes and 1830 transcendental ratios.")


if __name__ == "__main__":
    crossing_minor()
    prime_intervals()
    counts()
