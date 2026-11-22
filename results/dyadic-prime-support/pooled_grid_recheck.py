"""Independent interpolation, prime counts and analytic constants for the pooled theorem."""
from bisect import bisect_left,bisect_right
from collections import defaultdict
from csv import DictReader
from fractions import Fraction as F
from itertools import combinations_with_replacement,pairwise,product
from math import factorial,gcd,lcm
from pathlib import Path

from minimum_grids_recheck import HOSTS,evaluate,primality
from phase_rank import bernstein

EVIDENCE = Path(__file__).resolve().parent/"evidence"


def integral_row(row):
    scale = lcm(*(F(x).denominator for x in row))
    return [int(scale*x) for x in row]


def rank(input_rows):
    basis = {}
    for source in input_rows:
        row = integral_row(source)
        for k,prior in sorted(basis.items()):
            if row[k]:
                c,d = prior[k],row[k]
                row = [c*x-d*y for x,y in zip(row,prior)]
                divisor = gcd(*row)
                if divisor>1:
                    row = [x//divisor for x in row]
        pivot = next((j for j,x in enumerate(row) if x),None)
        if pivot is not None:
            basis[pivot] = row
    return len(basis)


def det(rows):
    # Keep the original column order and evaluate by successive fraction-free pivots.
    m = [list(map(F,row)) for row in rows]
    n = len(m)
    sign = prev = 1
    for k in range(n-1):
        p = next((i for i in range(k,n) if m[i][k]),None)
        if p is None:
            return F(0)
        if p!=k:
            m[p],m[k] = m[k],m[p]
            sign = -sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                m[i][j] = (m[k][k]*m[i][j]-m[i][k]*m[k][j])/prev
        prev = m[k][k]
    return sign*m[-1][-1]


def coefficient(norm,sign,z):
    return F(sign,norm)*(max(0,4*z-norm)-6*max(0,z-F(norm,2))+2*max(0,z-norm))


def classify():
    data = [{p:evaluate(p,h) for p in range(3,431,2) if primality(p) and evaluate(p,h)[0]<431} for h in (0,1)]
    primes = sorted(set(data[0])|set(data[1]))
    bounds = {F(1),F(2)}
    masses = {m for field in data for m,_ in field.values()}|{431}
    for m in masses:
        bounds.update(F(m,2**j) for j in range(10) if 2**j<=m<=2**(j+1))
    sides = [(l,l) for l in sorted(bounds)[:-1]]+list(pairwise(sorted(bounds)))

    def shape(cell):
        left,right = cell
        out = []
        for j in range(8):
            if 4*right*2**j>431:
                continue
            for field in data:
                pair = [[coefficient(*field[p],a*2**j) if p in field else 0 for p in primes] for a in cell]
                if any(x for row in pair for x in row):
                    out.append(pair)
        return out

    entries = [shape(cell) for cell in sides]
    cases = {}
    with (EVIDENCE/"pooled-grid-minors.tsv").open() as stream:
        for row in DictReader(stream,delimiter="\t"):
            key = tuple(int(row[k]) for k in ("cell_a","cell_b","target"))
            assert key not in cases
            cases[key] = ([primes.index(int(v)) for v in row["dependent_primes"].split(",")],
                          [int(v) for v in row["minor_rows"].split(",")])
    checked = 0
    for i,k in combinations_with_replacement(range(len(sides)),2):
        cells = (sides[i],sides[k])
        rows = [(s,pair) for s,n in enumerate((i,k)) for pair in entries[n]]
        endpoints = [integral_row(v) for s,pair in rows for v in pair]
        r = rank(endpoints)
        for target in (17,337,353,367):
            j = primes.index(target)
            if (i,k,target) not in cases:
                assert rank(endpoints+[[int(x==j) for x in range(len(primes))]])>r
                continue
            cols,chosen = cases.pop((i,k,target))
            rest = [c for c in cols if c!=j]
            n = len(rest)
            allowed = [(s,pair) for s,pair in rows if any(any(v[c] for c in cols) for v in pair)]
            dims = [0 if c[0]==c[1] else min(n+1,sum(v==s for v,_ in allowed)) for s,c in enumerate(cells)]
            if len(allowed)>n:
                samples = [[F(a+1,b+2) for a in range(b+1)] for b in dims]
                for point in product(*samples):
                    full = [[l[c]+point[s]*(r[c]-l[c]) for c in cols] for s,(l,r) in allowed]
                    assert rank(full)<=n
            minor = [rows[t] for t in chosen]
            degree = [0 if c[0]==c[1] else sum(v==s for v,_ in minor) for s,c in enumerate(cells)]
            nodes = [[F(t,d) for t in range(d+1)] if d else [F(0)] for d in degree]
            table = []
            for a in nodes[0]:
                line = []
                for b in nodes[1]:
                    point = (a,b)
                    line.append(det([[l[c]+point[s]*(r[c]-l[c]) for c in rest] for s,(l,r) in minor]))
                table.append(bernstein(line))
            weights = [x for column in zip(*table) for x in bernstein(column)]
            assert any(weights) and min(weights)*max(weights)>=0, (i,k,target)
            checked += 1
    assert not cases and checked==257
    a12 = lambda host,z,p: coefficient(*data[host][p],F(z)) if p in data[host] else 0
    for p in primes:
        signed = (-2*a12(0,12,p)+5*(a12(0,24,p)-a12(1,24,p))-6*(a12(0,28,p)-a12(1,28,p)))
        assert signed==F(90,97)*int(p==97)
    print("Independent Euler symbols, hinge weights and rational ranks, with the bivariate Bernstein helper from phase_rank.py, verify the entire pooled table.")


def real_octaves():
    directions = defaultdict(list)
    for p in range(3,26001,2):
        if p in set(HOSTS[0])|set(HOSTS[1]) or not primality(p):
            continue
        vector = tuple(s if m==p else 0 for m,s in (evaluate(p,0),evaluate(p,1)))
        if vector!=(0,0):
            directions[max(vector,tuple(-v for v in vector))].append(p)
    last = None
    for y in range(2,52001):
        # Here y=4t. Fractional and integral endpoints are tested by integer comparisons.
        counts = [bisect_left(ps,(y+1)//2)-bisect_right(ps,y//4) for ps in directions.values()]
        remaining = sorted(max(n-1,0) for n in counts)
        if sum(remaining[:-1])<2:
            last = y
    assert last==862
    print("Bisection in four independent prime lists gives the same sharp endpoint 431/2.")


def improved_envelope():
    # A separate majorant uses a smaller base and a different decay exponent.
    b,nu,Q = F(231,40),F(13,8),7_000_000
    log_hi = F(438,25)
    assert 4/F(693,1000)<b
    assert F(2719,1000)**35<(Q*b)**2
    taylor = lambda x,n: sum(x**k/factorial(k) for k in range(n+1))
    assert taylor(log_hi,110)>Q*b
    R = F(4823,500)
    assert taylor(F(1095)/(R*log_hi),80)>640
    assert F(1095)/(R*F(35,2)**2)<2-nu
    mass = (nu*b/((nu-1)*F(314,100))*(log_hi-F(14,5)+1/(nu-1))
            +F(247,1000)*(log_hi+1/nu)+F(6894,1000))
    assert mass/(2*640) + F(3,10**6) < F(1,14)
    assert (F(249,2000)-F(3,28))/6 > F(1,500)
    print("A second triangular-kernel majorant gives error < 1/14 per twist at log X >= 1095.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    classify()
    real_octaves()
    improved_envelope()
