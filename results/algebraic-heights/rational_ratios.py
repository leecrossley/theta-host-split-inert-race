"""Verify 72-place height brackets directly and exclude every small-denominator ratio."""
import csv
from bisect import bisect_right
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import isqrt
from pathlib import Path

from height_ratios import prime_data
from verify_rank import primes_through

HERE = Path(__file__).resolve().parent
UNIT = 10**72
GRID = 10**105
MAX_NORM = 8_000_000


def add(x,y):
    return x[0]+y[0],x[1]+y[1]


def scale(x,k):
    a,b = k*x[0],k*x[1]
    return (a,b) if k>=0 else (b,a)


def multiply(x,y):
    terms = [a*b for a in x for b in y]
    return min(terms),max(terms)


def root(x):
    z = isqrt(x.numerator*GRID**2//x.denominator)
    return Fraction(z,GRID),Fraction(z+1,GRID)


@cache
def log_bounds(p):
    # Decimal.ln is correctly rounded to nearest. One fixed-grid unit
    # exceeds its absolute rounding error throughout the prime range.
    with localcontext() as ctx:
        ctx.prec = 128
        a,b = Decimal(p).ln().as_integer_ratio()
    k = a*GRID//b
    return k-1,k+2


def norm_sums(primes,host):
    entries = []
    for p in primes:
        local = prime_data(host,p,True)
        if local is None or local[0]>MAX_NORM:
            continue
        n,sg = local
        a,b = scale(log_bounds(p),sg)
        entries.append((n,a,b,a//n,-((-b)//n)))
    entries.sort()
    norms = [row[0] for row in entries]
    assert len(set(norms))==len(norms)
    sums = [[0] for _ in range(4)]
    for row in entries:
        for column,x in zip(sums,row[1:]):
            column.append(column[-1]+x)
    def mass(left,right,offset):
        # Prefixes enclose sums of the same rounded summands. Subtract
        # like endpoints to retain precisely the sum over this index slice.
        return tuple(Fraction(sums[k][right]-sums[k][left],GRID) for k in (offset,offset+1))
    def tent(y):
        i,j,k = (bisect_right(norms,int(m*y)) for m in (1,2,4))
        inside = add(scale(mass(i,j,0),2),scale(mass(i,j,2),-2*y))
        outside = add(scale(mass(j,k,2),4*y),scale(mass(j,k,0),-1))
        return add(inside,outside)
    return tent


def read_table(name):
    with (HERE/"evidence"/name).open() as source:
        return list(csv.DictReader(source,delimiter="\t"))


def deep_brackets():
    catalog = read_table("exact-crossings.tsv")
    table = read_table("height-brackets-72.tsv")
    assert len(catalog)==len(table)==62
    ps = list(primes_through(MAX_NORM))
    tents = {host:norm_sums(ps,host) for host in "AB"}
    root_two = root(Fraction(2))
    inverse_root = (1/root_two[1],1/root_two[0])
    gap = add(root_two,(-1,-1))
    b = multiply(scale(root_two,2),multiply(gap,gap))
    log_two = tuple(Fraction(x,GRID) for x in log_bounds(2))
    b_log = multiply(b,log_two)
    found = []
    for reference,entry in zip(catalog,table):
        assert all(reference[k]==entry[k] for k in ("host","sign","left_eighths"))
        host = entry["host"]
        sg = int(entry["sign"])
        orientation = int(reference["orientation"])
        k = int(entry["height_floor_1e72"])
        assert k//(UNIT//10**9)==int(reference["lower_billionths"])
        for step,side in ((0,-orientation),(1,orientation)):
            y = Fraction(k+step,UNIT)
            tent = tents[host]
            delta = add(multiply(inverse_root,tent(2*y)),scale(tent(y),-1))
            residual = add(scale(delta,-8),scale(multiply(b_log,root(y)),-sg*Fraction(3947,2000)))
            assert residual[0]>0 if side>0 else residual[1]<0
        found.append((k,host,sg,int(entry["left_eighths"])))
    print("Two directed original-tent endpoint signs fix all 62 exact heights in disjoint intervals of width 10^-72.")
    ordered = sorted(found)
    assert all(a[0]+1<b[0] for a,b in zip(ordered,ordered[1:]))
    return ordered


def simplest(lo,hi):
    """A candidate inside an open rational interval, certified below by its Farey parents."""
    integer = lo.numerator//lo.denominator
    if integer+1 < hi:
        return Fraction(integer+1)
    if lo==integer:
        reciprocal = 1/(hi-integer)
        return integer+Fraction(1,reciprocal.numerator//reciprocal.denominator+1)
    return integer+1/simplest(1/(hi-integer),1/(lo-integer))


def minimum_denominator(heights):
    best = None
    intervals = []
    for early,late in combinations(heights,2):
        lo,hi = Fraction(late[0],early[0]+1),Fraction(late[0]+1,early[0])
        candidate = simplest(lo,hi)
        p,q = candidate.as_integer_ratio()
        ql = pow(p,-1,q)
        pl = (p*ql-1)//q
        pr,qr = p-pl,q-ql
        assert pr*ql-pl*qr==1 and ql+qr==q
        assert Fraction(pl,ql)<=lo<candidate<hi<=Fraction(pr,qr)
        record = (q,p,early,late)
        if best is None or record<best:
            best = record
        intervals.append((lo,hi))
    q,p,early,late = best
    for lo,hi in intervals:
        near = ((lo+hi)/2).limit_denominator(q-1)
        assert not lo<near<hi
    assert len(intervals)==1891
    print(f"All 1891 ratios in increasing-height order have rational denominator at least {q} if rational.")
    print(f"First candidate allowed by the brackets: {p}/{q}; early host/sign/left eighths {early[1:]}, late {late[1:]}.")


if __name__=="__main__":
    minimum_denominator(deep_brackets())
