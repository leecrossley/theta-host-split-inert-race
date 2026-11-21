"""Certify bounded product relations for 62 exact hits and the first eleven primes."""
import csv
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import chain
from math import isqrt
from pathlib import Path

from height_ratios import prime_data
from rational_ratios import deep_brackets

HERE = Path(__file__).resolve().parent
PRIMES = (2,3,5,7,11,13,17,19,23,29,31)
ROWS = ((156,1,32,11),(240,20,511,90),(360,1000,20606,697),(500,90000,2003443,72212))
LOG_DIGITS = 505


def log_bins():
    with (HERE/"evidence/height-brackets-512.tsv").open() as file:
        roots = [int(row["height_floor_1e512"]) for row in csv.DictReader(file,delimiter="\t")]
    assert len(roots)==62
    bins = []
    # Inputs have at most 519 digits. Each ln has absolute rounding error
    # below 10^-560, since the largest logarithm is below 14.
    with localcontext() as ctx:
        ctx.prec = 575
        epsilon = Decimal(1).scaleb(-545)
        intervals = ((Decimal(k).scaleb(-512),Decimal(k+1).scaleb(-512)) for k in roots)
        for lo,hi in chain(intervals,((Decimal(p),Decimal(p)) for p in PRIMES)):
            left = int((lo.ln()-epsilon).scaleb(LOG_DIGITS))
            right = int((hi.ln()+epsilon).scaleb(LOG_DIGITS))
            assert left==right and left>0
            bins.append(left)
    return bins


def leading_determinants(matrix):
    """Bareiss elimination. Row pivoting is only needed for a basis change."""
    a = [row[:] for row in matrix]
    n = len(a)
    previous = sign = 1
    out = []
    for k in range(n-1):
        if not a[k][k]:
            i = next(i for i in range(k+1,n) if a[i][k])
            a[i],a[k] = a[k],a[i]
            sign = -sign
        pivot = a[k][k]
        out.append(sign*pivot)
        for i in range(k+1,n):
            row = a[i]
            for j in range(k+1,n):
                row[j],remainder = divmod(pivot*row[j]-row[k]*a[k][j],previous)
                assert remainder==0
        previous = pivot
    return out+[sign*a[-1][-1]]


def check_lattice(bins,exponent,height,expected_gso,expected_margin):
    n = len(bins)
    denominator = 10**(LOG_DIGITS-exponent)
    rounded = []
    loss = 0
    for k in bins:
        z = (2*k+denominator)//(2*denominator)
        assert (2*(k+1)+denominator)//(2*denominator)==z
        rounded.append(z)
        loss += max(abs(k-z*denominator),abs(k+1-z*denominator))
    with (HERE/"evidence"/f"product-basis-{exponent}.tsv").open() as file:
        change = [list(map(int,row)) for row in csv.reader(file,delimiter="\t")]
    assert n==73 and len(change)==n and all(len(row)==n for row in change)
    assert abs(leading_determinants(change)[-1])==1
    lattice = [row+[sum(a*b for a,b in zip(row,rounded))] for row in change]
    gram = [[sum(a*b for a,b in zip(x,y)) for y in lattice] for x in lattice]
    minors = leading_determinants(gram)
    assert all(d>0 for d in minors)
    lower_square = min(Fraction(d,c) for d,c in zip(minors,[1]+minors[:-1]))
    g = isqrt(lower_square.numerator//lower_square.denominator)
    margin = isqrt(g*g-n*height**2)-(height*loss+denominator-1)//denominator
    assert g==expected_gso and margin==expected_margin>0
    print(f"Grid 10^{exponent}: Gram-Schmidt length floor {g}; for coefficient box {height}, every nonzero log sum has absolute value at least {margin}/10^{exponent}.")


def centres():
    for left,right,primes in ((62,79,(7,31)),(289,316,(3,7))):
        y = Fraction(left+right,16)
        norms = [prime_data("A",p,False)[0] for p in primes]
        assert len(set(norms))==2 and all(y<n<8*y for n in norms)
        assert all(not any(left<k*n<right for k in (1,2,4,8)) for n in norms)
    assert 90001**73>10**361 and 90001**62>10**307
    print("The two paired A cells each have at least two distinct active norms. Their opposite hits have the exact zero as geometric mean.")


def main():
    deep_brackets(512)
    bins = log_bins()
    for row in ROWS:
        check_lattice(bins,*row)
    centres()


if __name__=="__main__":
    main()
