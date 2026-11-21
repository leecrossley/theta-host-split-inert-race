"""Verify local Frobenius alternatives and three completely explicit indistinguishable prefixes."""
from argparse import ArgumentParser
from itertools import product
from math import gcd, prod
from pathlib import Path
from flint import fmpz

EXAMPLES=((3,7),(23,3),(31,3))
ROOT=Path(__file__).resolve().parent


def order(p):
    return next(d for d in (1,2,4) if pow(p,d,16) in (1,15))


def symbol(p,q):
    return int(fmpz(p).jacobi(q))


def prime_list(n):
    return [k for k in range(3,n+1,2) if fmpz(k).is_prime()]


def conductor_class(p,r,primes):
    # At quartic order 2 or 4 every good conductor pair gives the same term.
    m=16*prod(t for t in primes if t not in (p,r) and order(t)==1)
    a=p%m
    if order(p)==1:
        # Since p and the replacement are 3 modulo 4, reciprocity reverses the sign.
        b=next(k for k in range(1,p) if symbol(k,p)==symbol(p,r))
        a+=m*((b-a)*pow(m,-1,p)%p)
        m*=p
    assert gcd(a,m)==1
    return a,m


def mass_map(u,v,primes,limit):
    masses={}
    for t in primes:
        if t in (u,v):
            continue
        s,r=symbol(t,u),symbol(t,v)
        f=next(j for j in (1,2,4) if pow(t,j,16) in (1,15) and (s*r)**j==1)
        n=t**f
        if n<=limit:
            masses[n]=-8*s**f
    return masses


def witness(p,r,make,supplied):
    cap=p**max(2,order(p))
    assert r**max(2,order(r))<cap
    primes=prime_list(cap)
    a,m=conductor_class(p,r,primes)
    if make:
        q=a
        while q<=cap or not fmpz(q).is_prime():
            q+=m
    else:
        q=supplied
    assert q>cap and q%m==a and fmpz(q).is_prime()
    original=mass_map(p,r,primes,cap)
    changed=mass_map(q,r,primes,cap)
    assert sorted(n for n in original.keys()|changed.keys() if original.get(n,0)!=changed.get(n,0))==[cap]
    assert original.get(cap,0)==0 and changed[cap]==-8
    return q


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("--witnesses",action="store_true")
    opt=ap.parse_args()
    if opt.witnesses:
        for p,r in EXAMPLES:
            print(p,r,witness(p,r,True,None),sep="\t")
        return
    assert 547**4==89526025681 and 727**2==528529 and order(547)==4 and order(727)==2
    for a,s,t in product((1,3,5,7,9,11,13,15),(-1,1),(-1,1)):
        d=order(a)
        f=next(j for j in (1,2,4) if pow(a,j,16) in (1,15) and (s*t)**j==1)
        if d>1:
            assert f==d and s**f==1
        else:
            assert f==(1 if s==t else 2) and abs(s**f)==1
    rows=[tuple(map(int,row.split())) for row in (ROOT/"evidence"/"twins.tsv").read_text().splitlines()]
    assert [(p,r) for p,r,_ in rows]==list(EXAMPLES)
    for p,r,q in rows:
        witness(p,r,False,q)
    print("All 32 residue and sign combinations obey the claimed local detection rule.")
    print("Explicit prime-conductor alternatives agree at every smaller norm and first disagree at 81, 529 and 961.")
    print("The target thresholds are 528529/4 for 727 and 89526025681/4 for 547.")


if __name__=="__main__":
    main()
