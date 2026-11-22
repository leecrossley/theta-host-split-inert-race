"""Independently find all rank-changing phases using Smith invariants over Q[t]."""
from fractions import Fraction as F
from itertools import pairwise
from math import isqrt
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import invariant_factors

T=s.symbols("t")
DOMAIN=s.QQ.poly_ring(T)
PAIRS=((547,727),(571,599))
LIMIT=271


def norms(u,v):
    outcomes=[]
    for p in range(3,LIMIT+1,2):
        if any(p%d==0 for d in range(2,isqrt(p)+1)):
            continue
        a=1 if pow(p,(u-1)//2,u)==1 else -1
        b=1 if pow(p,(v-1)//2,v)==1 else -1
        # Every omitted odd prime on either field is greater than 271.
        d=next(d for d in (1,2,4) if pow(p,d,16) in (1,15) and (a*b)**d==1)
        if p**d<=LIMIT:
            outcomes.append((p,p**d,-a**d))
    return outcomes


def invariant(rows):
    if not rows or not rows[0]:
        return []
    matrix=DomainMatrix.from_Matrix(s.Matrix(rows)).convert_to(DOMAIN)
    return [s.Poly(p.as_expr(),T,domain=s.QQ).monic() for p in invariant_factors(matrix) if p]


def rowlist(t,hosts,fields):
    depth=5 if t*256<=271 else 4
    observed=[]
    primes=sorted({p for h in hosts for p,_,_ in fields[h]})
    for h in hosts:
        for j in range(1,depth+2):
            z=2**j
            entries={}
            for p,m,sign in fields[h]:
                if z*t<m<4*z*t:
                    factor=s.Rational(z,m)*T
                    entries[p]=sign*((2-2*factor) if m<2*z*t else 4*factor-1)
            if entries:
                observed.append([entries.get(p,s.S.Zero) for p in primes])
    return observed,primes


def expected(t,h):
    a=(t<=F(47,32) or t>F(31,16))
    b=(t<=F(79,64) or t>F(49,32))
    if h==(0,): return [31] if a else []
    if h==(1,): return [7] if b else []
    both_small=(t<=F(97,64) or t>F(31,16))
    return [p for p,v in ((7,b),(31,both_small),(47,both_small)) if v]


def main():
    if s.__version__!="1.14.0":
        raise RuntimeError("Use SymPy 1.14.0 for this independent replay.")
    fields=[norms(u,v) for u,v in PAIRS]
    split={F(1),F(2),F(271,256)}
    for field in fields:
        for _,m,_ in field:
            split.update(F(m,2**j) for j in range(9) if 2**j<m<2**(j+1))
    cells=0
    for lo,hi in pairwise(sorted(split)):
        for h in ((0,),(1,),(0,1)):
            mid=(lo+hi)/2
            row,ps=rowlist(mid,h,fields)
            sizes=[]
            for j in range(-1,len(ps)):
                reduced=row if j<0 else [r[:j]+r[j+1:] for r in row]
                inv=invariant(reduced)
                sizes.append(len(inv))
                if inv:
                    for f,_ in inv[-1].factor_list()[1]:
                        # A rational endpoint factor can vanish only at that endpoint.
                        if f.degree()==1 and (f.eval(lo)==0 or f.eval(hi)==0):
                            continue
                        assert f.count_roots(lo,hi)==0
            dep=5 if mid*256<=271 else 4
            zero=(2-int(mid>F(31,16)),3-int(mid>F(49,32)))
            assert sizes[0]==sum(dep+1-zero[k] for k in h)
            assert [p for p,r in zip(ps,sizes[1:]) if r<sizes[0]]==expected(mid,h)
            row,ps=rowlist(lo,h,fields)
            m=s.Matrix(row).subs(T,s.Rational(lo.numerator,lo.denominator))
            rank=m.rank()
            dep=5 if lo*256<=271 else 4
            zero=(2-int(lo>F(31,16)),3-int(lo>F(49,32)))
            assert rank==sum(dep+1-zero[k] for k in h)
            iso=[]
            for j,p in enumerate(ps):
                a=m.copy();a.col_del(j)
                if a.rank()<rank:
                    iso.append(p)
            assert iso==expected(lo,h)
            cells+=1
    assert cells==39
    print("SymPy 1.14.0 Smith factors and real-root counts certify all 39 field/cell cases, with direct ranks at every included endpoint.")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
