"""Check exact rational observation schedules by prime exponents."""
from fractions import Fraction
from itertools import combinations
from math import gcd, lcm
from horizon import order,prime_list,symbol


def prime_terms(left,right,primes,cap):
    terms={}
    for p in primes:
        if p in (left,right):
            continue
        a,b=symbol(p,left),symbol(p,right)
        d=order(p)
        f=d if d>1 or a==b else 2
        if p**f<=cap:
            terms[p]=(p**f,-a**f)
    return terms


def exponents(y,terms):
    result={}
    for p,(n,s) in terms.items():
        w=max(0,min(2-2*y/n,4*y/n-1))
        if w:
            result[p]=s*w
    return result


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    p,r=547,727
    n,m=p**order(p),r**order(r)
    assert m<n/4 and n==89526025681 and m==528529
    for t in (m,n):
        y=Fraction(t+1,4)
        assert Fraction(t,4)<y<t and (4*y/t-1)==Fraction(1,t)
    cap=240
    primes=prime_list(cap)
    pairs=list(combinations((3,7,11,19,23,31,43,47),2))
    target=prime_terms(3,7,primes,cap)
    # Both diagnostic bands for 3 and 7 contain the rational height 25.
    y=Fraction(25)
    original=exponents(y,target)
    for a,b in pairs:
        if (a,b)==(3,7):
            continue
        other=exponents(y,prime_terms(a,b,primes,cap))
        assert any(other.get(q,0)!=original.get(q,0) for q in (3,7))
    # Exponentiating a cleared log sum gives this rational product exactly.
    low=exponents(Fraction(25,3),target)
    scale=lcm(*(v.denominator for v in low.values()))
    product=Fraction(1)
    for t,v in low.items():
        product*=Fraction(t)**int(scale*v)
    for t,v in low.items():
        def valuation(num):
            count=0
            while num%t==0:
                num//=t
                count+=1
            return count
        assert valuation(product.numerator)-valuation(product.denominator)==int(scale*v)
    modulus=16*3*7*23
    cls=next(a for a in range(1,modulus,16) if gcd(a,modulus)==1 and
             (symbol(a,3),symbol(a,7),symbol(a,23))==(1,1,-1))
    for k in range(10):
        t=cls+k*modulus
        assert (symbol(t,3),symbol(t,7),symbol(t,23))==(1,1,-1)
    print("The two target diagnostic bands are disjoint and each final sample lies strictly in its own band.")
    print("A single rational observation distinguishes the 3/7 target from 27 comparison pairs at a point where both diagnostic bands overlap.")
    print("Exact integer factorization recovers a sample's prime coefficients. A reduced residue class separates the 3/7 and 3/23 high tails.")


if __name__=="__main__":
    main()
