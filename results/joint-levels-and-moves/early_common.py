"""Enclose eight early two-field moves by summing actual good prime ideals."""
from bisect import bisect_left
from itertools import product
from math import isqrt

from flint import arb, ctx

from deadline_arithmetic import BAD
from deadline_spectrum import PAIRS

GROUPS=((2053,4216,26759,6417),(6569842,10490706,3315107,20790365))


def small_primes(n):
    flags=bytearray([1])*(n+1)
    flags[:2]=bytes(2)
    for p in range(2,isqrt(n)+1):
        if flags[p]:
            flags[p*p::p]=bytes((n-p*p)//p+1)
    return [p for p in range(2,n+1) if flags[p]]


def degree_one_primes(cap,small):
    for residue in (1,15):
        length=(cap-residue)//16+1
        flags=bytearray([1])*length
        if residue==1:
            flags[0]=0
        for p in small:
            if p==2:
                continue
            lower=max(0,(p*p-residue+15)//16)
            start=lower+(-residue*pow(16,-1,p)-lower)%p
            if start<length:
                flags[start::p]=bytes((length-1-start)//p+1)
        pos=flags.find(b'\x01')
        while pos>=0:
            yield 16*pos+residue
            pos=flags.find(b'\x01',pos+1)


def symbols(p):
    t=[-1]*p
    for n in range(1,p):
        t[n*n%p]=1
    t[0]=0
    return t


def sums(edges):
    cap=edges[-1]
    small=small_primes(isqrt(cap))
    tables=[tuple(map(symbols,pair)) for pair in PAIRS]
    bins=[[[arb(0),arb(0)] for _ in edges] for _ in PAIRS]
    for p in degree_one_primes(cap,small):
        value=None
        for i,((u,v),(lt,rt)) in enumerate(zip(PAIRS,tables)):
            sg=lt[p%u]
            if p in BAD[i] or sg!=rt[p%v]:
                continue
            if value is None:
                value=arb(p).log()
            pair=bins[i][bisect_left(edges,p)]
            pair[0]-=sg*value
            pair[1]-=sg*value/p
    for p in small:
        for i,((u,v),(lt,rt)) in enumerate(zip(PAIRS,tables)):
            if p in BAD[i]:
                continue
            s,t=lt[p%u],rt[p%v]
            f=next(k for k in (1,2,4) if pow(p,k,16) in (1,15) and (s*t)**k==1)
            n=p**f
            if f>1 and n<=cap:
                val=-(s**f)*arb(p).log()
                pair=bins[i][bisect_left(edges,n)]
                pair[0]+=val
                pair[1]+=val/n
    out=[]
    for h in bins:
        a=b=arb(0)
        partial={}
        for e,(c,d) in zip(edges,h):
            a,b=a+c,b+d
            partial[e]=(a,b)
        out.append(partial)
    return out


def score(y,partial):
    a,b=partial[y]
    c,d=partial[2*y]
    e,f=partial[4*y]
    r=arb(2).sqrt()
    return 8*(2*(c-a)-2*y*(d-b)+4*y*(f-d)-(e-c))/(2*r*(r-1)**2*arb(y).sqrt())


def main():
    ctx.prec=144
    edges=sorted({k*y for group in GROUPS for y in group for k in (1,2,4,8)})
    results=sums(edges)
    assert min(GROUPS[0])>2048 and 2*max(GROUPS[0])<2**18
    assert min(GROUPS[1])>arb(15).exp() and 2*max(GROUPS[1])<2048*2**18
    for group in GROUPS:
        for signs,y in zip(product((-1,1),repeat=2),group):
            for sign,tab in zip(signs,results):
                value=(score(2*y,tab)-score(y,tab))/arb(2).log()
                assert sign*value>4,(y,sign,value)
    print(f'Eight direct prime-ideal witnesses through norm {edges[-1]} cover every H from 1 through exp(15) for the two-move square.')


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    main()
