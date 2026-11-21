"""Decode local prime factors from three actual finite width-two race samples."""
from bisect import bisect_right
from collections import Counter
from hashlib import sha256
from math import isqrt
from pathlib import Path
import sys
from flint import arb, ctx

LIMIT = 1_000_000
TYPE_PRIMES = (343057911569,343070636977)
OMITTED = {2,547,727,*TYPE_PRIMES}
EXAMPLES = {17,31,47,97,113,193,239,257,337,401,449,577,641,1009,65537,999983}
EVIDENCE = Path(__file__).resolve().parent/"evidence"


def primes_to(limit):
    prime = bytearray(b"\1")*(limit+1)
    prime[:2] = b"\0\0"
    for p in range(2,isqrt(limit)+1):
        if prime[p]:
            prime[p*p::p] = bytes((limit-p*p)//p+1)
    return [p for p in range(2,limit+1) if prime[p]]


def symbol(p,q):
    u = pow(p,(q-1)//2,q)
    assert u in (1,q-1)
    return 1 if u==1 else -1


class GoodPrimeScores:
    def __init__(self,primes,limit):
        terms = []
        for p in primes:
            if p in OMITTED:
                continue
            s,t = symbol(p,547),symbol(p,727)
            f = next(j for j in (1,2,4) if pow(p,j,16) in (1,15) and (s*t)**j==1)
            if p**f <= limit:
                terms.append((p**f,-8*s**f*arb(p).log()))
        self.ns = []
        self.sums = [arb(0)]
        self.ratios = [arb(0)]
        for n,w in sorted(terms):
            self.ns.append(n)
            self.sums.append(self.sums[-1]+w)
            self.ratios.append(self.ratios[-1]+w/n)

    def at_tick(self,k):
        # k=4y. The weight at all excluded outer endpoints is zero.
        l,m,u = [bisect_right(self.ns,i) for i in (k//4,k//2,k)]
        a,b = self.sums,self.ratios
        return 2*(a[m]-a[l])-arb(k)/2*(b[m]-b[l])+k*(b[u]-b[m])-(a[u]-a[m])


def factor(p,u):
    # Return f, sign, exponent for the denominator (1-sign*T^f)^exponent.
    if p == 2:
        return 2,1,1
    if p == 547:
        return 4,1,1
    if p == 727:
        return 2,1,2
    if p in TYPE_PRIMES:
        return 1,1,8
    order = next(k for k in (1,2,4) if pow(p,k,16) in (1,15))
    if order > 1:
        assert u == 0
        return order,1,8//order
    assert u in (-8,0,8)
    return (2,1,4) if u==0 else (1,u//8,8)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 192
    primes=primes_to(LIMIT)
    scores=GoodPrimeScores(primes,LIMIT)
    counts,digest=Counter(),sha256()
    examples=[]
    for p in primes:
        if p in OMITTED:
            f,s,d = factor(p,None)
        else:
            samples=[scores.at_tick(k) for k in (p-1,p,p+1)]
            u=-p*(samples[0]-2*samples[1]+samples[2])/arb(p).log()
            tr=u.unique_fmpz()
            assert tr is not None and abs(u-int(tr))<arb(10)**-37
            f,s,d=factor(p,int(tr))
            if p in EXAMPLES:
                midpoints=[]
                for v in samples:
                    m=(v*10**24+arb(1)/2).floor().unique_fmpz()
                    assert m is not None and abs(v*10**24-int(m)) < arb(1)/2
                    midpoints.append(str(m))
                examples.append("\t".join([str(p),*midpoints,str(tr)]))
        line=f"{p}\t{f}\t{s}\t{d}\n"
        counts[f,s,d]+=1
        digest.update(line.encode())
        if sys.argv[1:]==["--all"]:
            print(line,end="")
    assert len(primes)==78498
    sampledata="\n".join(examples)+"\n"
    if sys.argv[1:]==["--examples"]:
        print(sampledata,end="")
        return
    assert sampledata==(EVIDENCE/"finite-samples.tsv").read_text()
    assert digest.hexdigest()==(EVIDENCE/"factor-digest.txt").read_text().strip()
    if sys.argv[1:]==["--all"]:
        return
    assert not sys.argv[1:]
    print("All 78,498 primes through 1,000,000 decoded. 78,495 good odd primes were evaluated from their three finite scores.")
    print("Ordered (p, residue degree, sign, exponent) rows SHA-256: "+digest.hexdigest())
    print("Counts by (degree, sign, exponent): "+", ".join(f"{k}={v}" for k,v in sorted(counts.items())))
    print(f"{len(examples)} three-score examples reproduce the archived 24-decimal enclosures.")


if __name__ == "__main__":
    main()
