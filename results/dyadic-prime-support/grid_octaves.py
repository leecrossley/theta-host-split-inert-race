"""Build common prime-norm witnesses covering real octaves up to the analytic cutoff."""
from math import isqrt
from pathlib import Path
import hashlib
from support_recheck import BAD,characters,prime_data,small_primes

END=2**128
START=6554


def root(n,d):
    return isqrt(isqrt(n)) if d==4 else isqrt(n) if d==2 else n


def cover():
    divisors=small_primes(1<<17)
    fields=[characters(h) for h in (0,1)]
    omitted=set(BAD[0])|set(BAD[1])
    def good_prime(p):
        assert p < (1<<17)**2
        if p%2==0 or p in omitted:
            return False
        for d in divisors:
            if d*d>p:break
            if p%d==0:return False
        return True
    left=START
    records=[]
    while left<END:
        right=min(5*left//4,END)
        found=set()
        for d in (4,2,1):
            first=root(right,d)+1
            if first >= (1<<17)**2:
                continue
            last=min(root(2*left-1,d),first+10000,(1<<17)**2-1)
            for p in range(first+(first%2==0),last+1,2):
                if p in found or not good_prime(p):
                    continue
                m=p**d
                if not all(prime_data(p,h,fields[h])[0]==m for h in (0,1)):
                    continue
                assert right<m<2*left
                records.append((left,right,p,d))
                found.add(p)
                if len(found)==3:break
            if len(found)==3:break
        assert len(found)==3
        left=right
    report='left\tright\tprime\tdegree\n'+''.join('\t'.join(map(str,v))+'\n' for v in records)
    return records,report


def main():
    rows,report=cover()
    actual=Path(__file__).resolve().parent/'evidence/common-octaves.tsv'
    assert actual.read_text()==report
    print(f"{len(rows)//3} closed start intervals from 6554 to 2^128 are covered by {len(rows)} common prime-norm witnesses; largest prime {max(r[2] for r in rows)}.")
    print("Common-octave table SHA-256 "+hashlib.sha256(report.encode()).hexdigest()+".")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
