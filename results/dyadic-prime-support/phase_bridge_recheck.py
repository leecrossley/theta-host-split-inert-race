"""Recheck every exceptional-direction certificate by eliminating each side separately."""
from collections import Counter, defaultdict
from math import isqrt
from pathlib import Path

P = 1000000007
CONDUCTORS = (547,727,571,599)


def project(block):
    # Each row has four ground coefficients and two columns for the common
    # constant and the normalized quadratic term. Reduce the four ground columns.
    m=[list(row) for row in block]
    dimension=0
    for col in range(4):
        i=next((i for i in range(dimension,len(m)) if m[i][col]),None)
        if i is None:
            continue
        m[i],m[dimension]=m[dimension],m[i]
        inv=pow(m[dimension][col],-1,P)
        pivot=[v*inv%P for v in m[dimension]]
        for j in range(dimension+1,len(m)):
            q=m[j][col]
            if q:
                m[j]=[(a-q*b)%P for a,b in zip(m[j],pivot)]
        dimension+=1
    shapes=Counter()
    for a,_,b,_,_,_ in block:
        kind=(int(a!=0),int(b!=0),int((a+b)%P==0))
        shapes[kind]+=1
    upper=min(4,sum(min(2,n) for n in shapes.values()))
    assert dimension==upper
    assert all(not any(row[:4]) for row in m[dimension:])
    return [row[-2:] for row in m[dimension:]]


def run():
    assert all(P%d for d in range(2,isqrt(P)+1))
    path=Path(__file__).resolve().parent/'evidence/phase-proth-witnesses.tsv'
    blocks=defaultdict(list)
    lower=2**24
    now=0
    for line in path.read_text().splitlines()[1:]:
        j,code,offset,_=map(int,line.split())
        while now<j:
            lower=lower*13//10
            now+=1
        power=2**((lower.bit_length()-1)//2+1)
        prime=1+((lower-1)//power+offset)*power
        signs=[pow(prime,(q-1)//2,q) for q in CONDUCTORS]
        assert all(v in (1,q-1) for v,q in zip(signs,CONDUCTORS))
        u,v,s,t=[1 if v==1 else -1 for v in signs]
        assert code==sum(2**k for k,val in enumerate((u,v,s,t)) if val==1)
        a,b=-(u+v)//2,-(s+t)//2
        q=pow(prime,-1,P)
        z=int(u!=v)-int(s!=t)
        blocks[j].append((prime,[a%P,a*q%P,b%P,b*q%P,z%P,z*q*q%P]))
    total=0
    for group in blocks.values():
        rows=[v for _,v in sorted(group)]
        assert len(rows)==12
        cache={}
        def residual(ids):
            if ids not in cache:
                cache[ids]=project([rows[j] for j in ids])
            return cache[ids]
        for cut in range(13):
            for omit in range(12):
                left=tuple(j for j in range(cut) if j!=omit)
                right=tuple(j for j in range(cut,12) if j!=omit)
                terms=residual(left)+residual(right)
                first=next((x for x in terms if any(x)),None)
                assert first is not None
                u,v=first
                assert any((u*b-v*a)%P for a,b in terms)
                total+=1
    assert total==278616
    print("Independent elimination modulo 1,000,000,007 gives two independent residual equations in all 278,616 systems.")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    run()
