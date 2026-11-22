"""Check two-clock exclusions by exact rational collocation and verify every prime-norm bin independently."""
from bisect import bisect_left
from collections import defaultdict
from csv import DictReader
from fractions import Fraction as F
from itertools import combinations_with_replacement,pairwise
from math import isqrt
from pathlib import Path

HERE=Path(__file__).resolve().parent
HOSTS=((547,727,343057911569,343070636977),(571,599,413718278401,603339156001))
LISTS=({7,31,47,113,127,223,271},{3,7,79,113,127,239,241})


def primality(p):
    return p>=3 and p%2 and all(p%d for d in range(3,isqrt(p)+1,2))


def evaluate(p,h):
    u,v,*aux=HOSTS[h]
    assert p not in (2,u,v,*aux)
    a=1 if pow(p,(u-1)//2,u)==1 else -1
    b=1 if pow(p,(v-1)//2,v)==1 else -1
    f=next(f for f in (1,2,4) if pow(p,f,16) in (1,15) and (a*b)**f==1)
    return p**f,-a**f


def collocation(ms,knots):
    # Rows are free nodal values; the last knot and the full right tail are zero.
    mat=[[F(0)]*len(ms) for _ in knots[:-1]]
    for j,(m,_) in enumerate(ms):
        k=bisect_left(knots,m)
        if k==len(knots):continue
        if knots[k]==m:
            if k<len(mat):mat[k][j]=F(1)
        else:
            l,r=knots[k-1:k+1]
            u=(F(1,m)-1/r)/(1/l-1/r)
            mat[k-1][j]=u
            if k<len(mat):mat[k][j]=1-u
    pivots=[]
    for j in range(len(ms)):
        i=len(pivots)
        pivot=next((k for k in range(i,len(mat)) if mat[k][j]),None)
        if pivot is None:continue
        mat[i],mat[pivot]=mat[pivot],mat[i]
        lead=mat[i][j]
        mat[i]=[x/lead for x in mat[i]]
        for k in range(len(mat)):
            if k!=i and mat[k][j]:
                c=mat[k][j]
                mat[k]=[x-c*y for x,y in zip(mat[k],mat[i])]
        pivots.append(j)
    return {ms[j][1] for i,j in enumerate(pivots) if sum(x!=0 for x in mat[i])==1}


def check_orders():
    for h in range(2):
        data=sorted((evaluate(p,h)[0],p) for p in range(3,289,2) if primality(p) and evaluate(p,h)[0]<289)
        marks={F(1),F(2)}
        for m,p in data+[(289,17)]:
            d=1
            while 2*d<=m:d*=2
            marks.add(F(m,d))
        reps=marks-{F(2)}
        for l,r in pairwise(sorted(marks)):
            reps.update(((3*l+r)/4,(l+3*r)/4))
        union=set()
        for v,w in combinations_with_replacement(sorted(reps),2):
            knots=sorted({t*2**k for k in range(10) for t in (v,w) if t*2**k<=289})
            union.update(collocation(data,knots))
        assert union==LISTS[h]
    print("Independent rational collocation at alternative interior representatives gives the same two seven-prime lists.")


def common_octaves():
    groups=defaultdict(list)
    with (HERE/'evidence/common-octaves.tsv').open() as stream:
        for row in DictReader(stream,delimiter="\t"):
            l,r,p,d=(int(row[key]) for key in ("left","right","prime","degree"))
            assert primality(p) and all(evaluate(p,h)[0]==p**d for h in range(2))
            assert r<p**d<2*l
            groups[l,r].append(p)
    bounds=sorted(groups)
    assert bounds[0][0]==6554 and bounds[-1][1]==2**128
    assert all(l<r and len(set(groups[l,r]))==3 for l,r in bounds)
    assert all(r==v for (l,r),(v,w) in pairwise(bounds))
    for h in range(2):
        masses=sorted(evaluate(p,h)[0] for p in range(3,13109,2) if p not in HOSTS[h]
                      and primality(p) and evaluate(p,h)[0]<=13108)
        assert sum(F(289,2)<m<289 for m in masses)==2
        # For an arbitrary left endpoint, use the three norms after its predecessor.
        for left,third in zip(masses,masses[3:]):
            if left>6554:break
            if F(289,2)<left<=6554:
                assert third<2*left
        i=bisect_left(masses,F(289,2))
        assert masses[i+2]==289
    print("Trial division, Euler symbols and consecutive-norm triples reproduce every real-octave certificate.")


def formulas():
    coeff=defaultdict(list)
    with (HERE/'evidence/minimum-grid-formulas.tsv').open() as stream:
        for row in DictReader(stream,delimiter="\t"):
            key=int(row["host"]),int(row["prime"])
            coeff[key].append((F(row["phase"])*2**int(row["level_exponent"]),F(row["coefficient"])))
    for (h,target),terms in coeff.items():
        bound=int(max(4*z for z,c in terms))
        for p in range(3,bound+1,2):
            if p in HOSTS[h] or not primality(p):continue
            m,s=evaluate(p,h)
            # Hinge representation in the level argument z.
            def weight(z):
                return F(1,m)*(max(0,4*z-m)-6*max(0,z-F(m,2))+2*max(0,z-m))
            v=sum(c*s*weight(z) for z,c in terms)
            assert v==int(p==target)
    assert {p for h,p in coeff if h==0}==LISTS[0] and {p for h,p in coeff if h==1}==LISTS[1]
    print("An independent hinge expansion checks every rational prime coefficient in the fourteen formulas.")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_orders()
    common_octaves()
    formulas()
