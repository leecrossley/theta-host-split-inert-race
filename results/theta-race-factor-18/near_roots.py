"""Independent full-precision prime sums in the four boundary neighborhoods."""
from bisect import bisect_left
from itertools import compress
from math import isqrt
from flint import arb, fmpz

NEIGHBORHOODS = ((6288627,128),(20306533,128),(20077386,128),(61070400,128),
                 (562481520,4096),(1701283906,4096))
OMITTED = {2, 547, 727, 343057911569, 343070636977}


def legendre(p, q):
    return 1 if pow(p, (q-1)//2, q)==1 else -1


def higher_degrees(limit):
    out = {}
    for p in range(3, isqrt(limit)+1, 2):
        if p in OMITTED or not fmpz(p).is_prime():
            continue
        u, v = legendre(p,547), legendre(p,727)
        f = next(f for f in (1,2,4) if pow(p,f,16) in (1,15) and (u*v)**f==1)
        if f>1 and p**f<=limit:
            out[p**f]=p
    return out


def knots(higher):
    groups=[]
    for m,width in NEIGHBORHOODS:
        a,b=m-width,m+width
        keys={8*a,8*b}
        for j in (1,2,4,8):
            for n in range(j*a,j*b+1):
                degree_one = (n % 16 in (1,15) and n not in OMITTED and
                              fmpz(n).is_prime() and legendre(n,547)==legendre(n,727))
                if degree_one or n in higher:
                    keys.add(8*n//j)
        groups.append(sorted(keys))
    return groups


def signed_prefix(edges,higher):
    """Sieve each degree-one residue progression and bin direct Arb logs by norm."""
    limit=edges[-1]
    small=[n for n in range(3,isqrt(limit)+1,2) if fmpz(n).is_prime()]
    strikes={r:[(p,(-r*pow(16,-1,p))%p) for p in small] for r in (1,15)}
    sq1={i*i%547 for i in range(1,547)}
    sq2={i*i%727 for i in range(1,727)}
    ordinary=[arb(0) for _ in edges]
    reciprocal=[arb(0) for _ in edges]
    def add(n,p,sign):
        j=bisect_left(edges,n)
        z=arb(p).log()
        if sign<0:
            z=-z
        ordinary[j]+=z
        reciprocal[j]+=z/n
    for n,p in higher.items():
        add(n,p,-1)
    size=1<<18
    for r in (1,15):
        cap=(limit-r)//16+1
        for first in range(0,cap,size):
            last=min(cap,first+size)
            flag=bytearray([1])*(last-first)
            for p,rem in strikes[r]:
                if p*p>16*(last-1)+r:
                    break
                k=max(first,(p*p-r+15)//16)
                k+=(rem-k)%p
                if k<last:
                    flag[k-first::p]=bytes((last-1-k)//p+1)
            for k in compress(range(first,last),flag):
                p=16*k+r
                if p<2 or p in OMITTED:
                    continue
                u=p%547 in sq1
                if u==(p%727 in sq2):
                    add(p,p,-1 if u else 1)
    prefix={}
    x=y=arb(0)
    for j,edge in enumerate(edges):
        x+=ordinary[j];y+=reciprocal[j]
        prefix[edge]=(x,y)
    return prefix


def near_values(keys,prefix):
    root=arb(2).sqrt()
    scale=8/(2*root*(root-1)**2*arb(2).log())
    def level(k):
        y=arb(k)/8
        i,j,l=(prefix[(k*m)//8] for m in (1,2,4))
        return 2*(j[0]-i[0])-2*y*(j[1]-i[1])+4*y*(l[1]-j[1])-(l[0]-j[0])
    return {k:scale*(level(2*k)/root-level(k))/(arb(k)/8).sqrt() for k in keys}


def crossing(keys,values,sign,left):
    """Enclose the last exit on the left or first entry on the right of a strict excursion."""
    target=arb(19735)/10000
    point=[sign*values[k] for k in keys]
    assert all(z>target or z<target for z in point)
    assert (point[0]>target)==left and (point[-1]>target)!=left
    ix=max(i for i,z in enumerate(point) if z>target) if left else min(i for i,z in enumerate(point) if z>target)-1
    u,v=keys[ix:ix+2]
    start, end = arb(u)/8, arb(v)/8
    gl, gh = point[ix]*start.sqrt(), point[ix+1]*end.sqrt()
    slope = (gh-gl)/(end-start)
    derivative = slope - target/(2*start.union(end).sqrt())
    assert derivative<0 if left else derivative>0
    l, h = start, end
    for _ in range(72):
        mid=(l+h)/2
        weight=(mid-start)/(end-start)
        z=(1-weight)*gl+weight*gh-target*mid.sqrt()
        if z>0 or z<0:
            if (z>0)==left:l=mid
            else:h=mid
        else:
            # Mean-value bound if the current enclosure straddles the target.
            delta=(abs(z)/abs(derivative)).upper()
            bracket=(mid-delta).union(mid+delta)
            return u, v, bracket.intersection(l.union(h))
    return u, v, l.union(h)


def exact_roots():
    limit=8*sum(NEIGHBORHOODS[-1])
    higher=higher_degrees(limit)
    groups=knots(higher)
    allkeys=sorted({k for ks in groups for k in ks})
    edges=sorted({(k*m)//8 for k in allkeys for m in (1,2,4,8)})
    prefix=signed_prefix(edges,higher)
    scores=near_values(allkeys,prefix)
    assert max(z.rad() for z in scores.values())<arb("1e-27")
    roots=all_crossings(groups,scores)
    return allkeys,scores,roots


def all_crossings(groups, scores):
    return [crossing(g,scores,s,left) for g,s,left in
            zip(groups,(1,1,-1,-1,-1,-1),(True,False)*3,strict=True)]
