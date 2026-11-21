"""Independent Arb replay from prime-ideal tent prefix sums and exact rational starts."""
from bisect import bisect_right
from collections import Counter
from pathlib import Path
from flint import arb, ctx, fmpz

TOP = 8_000_000
FIELDS = ((547,727),(571,599))
HERE = Path(__file__).resolve().parent


def field_data(rational_primes, host):
    terms=[]
    u,v=FIELDS[host]
    for p in rational_primes:
        if p in (u,v):
            continue
        a,b=(1 if pow(p,(q-1)//2,q)==1 else -1 for q in (u,v))
        f=next(f for f in (1,2,4) if pow(p,f,16) in (1,15) and (a*b)**f==1)
        if p**f<=TOP:
            terms.append((p**f,(-a**f)*arb(p).log()))
    norms=[]
    total,divided=[arb(0)],[arb(0)]
    for n,mass in sorted(terms):
        norms.append(n)
        total.append(total[-1]+mass)
        divided.append(divided[-1]+mass/n)
    knots=sorted({8,TOP,*[k for n in norms for k in (n,2*n,4*n,8*n) if k<=TOP]})
    root=arb(2).sqrt()
    factor=8/(2*root*(root-1)**2*arb(2).log())

    def level(num,den):
        y=arb(num)/den
        i,j,k=(bisect_right(norms,m*num//den)for m in (1,2,4))
        first=2*(total[j]-total[i])-2*y*(divided[j]-divided[i])
        last=4*y*(divided[k]-divided[j])-(total[k]-total[j])
        return factor*(first+last)/y.sqrt()

    def move(num,den):
        return level(2*num,den)-level(num,den)

    return knots,move


def sign(v):
    if v>0:
        return 1
    assert v<0
    return -1


def main():
    ctx.prec=192
    ps=[p for p in range(3,TOP,2) if fmpz(p).is_prime()]
    rows=[l.split() for l in (HERE/'evidence/exact-crossings.tsv').read_text().splitlines()[1:]]
    expected={(h,int(s),int(l),int(r)):(int(ori),int(n)) for h,s,l,r,ori,n in rows}
    found={}
    critical=[]
    events=[]
    sample_cells={}
    sample_count=0
    sample=arb(2302529)/10**6
    knots_all=[]
    degrees={}
    target=arb(3947)/2000
    for i,name in enumerate(("A","B")):
        knots,move=field_data(ps,i)
        left_value=move(knots[0],8)
        prior_sign=0
        for l,r in zip(knots,knots[1:]):
            right_value=move(r,8)
            knots_all.append(((name,l,l),abs(left_value)))
            if l==8:
                assert left_value.contains(0) and right_value.contains(0)
                critical.append(((name,l,l),abs(left_value)))
                left_value=right_value
                continue
            yl,yr=arb(l)/8,arb(r)/8
            gl,gr=yl.sqrt()*left_value,yr.sqrt()*right_value
            a=(gr-gl)/(yr-yl)
            b=gl-a*yl
            dl,dr=sign(a*yl-b),sign(a*yr-b)
            if prior_sign*dl<0:
                critical.append(((name,l,l),abs(left_value)))
                ss=sign(left_value)
                events.append((abs(left_value),name,ss,l,l,2*ss*dl))
            stops=[left_value,right_value]
            if dl!=dr:
                y=b/a
                assert yl<y<yr
                at=(a*y+b)/y.sqrt()
                stops.insert(1,at)
                ss=sign(at)
                events.append((abs(at),name,ss,l,r,2*ss*dr))
                critical.append(((name,l,r),abs(at)))
            for s in (1,-1):
                for w,z in zip(stops,stops[1:]):
                    if sign(w-s*sample)!=sign(z-s*sample):
                        assert abs(w-s*sample)>arb(1)/10**10 and abs(z-s*sample)>arb(1)/10**10
                        sample_cells[name,l,r]=4*a*b
                        sample_count+=1
                    if sign(w-s*target)==sign(z-s*target):
                        continue
                    key=name,s,l,r
                    assert key in expected and key not in found
                    orientation,n=expected[key]
                    assert sign(move(n,10**9)-s*target)==-orientation
                    assert sign(move(n+1,10**9)-s*target)==orientation
                    assert dl==dr==orientation
                    found[key]=True
                    degrees[name,l,r]=4*a*b
            left_value=right_value
            prior_sign=dr
        ss=sign(left_value)
        events.append((abs(left_value),name,ss,TOP,TOP,-ss*prior_sign))
        critical.append(((name,TOP,TOP),abs(left_value)))
        knots_all.append(((name,TOP,TOP),abs(left_value)))
    sample_constants=sorted(sample_cells.values(),key=lambda x:float(x.mid()))
    assert sample_count==76 and len(sample_constants)==74
    assert all(not v.contains(0) for v in sample_constants) and all(x<y for x,y in zip(sample_constants,sample_constants[1:]))
    # The primary table gives exact levels by host and prime-knot tags. Enclose
    # each independently and check the endpoint convention for all 1,273 rows.
    by_tag={(h,sg,l,r):(v,delta) for v,h,sg,l,r,delta in events}
    assert len(by_tag)==len(events)==1273
    reference=(HERE/'evidence/critical-amplitudes.tsv').read_text().splitlines()[1:]
    assert len(reference)==len(events)
    before=-sum(e[-1] for e in events)
    previous=arb(0)
    for line in reference:
        h,sg,l,r,delta,floor,at,after=line.split()
        value,d=by_tag.pop((h,int(sg),int(l),int(r)))
        assert int(delta)==d and previous<value and int(floor)<value*10**15<int(floor)+1
        assert int(at)==before+(d//2 if abs(d)==2 else max(d,0))
        before+=d
        assert before==int(after)
        previous=value
    assert not by_tag and before==0
    assert set(found)==set(expected) and len(found)==62
    values=sorted(degrees.values(), key=lambda x:float(x.mid()))
    assert len(values)==60 and all(not v.contains(0)for v in values)
    assert all(a<b for a,b in zip(values,values[1:]))
    def check_nearest(data,kl,kh):
        points=dict(data)
        lo,hi=points[kl],points[kh]
        assert lo<target<hi
        assert all(k in (kl,kh) or v<lo or v>hi for k,v in data)
    check_nearest(critical,("B",1732,1732),("B",4382524,4382524))
    check_nearest(knots_all+[(k,v)for k,v in critical if k[1]!=k[2]],
                  ("B",1939489,1939489),("A",1622212,1622212))
    print("The independent 192-bit Arb prefix sums and deterministic prime tests verify all 62 disjoint rational root brackets, all 60 distinct branch constants, the 1,273 ordered amplitude changes, the four maximal intervals, and the 74 distinct branch constants at amplitude 2.302529.")


if __name__=="__main__":
    main()
