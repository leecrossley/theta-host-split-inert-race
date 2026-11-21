"""Count exact crossings by directed integer intervals on every prime cell through 10^6."""
from collections import Counter
from decimal import Decimal, localcontext
from functools import cache
from math import isqrt
from pathlib import Path

from verify_rank import FIELDS, primes_through

GRID = 1 << 116
LIMIT = 8_000_000          # knot index k = 8y
TARGET_NUM, TARGET_DEN = 3947, 2000


class Ball:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        assert lo <= hi
        self.lo, self.hi = lo, hi

    @classmethod
    def ratio(cls, a, b=1):
        assert b > 0
        return cls(a*GRID//b, -((-a*GRID)//b))

    def __add__(self, other):
        return Ball(self.lo+other.lo, self.hi+other.hi)

    def __neg__(self):
        return Ball(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        products = (self.lo*other.lo, self.lo*other.hi,
                    self.hi*other.lo, self.hi*other.hi)
        return Ball(min(products)//GRID, -((-max(products))//GRID))

    def scaled(self, a, b=1):
        assert b > 0
        l, h = sorted((self.lo*a, self.hi*a))
        return Ball(l//b, -((-h)//b))

    def inv(self):
        assert self.lo*self.hi > 0
        return Ball(GRID**2//self.hi, -((-GRID**2)//self.lo))

    def sqrt(self):
        assert self.lo >= 0
        return Ball(isqrt(self.lo*GRID), isqrt(self.hi*GRID)+1)

    def sign(self):
        if self.lo > 0:
            return 1
        if self.hi < 0:
            return -1
        assert self.lo == self.hi == 0
        return 0

    def abs(self):
        return self if self.lo >= 0 else -self if self.hi <= 0 else Ball(0,max(-self.lo,self.hi))

    def floor_at(self, scale):
        a, b = self.lo*scale//GRID, self.hi*scale//GRID
        assert a == b
        return a


@cache
def logarithm(p):
    # Decimal.ln is correctly rounded to nearest at the stated precision.
    # All inputs lie in [2,8*10^6], so one grid unit exceeds its roundoff.
    assert 2 <= p <= LIMIT
    with localcontext() as c:
        c.prec = 64
        a, b = Decimal(p).ln().as_integer_ratio()
    x = Ball.ratio(a,b)
    return Ball(x.lo-1,x.hi+1)


ONE = Ball.ratio(1)
ZERO = Ball.ratio(0)
SQRT2 = Ball.ratio(2).sqrt()
H = SQRT2.inv()
SCALE = (SQRT2.scaled(2)*(SQRT2-ONE)*(SQRT2-ONE)*logarithm(2)).inv().scaled(-8)
DB = (-H, H.scaled(3)+ONE, -H.scaled(2)-ONE.scaled(3), ONE.scaled(2))
DA = (H.scaled(8), -H.scaled(12)-ONE.scaled(4), H.scaled(4)+ONE.scaled(6), -ONE.scaled(2))
TARGET = Ball.ratio(TARGET_NUM,TARGET_DEN)
DEGREE_TARGET = Ball.ratio(2302529,10**6)
DEGREE_RADIUS = Ball.ratio(1,10**10)
HERE = Path(__file__).resolve().parent


def events(primes, field):
    u,v = FIELDS[field]
    squares = [{a*a % q for a in range(1,q)} for q in (u,v)]
    out = {}
    for p in primes:
        if p in (2,u,v):
            continue
        su,sv = (1 if p%q in s else -1 for q,s in zip((u,v),squares))
        if p%16 in (1,15) and su == sv:
            f = 1
        elif p*p%16 in (1,15):
            f = 2
        else:
            f = 4
        norm = p**f
        if norm > LIMIT:
            continue
        mass = SCALE*logarithm(p).scaled(su**f)
        for i in range(4):
            k = (1 << i)*norm
            if k <= LIMIT:
                assert k not in out
                out[k] = (mass*DA[i]).scaled(1,norm), mass*DB[i]
    return out


def score(a,b,y):
    return (a*y+b)*y.sqrt().inv()


def isolating_bin(a,b,kl,kr,target,orientation):
    # Bisection of the monotone physical branch, using an exact dyadic height.
    lo,hi,den = kl,kr,8
    for _ in range(60):
        mid=lo+hi
        lo,hi,den=2*lo,2*hi,2*den
        side=(score(a,b,Ball.ratio(mid,den))-target).sign()
        assert side
        if side == orientation:
            hi=mid
        else:
            lo=mid
    unit=10**9
    assert lo*unit//den == hi*unit//den
    return lo*unit//den


def one_field(primes, field):
    name = "AB"[field]
    changes=events(primes,field)
    keys=sorted({8, LIMIT, *changes})
    assert keys[0]==8
    a=b=ZERO
    rows=[]
    folds=[]
    vertices=[]
    amplitudes=[]
    degree_rows=[]
    prior=0
    for l,r in zip(keys,keys[1:]):
        if l in changes:
            da,db=changes[l]
            a,b=a+da,b+db
        yl,yr=Ball.ratio(l,8),Ball.ratio(r,8)
        vl,vr=score(a,b,yl),score(a,b,yr)
        dl,dr=(a*yl-b).sign(),(a*yr-b).sign()
        vertices.append(((name,l,l),vl.abs()))
        if prior*dl<0:
            amplitudes.append((vl.abs(),2*vl.sign()*dl,(name,l,l),vl.sign()))
        if l==8 or prior*dl<0:
            folds.append(((name,l,l),vl.abs()))
        if dl==0 or dr==0:
            assert l==8 and a.lo==a.hi==b.lo==b.hi==0
        heights=[vl,vr]
        if dl*dr<0:
            at=b*a.inv()
            assert yl.hi<at.lo and at.hi<yr.lo
            value=score(a,b,at)
            folds.append(((name,l,r),value.abs()))
            amplitudes.append((value.abs(),2*value.sign()*dr,(name,l,r),value.sign()))
            heights.insert(1,value)
        if name=="A" and (l,r) in ((62,79),(289,316)):
            # Each early cell contains both signed roots throughout the high-count ranges.
            assert min(vl.hi,vr.hi) < -Ball.ratio(247,100).hi
            assert max(vl.lo,vr.lo) > Ball.ratio(247,100).hi
        for sign in (1,-1):
            target=TARGET.scaled(sign)
            for left,right in zip(heights,heights[1:]):
                t=DEGREE_TARGET.scaled(sign)
                if (left-t).sign()!=(right-t).sign():
                    assert (left-t).abs().lo>DEGREE_RADIUS.hi and (right-t).abs().lo>DEGREE_RADIUS.hi
                    assert dl==dr and dl
                    degree_rows.append((name,sign,l,r,dl,0,(a*b).scaled(4)))
                sl,sr=(left-target).sign(),(right-target).sign()
                assert sl and sr
                if sl!=sr:
                    assert dl==dr==sr
                    lower=isolating_bin(a,b,l,r,target,sr)
                    # The same two coefficients serve both target signs in this cell.
                    rows.append((name,sign,l,r,sr,lower,(a*b).scaled(4)))
        prior=dr
    folds.append(((name,LIMIT,LIMIT),vr.abs()))
    vertices.append(((name,LIMIT,LIMIT),vr.abs()))
    amplitudes.append((vr.abs(),-vr.sign()*dr,(name,LIMIT,LIMIT),vr.sign()))
    return rows,folds,vertices,len(keys)-1,amplitudes,degree_rows


def neighbours(points, expected, grid):
    lo_tag,hi_tag = expected
    lookup=dict(points)
    low,high=lookup[lo_tag],lookup[hi_tag]
    assert low.hi<TARGET.lo and TARGET.hi<high.lo
    for tag,v in points:
        if tag not in expected:
            assert v.hi<low.lo or v.lo>high.hi
    return low.floor_at(grid),high.floor_at(grid)




def amplitude_map(events):
    """Integrate the changes in the four signed counts from the empty upper range."""
    ordered=sorted(events,key=lambda row:row[0].lo)
    assert len(ordered)==1273 and all(x[0].hi<y[0].lo for x,y in zip(ordered,ordered[1:]))
    counts=Counter()
    for _,delta,(h,_,_),s in ordered:
        counts[h,s]-=delta
    assert sum(counts.values())==60
    lines=["host\tsign\tleft_eighths\tright_eighths\tchange\tlower_quadrillionths\ttotal_at\ttotal_above"]
    maximum_spans=[]
    high_levels=[]
    for j,(value,delta,(h,l,r),s) in enumerate(ordered):
        before=sum(counts.values())
        at=before+(delta//2 if abs(delta)==2 else max(0,delta))
        counts[h,s]+=delta
        after=sum(counts.values())
        assert counts[h,s]>=0 and max(at,after)<=76
        endpoint=value.floor_at(10**15)
        lines.append(f"{h}\t{s}\t{l}\t{r}\t{delta}\t{endpoint}\t{at}\t{after}")
        if after==76:
            assert j+1<len(ordered)
            maximum_spans.append(((h,l,r),ordered[j+1][2]))
        if at>=75 or after>=75:
            assert Ball.ratio(23,10).hi < value.lo and value.hi < Ball.ratio(247,100).lo
            high_levels.append(j)
        assert at<76
    assert not any(counts.values()) and len(high_levels)>0
    expected=[(3737924,3737596),(4125553,4125551),(4125886,4126049),(4119601,4119586)]
    assert maximum_spans==[(("B",u,u),("B",v,v))for u,v in expected]
    return "\n".join(lines)+"\n"


def check_discriminants(rows, distinct):
    by_cell={}
    for h,_,l,r,_,_,value in rows:
        assert value.sign()
        by_cell[h,l,r]=value
    disjoint=sorted(by_cell.values(),key=lambda x:x.lo)
    assert len(disjoint)==distinct and all(a.hi<b.lo for a,b in zip(disjoint,disjoint[1:]))


def verify():
    primes=list(primes_through(LIMIT))
    data=[one_field(primes,i) for i in range(2)]
    rows=[r for entry in data for r in entry[0]]
    folds=[r for entry in data for r in entry[1]]
    vertices=[r for entry in data for r in entry[2]]
    event_table=amplitude_map([v for entry in data for v in entry[4]])
    assert event_table==(HERE/'evidence/critical-amplitudes.tsv').read_text()
    degree_rows=[v for entry in data for v in entry[5]]
    assert len(degree_rows)==76
    check_discriminants(degree_rows,74)
    assert [d[3] for d in data]==[131709,130860]
    assert Counter((h,s)for h,s,*_ in rows)=={("A",1):17,("A",-1):18,("B",1):13,("B",-1):14}
    assert len({r[5]for r in rows})==len(rows)==62
    critical=neighbours(folds,(("B",1732,1732),("B",4382524,4382524)),10**12)
    assert critical==(1822860194015,2033070577248)
    chamber=neighbours(vertices+[(k,v)for k,v in folds if k[1]!=k[2]],
                       (("B",1939489,1939489),("A",1622212,1622212)),10**12)
    assert chamber==(1973295623349,1973551727595)
    check_discriminants(rows,60)
    receipt="host\tsign\tleft_eighths\tright_eighths\torientation\tlower_billionths\n"
    receipt+=''.join('\t'.join(map(str,r[:-1]))+'\n' for r in rows)
    assert receipt==(HERE/'evidence/exact-crossings.tsv').read_text()
    print("All 262,569 prime cells through one million give 17 gain and 18 loss crossings on 547/727, and 13 gain and 14 loss crossings on 571/599. All 62 roots have distinct enclosures of width one billionth.")
    print("The full map has 1,273 distinct positive critical or terminal amplitudes. The maximum is 76 crossings on exactly four open intervals. Near 2.302529 they occupy 74 cells with distinct nonzero branch constants.")
    print("The count interval ends at the two B knots 1732/8 and 4382524/8, with magnitudes in (1.822860194015,1.822860194016) and (2.033070577248,2.033070577249).")
    print("The same 60 crossing cells persist between the B knot 1939489/8 and A knot 1622212/8. Their 60 nonzero discriminant constants lie in pairwise disjoint intervals.")


if __name__=="__main__":
    verify()
