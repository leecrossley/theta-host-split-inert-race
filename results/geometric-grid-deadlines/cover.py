"""Certify a two-sided prime-race filter on every real grid spacing in the cover."""
from pathlib import Path
import numpy as np
from flint import arb,ctx

ROOT=Path(__file__).resolve().parent
Q=1<<21
TURN=1<<24
WEIGHT=1<<25
KERNEL=1<<20
DEN=2*Q*TURN


def modes(rows,field):
    from moments import levels,closed
    pair=((547,727),(571,599))[field]
    a=arb(2).log()
    root=arb(2).sqrt()
    c=2*(root-1)**2
    mass,v1=levels(pair,2,closed)
    moment=(mass-v1)/2
    turns=[]
    lower=[]
    upper=[]
    for i,q,l,lo,hi in rows:
        if i!=field:
            continue
        t=arb(lo+hi,hi-lo)/(2*10**14)
        sn=(a*t/2).sin()
        g=(3-2*root*(a*t).cos())/(c*(t*t+arb(1)/4))
        mass-=2*g
        moment-=2*g*sn**2
        phase=t*a/(2*arb.pi())*TURN
        tick=int((phase.mid()+arb(1)/2).floor().unique_fmpz())
        assert abs(phase-tick)<1
        amp=4*g*abs(sn)/a*WEIGHT
        low=int(amp.lower().floor().unique_fmpz())
        high=int(amp.upper().ceil().unique_fmpz())
        assert low>=0 and low<=amp and amp<=high
        turns.append(tick)
        lower.append(low)
        upper.append(high)
    assert mass>0 and moment>0
    assert 2*(mass*moment).sqrt()/a<arb(163)/100
    assert sum(upper)*KERNEL<2**63
    assert 11*max(turns)*4*Q<2**62
    return np.array(turns,dtype=np.int64),lower,np.array(upper,dtype=np.int64)


def kernel(rate,uncertainty,lo,hi,n):
    # Distance to the nearest integer in units of 1/(2 Q TURN).
    angle=(lo+hi)*rate%DEN
    distance=np.minimum(angle,DEN-angle)-np.abs(rate)*(hi-lo)-4*Q*uncertainty
    scaled=2*n*np.maximum(distance,0)
    bins=np.maximum(1,scaled//(DEN//KERNEL))
    divisor=bins**2
    assert 2*n*(KERNEL//2)<=2**31
    return np.minimum(KERNEL,(KERNEL**3+divisor-1)//divisor)


def verify(cell,arrays):
    n,lo,hi,i,m,k,b=cell
    rates,low,high=arrays
    assert 1<=n<=2048 and Q<=lo<=hi<=2*Q
    assert 0<=i<len(rates) and 1<=m<=6
    assert (k==-1 and b==0) or (0<=k<len(rates) and k!=i and 1<=b<=4)
    second=rates[k] if b else 0
    reward=low[i]*m*(b+1)+(low[k]*b*(m+1) if b else 0)
    side=0
    mass=0
    for j in range(-m,m+1):
        for l in range(-b,b+1):
            weight=(m+1-abs(j))*(b+1-abs(l))
            tick=j*rates[i]+l*second
            error=abs(j)+abs(l)
            mass+=weight*int(kernel(tick,error,lo,hi,n))
            filtered=kernel(rates+tick,error+1,lo,hi,n)
            if (j,l)==(-1,0):
                filtered[i]=0
            if b and (j,l)==(0,-1):
                filtered[k]=0
            side+=weight*int(high@filtered)
    # 1.985 demanded from the full profile, 1.63 for the omitted spectrum.
    assert 200*(reward*KERNEL-side)>723*WEIGHT*mass,(lo,hi,i,k,n)


def certificates(field):
    file=ROOT/'evidence'/f'cells-{field}.tsv'
    return [tuple(map(int,row.split())) for row in file.read_text().splitlines()]


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=144
    rows=[tuple(map(int,row.split())) for row in (ROOT/'evidence/dyadic-zeros.tsv').read_text().splitlines()]
    for field in (0,1):
        arrays=modes(rows,field)
        cells=certificates(field)
        cover=[r for r in cells if r[1]<r[2]]
        assert cover[0][1]==Q and cover[-1][2]==2*Q
        assert all(left[2]==right[1] for left,right in zip(cover,cover[1:]))
        assert {r[0] for r in cover}=={2048}
        points=[r for r in cells if r[1]==r[2]]
        assert len(points)==1 and points[0][1]==Q
        for c in cells:
            verify(c,arrays)
        print(f"Field {field}: {len(cover):,} closed cells cover every log step from log 2 to log 4; the dyadic point uses N={points[0][0]}.")
    print("Integer distance certificates force both signs above 1.985 in the limiting profiles, uniformly in the offset.")


if __name__=="__main__":
    main()
