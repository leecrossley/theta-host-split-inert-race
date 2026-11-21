"""Replay the cover with endpoint arithmetic, sine tables and a cubic tail bound."""
from functools import lru_cache
from pathlib import Path
import numpy as np
from flint import arb,ctx

ROOT=Path(__file__).resolve().parent
Q=2**21
TURN=2**28
AMP=2**24
COUNT=2**18
UNIT=2**22
DEN=Q*TURN


def floor_ball(x):
    return int(x.lower().floor().unique_fmpz())


def ceil_ball(x):
    return int(x.upper().ceil().unique_fmpz())


def spectral_data(field,table):
    from moments import levels,series,cubic
    pair=((547,727),(571,599))[field]
    v0,v1,v2,v3=levels(pair,4,series)
    a=arb(2).log()
    r=arb(2).sqrt()
    m1=(v0-v1)/2
    m2=(3*v0-4*v1+v2)/8
    m3=(10*v0-15*v1+6*v2-v3)/32
    residual=2*(11*v0+75*m1-48*m2+16*m3)/(54*a)
    turns=[]
    amps=[]
    for host,_,_,left,right in table:
        if host!=field:
            continue
        t=arb(left+right,right-left)/(2*10**14)
        sine=(a*t/2).sin()
        g=(1+4*r*(sine/(r-1))**2)/(2*(t*t+arb(1)/4))
        residual-=4*g*cubic(sine*sine)/a
        wt=4*g*abs(sine)/a
        turns.append((floor_ball(arb(left)*a*TURN/(2*arb.pi()*10**14)),
                      ceil_ball(arb(right)*a*TURN/(2*arb.pi()*10**14))))
        amps.append((floor_ball(wt*AMP),ceil_ball(wt*AMP)))
    assert residual>0 and residual<arb(163)/100
    lt,ut=np.array(turns,dtype=np.int64).T
    lw,uw=np.array(amps,dtype=np.int64).T
    assert np.all(ut>=lt) and np.all(lt>=0) and np.all(uw>=lw) and np.all(lw>=0)
    assert int(max(ut))*11*2*Q < 2**62 and int(uw.sum())*UNIT<2**63
    return lt,ut,lw,uw


@lru_cache(maxsize=None)
def sine_envelopes(n):
    table=[UNIT]
    for j in range(1,COUNT//2+1):
        sn=(arb(j)/COUNT).sin_pi()
        table.append(min(UNIT,ceil_ball(UNIT/(n*sn)**2)))
    return np.array(table,dtype=np.int64)


def geometric_kernel(low,high,v,w,values):
    # The corners attaining the product endpoints depend on the sign.
    left=np.minimum(low*v,low*w)
    right=np.maximum(high*v,high*w)
    crosses=left//DEN!=right//DEN
    distance=np.where(crosses,0,np.minimum(left%DEN,(-right)%DEN))
    bins=distance//(DEN//COUNT)
    assert np.all((0<=bins)&(bins<=COUNT//2))
    return values[bins]


def evaluate(row,data):
    n,v,w,f,m,g,b=row
    lt,ut,lw,uw=data
    values=sine_envelopes(n)
    den=(m+1)*(b+1)
    main=int(lw[f])*m*(b+1)
    if b:
        main+=int(lw[g])*b*(m+1)
    penalty=0
    normal=0
    for j in range(-m,m+1):
        for k in range(-b,b+1):
            coeff=(m+1-abs(j))*(b+1-abs(k))
            jl,ju=sorted((int(lt[f])*j,int(ut[f])*j))
            kl,ku=sorted((int(lt[g])*k,int(ut[g])*k))
            low=jl+kl
            high=ju+ku
            normal+=coeff*int(geometric_kernel(low,high,v,w,values))
            envelope=geometric_kernel(lt+low,ut+high,v,w,values)
            if (j,k)==(-1,0):
                envelope[f]=0
            if b and (j,k)==(0,-1):
                envelope[g]=0
            penalty+=coeff*int(uw@envelope)
    assert den>0 and normal>=den*UNIT
    assert 200*(main*UNIT-penalty)>723*AMP*normal,(n,v,w)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=144
    table=[tuple(map(int,s.split())) for s in (ROOT/'evidence/dyadic-zeros.tsv').read_text().splitlines()]
    total=0
    for f in (0,1):
        data=spectral_data(f,table)
        rows=[tuple(map(int,s.split())) for s in (ROOT/'evidence'/f'cells-{f}.tsv').read_text().splitlines()]
        intervals=sorted((v,w) for _,v,w,*_ in rows if v!=w)
        assert intervals[0][0]==Q and intervals[-1][1]==2*Q
        assert all(x[0]<x[1] for x in intervals)
        assert all(x[1]==y[0] for x,y in zip(intervals,intervals[1:]))
        for row in rows:
            evaluate(row,data)
        total+=len(rows)
    print(f"Sine envelopes from directed arithmetic, endpoint products and cubic moments independently verify all {total:,} filter certificates.")


if __name__=="__main__":
    main()
