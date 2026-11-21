"""Two enclosures of the discrete triangular kernel at the dyadic step."""
from pathlib import Path
import sys
from flint import arb,ctx

ROOT=Path(__file__).resolve().parent


def quotient(x,n):
    den=x.sin_pi()
    if den.contains(0):
        return arb(1)
    value=((n*x).sin_pi()/(n*den))**2
    return value.upper().min(arb(1))


def factored(x,n):
    d=n
    while d%2==0:
        d//=2
    poly=arb(1)
    for k in range(1,(d+1)//2):
        poly+=2*(2*k*x).cos_pi()
    value=(poly/d)**2
    step=d
    while step<n:
        value*=(step*x).cos_pi()**2
        step*=2
    return value.upper().min(arb(1))


def check(field,record,table,method):
    n,i,k=record
    assert (field,n) in ((0,112),(1,160)) and 0<=i<len(table) and 0<=k<len(table) and i!=k
    a=arb(2).log()
    r=arb(2).sqrt()
    c=2*(r-1)**2
    series=[]
    for _,_,_,left,right in table:
        t=arb(left+right,right-left)/(2*10**14)
        angle=t*a/(2*arb.pi())
        sn=(a*t/2).sin()
        g=(3-2*r*(t*a).cos())/(c*(t*t+arb(1)/4))
        series.append((angle,4*g*abs(sn)/a))
    f,h=series[i][0],series[k][0]
    gain=(series[i][1]+series[k][1])/2
    error=arb(0)
    mass=arb(1)
    for j in (-1,0,1):
        for l in (-1,0,1):
            u=arb((2-abs(j))*(2-abs(l)))/4
            shift=j*f+l*h
            if j or l:
                mass+=u*method(shift,n)
            for z,(v,amp) in enumerate(series):
                if (z==i and (j,l)==(-1,0)) or (z==k and (j,l)==(0,-1)):
                    continue
                error+=u*amp*method(v+shift,n)
    assert gain-error>arb(723)/200*mass,(field,n,method.__name__)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=128
    table=[tuple(map(int,l.split())) for l in (ROOT/'evidence/dyadic-zeros.tsv').read_text().splitlines()]
    for field,record in enumerate(((112,5193,4381),(160,2,3391))):
        own=[v for v in table if v[0]==field]
        for method in (quotient,factored):
            check(field,record,own,method)
        print(f"Both evaluations certify the exact dyadic lattice in {2*record[0]-1} consecutive starts for field {field}.")


if __name__=="__main__":
    main()
