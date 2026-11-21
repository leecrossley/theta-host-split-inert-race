"""Replay every dyadic word with a shifted cosine grid and a cubic tail bound."""
from fractions import Fraction
from pathlib import Path
import numpy as np
from flint import arb,ctx
from independent import gamma_series,tent

ROOT=Path(__file__).resolve().parent
SIZE=1<<16
SCALE=1<<23
TICKS=10**14
PHASE=1<<64


def floor_mid(x):
    return int(x.mid().floor().unique_fmpz())


def majorant(u):
    return (11+u*(75+u*(16*u-48)))/54


def levels(pair):
    a=arb(2).log()
    r=arb(2).sqrt()
    c=2*(r-1)**2
    constant=sum((p*f/arb.pi()).log() for p in pair for f in (1,8,16,16))-8*(arb.const_euler()+a)
    values=[]
    for n in range(4):
        z=Fraction(2**n)
        arch=constant*tent(z)-8/c*(r*(gamma_series(2*z)+gamma_series(z/2))-3*gamma_series(z))
        # Up to norm 16 both fields have zero odd prime traces. The two
        # conductor characters cancel at 2 to odd order and add at even order.
        arithmetic=arb(0)
        for j in (2,4):
            m=2**j
            arithmetic+=2*a/arb(m).sqrt()*tent(z/m)
        values.append(arch-arithmetic)
    return values


def host_data(i,rows):
    pair=((547,727),(571,599))[i]
    a=arb(2).log()
    r=arb(2).sqrt()
    c=2*(r-1)**2
    v0,v1,v2,v3=levels(pair)
    moments=((v0-v1)/2,(3*v0-4*v1+v2)/8,(10*v0-15*v1+6*v2-v3)/32)
    allowance=2*(11*v0+75*moments[0]-48*moments[1]+16*moments[2])/(54*a)
    angles=[]
    weights=[]
    slope=arb(0)
    for host,_,_,lo,hi in rows:
        if host != i:
            continue
        t=arb(lo+hi,hi-lo)/(2*TICKS)
        sn=(a*t/2).sin()
        g=(1+4*r*sn**2/(r-1)**2)/(2*(t*t+arb(1)/4))
        allowance-=4*g*majorant(sn**2)/a
        w=-4*g*sn/a
        step=t*a/(4*arb.pi())
        tick=floor_mid(step*PHASE)
        weight=floor_mid(w*SCALE)
        angles.append(tick%PHASE)
        weights.append(weight)
        rounded=arb(weight)/SCALE
        allowance+=abs(w-rounded)+abs(rounded)*(arb.pi()/SIZE+arb(1)/SCALE)
        slope+=abs(rounded)*2*arb.pi()*abs(step-arb(tick)/PHASE)
    assert allowance>0
    assert sum(map(abs,weights))*SCALE < 2**63
    return np.array(angles,dtype=np.uint64),np.array(weights,dtype=np.int64),allowance,slope


def read(name):
    lines=(ROOT/'evidence'/name).read_text().splitlines()
    return [tuple(int(n) for n in line.split()) for line in lines]


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=144
    values=[]
    for j in range(SIZE):
        z=((2*j+1)*arb.pi()/SIZE).cos()*SCALE
        v=floor_mid(z)
        # Exact values near an integer can round down to its predecessor.
        assert abs(z-v)<1
        values.append(v)
    table=np.asarray(values,dtype=np.int64)
    quarter=np.uint64(1<<62)
    shift=np.uint64(48)
    rows=read('dyadic-zeros.tsv')
    cases=read('dyadic-starts.tsv')
    thresholds=dict(read('dyadic-margins.tsv'))
    assert sorted((n,w) for n,w,_ in cases)==[(n,c) for n in sorted(thresholds) for c in range(4**n)]
    for i in (0,1):
        t,w,base,slope=host_data(i,rows)
        for n,c,k in cases:
            # sin(x)=cos(x-pi/2). The selected cosine is at the centre
            # of the phase cell. Unsigned multiplication is modulo 2^64.
            err=base+slope*max(abs(2*k+1),abs(2*(k+n-1)+1))
            for j in range(n):
                m=2*(k+j)+1
                indices=((t*np.uint64(m%PHASE)-quarter)>>shift).astype(np.intp)
                total=int(np.dot(w,table[indices]))
                sign=1 if ((c>>((n-j)*2-1-i))%2) else -1
                assert sign*arb(total)/SCALE**2-err>arb(thresholds[n])/10000,(n,c,k,i,j)
    print(f"Shifted cosine cells and the cubic moment majorant independently certify every move in {len(cases):,} dyadic witnesses.")


if __name__=='__main__':
    main()
