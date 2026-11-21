"""Enclose the signed Fourier sums at integer dyadic exponents."""
from fractions import Fraction
from pathlib import Path
import subprocess
import sys
from flint import arb,ctx
from patterns import Host,HOSTS

ROOT=Path(__file__).resolve().parent
TICKS=10**14
AMP=1<<24
PHASE=1<<64
SINE_SIZE=1<<15


def integer(x):
    return int((x.mid()+arb(1)/2).floor().unique_fmpz())


def budgets(i,roots):
    h=Host(HOSTS[i],18,1)
    mass=h.profile(Fraction(1))
    moment=(mass-h.profile(Fraction(2)))/2
    plain=arb(0)
    phase_slope=arb(0)
    weight_error=arb(0)
    modes=[]
    for q,l,lo,hi in roots:
        t=arb(lo+hi,hi-lo)/(2*TICKS)
        sn=(h.a*t/2).sin()
        g=(3-2*h.root*(h.a*t).cos())/(h.c*(arb(1)/4+t*t))
        mass-=2*g
        moment-=2*g*sn**2
        weight=-4*g*sn/h.a
        theta=t*h.a/(4*arb.pi())
        tick=integer(theta*PHASE)
        amp=integer(weight*AMP)
        rounded=arb(amp)/AMP
        phase_slope+=abs(rounded)*abs(theta-arb(tick)/PHASE)*2*arb.pi()
        weight_error+=abs(weight-rounded)
        plain+=abs(rounded)
        modes.append((tick%PHASE,amp))
    assert mass>0 and moment>0
    tail=2*(mass*moment).sqrt()/h.a
    base=tail+weight_error+plain*(arb.pi()/SINE_SIZE+arb(1)/AMP)
    return modes,base,phase_slope


def read(name):
    return [tuple(map(int,line.split())) for line in (ROOT/'evidence'/name).read_text().splitlines()]


def evaluate(binary,host,roots,samples,table):
    modes,base,slope=budgets(host,roots)
    data=table+str(len(modes))+'\n'
    data+=''.join(f'{tick} {amp}\n' for tick,amp in modes)
    data+=''.join(f'{n} {code} {k}\n' for n,code,k in samples)
    result=subprocess.run([binary],input=data,text=True,check=True,stdout=subprocess.PIPE)
    out=[tuple(map(int,x.split())) for x in result.stdout.splitlines()]
    assert [x[:3] for x in out]==samples
    return [(row[:3],row[3:],base+slope*(2*(row[2]+row[0]-1)+1)) for row in out]


def positive_sample(case):
    n,code,k=case
    if k>=0:
        return case
    assert k<=-n
    # Evenness gives P(-x-a) = -P(x). Reverse the four-letter digits
    # and both signs to evaluate a negative start at a positive exponent.
    reverse=0
    for _ in range(n):
        reverse=4*reverse+3-code%4
        code//=4
    return n,reverse,-k-n


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=144
    roots=read('dyadic-zeros.tsv')
    samples=read('dyadic-starts.tsv')
    assert max(abs(k) for _,_,k in samples)==2143656793
    table_values=[]
    for j in range(SINE_SIZE):
        z=(2*arb.pi()*j/SINE_SIZE).sin()*AMP
        v=integer(z)
        assert abs(z-v)<1
        table_values.append(v)
    assert max(map(abs,table_values))<=AMP
    table=' '.join(map(str,table_values))+'\n'
    thresholds=dict(read('dyadic-margins.tsv'))
    assert len(thresholds)==len(read('dyadic-margins.tsv'))
    assert sorted((n,code) for n,code,_ in samples)==[(n,c) for n in sorted(thresholds) for c in range(4**n)]
    for i in (0,1):
        vals=evaluate(sys.argv[1],i,[x[1:] for x in roots if x[0]==i],
                      [positive_sample(x) for x in samples],table)
        for (n,code,k),row,error in vals:
            assert len(row)==n
            for j,v in enumerate(row):
                sign=1 if code & (1<<(2*(n-1-j)+1-i)) else -1
                assert sign*arb(v)/AMP**2-error > arb(thresholds[n])/10000,(n,code,k,i,j,float(sign*arb(v)/AMP**2-error))
    for n,margin in sorted(thresholds.items()):
        print(f"Length {n}: all {4**n:,} two-field words at dyadic starts have profile margin above {margin}/10000.")


if __name__=="__main__":
    main()
