"""Use a positive trigonometric symbol to verify every real coefficient vector."""
from collections import defaultdict
from pathlib import Path
from flint import arb,acb,ctx

ROOT=Path(__file__).resolve().parent


def rows(name):
    return [tuple(map(int,line.split())) for line in (ROOT/'evidence'/name).read_text().splitlines()]


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=160
    codes=rows('brackets.tsv')
    table=defaultdict(list)
    for f,kind,n,i,weight in rows('quadratures.tsv'):
        assert 0<=i<len(codes) and codes[i][0]==f and weight>0
        table[f,kind,n].append((codes[i],weight))
    target={4:('1.46','2.45'),8:('1.4','1.83'),16:('1.394','1.19'),
            32:('.15','.14'),64:('.021','.002'),128:('.0059','.00012'),256:('.00136','.00006')}
    assert set(table)=={(f,k,n) for f in (0,1) for k in (0,1) for n in target}
    a=arb(2).log()
    root=arb(2).sqrt()
    for (field,kind,n),entries in table.items():
        moments=[arb(0) for _ in range(n)]
        cost=arb(0)
        for (_,_,_,lo,hi),mass in entries:
            angle=a*arb(lo+hi)/(2*10**14)
            radius=a*(hi-lo)/(2*10**14)
            u=acb(0,angle).exp()
            v=acb(1)
            for j in range(1,n):
                v*=u
                moments[j]+=mass*(v.real+arb(0,j*radius))
            t=arb(lo+hi,hi-lo)/(2*10**14)
            sin=(t*a/2).sin()
            weight=(1+4*root*(sin/(root-1))**2)/(2*(t*t+arb(1)/4))
            if kind:
                weight*=2*abs(sin)/a
            assert weight>0
            cost+=mass/weight**2
            moments[0]+=mass
        derivative=2*sum(j*abs(moments[j]) for j in range(1,n))
        # S(theta) is even. Cover [0,pi] by 32 intervals, testing
        # midpoint values and a derivative bound on each closed interval.
        g=32
        cutoff=arb(target[n][kind])**2*cost
        for k in range(g):
            theta=arb.pi()*(2*k+1)/(2*g)
            symbol=moments[0]+2*sum(moments[j]*(j*theta).cos() for j in range(1,n))
            assert symbol-derivative*arb.pi()/(2*g)>cutoff,(field,kind,n,k)
    print("The trigonometric symbols stay above their required thresholds on 896 closed angular cells by an independent midpoint-and-derivative calculation.")


if __name__=='__main__':
    main()
