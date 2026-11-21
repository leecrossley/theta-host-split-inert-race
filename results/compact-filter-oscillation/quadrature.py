"""Two-sided finite prediction margins from positive quadratures at actual zeros."""
from collections import defaultdict
from pathlib import Path
from flint import arb,ctx

ROOT=Path(__file__).resolve().parent
# Common strict lower bounds, numerator in millionths, for both fields.
THRESHOLDS={4:(1460000,2450000),8:(1400000,1830000),16:(1394000,1190000),
            32:(150000,140000),64:(21000,2000),128:(5900,120),256:(1360,60)}


def read(name):
    return [tuple(map(int,row.split())) for row in (ROOT/'evidence'/name).read_text().splitlines()]


def verify_frames(retention):
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=144
    data=read('brackets.tsv')
    by=defaultdict(list)
    for f,kind,m,index,weight in read('quadratures.tsv'):
        assert 0<=index<len(data) and weight>0 and data[index][0]==f
        by[f,kind,m].append((index,weight))
    assert set(by)=={(f,k,n) for f in (0,1) for k in (0,1) for n in THRESHOLDS}
    a=arb(2).log()
    r=arb(2).sqrt()
    c=2*(r-1)**2
    for (f,kind,n),qs in sorted(by.items()):
        assert len(set(k for k,_ in qs))==len(qs)
        mass=sum(w for _,w in qs)
        moments=[arb(0) for _ in range(n)]
        reciprocal=arb(0)
        for z,w in qs:
            _,_,_,left,right=data[z]
            t=arb(left+right,right-left)/(2*10**14)
            omega=t*a
            g=(3-2*r*omega.cos())/(c*(t*t+arb(1)/4))
            strength=g*(2*abs((omega/2).sin())/a if kind else 1)
            assert strength>0
            reciprocal+=w/strength**2
            for j in range(1,n):
                moments[j]+=w*(j*omega).cos()
        floor=mass-2*sum(abs(b) for b in moments[1:])
        target=arb(THRESHOLDS[n][kind])/10**6
        assert floor*retention**2>target**2*reciprocal,(f,kind,n)


def main():
    verify_frames(arb(1))
    print("Positive quadratures at 710 genuine zero brackets certify all 28 width-two level and move radii by diagonal dominance.")
    for n,(level,move) in THRESHOLDS.items():
        print(f"{n} coordinates: strict two-sided level margin {level}/1000000 and normalized move margin {move}/1000000 in each field.")


if __name__=='__main__':
    main()
