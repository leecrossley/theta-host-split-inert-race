"""Compute both ellipse semiaxes directly and by the width-uniform coefficient bound."""
from flint import arb,ctx
from zeros import SCALE,targets

def semiaxes():
    ctx.prec=160
    a,r=arb(2).log(),arb(2).sqrt()
    c=2*(r-1)**2
    sums=[[arb(0) for _ in range(4)] for i in (0,1)]
    for p,q,label,lo,hi in targets():
        t=arb(lo+hi,hi-lo)/(2*SCALE)
        k=2*(a*t/2).sin()/a
        g=(3-2*r*(a*t).cos())/(c*(t*t+arb(1)/4))
        lower=1/(2*(t*t+arb(1)/4))
        i=0 if p in (547,727) else 1
        assert k>0 and g>=lower
        assert (g-(1+4*r*(a*t/2).sin()**2/(r-1)**2)*lower).contains(0)
        for j,z in enumerate((g,g*k,lower,lower*k)):
            sums[i][j]+=z
    return sums


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    sums=semiaxes()
    for host, bounds, readings in zip(("547/727","571/599"),
                                     ((99375,34686,66921,19735),(79731,21866,61027,15894)),sums):
        assert all(v>arb(k)/10000 for k,v in zip(bounds,readings))
        print(f"{host}: dyadic and all-width ellipse semiaxes exceed {bounds} / 10000.")
    assert sums[0][0]+sums[1][0]>arb(179106)/10000
    assert sums[0][1]+sums[1][1]>arb(56552)/10000
    print("At width two the level and move differences exceed 17.9106 and 5.6552 in both directions on every tail.")
