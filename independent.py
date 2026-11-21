"""Recompute the variance witnesses and the finite-width synthesis by different formulas."""
from fractions import Fraction
from flint import arb,ctx
from shell_energy import brackets,COUNTS,FLOORS,SHELLS,SCALE


def rational_check():
    # On [1,2] use the polynomial in ascending powers of x.
    p=[Fraction(1)]
    for factor in ([-1,1],[2,-1])*3:
        q=[Fraction(0)]*(len(p)+1)
        for i,a in enumerate(p):
            for j,b in enumerate(factor):
                q[i+j]+=a*b
        p=q
    d2=[(i+1)*(i+2)*p[i+2] if i+2 < len(p) else 0 for i in range(len(p))]
    gen=[a-b/4 for a,b in zip(d2,p)]
    # Two polynomial solutions to P'+uP=generator give exact exponential antiderivatives.
    def antiderivative(exponent):
        out=[Fraction(0)]*(len(p)+1)
        for i in reversed(range(len(p))):
            out[i]=(gen[i]-(i+1)*out[i+1])/exponent
        return out[:-1]
    high,low=antiderivative(Fraction(1,2)),antiderivative(Fraction(-1,2))
    def at(q,x):
        return sum(a*x**i for i,a in enumerate(q))
    for x in (1,2):
        assert at(high,x)==at(low,x)==0
    assert all(b-a==v for a,b,v in zip(high,low,p))


def main():
    ctx.prec=192
    totals=[arb(0) for _ in SHELLS]
    nrows=[0 for _ in SHELLS]
    for _,_,lo,hi in brackets():
        for i,T in enumerate(SHELLS):
            if SCALE*T < lo < hi < 2*SCALE*T:
                mid=arb(lo+hi)/(2*SCALE)
                t=mid+arb(0,1)*arb(hi-lo)/(2*SCALE)
                k=(t*arb(2).log()).cos()
                f=arb(12)/5*(2*k*k-1)-arb(8).sqrt()*k
                chord=4*(t/(2*T)).sin()**2
                totals[i]+=2*(f*chord/(t**2+arb(1)/4))**2
                nrows[i]+=1
    assert tuple(nrows)==COUNTS
    for i,T in enumerate(SHELLS):
        assert (totals[i]*T**3/arb(2*T).log()).floor().unique_fmpz()==FLOORS[i]
    assert Fraction(8,2)+Fraction(144,50)==Fraction(172,25)
    rational_check()
    print("Independent half-angle and double-angle evaluations reproduce all six shell certificates.")
    print("The constant-free contrast has exact mean numerator square 172/25.")
    print("Exact rational polynomial antiderivatives give the same compact kernel on the full interval.")


if __name__ == "__main__":
    main()
