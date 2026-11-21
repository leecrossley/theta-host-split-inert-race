"""Check the gamma integral anew with an exponential series and a uniform tail bound."""
from flint import arb, ctx

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")


def main():
    ctx.prec=112
    a=arb(2).log()
    b=arb(2).sqrt()
    c=2*(b-1)**2
    B=b*c
    conductors=(547,4376,8752,8752,727,5816,11632,11632)
    constant=sum((q/arb.pi()).log() for q in conductors)+8*(arb(3)/4).digamma()
    powers=((4,2,-2),(16,2,-2),(31,31,8),(47,47,-8),(49,7,-8),(64,2,-2))
    prime_terms=[(arb(n).log(),m*arb(p).log()/arb(n).sqrt()) for n,p,m in powers]

    def g(s):
        # Direct weight in centred coordinates, independent of the hyperbolic-sine implementation.
        v = s.exp()/2
        return b/B*(-s/2).exp()*(2*(1-v)).min(4*v-1).max(0)

    def j(s):
        return (g(s+a/2)-g(s-a/2))/a

    count=256
    us=[arb(3)/2+2*k for k in range(count)]
    denoms=[u*u-arb(1)/4 for u in us]
    start=arb(3)/2+2*count
    tail=8*28*(1/start**2+1/(2*start))
    assert tail<arb(22)/100
    for k in range(1014,3214):
        x=arb(2*k+1,1)/2000
        jx=j(x)
        prime=sum(coeff*j(x-logn) for logn,coeff in prime_terms)
        # The four knots of the odd kernel's second derivative.
        distances=[abs(x+3*a/2),abs(x+a/2),abs(x-a/2),abs(x-3*a/2)]
        integral=arb(0)
        for u,den in zip(us,denoms):
            # D(x+a/2)-D(x-a/2) = b(e_{+3a/2}-e_{-3a/2})
            #                         -(b+3)(e_{+a/2}-e_{-a/2}).
            e0,e1,e2,e3=[(-u*v).exp() for v in distances]
            jump=b*(e0-e3)-(b+3)*(e1-e2)
            integral += (jx/(2*u)+jump/(a*c))/den
        lower=prime+constant*jx-8*integral-tail
        assert lower > -arb(9)/5, (k,lower)
    print("independent exponential series: 2200 cells, 256 terms, uniform tail < 0.22, profile > -9/5")


if __name__=="__main__":
    main()
