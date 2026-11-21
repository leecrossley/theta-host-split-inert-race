"""Exhibit two two-point spectra with equal dyadic moments and unequal signed moves."""
from flint import arb, ctx


def verify():
    with ctx.workprec(192):
        a = arb(2).log()
        p = 4*arb.pi()/a
        scale = arb(2).sqrt()
        def denominator(t):
            return 1/(arb(1)/4+t*t)
        def balance(t):
            return (denominator(t)+denominator(t+5*p)
                    -denominator(t+p)-denominator(t+2*p))
        low,high = arb(1),arb(1000)
        assert balance(low)>0 and balance(high)<0
        for _ in range(140):
            mid = (low+high)/2
            if balance(mid)>0:
                low=mid
            elif balance(mid)<0:
                high=mid
            else:
                raise ArithmeticError("Unresolved sign in the equal-moment bracket")
        t = low.union(high)
        assert 108<t<109
        for x,limit in ((arb(1),arb(7)/10**6),(arb(43380915948)/1000,arb(24)/10**6)):
            def signed(j):
                v=t+j*p
                return (2*(3-2*scale*(a*v).cos())*(a*v/2).sin()*(v*x).sin()
                        *denominator(v)/(a*(scale-1)**2))
            assert signed(0)+signed(5)-signed(1)-signed(2) < -limit
    print("Two simple two-ordinate spectra with the natural weights have exactly identical dyadic moments. Their signed moves differ at log centres 1 and 43,380,915.948.")
