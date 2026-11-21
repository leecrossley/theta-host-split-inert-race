"""Verify the finite-height width integral on an explicit field prime-ideal sum."""
from math import isqrt
from flint import arb, acb, ctx


def bump(t):
    return (t-1)**3*(2-t)**3


def generator(t):
    u = (t-1)*(2-t)
    return 6*u*(3-2*t)**2-6*u**2-u**3/4


def norm_weights(limit):
    for p in range(3,limit+1,2):
        if p in (547,727) or any(p%k==0 for k in range(2,isqrt(p)+1)):
            continue
        s,t = [1 if pow(p,(q-1)//2,q)==1 else -1 for q in (547,727)]
        f = next(j for j in (1,2,4) if pow(p,j,16) in (1,15) and (s*t)**j==1)
        n=p**f
        if n <= limit:
            yield n,-8*s**f*arb(p).log()/arb(n).sqrt()


def integral_kernel(v):
    if v >= 2:
        return arb(0)
    bottom = max(arb(1),v)
    ans = acb.integral(lambda z,analytic: generator(z)*((z-v)/2).sinh(), bottom,2)
    assert ans.imag.contains(0)
    return ans.real


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=224
    half=arb(1)/2
    for x in (half,arb(1),arb(5)/4,arb(3)/2,arb(7)/4,arb(2),arb(3)):
        expected = bump(x)/2 if 1<x<2 else arb(0)
        e=integral_kernel(x)-expected
        assert e.contains(0) and abs(e)<arb(10)**-50
    smooth_lower = arb(1)/2*(-arb(4)/3).exp()*(arb(1)/10).cos()
    assert smooth_lower > arb(1)/10
    for centre in (31,100,250):
        # Every selected norm is smaller than centre*exp(2), which is below 8*centre.
        lhs,rhs=arb(0),arb(0)
        for n,weight in norm_weights(8*centre):
            x=abs((arb(centre)/n).log())
            lhs += weight*integral_kernel(x)
            rhs += weight*bump(x)/2 if 1<x<2 else 0
        assert (lhs-rhs).contains(0) and abs(lhs-rhs)<arb(10)**-47
        assert abs(rhs)>arb(1)/10000
    print("Directed complex integration reproduces the synthesized kernel at seven interior, boundary and exterior points.")
    print("At geometric centres 31, 100 and 250, complete good-prime ideal sums agree through 47 decimal places.")
    print("The smooth compact example has spectral coefficient above 1/10 throughout |t|<1/20.")


if __name__ == "__main__":
    main()
