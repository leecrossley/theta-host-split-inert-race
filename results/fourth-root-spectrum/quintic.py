"""Check the compact spectral weights, their real Fourier integrals and both finite zero tests."""
from fractions import Fraction
from math import comb, factorial
from flint import acb, arb, ctx

from spectrum import arithmetic_cost, zero_tail


def spline(w, order=0):
    """Rational quintic multiplier or its first two derivatives."""
    w = Fraction(w)
    power = 5-order
    return sum((-1)**j*comb(6,j)*max(w+6-2*j,0)**power
               for j in range(7))*Fraction(factorial(5),2112*factorial(power))


def ball(q):
    return arb(q.numerator)/q.denominator


def check_fourier():
    assert [spline(0,j) for j in range(3)] == [1,0,Fraction(-5,11)]
    # The inverse transform is (20/(11 pi)) sinc(v)^6. Each finite
    # integral is enclosed by Arb. The omitted real tails use |sin v| <= 1.
    norm = 20/(11*arb.pi())
    limit = 24
    samples = (Fraction(0),Fraction(1,2),Fraction(2),Fraction(13,4),
               Fraction(4),Fraction(11,2),Fraction(6),Fraction(33,5))
    checked = 0
    for point in samples:
        w = ball(point)
        for j in range(3):
            def integrand(v, analytic):
                return (acb(0,1)*v)**j * v.sinc()**6 * (acb(0,1)*v*w).exp()
            integral = sum((acb.integral(integrand,k,k+1) for k in range(limit)), acb(0))
            finite = 2*norm*integral.real
            tail = 2*norm/((5-j)*arb(limit)**(5-j))
            assert finite.is_finite() and finite.rad() < arb(1)/10**18
            assert abs(finite-ball(spline(point,j))) < tail
            assert tail < arb(3)/10**5
            assert spline(-point,j)==(-1)**j*spline(point,j)
            checked += 1
    print(f"The quintic and its first two derivatives pass {checked} outward real-Fourier checks at interior points, knots, the band edge and beyond it. The normalized variance is exactly 5/11.")


def check_costs():
    x = arb(12494)
    scale = arb(5*10**6)/6
    inv = (-x/4).exp()/scale
    variance = arb(5)/11
    prime = 10**5*(variance/scale**2+(x+3)*variance.sqrt()*inv+(x+3)*(-x/6).exp())
    norm = 20/(11*arb.pi())
    # (|u| lambda)^-6 dominates the outside of [-1/2,1/2].
    outer = 2*norm*(arb(32)/5+2+arb(1)/3)
    third = 2*norm*(arb(1)/4+arb(1)/2)
    assert outer < 12 and third < 1
    # The constant 14 covers the tail of the quadratic and the real
    # exponential Taylor remainder, both after the local truncation.
    assert 12 + (arb(1)/4).exp()/48 < 14
    assert 2*norm*arb(32)/5*inv**5 < arb(1)/1000
    strip = 14*10**6*(-x/4).exp()/scale**3
    normalization = 14*(10**6+1)*(-x/4).exp()/scale**3
    assert prime < arb(7)/10**8
    assert strip+normalization < arb(1)/10**10
    assert prime + strip + normalization + arithmetic_cost() < arb(8)/10000
    b = arb(100).exp()
    projected = 10**6*(b*variance.sqrt()*inv+variance*inv**2+12*inv**5) + 100*zero_tail(b,11632)
    assert projected < arb(1)/10**10
    cut = arb(10)**12
    mass = zero_tail(arb(8590),11632)
    a = (2*11632*cut).log()+1
    for reference,ceiling,units,budget,top in ((1620548,1620550,4,22,405153),
            (1620548,1177000,4,12,405153),(3895692,3895692,2,30,973939)):
        upper = arb(ceiling+1)
        eps = arb(units)/10**8
        p = arb(49)/50*upper/reference
        assert 0 < p < 1 and arb(ceiling)/4+arb(5000000).log() < top
        near = (2000*eps*(eps*upper).exp()+50*((eps*upper).cosh()-1))*mass
        far = (1000*(eps*upper).exp()+50)*160/cut*(a/(1-p)+p/(1-p)**2)
        rest = 50*zero_tail(arb(top).exp(),11632)
        assert near+far+rest+strip+normalization+projected+prime < arb(budget)/1000
    print("The normalized prime correction is below 0.00000007. The strip Taylor and normalization errors together cost less than 10^-10. The strict 0.0008 reconstruction bound holds from log centre 12,494 with cutoff 5,000,000 exp(x/4).")
    print("The first-400H ordinate region now ends at exp(405153). The separate later interval uses exp(973939). All three projected margins survive.")


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 120
    # The weights vanish exactly at and beyond the cutoff.
    assert all(spline(w, j) == 0 for w in (6, Fraction(33, 5)) for j in range(3))
    check_fourier()
    check_costs()


if __name__ == "__main__":
    main()
