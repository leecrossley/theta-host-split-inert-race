"""Check the stop band and support constants for the positive Logan window."""
from fractions import Fraction
from flint import arb,ctx


def main():
    ctx.prec = 160
    c,v = arb(8550),arb(8586)
    gap = (v*v-c*c-arb(1)/4).sqrt()
    assert 8590*11632 < 10**8 and gap>780 and v/(2*gap)<6
    assert 3*c/780 < 34 and c-6>8500
    assert 2/(1-(-arb(2)).exp()) < 3
    assert 34*10**6*(-arb(500)).exp() < arb(1)/1000
    assert 4*10**6/arb(2**40) < arb(1)/1000
    assert 161000+Fraction(324000,16999)<161020
    assert Fraction(17001,16999)<Fraction(3001,2999)
    assert arb(3342843)/10**6-arb(2)/100 > arb(332)/100*(1+1/arb(2**40))
    print('The tilted sphere kernel has stop-band multiplier below 34 exp(-8544 L). The verified finite-height range is 8590.')
    print('The improved endpoint is below exp(161020) times H^(17001/16999). The two signed moment tests also clear 3.32 log 2.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    main()
