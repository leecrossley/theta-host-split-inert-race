"""Certify positive modulated averages for every starting height."""
from fractions import Fraction

from flint import arb, ctx

from finite_deadlines import OMITTED, arithmetic, prepare, spectra, uniform, witnesses
from finite_spectrum import brackets, ordinate

CARRIER = (8752, 1093, 251203741273248, 251203741273252)
ORDER = 13
BOX = 54


def square_and_higher_powers():
    a, r = arb(2).log(), arb(2).sqrt()
    B = 2*r*(r-1)**2
    area = 6*r-2*(1+r)*(2*r).sqrt()
    assert 0 < area < 1
    # The two degree-one classes modulo 16; Theorem 1, Table 1 at 10^10.
    eps = arb(9)/1000
    square = (2*area+2*eps*(2+r))/B
    ell = arb(100)
    Z = (ell/2).exp()
    assert Z > 10**10
    square += arb(816)/100/B*(2/Z).sqrt()
    l4 = ell+arb(4).log()
    high = arb(816)/100/B*(sum((l4/k-ell/2).exp() for k in (3,4))+
                            l4/a*(l4/5-ell/2).exp())
    bad = 8*(3*a+sum(arb(p).log() for p in OMITTED[1:]))/(B*Z)
    cost = (square+2*high+2*bad)/a
    assert cost < arb(2353)/1000
    return cost


def transfers(method):
    ctx.prec = 144 if method=='sinc' else 160
    a, r = arb(2).log(), arb(2).sqrt()
    g = ordinate(CARRIER)
    data = brackets()
    assert data.count(CARRIER)==1 and all(
        t==CARRIER or abs(ordinate(t)-g)>0 for t in data)
    assert ORDER*g < 33
    main, leakage = arb(0),arb(0)
    ks = range(-ORDER,ORDER+1)
    weights = {k:arb(ORDER+1-abs(k))/(ORDER+1) for k in ks}

    def box_pair(v):
        if method=='sinc':
            return (BOX*v).sinc()**2
        if v.is_zero():
            return arb(1)
        return (1-(2*BOX*v).cos())/(2*(BOX*v)**2)

    for row in data:
        t=ordinate(row)
        u=a*t/2
        numerator = (3-2*r*(2*u).cos() if method=='sinc' else
                     (r-1)**2+4*r*u.sin()**2)
        coeff=abs(numerator*2*u.sin()/(2*(r-1)**2*a*(t*t+arb(1)/4)))
        for k in ks:
            if row==CARRIER and k==1:
                main=2*weights[k]*coeff
            else:
                leakage+=2*weights[k]*coeff*abs(box_pair(t-k*g))
    norm=1+2*sum(weights[k]*abs(box_pair(k*g)) for k in range(1,ORDER+1))
    tail=50*arb(127)/1000*(ORDER+1)/(BOX*(48-ORDER*g))**2
    assert main>arb(633)/100 and leakage<arb(62)/100 and norm<arb(100007)/100000
    assert tail<arb(14)/100000
    # All large off-line zeros and all trivial zeros together cost below 0.001.
    assert (main-leakage-tail-arb(1)/1000)/norm-square_and_higher_powers()>arb(332)/100
    print(f'{method}: 6,000 positive ordinates, 27 harmonics, and the unlisted tail certify both 3.32 log 2 moves in the shorter kernel.')


def stop_band():
    ctx.prec=144
    c,v=arb(8520),arb(8557)
    gap=(v*v-c*c-arb(1)/4).sqrt()
    assert 8590*11632 < 10**8 and gap>700 and v/(2*gap)<7
    assert 3*c/700 < 37 and c-7>8500
    assert 14*37*10**6*(-arb(41)).exp() < arb(1)/10000
    assert 14*48*(-arb(150)).exp()/(1-(-arb(200)).exp()) < arb(1)/10000
    assert 218+Fraction(600,16999)<219
    assert 250000*(Fraction(1001,1000)-Fraction(17001,16999))>219
    print('With the tilted sphere, every completion is before exp(219) times H^(17001/16999).')


def early_witnesses():
    # This single pair is after every H up to exp(100) and ends before exp(219).
    small = [(10320,1,80),(10348,-1,80),(10688,1,80),(10768,-1,80),
             (10996,1,80),(11184,-1,80),(11416,1,80),(11520,-1,80)]
    for method in ('direct','split'):
        prepare(1000,method)
        radius,rest=uniform(1000)
        found=set()
        for row in small:
            assert row in witnesses(1000)
            m,sign,tier=row
            x=arb(m)/100
            assert x-radius-3*arb(2).log()/2>100 and x+radius<219
            if method=='direct':
                value=sum(c*(t*x).sin() for t,c in spectra[1000,method])
            else:
                left,right=divmod(m,100)
                value=arb(0)
                for t,c in spectra[1000,method]:
                    sa,ca=(t*left).sin_cos()
                    sb,cb=(t*right/100).sin_cos()
                    value+=c*(sa*cb+ca*sb)
            if sign*value-rest-arithmetic(tier)>arb(332)/100:
                found.add(sign)
        assert found=={-1,1}
    print('Both finite implementations find stronger witnesses after exp(100) and before exp(219).')


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    stop_band()
    for method in ('sinc','cosine'):
        transfers(method)
    early_witnesses()
