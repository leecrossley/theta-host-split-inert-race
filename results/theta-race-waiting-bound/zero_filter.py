"""Certified constants for the unconditional zero-cancelling positive averages."""
from fractions import Fraction as Q
from math import factorial
from flint import arb, ctx

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")

FIRST = 8590
SPLIT = 200000
BASE = 1000000
G = Q(2513,1000)


def scalar_bounds(pi, log_first, log_split, log_base):
    # Use a single upper/lower tuple for either rational or ball arithmetic.
    pi_lo, pi_hi = pi
    l0, u0 = log_first
    l1, u1 = log_split
    lb, ub = log_base
    # Corollary 1.2 of Bennett, Martin, O'Bryant and Rechnitzer.
    assert 8*(1/pi_lo+(Q(247,1000)+Q(6894,1000)/l0)/FIRST)<3
    # For T >= BASE, pi*sum(1/gamma) <= 1.3*log(12000T)^2 minus the lower-end discount.
    assert 3*pi_hi/4*(Q(1,2)+1/lb)<Q(13,10)
    pre=27*pi_hi*Q(FIRST,FIRST-33)*Q(3,4)*(u1+(u1**2-l0**2)/2)
    assert pre<5406
    for start, upper, target in ((FIRST,u0,Q(302,100)), (SPLIT,u1,Q(333,100))):
        second_moment=Q(3,2)*(upper+1)/start
        loss=(pi_hi*G/2)**2*second_moment
        margin=((Q(633,100)*(1-loss)-Q(62,100)-Q(14,100000)-Q(1,1000))/
                Q(100007,100000)-Q(2353,1000))
        assert margin>target
    assert 221-3*pi_lo*l0**2/4<0
    assert 221+2*pre-3*pi_lo*l1**2/4<10000
    # The derivative in h of 1.3*log(12000*(BASE+4h))^2 decreases from this value.
    assert 8*Q(13,10)*ub/BASE<Q(1,2)
    radius=108+1+Q(4,BASE)+pre+Q(13,10)*ub**2
    assert 1+radius-Q(BASE,4)<-100
    assert 2+2*(108+1+Q(4,BASE))<221


def rational_check():
    # Locate ln(n) with rational arithmetic after range reduction to [1,2].
    def log_interval(n):
        def local(x):
            z=(x-1)/(x+1)
            acc=Q(0)
            power=z
            for j in range(48):
                acc+=2*power/(2*j+1)
                power*=z*z
            return acc, acc+2*power/(97*(1-z*z))
        k=n.bit_length()-1
        lo,hi=local(Q(n,1<<k))
        a,b=local(Q(2))
        return lo+k*a, hi+k*b
    def atan_interval(z):
        n=16
        v=sum(((-1)**j*z**(2*j+1)/(2*j+1) for j in range(n)),Q(0))
        return v,v+z**(2*n+1)/(2*n+1)
    a,b=atan_interval(Q(1,5))
    c,d=atan_interval(Q(1,239))
    pi=(16*a-4*d,16*b-4*c)
    scalar_bounds(pi,*(log_interval(12000*n) for n in (FIRST,SPLIT,BASE)))
    # For the unbounded high tail the per-box attenuation is strictly stronger than exp(-1).
    assert 4*(1-Q(33,BASE))>3
    assert sum((Q(100)**k / factorial(k) for k in range(15)),Q(0))>14*10**10


def ball_check():
    ctx.prec=160
    unit=arb(1)
    def rational_ball(x):
        return arb(x.numerator)/x.denominator
    # Mixed arithmetic in scalar_bounds uses Fractions. Supply rational outer bounds.
    def enclose(v):
        den=10**25
        i=(v*den).floor().unique_fmpz()
        assert i is not None
        return Q(int(i),den), Q(int(i)+1,den)
    pi=arb.pi()
    scalar_bounds(enclose(pi),*(enclose(arb(12000*n).log()) for n in (FIRST,SPLIT,BASE)))
    assert 4*(1-33*unit/BASE)>unit.exp()
    assert 14*10**6*(-arb(100)).exp()<arb(1)/10000
    # Deleting any subset of conjugate-reflection quartets below a cutoff only improves these costs.
    for d in (arb(1)/100, arb(1)/10000):
        for delta in (arb(1)/1000, arb(1)/3):
            g=rational_ball(G)
            response=1-(d*g).sin()**2/(d*delta).cosh()**2
            assert response>1-(d*g)**2 and response<1
    for h,excess in ((2*10**9,Q(1,100000)),(2*10**10,Q(1,1000000))):
        t=BASE+4*h
        logt=arb(12000*t).log()
        assert 10000+arb(13)/5*logt**2 < rational_ball(excess)*h
        assert (arb(104)/5)*logt/t < rational_ball(excess)


if __name__=='__main__':
    rational_check()
    ball_check()
    print("two independent constant checks: finite zero cancellation preserves both 3.02 and 3.33 log 2 moves")
    print("short-box tail below 0.0001; low shifted filter radius below 5406; explicit log-squared deadlines")
    print("the stronger moves finish before H^1.00001 at log H >= 2000000000 and H^1.000001 at log H >= 20000000000")
