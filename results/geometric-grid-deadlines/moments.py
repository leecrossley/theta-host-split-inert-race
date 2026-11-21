"""Compact-test spectral moments in the two specified fields."""
from fractions import Fraction
from flint import arb


def real(x):
    return arb(x.numerator)/x.denominator


def tent(z):
    s=min(z,1/z)
    if s<=Fraction(1,2):
        return arb(0)
    v=(2*real(s)).sqrt()
    root=arb(2).sqrt()
    return (v-1/v)/(root*(root-1)**2)


def closed(z):
    s=min(z,1/z)
    if s==1:
        return arb(2).log()
    x=real(s)
    return ((1+x)*(1+x).log()+(1-x)*(1-x).log())/(2*x.sqrt())


def series(z):
    s=min(z,1/z)
    if s==1:
        return arb(2).log()
    x=real(s)
    terms=sum((x**(2*k+1)*x.sqrt()/((2*k+1)*(2*k+2)) for k in range(128)),arb(0))
    error=x**257*x.sqrt()/(257*258*(1-x*x))
    return terms+arb(0,error)


def levels(pair,count,arch):
    a=arb(2).log()
    r=arb(2).sqrt()
    c=2*(r-1)**2
    constant=sum((p*f/arb.pi()).log() for p in pair for f in (1,8,16,16))-8*(arb.const_euler()+a)
    out=[]
    for j in range(count):
        z=Fraction(2**j)
        gamma=r*(arch(2*z)+arch(z/2))-3*arch(z)
        arithmetic=arb(0)
        for k in (2,4):
            n=2**k
            arithmetic+=2*a*tent(z/n)/arb(n).sqrt()
        out.append(constant*tent(z)-8*gamma/c-arithmetic)
    return out


def cubic(u):
    return (11+u*(75+u*(16*u-48)))/54
