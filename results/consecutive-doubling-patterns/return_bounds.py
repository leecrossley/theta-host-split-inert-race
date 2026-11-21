"""Verify the tail, dimension, recurrence and word-count budgets."""
from fractions import Fraction
from flint import arb, ctx
from patterns import HOSTS, Host

MAX_CONDUCTOR=11632


def budget(T, bits):
    pi=arb.pi()
    a=arb(2).log()
    t=arb(T)
    q=arb(MAX_CONDUCTOR)
    log_qt=(q*t).log()
    # Corollary 1.2 of Bennett, Martin, O'Bryant and Rechnitzer.
    # The difference between the main terms dominates the error at T
    # and increases for t >= T, uniformly in the conductor.
    assert t/pi*(2*pi*arb(1).exp()).log() > arb(247)/1000*log_qt+arb(6894)/1000
    assert (2*pi*arb(1).exp()).log()/pi > arb(247)/(1000*t)
    r=arb(2).sqrt()
    assert (3+2*r)/(2*(r-1)**2) < 17
    # Eight factors per host; both signs are counted in N.
    tail=8*34*(log_qt+1)/(pi*t)
    # Sixteen factors, but one phase coordinate for a conjugate pair.
    d=8*t*log_qt/pi
    error=2/a*(180*2*pi/(2**bits)+2*tail)
    return bits*d,error


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=144
    exponent,err=budget(250000,18)
    assert exponent < 252000000 and err < arb(6)/100
    exponent2,err2=budget(600,10)
    assert exponent2 < 241000 and err2 < arb(179)/10
    a=arb(2).log()
    for pair,n in zip(HOSTS,(57,56)):
        h=Host(pair,3,1)
        mass=h.profile(Fraction(1))
        assert 2*mass < n*9*a
    print("Common returns: exponent below 252,000,000 and profile error below 0.06.")
    print("Margin-nine words: at most 2^241000 per length on a tail. Monotone runs of 57 and 56 moves are excluded in the two fields.")


if __name__ == "__main__":
    main()
