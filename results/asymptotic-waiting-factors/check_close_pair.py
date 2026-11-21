"""Certify a close pair of opposite signs by closed gamma terms and a positive series."""
from flint import arb, ctx


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 160
    a = arb(2).log()
    root = arb(2).sqrt()
    c = 2*(root-1)**2
    conductors = (547, 4376, 8752, 8752, 727, 5816, 11632, 11632)
    constant = sum((q/arb.pi()).log() for q in conductors)-8*(arb.const_euler()+a)

    def closed(t):
        z = (-abs(t)).exp()
        return ((1+z)*(1+z).log()+(1-z)*(1-z).log())/(2*z.sqrt())

    def series(t):
        # Every evaluation stays at least 0.34 from a cusp. Bound the remaining
        # positive terms by a geometric progression with the first denominator.
        radius = abs(t)
        assert radius > arb(34)/100
        n = 128
        u = [2*k+arb(3)/2 for k in range(n+1)]
        value = sum(((-v*radius).exp()/(v*v-arb(1)/4) for v in u[:-1]), arb(0))
        bound = (-u[n]*radius).exp()/((u[n]**2-arb(1)/4)*(1-(-2*radius).exp()))
        return value+arb(0).union(bound.upper())

    def even(x, kernel):
        assert 0 < x < a and (x+a).exp() < 3
        tent = 4*((a-x)/2).sinh()/(root*c)
        return constant*tent-8*(root*(kernel(x+a)+kernel(x-a))-3*kernel(x))/c

    s = arb(273)/100000
    u = (even(a/2+s, closed)-even(a/2-s, closed))/a
    v = (even(a/2+s, series)-even(a/2-s, series))/a
    assert u.overlaps(v)
    for value in (u,v):
        assert -arb(1990)/1000 < value < -arb(1989)/1000
    assert arb(1601)/800 < 2*(2*s).exp() < arb(2011)/1000
    print("Closed gamma terms and 128 positive series terms agree: -1.990 < U(0.00273) < -1.989.")
    print("The opposite point follows by oddness and both fit strictly before factor 2.011.")


if __name__ == "__main__":
    main()
