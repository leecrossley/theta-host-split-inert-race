"""Enclose the sharp balanced amplitude and the two near-four deadlines."""
from flint import acb, arb, ctx


def kernel_closed(t):
    if t == 0:
        return arb(2).log()
    assert t > 0
    z = (-t).exp()
    m = -(-t).expm1()
    return ((1+z)*(1+z).log()+m*m.log())/(2*z.sqrt())


def kernel_integral(t):
    # Integrate the geometric expansion of the positive gamma series.
    # r = 1-u resolves the removable endpoint at zero displacement.
    if t == 0:
        return acb.integral(lambda r, _: 1/(2-r), 0, 1).real
    assert t > 0
    z = (-t).exp()
    m = -(-t).expm1()
    def integrand(r, analytic):
        return r/((m+z*r)*(1+z-z*r))
    knots = [arb(0)] + [arb(2)**j for j in range(-40,1)]
    total = sum((acb.integral(integrand, lo, hi) for lo,hi in zip(knots,knots[1:])),acb(0))
    assert total.is_finite() and total.imag.contains(0)
    return (-3*t/2).exp()*total.real


def evaluate(kernel):
    a = arb(2).log()
    s = arb(2).sqrt()
    c = 2*(s-1)**2
    conductors = (547,4376,8752,8752,727,5816,11632,11632)
    k = sum((q/arb.pi()).log() for q in conductors)-8*(arb.const_euler()+a)
    assert k > 47
    t0,ta,t2a = (kernel(t) for t in (arb(0),a,2*a))
    v0 = k*4*(a/2).sinh()/(s*c)-8*(2*s*ta-3*t0)/c
    va = -8*(s*(t2a+t0)-3*ta)/c
    def interior(x):
        assert 0 < x < a
        tent = 4*((a-x)/2).sinh()/(s*c)
        return k*tent-8*(s*(kernel(x+a)+kernel(a-x))-3*kernel(x))/c
    high_w = ((4-arb(3)/10**9)/2).log()
    low_w = ((4-arb(11)/10**9)/2).log()
    lower = (interior((a-high_w)/2)-interior((a+high_w)/2))/a
    upper = ((v0-va)*(v0-interior(low_w))).sqrt()/a
    return (v0-va)/a,lower,upper


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 160
    primary = evaluate(kernel_closed)
    ctx.prec = 144
    secondary = evaluate(kernel_integral)
    for x,y in zip(primary,secondary):
        assert x.overlaps(y)
    threshold = arb(269364)/1000
    for endpoint,lower,upper in (primary,secondary):
        assert arb(2693640010193614394)/10**16 < endpoint < arb(2693640010193614395)/10**16
        assert lower > threshold + arb(15)/10**8
        assert upper < threshold - arb(8)/10**8
    print("Closed gamma terms and certified rational-integral quadrature of the gamma kernel, within one explicit formula, enclose the factor-four endpoint between 269.3640010193614394 and 269.3640010193614395.")
    print("At 3.999999997 the symmetric pair beats 269.364 by more than 0.00000015; at 3.999999989 the critical-line upper bound is below 269.364 by more than 0.00000008.")


if __name__ == "__main__":
    main()
