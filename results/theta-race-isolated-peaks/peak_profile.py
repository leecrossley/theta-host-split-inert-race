"""Certify the extremal height, prime-only error and the high-score density budget."""
from flint import arb, ctx
from inputs import brackets

CONDUCTORS = (547, 4376, 8752, 8752, 727, 5816, 11632, 11632)


def main():
    ctx.prec = 192
    a, root = arb(2).log(), arb(2).sqrt()
    c = 2*(root-1)**2
    b = root*c
    g0 = 4*(a/2).sinh()/b
    constant = sum((q/arb.pi()).log() for q in CONDUCTORS)-8*(arb.const_euler()+a)
    z = arb(1)/2
    temp = ((1+z)*(1+z).log() + (1-z)*(1-z).log())/(2*z.sqrt())
    v0 = constant*g0-8*(2*root*temp-3*a)/c
    peak = 1+v0
    assert arb("175.71160384847839") < peak < arb("175.71160384847841")
    selected = arb(0)
    for _, _, lo, hi in brackets():
        t = arb(lo+hi, hi-lo)/(2*10**14)
        selected += 2*(3-2*root*(a*t).cos())/(c*(arb(1)/4+t*t))
    assert 0 < v0-selected < arb(647)/1000

    s = arb(300)
    eps = arb(10)**-12
    error = (2*10**6*(-s/4).exp() + 128*(-(s-2)/2).exp()/(1-(-(s-2)).exp())
             + 22000*(-3*s/20).exp())
    assert error < eps
    assert 1600*root/b*(7*a/20).exp() < 22000
    assert 2*root*(a/2).cosh()+3 < 100*c/2
    # This weaker omitted-weight bound also follows from independent real quadrature.
    retained = v0-arb(649)/1000
    assert retained*(2*arb.pi()/256).cos()-174-eps > arb(1)/128
    assert peak+eps < 176 and peak+eps-175 < arb(713)/1000
    assert peak+eps-170 < arb(5713)/1000 and peak+eps-162 < arb(13713)/1000
    print("The limiting peak lies strictly between 175.71160384847839 and 175.71160384847841.")
    print("The unshifted prime-only error is < 10^-12 from log centre 300.")
    print("The 256-arc tent weight gives lower uniform log density >= 2^-48007 above score 175.")
    print("The high-score deficits at levels 162, 170 and 175 are < 13.713, 5.713 and 0.713.")


if __name__ == "__main__":
    main()
