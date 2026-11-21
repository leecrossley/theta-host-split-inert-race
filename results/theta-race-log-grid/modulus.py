"""Bound the uniform modulus of the dyadic zero profile from its positive even companion."""
from flint import arb, ctx

CONDUCTORS = (547, 4376, 8752, 8752, 727, 5816, 11632, 11632)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 160
    a, root = arb(2).log(), arb(2).sqrt()
    c = 2*(root-1)**2
    b = root*c
    lip = 2*(a/2).cosh()/b
    conductors = sum((q/arb.pi()).log() for q in CONDUCTORS)
    arch = conductors + 8*(arb(3)/4).digamma()
    tail_integral = (-arb(3)/2).exp()/((arb(3)/2)*(1-(-arb(2)).exp()))
    assert lip*(abs(arch)+16+32*tail_integral) < 309 and 16*lip < 70
    denominator = 2000000
    h = arb(1)/denominator
    cut = 500000
    assert 4*cut == denominator and a+arb(1)/4 < arb(3).log()
    half = arb(1)/2
    test = ((1+half)*(1+half).log()+(1-half)*(1-half).log())/(2*half.sqrt())
    arch_star = conductors-8*(arb.const_euler()+a)
    full = arch_star*4*(a/2).sinh()/b - 8*(2*root*test-3*a)/c
    assert 0 < full < 175
    logh, logcut = (1/h).log(), arb(cut).log()
    small = h*((309+70*logh)/3
               + (arb(4)/15)*((309+70*logh)*logcut-35*logcut**2))
    remaining = arb(350)/(2*(2*cut+1))
    error = (2/a)*(8/arb.pi())*(small+remaining)
    assert error < arb(141)/10000
    original_tail = (arb(144)/100)/(48*arb(3)/32)**2
    assert arb(4917)/1000-arb(285)/100-original_tail-arb(7)/1000-error > arb(19735)/10000
    print("The gamma-integral bound gives V(0)-V(s) <= s*(309+70*log(1/s)) through s=1/4.")
    print("Absolute-sine Fourier summation bounds the profile change by 0.0141 for every |h| <= 1/2000000.")


if __name__ == "__main__":
    main()
