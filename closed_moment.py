"""Independently bound the same modulus by the closed positive gamma series in geometric blocks."""
from flint import arb, ctx


def gamma_series(x):
    z = (-x).exp()
    return ((1+z)*(1+z).log()+(1-z)*(1-z).log())/(2*z.sqrt())


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 192
    a, root = arb(2).log(), arb(2).sqrt()
    c = 2*(root-1)**2
    b = root*c
    lip = 2*(a/2).cosh()/b
    conductors = (547,4376,8752,8752,727,5816,11632,11632)
    ach = sum((q/arb.pi()).log() for q in conductors)-8*(arb.const_euler()+a)
    assert ach > 0
    Ta = gamma_series(a)
    full = ach*4*(a/2).sinh()/b - 8*(2*root*Ta-3*a)/c
    h = arb(1)/2000000
    cut = 500000
    number = 0
    first = 1
    bound = arb(0)
    while first <= cut:
        last = min(cut,max(first,103*first//100))
        s = last*h
        increment = ach*lip*s + 8/c*(root*(gamma_series(a-s)-Ta) + 3*(a-gamma_series(s)))
        reciprocal_mass = (1/arb(2*first-1)-1/arb(2*last+1))/2
        bound += reciprocal_mass*increment
        first = last+1
        number += 1
    bound += 2*full/(2*(2*cut+1))
    error = (2/a)*(8/arb.pi())*bound
    assert error < arb(85)/10000
    independent_tail = (arb(152)/100)/(48*arb(3)/32)**2
    assert arb(4917)/1000-arb(285)/100-independent_tail-arb(7)/1000-error > arb(19735)/10000
    print(f"Closed gamma-series calculation in {number} blocks independently bounds the profile change by 0.0085.")
    print("The independent 1.52 spectral tail still leaves the strict prime-race margin on the full interval.")


if __name__ == "__main__":
    main()
