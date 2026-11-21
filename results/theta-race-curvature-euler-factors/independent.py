"""Cross-check the sample inversion with convergent gamma series and primitive character values."""
from fractions import Fraction
from functools import lru_cache
from math import gcd
from flint import arb, ctx, dirichlet_char
from curvature import CASES, CONDUCTORS, ball, check_inputs, locations, read_samples

LABELS = ((547,546),(4376,1093),(8752,1093),(8752,5469),
          (727,726),(5816,1453),(11632,7269),(11632,1453))


@lru_cache(maxsize=None)
def series(ratio):
    ratio = min(ratio, 1/ratio)
    if ratio == 1:
        return arb(2).log()
    z = ball(ratio)
    term = z*z.sqrt()
    z2 = z*z
    total = arb(0)
    for j in range(640):
        # (2j+3/2)^2-1/4 = (2j+1)(2j+2).
        total += term/((2*j+1)*(2*j+2))
        term *= z2
    tail = term/((1281*1282)*(1-z2))
    assert tail < arb(10)**-48
    return total+arb(0,1)*tail


def correction(x):
    r, a = arb(2).sqrt(), arb(2).log()
    c = 2*(r-1)**2
    g = 4*((a-abs(ball(x).log())).max(0)/2).sinh()/(r*c)
    constant = sum((q/arb.pi()).log() for q in CONDUCTORS)-8*(arb.const_euler()+a)
    return constant*g-8/c*(r*series(2*x)+r*series(x/2)-3*series(x))


def characters(p):
    indices = []
    for q, label in LABELS:
        chi = dirichlet_char(q, label)
        assert chi.is_primitive() and chi.conductor() == q
        if gcd(p,q) != 1:
            continue
        n = int(chi.group().exponent())
        a = 4*int(chi.chi_exponent(p))
        assert a % n == 0
        indices.append((a//n)%4)
    assert indices.count(1) == indices.count(3)
    return [indices.count(k) for k in (0,2,1)]


def main():
    check_inputs()
    ctx.prec = 288
    table = read_samples()
    origin = correction(Fraction(1))
    c = 2*(arb(2).sqrt()-1)**2
    for (p, m), values in table.items():
        n = p**m
        xs = locations(n)
        f = [e/2-origin+correction(x) for e,x in zip(values,xs)]
        half = (ball(xs[2])/n).log()/2
        der_left = (f[1]*half.cosh()-f[0])/(2*half.sinh())
        der_right = (f[2]-f[1]*half.cosh())/(2*half.sinh())
        encoded = (der_right-der_left)*c*arb(n).sqrt()/arb(p).log()
        plus, minus, pairs = characters(p)
        def power(j):
            return plus+(-1)**j*minus+(0 if j%2 else 2*(-1)**(j//2)*pairs)
        if p == 2:
            exact = power(m+1)-3*power(m)+2*(0 if m==1 else power(m-1))
        else:
            exact = -3*power(m)
        assert encoded.contains(exact) and abs(encoded-exact)<arb(10)**-15
    # Exact recovery for each possible conjugation-stable fourth-root multiset of degree <= 8.
    cases=0
    for plus in range(9):
        for minus in range(9-plus):
            for pairs in range((8-plus-minus)//2+1):
                u,v,z = plus-minus,plus+minus-2*pairs,plus+minus+2*pairs
                assert (z+v+2*u,z+v-2*u,z-v) == (4*plus,4*minus,4*pairs)
                cases+=1
    assert set(p for p,m in table) == set(CASES)
    print("Independent gamma-series corrections and two one-sided derivatives recover every sample jump.")
    print("All ten example primes match direct primitive Conrey characters, including 2, 547 and 727.")
    print(f"Exact multiplicity inversion also holds on all {cases} fourth-root multisets of degree at most eight.")


if __name__ == "__main__":
    main()
