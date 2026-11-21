"""Certify all four joint doubling sign patterns by two prime-ideal counts."""
from fractions import Fraction
from math import lcm
from flint import arb, ctx, fmpz, nmod_poly

STARTS = {17: (-1, 1), 97: (-1, -1), 162: (1, -1), 658: (1, 1)}
HOSTS = ((547, 727), (571, 599))


def coefficient(p, u, v):
    split = (int(fmpz(p).jacobi(u)), int(fmpz(p).jacobi(v)))
    f = next(j for j in (1, 2, 4)
             if pow(p, j, 16) in (1, 15) and (split[0] * split[1])**j == 1)
    return p**f, -8 * split[0]**f


def degrees(poly, p):
    _, factors = nmod_poly(poly, p).factor()
    assert all(m == 1 for _, m in factors)
    return [int(z.degree()) for z, _ in factors]


def polynomial_coefficient(p, u, v):
    # The three linearly disjoint defining fields are generated in this order.
    quartic = degrees([2, 0, -4, 0, 1], p)
    real = degrees([-u*v, 0, 1], p)
    imaginary = degrees([u, 0, 1], p)
    f, = {lcm(d, e) for d in quartic for e in real}
    g, = {lcm(f, d) for d in imaginary}
    return p**f, -8 if g == f else 8


def ball(r):
    return arb(r.numerator) / r.denominator


def weight(n, y):
    return max(Fraction(0), min(2 - Fraction(2*y, n), Fraction(4*y, n) - 1))


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 144
    b = 2*arb(2).sqrt()*(arb(2).sqrt()-1)**2
    bounds = {}
    for i, (u, v) in enumerate(HOSTS):
        terms = []
        for p in range(3, 8*max(STARTS)+1, 2):
            if p in (u, v) or not fmpz(p).is_prime():
                continue
            n, c = coefficient(p, u, v)
            assert (n, c) == polynomial_coefficient(p, u, v)
            if n <= 8*max(STARTS):
                terms.append((n, c*arb(p).log()))
        for y, signs in STARTS.items():
            raw = [sum((logweight * ball(weight(n, s)) for n, logweight in terms), arb(0))
                   for s in (y, 2*y)]
            delta = (raw[1]/arb(2*y).sqrt()-raw[0]/arb(y).sqrt())/(b*arb(2).log())
            assert delta * signs[i] > 10
            # Display a pair of rational bounds, verified before printing.
            lo = int((1000*delta).floor().unique_fmpz())
            assert arb(lo)/1000 < delta < arb(lo+1)/1000
            bounds[i, y] = (lo, lo+1)
    for y in STARTS:
        print(f"Start {y}: normalized moves in {bounds[0,y]}/1000 and {bounds[1,y]}/1000.")
    print("All four sign pairs exceed magnitude 10 in both races. The Frobenius and polynomial counts agree.")


if __name__ == "__main__":
    main()
