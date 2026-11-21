"""Check the spectral-count constants and a compact differential test."""
from fractions import Fraction as Q
from math import factorial

from flint import acb, acb_poly, arb, ctx

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")


def log_bounds(x):
    z = (x-1)/(x+1)
    total = Q(0)
    power = z
    for j in range(50):
        total += 2*power/(2*j+1)
        power *= z*z
    return total, total+2*power/(101*(1-z*z))


def atan_bounds(x):
    total = sum(((-1)**j*x**(2*j+1)/(2*j+1) for j in range(24)), Q(0))
    return total, total+x**49/49


def comparison():
    # 200 = 2^7 * 25/16. Machin's formula and two atanh series suffice.
    a, b = atan_bounds(Q(1, 5))
    c, d = atan_bounds(Q(1, 239))
    pi0, pi1 = 16*a-4*d, 16*b-4*c
    lo, hi = log_bounds(Q(2))
    x, y = log_bounds(Q(25, 16))
    log0, log1 = 7*lo+x, 7*hi+y
    assert 14487 < 8590*log0/pi1 < 8590*log1/pi0 < 14488
    assert 2/(1-pi1/8590) < Q(2000732, 1000000)
    assert 11632*8590 < 10**8
    count = 8590*arb(200).log()/arb.pi()
    assert 14487 < count < 14488
    assert 2*(arb.pi()/8590).exp() < arb(2000732)/10**6


def integrate_polynomial(poly, z):
    # Independent integral of exp(zx) times a polynomial on [-1,1].
    # Truncate exp at degree 159 and bound its uniform Taylor remainder.
    degree = 160
    coeffs = poly.coeffs()
    powers = [acb(1)]
    for k in range(1, degree):
        powers.append(powers[-1]*z/k)
    ans = acb(0)
    for j, coefficient in enumerate(coeffs):
        ans += 2*coefficient*sum((powers[k]/(j+k+1)
                                 for k in range(j % 2, degree, 2)), acb(0))
    norm = sum((abs(c) for c in coeffs), arb(0))
    error = 2*norm*abs(z).exp()*abs(z)**degree/factorial(degree)
    return ans+acb(arb(0, error), arb(0, error))


def signed_test():
    target = acb(0, 4)
    excluded = [acb(0, -4), acb(arb(1)/8, 4), acb(arb(1)/8, -4),
                acb(arb(1)/5, 9), acb(arb(1)/5, -9)]
    polynomial = acb_poly([1])
    for point in excluded:
        polynomial *= acb_poly([-point, 1])*(1/(target-point))
    # The compact polynomial bump and its first seven derivatives vanish at both ends.
    bump = acb_poly([1, 0, -1])**8
    derivative = bump
    probe = acb_poly([])
    for order, weight in enumerate(polynomial.coeffs()):
        probe += (-1)**order*weight*derivative
        derivative = derivative.derivative()
    for z in [target]+excluded:
        integrated = integrate_polynomial(probe, z)
        factored = polynomial(z)*integrate_polynomial(bump, z)
        assert integrated.overlaps(factored)
        if z == target:
            assert abs(integrated) > arb(1)/10
        else:
            assert integrated.contains(0) and abs(integrated) < arb(10)**(-35)
    # Check the individual growing zero cancellation directly from its two weights.
    for z in excluded[1:]:
        u = z.real
        d = arb.pi()/(2*abs(z.imag))
        w = (u*d).exp()
        raw = (w*(-d*z).exp()+(d*z).exp()/w)/(w+1/w)
        assert raw.contains(0) and abs(raw) < arb(10)**(-50)


if __name__ == "__main__":
    ctx.prec = 256
    comparison()
    signed_test()
    print("rational and ball calculations: 14487 < 8590 log(200)/pi < 14488")
    print("one positive rightmost ordinate gives an infimum strictly below 2.000732")
    print("a direct compact integral and a differential multiplier agree at six complex exponents")
    print("the target survives; all five other moments, including a larger real part at the same ordinate, vanish")
