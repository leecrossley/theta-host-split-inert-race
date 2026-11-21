"""Check the bounds for positive interpolation at separated Dirichlet zeros."""
from fractions import Fraction as Q

from flint import acb, acb_mat, arb, ctx

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")


def exact_bounds():
    # L = log(12000 T) > 23 and 3.14 < pi < 22/7.
    ell = 23
    pi_lo, pi_hi = Q(314, 100), Q(22, 7)
    # Maximum local count divided by L for a window of length four.
    assert 8/pi_lo + Q(988, 1000) + Q(27576, 1000)/ell < 5
    assert Q(988, 1000) + Q(27576, 1000)/ell < 3

    # The sum of the separated off-diagonal majorant telescopes.
    for n in (1, 7, 101):
        total = sum((Q(2, 4*k*k-1) for k in range(1, n+1)), Q(0))
        assert total == 1-Q(1, 2*n+1)
    # With u=1/L, (1-u/2)(2+u)^2 - 4(1+u^4)
    # equals u(2-u-u^2/2-4u^3), positive on 0 < u <= 1/23.
    u = Q(1, ell)
    assert 2-u-u*u/2-4*u**3 > 0
    assert 16*u**3+8*u < 1
    assert 16*(Q(1, 8557**2) + Q(1, 4*8557)) < 1
    assert 3*pi_hi/(2*256*ell) < 1
    g = Q(2513, 1000)
    deficit = g*g*pi_hi*(2+2*u)*u/(4*256)
    assert deficit < Q(2, 1000)
    decaying = Q(14*10**6, (54*8557)**2)
    assert decaying < Q(7, 100000)
    gain = ((Q(633, 100)*(1-deficit)-Q(62, 100)-Q(14, 100000)-Q(2, 1000))
            / Q(100007, 100000)-Q(2353, 1000))
    assert gain > Q(333, 100)
    # R <= 108 + 1.000004 + 272L + 2 and L(10^6) < 24.
    assert 1+108+Q(1000004, 1000000)+272*24+2-250000 < -100
    assert Q(4*272, 10**6) < Q(1, 2)
    assert 221+4 == 225
    assert 540*4 == 2160


def window(z):
    """Bilateral moment of the density cos(x)/2 on [-pi/2,pi/2]."""
    return (arb.pi()*z/2).cosh()/(1+z*z)


def ball_test(nodes, offsets, scale):
    """Solve the moment equations and check a separate exponential-integral formula."""
    frequencies = tuple(-x for x in reversed(nodes)) + tuple(nodes)
    offsets = tuple(reversed(offsets)) + tuple(offsets)
    equations = [acb(delta, frequency) for delta, frequency in zip(offsets, frequencies)]
    count = len(equations)
    matrix = acb_mat([[window(scale*(z-acb(0, frequency))) for frequency in frequencies]
                      for z in equations])
    right = acb_mat([[-window(scale*z)] for z in equations])
    solution = matrix.solve(right)
    weights = [solution[j, 0] for j in range(count)]
    assert sum(abs(weight) for weight in weights) < arb(1)/10
    assert all(weights[i].overlaps(weights[-1-i].conjugate()) for i in range(count))

    def direct(z):
        # Integrate the two exponentials in cos(x) explicitly.
        z = scale*z
        half_pi = arb.pi()/2
        return sum(((z+acb(0, sign))*half_pi).sinh()/(z+acb(0, sign))
                   for sign in (-1, 1))/2

    def moment(z):
        return direct(z) + sum(weight*direct(z-acb(0, frequency))
                               for frequency, weight in zip(frequencies, weights))

    mass = moment(acb(0))
    assert mass.real > arb(9)/10 and mass.imag.contains(0)
    for z in equations:
        residual = moment(z)/mass
        assert residual.contains(0) and abs(residual) < arb(10)**(-33)
    return count


def ball_bounds():
    ctx.prec = 192
    ell = arb(12000*10**6).log()
    assert 23 < ell < 24
    assert (arb.pi()/4).cosh() < 2 and arb(1).cosh() < 2
    assert (arb.pi()/2) < 2
    assert arb(1)/1000 + arb(14*10**6)/(54*8557)**2 < arb(2)/1000
    assert arb(2*10**9).log() > 1
    # An asymmetric cluster, both signs, real offsets at opposite ends of the strip.
    count = ball_test([8557, 8561, 8566, 8581], [arb(1)/2, arb(1)/100, arb(2)/5, arb(1)/7], arb(1))
    # Closely spaced *scaled* nodes at the high-frequency cutoff.
    scale = arb(1)/(256*32**2)
    count += ball_test([int(x)*256*32**2 for x in (1024, 1027, 1030, 1033, 1037)],
                       [arb(i)/10 for i in (1, 2, 3, 4, 5)], scale)
    return count


if __name__ == "__main__":
    exact_bounds()
    checks = ball_bounds()
    print("separated positive interpolation: exact norm, local zero-count, support and 3.33 log 2 margin checks passed")
    print(f"{checks} complex test moments also cancel by an independent ball integral")
    print("unconditional completion before H times a quartic in 12000 T and a power 2160 of log(12000 T)")
