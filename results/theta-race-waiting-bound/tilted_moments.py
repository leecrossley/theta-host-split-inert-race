"""Check the spectral tilt and the exceptional-zero case on a finite test spectrum."""
from fractions import Fraction as Q
from itertools import product

from flint import acb, arb, ctx

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")


def disk(q):
    return arb(q.numerator)/q.denominator


def conditions():
    # For 1/2 < sigma < 1 the density exponent has a strict power saving.
    for sigma in (Q(1, 2)+Q(1, 10**12), Q(5, 8), Q(3, 4), Q(9999, 10000)):
        exponent = 3*(1-sigma)/(2-sigma)
        assert 1-exponent == (2*sigma-1)/(2-sigma) > 0
    assert ((Q(633, 100)-Q(62, 100)-Q(14, 100000)-Q(1, 1000))
            /Q(100007, 100000)-Q(2353, 1000)) > Q(333, 100)


def direct_and_factored():
    """Enumerate the atoms independently of the cosh transfer formula."""
    ctx.prec = 224
    g = 4
    # A double selected pole, a smaller real part at the same ordinate,
    # and one farther-right pole at a different ordinate.
    sigma, forward = Q(5, 8), Q(3, 4)
    factors = [(Q(0), g), (Q(0), 2*g)]
    factors.extend((forward-sigma, abs(9+k*g)) for k in (-1, 0, 1))
    atom_lists = []
    for offset, frequency in factors:
        radius = arb.pi()/(2*frequency)
        scale = (disk(offset)*radius).exp()
        total = scale+1/scale
        # The correction measure has the law of X-X', with X,X' independent.
        atom_lists.append([(-2*radius, 1/total**2), (arb(0), (scale**2+1/scale**2)/total**2),
                           (2*radius, 1/total**2)])
    atoms = [(sum(v[0] for v in row), product_mass(row)) for row in product(*atom_lists)]
    assert sum(mass for _, mass in atoms).contains(1)
    assert all(mass>0 for _,mass in atoms)

    def moment(z):
        direct = sum(mass*(z*x).exp() for x,mass in atoms)
        factored = acb(1)
        for offset, frequency in factors:
            d = arb.pi()/(2*frequency)
            v = disk(offset)
            factored *= ((z-v)*d).cosh()*((z+v)*d).cosh()/(v*d).cosh()**2
        assert direct.overlaps(factored)
        return direct

    for sign in (-1, 1):
        for harmonic in (-1, 0, 1):
            selected = moment(acb(0, sign*g+harmonic*g))
            if sign+harmonic == 0:
                assert selected.contains(1)
            else:
                assert selected.contains(0) and abs(selected)<arb(10)**(-55)
            other = moment(acb(disk(forward-sigma), sign*9+harmonic*g))
            assert other.contains(0) and abs(other)<arb(10)**(-55)
    # The absolute translated contribution from a smaller beta decays even at a matched ordinate.
    difference = disk(Q(9, 16)-sigma)
    assert (difference*1000).exp()<arb(10)**(-25)
    return len(atoms)


def product_mass(row):
    mass = arb(1)
    for _, weight in row:
        mass *= weight
    return mass


if __name__ == "__main__":
    conditions()
    size = direct_and_factored()
    print(f"tilted positive average: {size} enumerated atoms agree with the factored moments at all 12 signed harmonics")
    print("the matched pair survives, the growing terms vanish and the critical-line case exceeds 3.33 log 2")
