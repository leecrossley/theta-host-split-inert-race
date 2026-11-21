"""Check the signed prime-band secants and the order conditions in the mixed-field proof."""
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product

from check_jump import DA, DB


@dataclass(frozen=True)
class Quadratic:
    """An exact a + b/sqrt(2)."""

    a: Fraction = Fraction()
    b: Fraction = Fraction()

    def __neg__(self):
        return Quadratic(-self.a, -self.b)

    def __sub__(self, other):
        return Quadratic(self.a - other.a, self.b - other.b)

    def __mul__(self, other):
        return Quadratic(self.a * other.a + Fraction(self.b * other.b, 2),
                         self.a * other.b + self.b * other.a)

    def __truediv__(self, other):
        denominator = other.a**2 - Fraction(other.b**2, 2)
        assert denominator
        return self * Quadratic(Fraction(other.a, denominator),
                                Fraction(-other.b, denominator))

    def positive(self):
        if self.a >= 0 and self.b >= 0:
            return self.a > 0 or self.b > 0
        if self.a <= 0 and self.b <= 0:
            return False
        return (2 * self.a**2 > self.b**2) == (self.a > 0)

    def absolute(self):
        return self if self.positive() else -self


ZERO = Quadratic()
ROOT_TWO = Quadratic(0, 2)


def band_rows():
    """Recover (n U, V) from the four independently replayed knot jumps."""
    return [
        tuple(Quadratic(sum(row[0] for row in jumps[:k]),
                        sum(row[1] for row in jumps[:k]))
              for jumps in (DA, DB))
        for k in range(5)
    ]


def signed_types(left, right):
    return {(ZERO, ZERO), left, right, tuple(-v for v in left),
            tuple(-v for v in right)}


def secants():
    off, outer, middle, inner, right_off = band_rows()
    assert off == right_off == (ZERO, ZERO)
    ratios = [
        (outer, middle, {Quadratic(Fraction(1,8)), Quadratic(Fraction(1,4)),
                         Quadratic(0,Fraction(1,2)), Quadratic(Fraction(3,4),1)}),
        (middle, inner, {Quadratic(Fraction(3,2),-2), Quadratic(0,Fraction(1,2)),
                        Quadratic(Fraction(1,2)), Quadratic(1)}),
    ]
    for first, second, expected in ratios:
        types = signed_types(first,second)
        slopes = set()
        for (u,v),(w,z) in combinations(types,2):
            assert u != w and v != z
            slopes.add(((v-z)/(u-w)).absolute())
        assert slopes == expected
        for a,b in combinations(slopes,2):
            low, high = (a,b) if (b-a).positive() else (b,a)
            assert not (ROOT_TWO*low-high).positive()
        assert ((second[1]/second[0])/(first[1]/first[0])) == Quadratic(0,4)
    print("Both adjacent-band secant tables agree with the original knot jumps. Distinct absolute secants differ by a factor of at least sqrt(2). Each late-to-early within-band factor is 2 sqrt(2).")


def ordering():
    for signs in product((-1,1), repeat=4):
        # Only order matters. |x_1| < ... < |x_4| makes this permutation
        # identical to that for any positive strictly increasing scales.
        nodes = [sign*i for i,sign in enumerate(signs,1)]
        order = sorted(range(4),key=nodes.__getitem__)
        positions = [order.index(i) for i in range(4)]
        assert abs(positions[0]-positions[1]) == 1
        assert (positions[2]-positions[3]) % 2 == 1
    print("For all 16 target-sign orders both chronological pairs have opposite barycentric weight signs.")


def windows_and_counts():
    f = Fraction(101,100)
    r = Fraction(11,10)
    assert f**8 < r and r**6 < 2
    assert (f-1)/8 - (f+1)/3200 > Fraction(1,2000)
    assert 8*10**9 // 2000 > 300
    assert Fraction(19,16)**4 < 2
    assert (8*10**9)**4 * 2 < 2**1582
    assert (9848+2)//3 == 3283
    assert (4923+2)//3 == 1641 and (4925+2)//3 == 1642
    assert (1072+2)//3 == 358
    print("Two prime fourth powers suffice in each large factor-11/10 interval. The combined certified crossing counts yield 3,283 classes overall, 1,641 for gains and 1,642 for losses.")


if __name__ == "__main__":
    secants()
    ordering()
    windows_and_counts()
