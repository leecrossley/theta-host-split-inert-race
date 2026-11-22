"""Rebuild local affine coefficients from prime logs and verify rational endpoint signs."""

from bisect import bisect_right
from fractions import Fraction
from math import isqrt

from flint import arb, ctx, fmpz

LIMIT = 16 * 46696
SCALE = 10**18


def real(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


class LocalWeights:
    def __init__(self):
        rows = []
        for p in range(3, LIMIT + 1, 2):
            if p in (571, 599) or not fmpz(p).is_prime():
                continue
            signs = [1 if pow(p, (m - 1) // 2, m) == 1 else -1 for m in (571, 599)]
            deg = next(
                d
                for d in (1, 2, 4)
                if pow(p, d, 16) in (1, 15) and (signs[0] * signs[1]) ** d == 1
            )
            if p**deg <= LIMIT:
                rows.append((p**deg, -(signs[0] ** deg) * arb(p).log()))
        self.norms = []
        self.weights = [arb(0)]
        self.reciprocals = [arb(0)]
        for norm, weight in sorted(rows):
            self.norms.append(norm)
            self.weights.append(self.weights[-1] + weight)
            self.reciprocals.append(self.reciprocals[-1] + weight / norm)
        self.log_two = arb(2).log()
        root = arb(2).sqrt()
        self.factor = 8 / (2 * root * (root - 1) ** 2)

    def level_coefficients(self, y, k):
        lo, mid, hi = (bisect_right(self.norms, int(s * k * y)) for s in (1, 2, 4))
        first = self.weights[mid] - self.weights[lo]
        last = self.weights[hi] - self.weights[mid]
        reciprocal = self.reciprocals
        u = 4 * (reciprocal[hi] - reciprocal[mid]) - 2 * (reciprocal[mid] - reciprocal[lo])
        return (self.factor * arb(k).sqrt() * u, self.factor * (2 * first - last) / arb(k).sqrt())

    def move_coefficients(self, y, k):
        low = self.level_coefficients(y, k)
        high = self.level_coefficients(y, 2 * k)
        return tuple((b - a) / self.log_two for a, b in zip(low, high))

    def move(self, y, k):
        u, v = self.move_coefficients(y, k)
        return (u * real(y) + v) / real(y).sqrt()


def enclosure(race, k, cutoff, target):
    x, z = Fraction(cutoff, SCALE), Fraction(cutoff + 1, SCALE)
    sign = 1 if k == 1 else -1
    assert sign * (race.move(x, k) - target) < 0 < sign * (race.move(z, k) - target)
    return real(x).union(real(z))


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 238
    race = LocalWeights()
    centre = Fraction("46695.215")
    u, v = race.move_coefficients(centre, 1)
    s, t = race.move_coefficients(centre, 2)
    y = (t - v) / (u - s)
    height = y.sqrt()
    margin = (u * y + v) / height
    ds = (u * y - v) / (2 * y * height)
    dt = (s * y - t) / (2 * y * height)
    widening = 1 / ds - 1 / dt
    assert arb("1472.01872576754168") < widening < arb("1472.01872576754169")
    log_constant = y * race.log_two / widening
    assert arb("21.98793864880736") < log_constant < arb("21.98793864880737")
    for left, right, q, counts in (
        (None, 46695333262779289365032, race.move(Fraction(46695), 1), (4, 97121)),
        (
            46695215212591948194411,
            46695216684610667161768,
            margin - arb(1) / 10**6,
            (680, 21987939),
        ),
    ):
        lo = real(46695) if left is None else enclosure(race, 1, left, q)
        hi = enclosure(race, 2, right, q)
        for width, count in zip((hi - lo, (hi / lo).log() / race.log_two), counts, strict=True):
            assert (count - 1) * width < 1 < count * width
    p = 380207
    assert all(p % d for d in range(2, isqrt(p) + 1))
    assert p % 16 == 15 and all(pow(p, (m - 1) // 2, m) == 1 for m in (571, 599))
    assert 8 * 46696 < p < 16 * 46695
    print(
        "Independent affine coefficients and endpoint signs verify all four optimal phase counts."
    )
    print("The octave sampling constant lies in (21.98793864880736, 21.98793864880737).")


if __name__ == "__main__":
    main()
