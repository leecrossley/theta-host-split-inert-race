"""Certify dyadic moves and cover a continuous range of starting heights."""
from bisect import bisect_right
from fractions import Fraction
from math import isqrt

from flint import arb, ctx

if not __debug__:
    raise RuntimeError("Run these checks with assertions enabled.")

MAX_HEIGHT = 1_000_000
MAX_START = 4_000_000


def quadratic_signs(p):
    residue = [-1]*p
    for a in range(1, p):
        residue[a*a % p] = 1
    return residue


def prime_ideals(limit):
    mask = bytearray([1]) * (limit + 1)
    for p in range(2, isqrt(limit) + 1):
        if mask[p]:
            mask[p*p::p] = bytes((limit - p*p)//p + 1)
    t1, t2 = quadratic_signs(547), quadratic_signs(727)
    ideals = []
    for p in range(3, limit + 1, 2):
        if not mask[p] or p in (547, 727):
            continue
        a, b = t1[p % 547], t2[p % 727]
        # Residue degree in Q(zeta_16)^+ and in Q(sqrt(547*727)).
        if p % 16 in (1, 15) and a == b:
            degree = 1
        elif (p*p) % 16 in (1, 15):
            degree = 2
        else:
            degree = 4
        norm = p**degree
        if norm <= limit:
            ideals.append((norm, p, -a**degree))
    return sorted(ideals)


def locate_moves():
    grid = [1]
    while grid[-1] < MAX_START:
        grid.append((9*grid[-1] + 7)//8)
    sample = sorted(set(grid + [2*y for y in grid]))
    edges = sorted({v*k for v in sample for k in (1, 2, 4)})
    ideals = iter(prime_ideals(edges[-1]))
    row = next(ideals, None)
    prefix = {}
    coeff = arb(0)
    coeff_over_n = arb(0)
    for edge in edges:
        while row and row[0] <= edge:
            n, p, sign = row
            c = 8*sign*arb(p).log()
            coeff += c
            coeff_over_n += c/n
            row = next(ideals, None)
        prefix[edge] = coeff, coeff_over_n

    b = arb(2).sqrt()
    scores = {}
    for y in sample:
        a0, b0 = prefix[y]
        a1, b1 = prefix[2*y]
        a2, b2 = prefix[4*y]
        total = 2*(a1-a0) - 2*y*(b1-b0) + 4*y*(b2-b1) - (a2-a1)
        scores[y] = total/(2*b*(b-1)**2*arb(y).sqrt())
    threshold = (arb(19735)/10000) * arb(2).log()
    gains, losses = [], []
    for y in grid:
        change = scores[2*y] - scores[y]
        if change > threshold:
            gains.append(y)
        if change < -threshold:
            losses.append(y)
    return gains, losses


def cover(low, high, gains, losses):
    # On each half-open interval the next strictly later witnesses are fixed.
    checks = sorted({low, high, *(y for y in gains + losses if low <= y <= high)})
    worst = Fraction(0)
    for h in checks:
        g, l = gains[bisect_right(gains, h)], losses[bisect_right(losses, h)]
        worst = max(worst, Fraction(2*max(g, l), h))
    return worst


def main():
    ctx.prec = 112
    gains, losses = locate_moves()
    for h in (1, 100, 1000, MAX_HEIGHT):
        g, l = gains[bisect_right(gains, h)], losses[bisect_right(losses, h)]
        print(f"H={h}: gain {g} to {2*g}, loss {l} to {2*l}")
    for h, bound in ((1, 18), (100, Fraction(6851,657)), (1000, Fraction(6851,657))):
        w = cover(h, MAX_HEIGHT, gains, losses)
        assert w <= bound
        print(f"starting heights [{h},{MAX_HEIGHT}]: verified endpoint <= ({w}) H")


if __name__ == "__main__":
    main()
