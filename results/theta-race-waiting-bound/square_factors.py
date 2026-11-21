"""Certify two complex zero-free disks for the factors of the base-field zeta."""
from math import isqrt

from flint import acb, arb, ctx

if not __debug__:
    raise RuntimeError("Run these checks with assertions enabled.")

P, Q = 547, 727
CUTOFF = 2**18
TWICE_CENTRE = 6574356596173
SCALE = 10**12
# xi(5)=i. Values are exponents of i, with 9 marking a nonunit.
QUARTIC = (9, 0, 9, 3, 9, 1, 9, 2, 9, 2, 9, 1, 9, 3, 9, 0)


def quadratic(p):
    table = bytearray([2]) * p
    table[0] = 9
    for a in range(1, p):
        table[a*a % p] = 0
    return table


def period(twist, j, quad):
    conductor = (P*Q if twist else 1) * (1, 16, 8, 16)[j]
    table = bytearray(conductor)
    re = im = max_norm = 0
    for n in range(1, conductor + 1):
        k = quad[n % len(quad)] if twist else 0
        if k == 9 or (j and n % 2 == 0):
            table[n % conductor] = 9
            continue
        k = (k + j*QUARTIC[n % 16]) % 4
        table[n % conductor] = k
        re += (1, 0, -1, 0)[k]
        im += (0, 1, 0, -1)[k]
        max_norm = max(max_norm, re*re + im*im)
    assert (re, im) == (0, 0)
    return table, isqrt(max_norm) + 1


def main():
    ctx.prec = 112
    t = arb(TWICE_CENTRE) / SCALE
    radius = arb(2001) / 10**6  # covers disks of radius 2/1000 at 1 and 1+2i gamma
    p1, p2 = quadratic(P), quadratic(Q)
    quad = bytes(9 if n % P == 0 or n % Q == 0 else
                 (p1[n % P] + p2[n % Q]) % 4 for n in range(P*Q))
    keys = [(tw, j) for tw in (0, 1) for j in range(4) if tw or j]
    periods = [period(*key, quad) for key in keys]
    roots = [acb(1), acb(0, 1), acb(-1), acb(0, -1)]
    vals = [[acb(0) for _ in range(4)] for _ in keys]
    for n in range(1, CUTOFF + 1):
        nn = arb(n)
        logn = nn.log()
        sn, cs = (-t * logn).sin_cos()
        base = 1 / nn
        c = acb(cs, sn) * base
        for row, (table, _) in zip(vals, periods):
            k = table[n % len(table)]
            if k == 9:
                continue
            a, b = roots[k] * base, roots[k] * c
            row[0] += a
            row[1] -= a * logn
            row[2] += b
            row[3] -= b * logn

    # Uniform second derivative, from |C(x)| <= min(x,q), q <= exp(16).
    sigma = arb(99)/100
    scale = (arb(16)/100).exp()
    ell = arb(16)
    i1 = scale * (ell**2/2 + ell/sigma + 1/sigma**2)
    i2 = scale * (ell**3/3 + ell**2/sigma + 2*ell/sigma**2 + 2/sigma**3)
    assert 2*i1 + 8*i2 < 30000
    logn = arb(CUTOFF).log()
    for key, (table, M), row in zip(keys, periods, vals):
        for j, height in enumerate((arb(0), t)):
            s = acb(1, height)
            err = M * (1 + abs(s)) / CUTOFF
            derr = M * (logn + 1 + abs(s)*(logn + 1)) / CUTOFF
            margin = (abs(row[2*j]) - err - radius*(abs(row[2*j + 1]) + derr)
                      - 15000*radius**2)
            assert margin > arb(1)/100, (key, j, margin)
        print(f"base twist={key[0]} quartic={key[1]} conductor={len(table)} "
              f"partial-sum bound={M}: both disk margins > 1/100")

    # zeta(s)=s/(s-1)-s integral {u} u^(-s-1) du.
    assert 1-radius - (1+radius)*radius/sigma > arb(9)/10
    assert abs(acb(1,t).zeta()) - 12*radius > arb(1)/100
    print("Riemann zeta: regularized central disk and nonreal disk exclude zero")


if __name__ == "__main__":
    main()
