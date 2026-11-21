"""Independent series for the selected zero and direct prime-by-prime window sums."""
from math import comb, isqrt

from flint import acb, arb, ctx, fmpz

if not __debug__:
    raise RuntimeError("Run these checks with assertions enabled.")


def completed(t):
    # Split the Dirichlet series after one period and expand each Hurwitz tail at 3/2.
    q, kmax = 547, 70
    s = acb(arb(1)/2, t)
    assert abs(s) < 4 and fmpz(q).is_prime()
    first = acb(0)
    moments = [0] * (kmax + 1)
    for n in range(1, q):
        chi = 1 if pow(n, (q-1)//2, q) == 1 else -1
        first += chi * acb(n)**(-s)
        d, power = 2*n-q, 1
        for k in range(kmax + 1):
            moments[k] += chi*power
            power *= d
    assert moments[0] == 0
    tail = acb(0)
    coeff = acb(1)
    for k in range(1, kmax + 1):
        coeff *= -(s + k - 1)/k
        tail += coeff * (arb(moments[k])/(2*q)**k) * (s+k).zeta(arb(3)/2)
    err = 8*arb(q).sqrt()*comb(kmax+4,3) / 3**(kmax+1)
    value = first + acb(q)**(-s)*tail + acb(arb(0,err),arb(0,err))
    u = (s+1)/2
    return (q/arb.pi())**u * u.gamma() * value


def direct():
    # Compute the original weights separately in each witness window, with Euler's criterion.
    signs = {5:1, 9:-1, 2474870:1, 1220778:-1}
    evaluations = sorted(set(signs) | {2*y for y in signs})
    limit = 4*max(evaluations)
    sieve = bytearray([1])*(limit+1)
    for p in range(2,isqrt(limit)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit-p*p)//p+1)
    quadratic = {q:[pow(a,(q-1)//2,q) for a in range(q)] for q in (547,727)}
    roots = {p: arb(p).sqrt() for p in evaluations}
    weights = {y:arb(0) for y in evaluations}
    for p in range(3,limit+1,2):
        if not sieve[p] or p in (547,727):
            continue
        u = 1 if quadratic[547][p%547] == 1 else -1
        v = 1 if quadratic[727][p%727] == 1 else -1
        f = next(k for k in (1,2,4) if (u*v)**k == 1 and pow(p,k,16) in (1,15))
        norm = p**f
        if norm > limit:
            continue
        ns = [y for y in evaluations if y < norm < 4*y]
        if ns:
            coeff = -(8//f)*u**f*arb(norm).log()
            for y in ns:
                t = arb(y)/norm
                weights[y] += coeff * (2*(1-t) if norm <= 2*y else 4*t-1)
    b = arb(2).sqrt()
    denominator = 2*b*(b-1)**2*arb(2).log()
    for y, sign in signs.items():
        slope = (weights[2*y]/roots[2*y]-weights[y]/roots[y])/denominator
        assert slope*sign > arb(19735)/10000
        lo = int(slope.floor().unique_fmpz())
        print(f"direct slope at y={y} lies in ({lo}, {lo+1})")


def main():
    ctx.prec = 144
    left, right = [completed(arb(n)/10**12) for n in (3287178298086,3287178298087)]
    assert left.imag.contains(0) and right.imag.contains(0)
    assert left.real*right.real < 0
    print("independent Hurwitz series confirms the crossing of the conductor-547 factor")
    direct()


if __name__ == "__main__":
    main()
