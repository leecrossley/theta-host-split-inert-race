"""Two independent finite evaluations of the normalized dyadic move at its prime-norm knots."""
from bisect import bisect_right
from math import isqrt
from flint import arb, fmpz

MAX_START = 214710
OMITTED = {2, 547, 727, 343057911569, 343070636977}


def all_primes(limit):
    flag = bytearray([1]) * (limit + 1)
    flag[:2] = bytes(2)
    for p in range(2, isqrt(limit) + 1):
        if flag[p]:
            flag[p*p::p] = bytes((limit-p*p)//p + 1)
    return (p for p in range(3, limit+1, 2) if flag[p])


def event_sweep():
    """At k = 8y accumulate the changes of intercept and slope in sqrt(y)D(y)."""
    limit = 8*MAX_START
    root = arb(2).sqrt()
    scale = 2*root*(root-1)**2*arb(2).log()
    c0 = (-1/root, 3/root+1, -2/root-3, arb(2))
    c1 = (8/root, -12/root-4, 4/root+6, arb(-2))
    squares = [{i*i % q for i in range(1, q)} for q in (547, 727)]
    changes = {}
    for p in all_primes(limit):
        if p in OMITTED:
            continue
        u = 1 if p % 547 in squares[0] else -1
        v = 1 if p % 727 in squares[1] else -1
        if p % 16 in (1, 15) and u == v:
            degree = 1
        elif p*p % 16 in (1, 15):
            degree = 2
        else:
            degree = 4
        n = p**degree
        if n > limit:
            continue
        mass = -8*u**degree*arb(p).log()/scale
        for j, da, db in zip((1, 2, 4, 8), c0, c1, strict=True):
            k = j*n
            if k <= limit:
                a, b = changes.get(k, (arb(0), arb(0)))
                changes[k] = a+mass*da, b+mass*db/n
    a = b = arb(0)
    pieces = {}
    for k in sorted({8, limit, *changes}):
        da, db = changes.get(k, (arb(0), arb(0)))
        a, b = a+da, b+db
        y = arb(k)/8
        pieces[k] = (a, b, (a+b*y)/y.sqrt())
    return pieces


def direct_prefix():
    """Independent primality tests, Euler criteria and prefix sums of the original tent."""
    limit = 8*MAX_START
    terms = []
    for p in range(3, limit+1, 2):
        if p in OMITTED or not fmpz(p).is_prime():
            continue
        u, v = (int(pow(p, (q-1)//2, q) == 1)*2-1 for q in (547, 727))
        degree = next(f for f in (1, 2, 4)
                      if pow(p, f, 16) in (1, 15) and (u*v)**f == 1)
        n = p**degree
        if n <= limit:
            terms.append((n, -(u**degree)*arb(p).log()))
    norms, prefix, divided = [], [arb(0)], [arb(0)]
    for n, mass in sorted(terms):
        norms.append(n)
        prefix.append(prefix[-1]+mass)
        divided.append(divided[-1]+mass/n)
    keys = {8, limit}
    for n in norms:
        keys.update(k for k in (8*n, 4*n, 2*n, n) if k <= limit)
    root = arb(2).sqrt()
    normal = 2*root*(root-1)**2

    def level(k):
        # W(n/y) = 2-2y/n on (y,2y], then 4y/n-1 on (2y,4y].
        i, j, l = (bisect_right(norms, (k*m)//8) for m in (1, 2, 4))
        y = arb(k)/8
        inner = 2*(prefix[j]-prefix[i])-2*y*(divided[j]-divided[i])
        outer = 4*y*(divided[l]-divided[j])-(prefix[l]-prefix[j])
        return 8*(inner+outer)/(normal*y.sqrt())

    logtwo = arb(2).log()
    return {k: (level(2*k)-level(k))/logtwo for k in sorted(keys)}
