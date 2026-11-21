"""Validate extra Hardy brackets and count the omitted spectrum character by character."""
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from math import isqrt
from os import environ

from flint import acb, arb, ctx, dirichlet_char

from finite_spectrum import ROOT, brackets, hardy_check, ordinate
from verify import odd_factors

LIMIT = 65536
CENTRES = range(0, 295, 3)


def extra_brackets():
    extra = [tuple(map(int, line.split())) for line in
             (ROOT/"evidence/extra-zeros.tsv").read_text().splitlines()]
    assert len(extra) == 702 and extra == sorted(set(extra))
    assert all(0 < lo < hi <= 300*10**14 for _, _, lo, hi in extra)
    return extra


def all_brackets():
    family = {(q, k) for q, k, _ in odd_factors()}
    table = sorted(brackets() + extra_brackets())
    last = {}
    for q, k, lo, hi in table:
        assert (q, k) in family and last.get((q, k), 0) < lo < hi
        assert not (lo <= 294*10**14 <= hi)
        last[q, k] = hi
    assert len(table) == 6702 and set(last) == family
    return table


def primes(n):
    flags = bytearray([1])*(n+1)
    flags[:2] = bytes(2)
    for p in range(2, isqrt(n)+1):
        if flags[p]:
            flags[p*p::p] = bytes(1+(n-p*p)//p)
    return [p for p in range(2, n+1) if flags[p]]


def exterior(job):
    q, k, positive, reflected = job
    ctx.prec = 128
    chi = dirichlet_char(q, k)
    assert chi.is_primitive() and chi.parity() == 1
    assert (chi*dirichlet_char(q, pow(k, -1, q))).is_principal()
    sums = [arb(0) for _ in CENTRES]
    for p in primes(LIMIT):
        character = chi(p)
        if character.is_zero():
            continue
        n = p
        logp = arb(p).log()
        base = character/(p*arb(p).sqrt())
        rotation = acb(0, -3*logp).exp()
        coefficient, step = base, rotation
        while n <= LIMIT:
            value = logp*coefficient
            for i in range(len(sums)):
                sums[i] += value.real
                value *= step
            n *= p
            coefficient *= base
            step *= rotation
    tail = 2*(arb(LIMIT).log()+2)/arb(LIMIT).sqrt()
    plus = list(map(ordinate, positive))
    minus = list(map(ordinate, reflected))
    for v, part in zip(CENTRES, sums, strict=True):
        s = acb(arb(3)/2, v)
        derivative = ((q/arb.pi()).log()+((s+1)/2).digamma().real)/2-part
        known = sum((1/(1+(t-v)**2) for t in plus), arb(0))
        known += sum((1/(1+(t+v)**2) for t in minus), arb(0))
        assert derivative+tail-known < arb(4)/13, (q, k, v)
    return len(sums)


def check():
    table = all_brackets()
    by_character = defaultdict(list)
    for row in table:
        by_character[row[:2]].append(row)
    jobs = [(q, k, rows, by_character[q, pow(k, -1, q)])
            for (q, k), rows in by_character.items()]
    assert len(jobs) == 8 and all(reflected for _, _, _, reflected in jobs)
    workers = max(1, min(12, int(environ.get("THETA_WORKERS", "6"))))
    extra = extra_brackets()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        checked = sum(pool.map(hardy_check, [extra[i:i+12] for i in range(0,len(extra),12)]))
        counts = sum(pool.map(exterior, jobs))
    assert checked == 702 and counts == 792
    count = sum(hi < 294*10**14 for _, _, _, hi in table)
    print(f"1,404 additional Hardy signs and 792 exterior Euler sums certify all {count} positive zeros through ordinate 294, with multiplicity.")
    print("The eight primitive characters use the conjugate character for each reflected bracket. The retained list has 6,702 positive brackets.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Use the verifier with assertions enabled.")
    check()
