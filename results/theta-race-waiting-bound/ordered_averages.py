"""Enclose every selected positive spherical average of the actual signed prime move."""
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from os import environ
import lzma
import sys

from flint import arb, ctx

from carrier_window import square_and_higher_powers
from deadline_spectrum import all_brackets
from finite_deadlines import arithmetic
from finite_spectrum import ROOT, ordinate
from finite_sphere import CARRIER, DENOMINATORS, parameters, sphere

D = 294
KEEP = 270


def supports():
    table = lzma.decompress((ROOT/"evidence/ordered-supports.tsv.xz").read_bytes())
    rows = [tuple(map(int, line.split())) for line in table.splitlines()]
    assert len(set(rows)) == len(rows) == 49868
    for top, sign, m, k in rows:
        assert sign in (-1, 1) and k in DENOMINATORS
        assert 20*64 <= m < top*64 and top <= 12548
    return rows


@lru_cache(maxsize=2)
def spectrum(method):
    a, r = arb(2).log(), arb(2).sqrt()
    modes = []
    for row in all_brackets():
        t = ordinate(row)
        assert 0 < t < CARRIER and not t.contains(KEEP)
        u = a*t/2
        numerator = 3-2*r*(2*u).cos() if method == "direct" else (r-1)**2+4*r*u.sin()**2
        beta = -2*numerator*u.sin()/(a*(r-1)**2*(t*t+arb(1)/4))
        modes.append((t, beta, (CARRIER**2-t*t).sqrt()))
    return modes


@lru_cache(maxsize=None)
def prime_cost(ell):
    assert ell >= 21
    value = arithmetic(ell)
    if ell >= 100:
        value = value.min(square_and_higher_powers())
    if ell < 500:
        return value
    x = arb(ell)
    q = arb(6362704)
    cuts = [arb(91)/10] + [arb(v) for v in (12,15,18,21,24,27,30,36)]
    def mass(h):
        t = h.exp()
        assert t > 8590
        return 80*((2*q*t).log()+1)/(t-2)
    base = 132*208*(-x/(20*(q.log()+cuts[0]))).exp()
    for low, high in zip(cuts, cuts[1:]):
        base += 132*mass(low)*(-x/(20*(q.log()+high))).exp()
    base += 132*mass(cuts[-1])+20/(arb(1).exp()*x)
    base += 10**6*((-x/10).exp()+(-x/3).exp())
    base += 128*(-(x-2)/2).exp()/(1-(-2*(x-2)).exp())
    return value.min(base)


def kernel_constants():
    a, r = arb(2).log(), arb(2).sqrt()
    # Small imaginary arguments use the total variation of the signed kernel.
    assert 5*(4/a)*(3*a/4).exp() < 66
    # Above 2 the rational denominator costs (1+t^2)/t^2 <= 5/4.
    assert (arb(5)/4)*(2*r*(a/2).cosh()+3)*(a/4).cosh()/((r-1)**2*a) < 66
    assert 2*8*13 == 208
    assert (arb(91)/10).exp() > 8590


def check_rows(job):
    top, method, rows = job
    ctx.prec = 128 if method == "direct" else 160
    a, offset = arb(2).log(), 3*arb(2).log()/2
    length, _, _ = parameters(top, 16)
    weights = [(t, b*sphere(p,t,length,method)) for t,b,p in spectrum(method)]
    p0 = (arb(CARRIER)**2-D**2).sqrt()
    line = 50*(arb(127)/1000)*sphere(p0,arb(D),length,method)
    barrier = (arb(8590)**2-CARRIER**2-arb(1)/4).sqrt()
    assert barrier > 1200 and 8590/(2*barrier) < 8
    high = 30*10**6*((arb(top)+7)/2+arb(1)/16-8492*length).exp()
    by_den = defaultdict(list)
    for sign,m,k in rows:
        x = arb(m)/64
        _, _, radius = parameters(top,k)
        earliest = x-radius
        assert 20 < earliest and x < top+7
        trivial = 48*(-3*(earliest-2)/2).exp()/(1-(-2*(earliest-2)).exp())
        assert high+trivial < arb(1)/1000
        ell = int((earliest-offset).floor().unique_fmpz())
        by_den[k].append((sign, m, ell))
    base_m = min(m for _,m,_ in rows)
    for k, centres in by_den.items():
        d = 1/arb(k)
        base_cost = arb(19735)/10000 + arb(1)/1000 + line*arb(1).min(1/(D*d)**2)
        use = []
        for t, amplitude in weights:
            box = (d*t).sinc()**2 if method == "direct" else (1-(2*d*t).cos())/(2*(d*t)**2)
            coefficient = amplitude*box
            if t > KEEP:
                base_cost += abs(coefficient)
            elif method == "direct":
                use.append((t,coefficient))
            else:
                sp, cp = (t*base_m/64).sin_cos()
                use.append((t,coefficient*sp,coefficient*cp))
        for sign,m,ell in centres:
            if method == "direct":
                x = arb(m)/64
                total = sum((b*(t*x).sin() for t,b in use), arb(0))
            else:
                delta = arb(m-base_m)/64
                total = arb(0)
                for t,sine,cosine in use:
                    sd,cd = (t*delta).sin_cos()
                    total += sine*cd+cosine*sd
            assert sign*total > base_cost+prime_cost(ell), (top,k,sign,m,ell)
    return len(rows)


def check(method):
    if not __debug__ or method not in ("direct","split"):
        raise ValueError("Use direct or split, with assertions enabled.")
    ctx.prec = 128
    kernel_constants()
    groups = defaultdict(list)
    rows = supports()
    for top, sign, m, k in rows:
        groups[top].append((sign,m,k))
    workers = max(1,min(12,int(environ.get("THETA_WORKERS","6"))))
    jobs = [(t,method,g) for t,g in sorted(groups.items(),key=lambda pair:-len(pair[1]))]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        assert sum(pool.map(check_rows,jobs)) == len(rows)
    print(f"{method}: all 49,868 signed spherical averages exceed the threshold after the odd-zero, base-field and prime-power errors.")


if __name__ == "__main__":
    check(sys.argv[1])
