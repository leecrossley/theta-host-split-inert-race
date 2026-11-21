"""Verify a continuous starting-height cover by positive spherical averages."""
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from os import environ
from sys import argv

from flint import acb, acb_series, arb, ctx

from finite_spectrum import ROOT, brackets, ordinate
from verify import odd_factors

FIRST = 2000
LAST = 12500
STEP = 128
D = 63
C = 400
CARRIER = 8500
DENOMINATORS = (16, 24, 32, 48, 64)


def blocks():
    return [(lo, min(lo+STEP, LAST)) for lo in range(FIRST, LAST, STEP)]


@lru_cache(maxsize=2)
def coefficients(method):
    a = arb(2).log()
    r = arb(2).sqrt()
    data = []
    for row in brackets():
        t = ordinate(row)
        assert 0 < t < CARRIER and not t.contains(190)
        u = a*t/2
        if method == "direct":
            numerator = 3-2*r*(2*u).cos()
        else:
            numerator = (r-1)**2+4*r*u.sin()**2
        beta = -2*numerator*u.sin()/(a*(r-1)**2*(t*t+arb(1)/4))
        data.append((t,beta,(CARRIER**2-t*t).sqrt()))
    return data


def sphere(p, t, radius, method):
    c = arb(CARRIER)
    if method == "direct":
        return (c/p)*(radius*p).sinh()/(radius*c).sinh()
    return (c/p)*(-radius*t*t/(c+p)).exp()*(-(-2*radius*p).expm1())/(-(-2*radius*c).expm1())


def parameters(top, den):
    length = arb(top+82)/(2*8492)
    d = 1/arb(den)
    return length, d, length+2*d


def extend_low_count():
    # The predecessor covers ordinates through 48. Five more exterior
    # logarithmic derivatives cover the rest through 63.
    ctx.prec = 144
    ts = [ordinate(row) for row in brackets()]
    factors = list(odd_factors())
    for v in (51,54,57,60,63):
        z = acb(arb(3)/2,v)
        total = arb(0)
        for q, _, character in factors:
            jet = acb_series([z,1],2).dirichlet_l(character)
            total += (((q/arb.pi()).log()+((z+1)/2).digamma())/2+jet[1]/jet[0]).real
        used = sum((1/(1+(g-v)**2)+1/(1+(g+v)**2) for g in ts),arb(0))
        assert 0 < total-used < arb(1)/4 < arb(4)/13


def check_high(top):
    length, _, _ = parameters(top, 16)
    g = arb(8590)
    q = (g*g-CARRIER**2-arb(1)/4).sqrt()
    assert q > 1200 and g/(2*q) < 8
    # The coefficient sum is under 10^6. The two long boxes cost at most exp(d).
    rest = 10**6*30*((arb(top+7)/2+arb(1)/16)-8492*length).exp()
    trivial = 48*(-arb(3)*(FIRST-2)/2).exp()/(1-(-2*arb(FIRST-2)).exp())
    assert rest+trivial < arb(1)/1000


def verify_block(job):
    lo, top, rows, method = job
    ctx.prec = 128 if method == "direct" else 160
    a = arb(2).log()
    offset = 3*a/2
    width = arb(C//2).log()
    check_high(top)
    length, _, _ = parameters(top,16)
    data = coefficients(method)
    tf = [(t, beta*sphere(p,t,length,method)) for t,beta,p in data]
    by_den = defaultdict(list)
    windows = {-1:[],1:[]}
    for sign, m, den in rows:
        assert sign in windows and den in DENOMINATORS
        x = arb(m)/64
        assert lo-2 <= x < top+7
        _, _, radius = parameters(top,den)
        assert x-radius > FIRST and x-radius-offset > 100
        by_den[den].append((sign,m))
        windows[sign].append((x+radius-width,x-radius))

    for sign, intervals in windows.items():
        if method == "direct":
            frontier = arb(lo)+offset
            for start, stop in intervals:
                assert start < frontier < stop
                frontier = stop
            assert frontier > top+offset
        else:
            # Merge the open intervals after an independent sort by left endpoint.
            ordered = sorted(intervals,key=lambda pair: float(pair[0].mid()))
            frontier = arb(lo)+offset
            for start,stop in ordered:
                assert start < frontier
                frontier = frontier.max(stop)
            assert frontier > top+offset

    for den, points in by_den.items():
        _, d, _ = parameters(top,den)
        p0 = (arb(CARRIER)**2-D**2).sqrt()
        line = 50*(arb(127)/1000)*sphere(p0,arb(D),length,method)*min(arb(1),1/(D*d)**2)
        cost = arb(2353)/1000+arb(1)/1000+line+arb(19735)/10000
        retained = []
        for t, value in tf:
            if method == "direct":
                box = (d*t).sinc()**2
            else:
                box = (1-(2*d*t).cos())/(2*(d*t)**2)
            b = value*box
            if t > 190:
                cost += abs(b)
            elif method == "direct":
                retained.append((t,b))
            else:
                sn, cs = (t*lo).sin_cos()
                retained.append((t,b*sn,b*cs))
        for sign, m in points:
            if method == "direct":
                x = arb(m)/64
                total = sum((b*(t*x).sin() for t,b in retained),arb(0))
            else:
                delta = arb(m-64*lo)/64
                total = arb(0)
                for t, sine, cosine in retained:
                    sd,cd = (t*delta).sin_cos()
                    total += sine*cd+cosine*sd
            assert sign*total > cost, (top,den,m,sign,sign*total,cost)
    return len(rows)


def main(method):
    if not __debug__ or method not in ("direct","split"):
        raise ValueError("Use direct or split, with assertions enabled.")
    if method == "direct":
        extend_low_count()
        print("Five exterior logarithmic derivatives extend low-ordinate completeness through 63.",flush=True)
    groups = defaultdict(list)
    table = ROOT/'evidence'/'sphere-cover.tsv'
    for line in table.read_text().splitlines():
        top,sign,m,den = map(int,line.split())
        groups[top].append((sign,m,den))
    ranges = blocks()
    assert set(groups) == {hi for _,hi in ranges}
    count = sum(map(len,groups.values()))
    jobs = [(lo,hi,groups[hi],method) for lo,hi in ranges]
    workers = max(1,min(12,int(environ.get("THETA_WORKERS","6"))))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        assert sum(pool.map(verify_block, jobs)) == count
    print(f"{method}: {count} signed spherical averages cover both moves at every real log height from {FIRST} through {LAST}.")
    print(f"{method}: all high zeros, omitted line ordinates, arithmetic corrections and strict endpoint inequalities passed.")


if __name__ == "__main__":
    main(argv[1])
