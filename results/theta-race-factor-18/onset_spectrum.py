"""Verify the finite zero list and positive averages used for the first factor-18 miss."""

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from os import environ
from pathlib import Path
import lzma
import sys

from flint import acb, acb_series, arb, ctx, dirichlet_char
from splitting_square_bound import CAP, correction

ROOT = Path(__file__).resolve().parent
DEN = 10**14
FAMILY = {(547,546), (4376,1093), (8752,1093), (8752,5469),
          (727,726), (5816,1453), (11632,7269), (11632,1453)}
TABLE_HASH = "9a26df7437ee2b66c5d89fc458c0327c98f54b171e38bc575004f2b9a7c4eb6c"


@lru_cache(maxsize=1)
def table():
    data = lzma.decompress((ROOT/'evidence/onset-zeros.tsv.xz').read_bytes())
    assert sha256(data).hexdigest() == TABLE_HASH
    rows = [tuple(map(int,line.split())) for line in data.splitlines()]
    assert len(rows) == 13946 and rows == sorted(rows) and {row[:2] for row in rows} == FAMILY
    previous = {}
    for q,k,lo,hi in rows:
        assert previous.get((q,k),0) < lo < hi < 1365*DEN
        assert not lo <= 294*DEN <= hi
        previous[q,k] = hi
    return rows


def ordinate(row):
    _,_,lo,hi = row
    return arb(lo+hi,hi-lo)/(2*DEN)


def hardy(batch):
    ctx.prec = 112
    chars = {}
    for q,k,lo,hi in batch:
        if (q,k) not in chars:
            chars[q,k] = dirichlet_char(q,k)
        chi = chars[q,k]
        a,b = (chi.hardy_z(arb(n)/DEN) for n in (lo,hi))
        assert a.imag.contains(0) and b.imag.contains(0) and a.real*b.real < 0
    return len(batch)


def check_hardy():
    rows = table()
    jobs = [rows[i:i+24] for i in range(0,len(rows),24)]
    workers = max(1,min(12,int(environ.get("THETA_WORKERS","6"))))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        assert sum(pool.map(hardy,jobs)) == 13946
    print("The 27,892 Hardy endpoint signs certify the 13,946 disjoint positive brackets on the eight primitive characters.")


def check_exterior(job):
    q,k,plus,reflected = job
    ctx.prec = 128
    chi = dirichlet_char(q,k)
    assert chi.is_primitive() and chi.parity() == 1
    assert (chi*dirichlet_char(q,pow(k,-1,q))).is_principal()
    limit = 16384
    mask = bytearray([1])*(limit+1)
    mask[:2] = bytes(2)
    for p in range(2,isqrt(limit)+1):
        if mask[p]:
            mask[p*p::p] = bytes(1+(limit-p*p)//p)
    centres = range(0,295,3)
    sums = [arb(0) for _ in centres]
    for p in range(2,limit+1):
        if not mask[p]:
            continue
        base = chi(p)/(p*arb(p).sqrt())
        if base.is_zero():
            continue
        logp = arb(p).log()
        phase = acb(0,-3*logp).exp()
        n, value, rotation = p, logp*base, phase
        while n <= limit:
            sample = value
            for j in range(len(sums)):
                sums[j] += sample.real
                sample *= rotation
            n *= p
            value *= base
            rotation *= phase
    plus, minus = list(map(ordinate,plus)), list(map(ordinate,reflected))
    error = 2*(arb(limit).log()+2)/arb(limit).sqrt()
    for v, total in zip(centres,sums):
        z = acb(arb(3)/2,v)
        derivative = ((q/arb.pi()).log()+((z+1)/2).digamma().real)/2-total
        accounted = sum((1/(1+(t-v)**2) for t in plus),arb(0))
        accounted += sum((1/(1+(t+v)**2) for t in minus),arb(0))
        assert derivative + error - accounted < arb(4)/13
    z = acb(arb(3)/2)
    jet = acb_series([z,1],2).dirichlet_l(chi)
    total = (((q/arb.pi()).log()+((z+1)/2).digamma())/2+jet[1]/jet[0]).real
    return len(sums), total-sum((1/(1+t*t) for t in plus+minus),arb(0))


def exterior():
    ctx.prec = 128
    groups = defaultdict(list)
    for row in table():
        groups[row[:2]].append(row)
    jobs = [(q,k,rs,groups[q,pow(k,-1,q)]) for (q,k),rs in groups.items()]
    # Arb objects stay local to a process to avoid serializing their ball bounds.
    values = [check_exterior(job) for job in jobs]
    assert sum(n for n,_ in values) == 792
    omitted = sum((m for _,m in values),arb(0))
    assert 0 < omitted < arb(47)/1000
    assert sum(hi < 294*DEN for _,_,_,hi in table()) == 4194
    print("792 exterior Euler sums prove completeness through 294, with 4,194 positive zeros. The total unlisted line mass at zero is below 0.047.")


def terms(method):
    a,r = arb(2).log(),arb(2).sqrt()
    for row in table():
        t = ordinate(row)
        u = a*t/2
        odd = 3-2*r*(2*u).cos() if method == 'direct' else (r-1)**2+4*r*u.sin()**2
        yield t,-2*odd*u.sin()/(a*(r-1)**2*(arb(1)/4+t*t))


def check_positive(method):
    assert method in ('direct','split')
    ctx.prec = 128 if method == 'direct' else 160
    a=arb(2).log()
    carrier, length = arb(8500),1/arb(128)
    barrier = (arb(8590)**2-carrier**2-arb(1)/4).sqrt()
    assert barrier>1200 and 8590/(2*barrier)<8
    high = 30*10**6*(arb(43)/2+arb(1)/48-8492*length).exp()
    trivial = 48*(-3*(arb(21)-2)/2).exp()/(1-(-2*(arb(21)-2)).exp())
    assert high+trivial < arb(1)/1000
    assert 11632*8590 < 10**8
    samples = [tuple(map(int,line.split())) for line in
               (ROOT/'evidence/onset-supports.tsv').read_text().splitlines()]
    assert len(samples) == len(set(samples)) == 22
    modes=list(terms(method))
    windows = {-1:[],1:[]}
    def strict_support(m,den,sign,k):
        assert sign in (-1,1) and k in (96,128,144)
        x,d = arb(m)/den,1/arb(k)
        radius = length+2*d
        assert 21 < x-radius and x+radius < 43
        ylo = x-radius-3*a/2
        yhi = x+radius-3*a/2
        ell = int(ylo.floor().unique_fmpz())
        assert yhi < 42 and 2*arb(2).sqrt()*(yhi/2).exp()<CAP
        threshold = arb(19735)/10000 + correction(ell) + arb(1)/1000
        threshold += 50*(arb(47)/1000)/(294*d)**2
        total = arb(0)
        for t,b in modes:
            p = (carrier**2-t*t).sqrt()
            if method == 'direct':
                sphere = carrier*(length*p).sinh()/(p*(length*carrier).sinh())
                box = (d*t).sinc()**2
                sine = (t*x).sin()
            else:
                sphere = carrier/p*(-length*t*t/(carrier+p)).exp()
                sphere *= (-(-2*length*p).expm1())/(-(-2*length*carrier).expm1())
                box = (1-(2*d*t).cos())/(2*(d*t)**2)
                lo = m//den
                sn,cs = (lo*t).sin_cos()
                sp,cp = ((x-lo)*t).sin_cos()
                sine = sn*cp+cs*sp
            total += b*sphere*box*sine
        assert sign*total > threshold, (m,sign,k,total,threshold)
        return x-radius,x+radius

    for m,sign,k in samples:
        left,right = strict_support(m,512,sign,k)
        windows[sign].append((right-(arb(33)/4).log(),left))
    first = 9*arb(10).log()+3*a/2
    for sign,intervals in windows.items():
        goal = arb(40041 if sign<0 else 41330)/1000
        if method == 'direct':
            frontier = first
            for lo,hi in intervals:
                assert lo < frontier < hi
                frontier = hi
            assert frontier > goal
        else:
            target = goal
            for lo,hi in reversed(intervals):
                if lo < target < hi:
                    target = lo
            assert target < first
    print(f"{method}: the 22 positive averages cover every real start from one billion through centred log height 40.041 at factor 16.5; gains continue through 41.33.")

    # Two more strict loss witnesses bracket the return following the first miss.
    # Evaluate both at their exact centres, then compare the open eligible ranges.
    first = strict_support(164137,4096,-1,144)
    last = strict_support(173371,4096,-1,128)
    narrow=arb(400507)/10000
    assert first[1]-(arb(33)/4).log() < arb(40041)/1000 < narrow < first[0]
    assert last[1]-(arb(399)/40).log()<narrow
    assert arb(4133)/100 < last[0]
    print(f"{method}: two further strict loss supports extend the sharp-factor range through 40.0507 and prove that 19.95H works for both directions through centred log height 41.33.")


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    task = sys.argv[1]
    if task == 'hardy':
        check_hardy()
    elif task == 'exterior':
        exterior()
    else:
        check_positive(task)
