"""Certify the entire unit-height band of the sixteen primitive factors."""
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from math import gcd
from os import environ
from pathlib import Path

from flint import acb, acb_series, arb, ctx, dirichlet_char

ROOT = Path(__file__).resolve().parent
TICKS = 10**14
QUARTIC = {sign*pow(5,k,16) % 16:k for sign in (1,-1) for k in range(4)}


def character_family():
    family = {}
    for field, pair in enumerate(((547,727),(571,599))):
        for p in pair:
            for mult, residue in ((1,1),(8,5),(16,5),(16,13)):
                q = p*mult
                label = next(n for n in range(p-1,q,p) if n % mult == residue % mult)
                chi = dirichlet_char(q,label)
                assert chi.is_primitive() and chi.parity() == 1 and chi.conductor() == q
                exponent = int(chi.group().exponent())
                local_five = next(n for n in range(1,16*p,p) if n % 16 == 5)
                num = 4*int(chi.chi_exponent(local_five))
                assert num % exponent == 0
                power = (num//exponent) % 4
                assert power in {1:(0,),8:(2,),16:(1,3)}[mult]
                for n in range(1,16*p,2):
                    if n % p == 0:
                        continue
                    legendre = 0 if pow(n,(p-1)//2,p) == 1 else 2
                    expected = legendre + power*QUARTIC[n % 16]
                    assert (4*int(chi.chi_exponent(n))-exponent*expected) % (4*exponent) == 0
                family[q,label] = field,p,power
    return family


def rows():
    raw = [tuple(map(int,line.split())) for line in
           (ROOT/'evidence/band-brackets.tsv').read_text().splitlines()]
    assert raw == sorted(raw) and len(raw) == 810
    by = defaultdict(list)
    for f,q,l,a,b in raw:
        assert f in (0,1) and 0 < a < b < 35*TICKS
        assert not by[q,l] or by[q,l][-1][1] < a
        by[q,l].append((a,b))
    return raw,by


def ball(pair):
    a,b = pair
    return arb(a+b,b-a)/(2*TICKS)


def signs(job):
    q,label,pairs = job
    ctx.prec = 112
    chi = dirichlet_char(q,label)
    for a,b in pairs:
        left,right = (chi.hardy_z(arb(n)/TICKS) for n in (a,b))
        assert left.imag.contains(0) and right.imag.contains(0) and left.real*right.real < 0
    return len(pairs)


def check():
    family = character_family()
    raw,by = rows()
    assert set(by) == set(family)
    assert all(f == family[q,l][0] for f,q,l,_,_ in raw)
    jobs = [(q,l,items[i:i+12]) for (q,l),items in by.items()
            for i in range(0,len(items),12)]
    with ProcessPoolExecutor(max_workers=max(1,min(8,int(environ.get('THETA_WORKERS','4'))))) as pool:
        assert sum(pool.map(signs,jobs)) == len(raw)
    ctx.prec = 160
    z = acb(arb(3)/2)
    for q,label in family:
        chi = dirichlet_char(q,label)
        jet = acb_series([z,1],2).dirichlet_l(chi)
        total = (((q/arb.pi()).log()+((z+1)/2).digamma())/2 + jet[1]/jet[0]).real
        used = sum(1/(1+ball(pair)**2) for n in (label,pow(label,-1,q)) for pair in by[q,n])
        remainder = total-used
        assert 0 < remainder < arb(1)/5
    band = sorted((a,b,f) for f,_,_,a,b in raw if b <= TICKS)
    assert all(a>TICKS or b<TICKS for _,_,_,a,b in raw)
    assert all(first[1] < second[0] for first,second in zip(band,band[1:]))
    assert [sum(f==j for _,_,f in band) for j in (0,1)] == [9,7]
    print('Verified 810 positive Hardy brackets (1,620 strict endpoint signs) and the sixteen primitive field characters.')
    print('All sixteen completed logarithmic derivatives at 3/2 leave positive remainders below 1/5. The closed unit-height band has nine and seven simple positive zeros, mutually disjoint, their negative partners, and no real zeros.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    check()
