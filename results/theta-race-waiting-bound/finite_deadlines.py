"""Verify signed averages and cover every real height in each claimed interval."""
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from os import environ
from sys import argv

from flint import arb, ctx
from finite_spectrum import ROOT, brackets, ordinate

PARAMS = {1000:(560,16,256,526),2000:(1070,20,400,1017)}
STRATA = (0,23,30,40,55,80,120,200)
OMITTED = (2,547,727,343057911569,343070636977)
spectra = {}


def witnesses(top):
    data = [tuple(map(int,line.split())) for line in
            (ROOT/f'evidence/finite-{top}.tsv').read_text().splitlines()]
    assert data==sorted(set(data)) and len(data)==PARAMS[top][3]
    assert all(0<m<100*(top+9) and sign in (-1,1) and tier in STRATA for m,sign,tier in data)
    return data


def logfloor(tier):
    return 9*arb(10).log() if tier==0 else arb(tier)


@lru_cache(maxsize=None)
def arithmetic(tier):
    ell = logfloor(tier)
    a,r = arb(2).log(),arb(2).sqrt()
    B = 2*r*(r-1)**2
    Z = (ell/2).exp()
    if tier<50:
        widths = []
        for j in range(32):
            v = arb(j)/32
            u = 2/(1+v).sqrt()-(2/(2-v)).sqrt()
            assert u*Z>32
            widths.append(u*(2*Z).log()/(u*Z/16).log())
        square = 4*sum(widths)/32/B
    else:
        assert Z>10**10
        square = 2+8*arb(24)/1000/B
    square += arb(816)/100*(2/Z).sqrt()/B
    l4 = arb(4).log()+ell
    high = arb(816)/100/B * (sum((l4/k-ell/2).exp() for k in (3,4))+
                               l4/a*(l4/5-ell/2).exp())
    ram = 8*(3*a+sum(arb(p).log() for p in OMITTED[1:]))/B/Z
    return (square+2*high+2*ram)/a


def uniform(top):
    n,den,C,_ = PARAMS[top]
    h,d = 1/arb(3000),1/arb(den)
    R = n*h+2*d
    line = 50*arb(127)/1000/(48*d)**2
    a = arb(2).log()
    high = 10**6*((top+arb(9))/2+d).exp()*((h/2).cosh()/(8000*h))**n
    xmin = 9*arb(10).log()+3*a/2
    trivial = 48*(-3*(xmin-2)/2).exp()/(1-(-2*(xmin-2)).exp())
    assert high+trivial < arb(2)/1000
    return R,line+arb(2)/1000


def geometry(top,data):
    R,_ = uniform(top)
    offset = 3*arb(2).log()/2
    size = (arb(PARAMS[top][2])/2).log()
    for sign in (-1,1):
        frontier = 9*arb(10).log()+offset
        for m,s,tier in data:
            if s!=sign:
                continue
            x = arb(m)/100
            assert x-R-offset>logfloor(tier)
            assert x+R-size < frontier < x-R
            frontier = x-R
        assert frontier>top+offset


def prepare(top,method):
    ctx.prec = 128 if method=='direct' else 144
    n,den,_,_ = PARAMS[top]
    a,r = arb(2).log(),arb(2).sqrt()
    c = 2*(r-1)**2
    vals = []
    for row in brackets():
        t = ordinate(row)
        u = a*t/2
        if method=='direct':
            G = (3-2*r*(2*u).cos())/(c*(t*t+arb(1)/4))
            W = (t/3000).sinc()**n*(t/den).sinc()**2
        else:
            # n is even. A two-box response is (1-cos(2v))/(2v^2).
            G = ((r-1)**2+4*r*u.sin()**2)/(c*(t*t+arb(1)/4))
            h,d = t/3000,t/den
            W = ((1-(2*h).cos())/(2*h*h))**(n//2)*(1-(2*d).cos())/(2*d*d)
        vals.append((t,-4*G*u.sin()*W/a))
    spectra[top,method] = vals


def signed_values(job):
    top,method,rows = job
    if (top,method) not in spectra:
        prepare(top,method)
    _,rest = uniform(top)
    vals = spectra[top,method]
    for m,sign,tier in rows:
        x = arb(m)/100
        if method=='direct':
            result = sum(b*(t*x).sin() for t,b in vals)
        else:
            # Split the phase before multiplication, independently of the direct sine.
            hi,lo = divmod(m,100)
            a,b = arb(hi),arb(lo)/100
            result = arb(0)
            for t,v in vals:
                sa,ca = (t*a).sin_cos()
                sb,cb = (t*b).sin_cos()
                result += v*(sa*cb+ca*sb)
        assert sign*result-rest-arithmetic(tier)>arb(19735)/10000,(top,method,m,sign,result)
    return len(rows)


def main(method):
    assert method in ('direct','split')
    ctx.prec = 144
    for top,(_,_,C,count) in PARAMS.items():
        data = witnesses(top)
        geometry(top,data)
        jobs = [(top,method,data[i:i+16]) for i in range(0,count,16)]
        with ProcessPoolExecutor(max_workers=max(1,min(12,int(environ.get('THETA_WORKERS','6'))))) as pool:
            assert sum(pool.map(signed_values,jobs))==count
        print(f'{method}: {count} signed positive averages cover the full real range from 10^9 through exp({top}) with completion factor {C}.',flush=True)


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    main(argv[1])
