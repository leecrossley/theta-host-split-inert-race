"""Certify common positive averages and continuous coverage of all real starts."""
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from os import environ
from sys import argv

from flint import arb, ctx
from deadline_arithmetic import TIERS, costs, logfloor
from deadline_spectrum import ROOT, ordinate, tables

ORDERS={2:(18,860,15),4:(80,695,0)}
spectra={}


def geometry(dim):
    bits,count,begin=ORDERS[dim]
    data=[tuple(map(int,l.split())) for l in (ROOT/f'evidence/common-{dim}.tsv').read_text().splitlines()]
    assert len(data)==count and data==sorted(set(data))
    a=arb(2).log()
    radius=arb(1150)/3000+arb(2)/8
    size=(bits-1)*a
    for code in range(2**dim):
        frontier=arb(begin)+3*a/2
        for m,s,tier in data:
            assert 0<m<=208000 and 0<=s<2**dim and tier in TIERS
            if s!=code:
                continue
            x=arb(m)/100
            assert x-radius-3*a/2>logfloor(tier)
            assert x+radius-size<frontier<x-radius
            frontier=x-radius
        assert frontier>2000+3*a/2
    return data


def errors():
    a,r=arb(2).log(),arb(2).sqrt()
    h=1/arb(3000)
    high=10**6*(1040+arb(1)/8).exp()*((h/2).cosh()/(8000*h))**1150
    xmin=9*arb(10).log()+3*a/2
    trivial=48*(-3*(xmin-2)/2).exp()/(1-(-2*(xmin-2)).exp())
    assert high+trivial<arb(2)/1000
    assert (3+2*r)/(2*(r-1)**2)*(1+arb(48)**2)/(arb(1)/4+48**2)<18
    assert 2/a*(3+2*r)/(2*(r-1)**2)*(1+arb(48)**2)/(arb(1)/4+48**2)<50
    return tuple(c*arb(3)/10/(48/arb(8))**2+arb(2)/1000 for c in (18,50))


def prepare(method):
    ctx.prec=128 if method=='direct' else 144
    a,r=arb(2).log(),arb(2).sqrt()
    modes=[[],[]]
    for row in tables():
        t=ordinate(row)
        u=a*t/2
        if method=='direct':
            g=(3-2*r*(2*u).cos())/(2*(r-1)**2*(arb(1)/4+t*t))
            window=(t/3000).sinc()**1150*(t/8).sinc()**2
        else:
            g=((r-1)**2+4*r*u.sin()**2)/(2*(r-1)**2*(arb(1)/4+t*t))
            h,d=t/3000,t/8
            window=((1-(2*h).cos())/(2*h*h))**575*(1-(2*d).cos())/(2*d*d)
        modes[row[0]].append((t,2*g*window,-4*g*u.sin()*window/a))
    spectra[method]=modes


def checked(job):
    dim,method,rows=job
    if method not in spectra:
        prepare(method)
    residual=errors()
    signs=list(product((-1,1),repeat=dim))
    for m,code,tier in rows:
        x=arb(m)/100
        left,right=divmod(m,100)
        for field,modes in enumerate(spectra[method]):
            level=move=arb(0)
            for t,l,v in modes:
                if method=='direct':
                    sn,cs=(t*x).sin_cos()
                else:
                    sa,ca=(t*left).sin_cos()
                    sb,cb=(t*right/100).sin_cos()
                    sn,cs=sa*cb+ca*sb,ca*cb-sa*sb
                if dim==4:
                    level+=l*cs
                move+=v*sn
            el,ed=costs(field,tier)
            if dim==4:
                assert signs[code][2*field]*level-el-residual[0]>2,(dim,m,field,level)
            direction=signs[code][field if dim==2 else 2*field+1]
            assert direction*move-ed-residual[1]>4,(dim,m,field,move)
    return len(rows)


def main(method):
    assert method in ('direct','split')
    ctx.prec=144
    for dim,(bits,count,_) in ORDERS.items():
        data=geometry(dim)
        jobs=[(dim,method,data[i:i+16]) for i in range(0,count,16)]
        with ProcessPoolExecutor(max_workers=max(1,min(12,int(environ.get('THETA_WORKERS','6'))))) as pool:
            assert sum(pool.map(checked,jobs))==count
        print(f'{method}: {count} common positive averages cover all {2**dim} sign choices through exp(2000) with factor 2^{bits}.',flush=True)


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    main(argv[1])
