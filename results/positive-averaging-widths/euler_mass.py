"""Independently count the low zeros from a finite Euler sum at real part two."""
from math import gcd, isqrt
from flint import arb,ctx
from band_zeros import QUARTIC, ball, character_family, rows


def primes_through(n):
    keep = bytearray([1])*(n+1)
    for p in range(2,isqrt(n)+1):
        if keep[p]:
            keep[p*p::p] = bytes(len(keep[p*p::p]))
    return [p for p in range(2,n+1) if keep[p]]


def check():
    ctx.prec = 144
    stop = 4096
    primes = primes_through(stop)
    logs = {p:arb(p).log() for p in primes}
    tail = (arb(stop).log()+1)/stop
    disp = arb(3)/2
    _,by = rows()
    for (q,label),(_,p,j) in character_family().items():
        total = ((q/arb.pi()).log()+2-arb.const_euler()-2*arb(2).log())/2
        for v in primes:
            if gcd(v,q)>1:
                continue
            legendre = 0 if pow(v,(p-1)//2,p)==1 else 2
            phase = legendre + (j*QUARTIC[v % 16] if j else 0)
            power,n = 1,v
            while n <= stop:
                total -= (1,0,-1,0)[(power*phase) % 4]*logs[v]/(n*n)
                n *= v
                power += 1
        retained = sum(disp/(disp**2+ball(pair)**2)
                       for n in (label,pow(label,-1,q)) for pair in by[q,n])
        assert 0 < total-tail-retained and total+tail-retained < arb(1)/4
    print('At real part two, independent finite Euler sums through 4,096, including an absolute tail, place each omitted positive logarithmic-derivative mass below 1/4.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    check()
