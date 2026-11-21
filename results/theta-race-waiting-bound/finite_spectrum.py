"""Count the low ordinates and bound the retained and unretained zero masses."""
from concurrent.futures import ProcessPoolExecutor
from os import environ
from pathlib import Path

from flint import acb, acb_series, arb, ctx, dirichlet_char
from verify import odd_factors

ROOT = Path(__file__).resolve().parent
SCALE = 10**14


def brackets():
    table = [tuple(map(int,l.split())) for l in
             (ROOT/'evidence/finite-zeros.tsv').read_text().splitlines()]
    assert len(table)==6000 and table==sorted(table)
    family = {(q,k) for q,k,_ in odd_factors()}
    assert {row[:2] for row in table}==family
    last = {}
    for q,k,lo,hi in table:
        assert last.get((q,k),0)<lo<hi
        last[q,k] = hi
    return table


def ordinate(row):
    _,_,lo,hi = row
    return arb(lo+hi,hi-lo)/(2*SCALE)


def hardy_check(batch):
    ctx.prec = 112
    chars = {}
    for q,k,lo,hi in batch:
        if (q,k) not in chars:
            chars[q,k] = dirichlet_char(q,k)
        chi = chars[q,k]
        left,right = (chi.hardy_z(arb(e)/SCALE) for e in (lo,hi))
        assert left.imag.contains(0) and right.imag.contains(0) and left.real*right.real<0
    return len(batch)


def check():
    ctx.prec = 144
    data = brackets()
    jobs=[data[i:i+24] for i in range(0,len(data),24)]
    with ProcessPoolExecutor(max_workers=max(1,min(12,int(environ.get('THETA_WORKERS','6'))))) as pool:
        assert sum(pool.map(hardy_check,jobs))==len(data)
    ts = list(map(ordinate,data))
    for v in range(0,49,3):
        z = acb(arb(3)/2,v)
        total = arb(0)
        for q,_,chi in odd_factors():
            jet = acb_series([z,1],2).dirichlet_l(chi)
            total += (((q/arb.pi()).log()+((z+1)/2).digamma())/2+jet[1]/jet[0]).real
        used = sum(1/(1+(t-v)**2)+1/(1+(t+v)**2) for t in ts)
        assert 0<total-used<arb(1)/5
        if v==0:
            assert total-used<arb(127)/1000
    a,r = arb(2).log(),arb(2).sqrt()
    assert arb(1)/5 < 1/(1+(arb(3)/2)**2)
    assert 2*(3+2*r)/(a*2*(r-1)**2)*(1+arb(48)**2)/(arb(1)/4+48**2) < 50
    assert 11632*8000<10**8
    print('Verified 6,000 true positive zero brackets and 12,000 strict Hardy signs on the eight primitive odd factors.')
    print('Seventeen exterior logarithmic derivatives have residuals below 1/5. The line is complete through 48 and the unlisted line mass at the origin is below 127/1000.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    check()
