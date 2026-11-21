"""Check two actual low spectra and bound every omitted critical-line zero."""
from concurrent.futures import ProcessPoolExecutor
from os import environ
from pathlib import Path

from flint import acb, acb_series, arb, ctx, dirichlet_char

ROOT=Path(__file__).resolve().parent
SCALE=10**14
PAIRS=((547,727),(571,599))


def tables():
    rows=[tuple(map(int,line.split())) for line in (ROOT/'evidence/common-deadline-zeros.tsv').read_text().splitlines()]
    assert rows==sorted(rows) and len(rows)==len(set(rows))
    expected={}
    for i,primes in enumerate(PAIRS):
        for p in primes:
            for m,idx in ((1,1),(8,5),(16,5),(16,13)):
                q=m*p
                label=next(n for n in range(p-1,q,p) if n%m==idx%m)
                expected[q,label]=i
    assert {row[1:3] for row in rows}==set(expected)
    last={}
    for i,q,k,lo,hi in rows:
        assert expected[q,k]==i and last.get((q,k),0)<lo<hi<125*SCALE
        last[q,k]=hi
    return rows


def ordinate(row):
    a,b=row[-2:]
    return arb(a+b,b-a)/(2*SCALE)


def signs(rows):
    ctx.prec=112
    chars={}
    for _,q,k,a,b in rows:
        if (q,k) not in chars:
            chi=dirichlet_char(q,k)
            assert chi.is_primitive() and chi.parity()==1 and chi.conductor()==q
            chars[q,k]=chi
        left,right=(chars[q,k].hardy_z(arb(n)/SCALE) for n in (a,b))
        assert left.imag.contains(0) and right.imag.contains(0) and left.real*right.real < 0
    return len(rows)


def masses(rows):
    ctx.prec=144
    for i in (0,1):
        part=[row for row in rows if row[0]==i]
        zeros=list(map(ordinate,part))
        chars=sorted({row[1:3] for row in part})
        caps=[]
        for v in range(0,49,2):
            s=acb(arb(3)/2,v)
            total=arb(0)
            for q,k in chars:
                jet=acb_series([s,1],2).dirichlet_l(dirichlet_char(q,k))
                total+=(((q/arb.pi()).log()+((s+1)/2).digamma())/2+jet[1]/jet[0]).real
            omitted=total-sum(1/(1+(v-t)**2)+1/(1+(v+t)**2) for t in zeros)
            assert omitted>0 and omitted<arb(49)/100,(i,v,omitted)
            caps.append(omitted)
        assert caps[0]<arb(30)/100, (i,caps[0])
        print(f'{PAIRS[i][0]}/{PAIRS[i][1]}: {len(part)} positive brackets. The line is complete through 48 and its omitted mass is below 0.30 at zero.')


def main():
    rows=tables()
    with ProcessPoolExecutor(max_workers=max(1,min(12,int(environ.get('THETA_WORKERS','6'))))) as pool:
        assert sum(pool.map(signs,(rows[i:i+24] for i in range(0,len(rows),24))))==len(rows)
    print(f'{len(rows)} genuine brackets, {2*len(rows)} strict Hardy signs, sixteen primitive characters.')
    masses(rows)


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    main()
