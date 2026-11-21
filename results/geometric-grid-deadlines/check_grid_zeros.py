"""Verify every sign-changing bracket in the two positive-ordinate tables."""
from collections import defaultdict
from multiprocessing import Pool
from os import environ
from pathlib import Path
from flint import arb,ctx,dirichlet_char

ROOT=Path(__file__).resolve().parent
TICKS=10**14


def verify(batch):
    ctx.prec=100
    q,l,rows=batch
    ch=dirichlet_char(q,l)
    assert ch.is_primitive() and ch.parity()==1 and ch.conductor()==q
    for a,b in rows:
        x,y=(ch.hardy_z(arb(t)/TICKS) for t in (a,b))
        assert x.imag.contains(0) and y.imag.contains(0) and x.real*y.real<0,(q,l,a,b)
    return len(rows)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    expected={}
    for i,pair in enumerate(((547,727),(571,599))):
        for p in pair:
            for d,x in ((1,1),(8,5),(16,5),(16,13)):
                q=p*d
                l=next(n for n in range(p-1,q,p) if n%d==x%d)
                expected[q,l]=i
    groups=defaultdict(list)
    for line in (ROOT/'evidence/dyadic-zeros.tsv').read_text().splitlines():
        i,q,l,a,b=map(int,line.split())
        assert i==expected[q,l] and 0<a<b
        groups[q,l].append((a,b))
    assert groups.keys()==expected.keys()
    jobs=[]
    for (q,l),rows in sorted(groups.items()):
        rows.sort()
        assert all(x[1]<y[0] for x,y in zip(rows,rows[1:]))
        jobs.extend((q,l,rows[n:n+24]) for n in range(0,len(rows),24))
    with Pool(max(1,min(8,int(environ.get('THETA_WORKERS','4'))))) as pool:
        count=sum(pool.map(verify,jobs))
    print(f"{count:,} positive disjoint-within-factor zero brackets verified by {2*count:,} strict Hardy endpoint signs on all sixteen primitive factors.")


if __name__=='__main__':
    main()
