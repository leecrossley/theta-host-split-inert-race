"""Exclude mixed frequencies using independently computed samples and Chebyshev interpolation."""
from fractions import Fraction
from itertools import product
from multiprocessing import Pool
from os import environ
from pathlib import Path
from flint import acb,arb,ctx
from residues import Evaluator,ball

TICKS=10**14
ORDER=26
ROOT=Path(__file__).resolve().parent

def questions(rows,family):
    spans=[(Fraction(l,TICKS),Fraction(h,TICKS)) for _,_,_,l,h in rows]
    for k,(l,h) in enumerate(spans):
        if k in (range(5,9) if family==0 else range(5)):
            yield l,h
    # Start the first nonzero signed coefficient at +1; vary all later ones.
    for k,(l,h) in enumerate(spans):
        for tail in product((0,-1,1),repeat=len(rows)-k-1):
            if not any(tail):
                continue
            a,b=l,h
            for sign,(lo,hi) in zip(tail,spans[k+1:]):
                if sign==1:
                    a,b=a+lo,b+hi
                elif sign==-1:
                    a,b=a-hi,b-lo
            assert a*b>0
            yield (a,b) if a>0 else (-b,-a)

def one(key):
    p,j=key
    ctx.prec=144
    ev=Evaluator(p,j)
    family=0 if p in (547,727) else 1
    rows=[tuple(map(int,line.split())) for line in (ROOT/"evidence"/"brackets.tsv").read_text().splitlines()]
    central=ev.completed(arb(0)).conjugate()
    assert abs(central)>0
    indices=(1,1,2,0,1,1,2,2,1)
    for (cp,_,_,l,h),slot in zip(rows,indices,strict=True):
        if (p,j)==(cp,slot):
            a=ev.completed(arb(l)/TICKS)*central
            b=ev.completed(arb(h)/TICKS)*central
            assert a.imag.contains(0) and b.imag.contains(0) and a.real*b.real<0
    # On each cell of length 1/4, approximate at the 26 Chebyshev roots.
    cosine=[[((2*m+1)*k*arb.pi()/(2*ORDER)).cos() for m in range(ORDER)] for k in range(ORDER)]
    polys={}
    radius=arb(1)/8
    # Lagrange interpolation error componentwise, times 2 for the complex norm.
    error=96*ev.q*arb(4)**(-ORDER)
    checked=0
    for l,h in questions(rows,family):
        cell=int(2*(l+h))
        center=(2*cell+1)*radius
        if cell not in polys:
            vals=[ev.value(center+radius*x) for x in cosine[1]]
            coeff=[sum((z*c for z,c in zip(vals,row)),acb(0))*(1 if k==0 else 2)/ORDER
                   for k,row in enumerate(cosine)]
            polys[cell]=coeff
        # Intersections with a cell edge require an interpolation estimate on
        # its own cell. All measured intervals fit in the selected closed cell.
        assert Fraction(cell,4)<=l<=h<=Fraction(cell+1,4)
        z=(arb(ball((l+h)/2),ball((h-l)/2))-center)/radius
        b1=b2=acb(0)
        coeff=polys[cell]
        for a in reversed(coeff[1:]):
            b1,b2=2*z*b1-b2+a,b1
        val=z*b1-b2+coeff[0]
        assert abs(val)>error,(p,j,l,h)
        checked+=1
    return checked

def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    with Pool(min(8,int(environ.get("THETA_WORKERS","4")))) as pool:
        total=sum(pool.map(one,[(p,j) for p in (547,727,571,599) for j in range(4)]))
    assert total==157384
    print("Residue sums and Chebyshev interpolation independently verify all 157,384 exclusions and nine sign changes.")

if __name__=="__main__":
    main()
