"""Certify the complete interiors of the two long excursion gaps by linear interpolation."""
from pathlib import Path
import sys
from flint import arb, ctx
from band_values import values

# Ends lie inside the crossing neighborhoods checked directly from the prime knots.
SLABS = ((1, 6288755, 20306405, 4096), (-1, 20077514, 61070272, 8192),
         (-1, 562485616, 1701279810, 131072))


def gap_error(a, b):
    """Uniform secant error for sqrt(y) times either signed move on [a,b]."""
    h = b-a
    root = arb(2).sqrt()
    norm = 8/(2*root*(root-1)**2*arb(2).log())
    slope_counts = sum((M*(j*h+1)/(j*arb(a)) for j,M in
                        ((8,8/root), (4,12/root+4), (2,4/root+6), (1,arb(2)))),arb(0))
    return h*norm*arb(8*b).log()*slope_counts/4


def base_cells():
    return [(s,x,min(x+step,h)) for s,l,h,step in SLABS for x in range(l,h,step)]


def certified(s,a,b,scores):
    threshold=arb(19735)/10000
    left=(threshold-s*scores[a])*arb(a).sqrt()
    right=(threshold-s*scores[b])*arb(b).sqrt()
    assert left>0 and right>0, (s,a,b,left,right)
    error=gap_error(a,b)
    return left>error and right>error


def discover(binary, path):
    cells=base_cells()
    scores={}
    splits=[]
    while cells:
        fresh={k for _,a,b in cells for k in (a,b) if k not in scores}
        scores.update(values(binary,fresh,1))
        work=[]
        for s,a,b in cells:
            if not certified(s,a,b,scores):
                m=(a+b)//2
                assert a<m<b
                splits.append(m)
                work.extend(((s,a,m),(s,m,b)))
        cells=work
    path.write_text("".join(str(k)+"\n" for k in sorted(set(splits))))


def verify(binary,path):
    splits=list(map(int,path.read_text().splitlines()))
    assert splits==sorted(set(splits))
    cells=[]
    for s,l,h,step in SLABS:
        nodes=sorted({l,h,*range(l,h,step),*(k for k in splits if l<k<h)})
        cells.extend((s,a,b) for a,b in zip(nodes,nodes[1:]))
    assert all(any(l<k<h for _,l,h,_ in SLABS) for k in splits)
    scores=values(binary,{k for _,a,b in cells for k in (a,b)},1)
    for s,a,b in cells:
        assert certified(s,a,b,scores), (s,a,b)
    for sign,l,h,_ in SLABS:
        count=sum(s==sign and l<=a<b<=h for s,a,b in cells)
        label="gain" if sign==1 else "loss"
        print(f"{count} interpolation intervals exclude a strict {label} from {l} through {h}")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=144
    binary=Path(sys.argv[1]).resolve()
    path=Path(__file__).parent/"evidence"/"late-subdivisions.txt"
    if len(sys.argv)==3 and sys.argv[2]=="--write-mesh":
        discover(binary,path)
    elif len(sys.argv)!=2:
        raise SystemExit("Usage: python late_slabs.py SIEVE [--write-mesh]")
    verify(binary,path)
