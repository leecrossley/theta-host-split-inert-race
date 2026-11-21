"""Certify both dyadic moves for each real starting height up to one billion."""
from bisect import bisect_right
from fractions import Fraction
from pathlib import Path
import argparse
import subprocess

from flint import arb, ctx, fmpz

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")

def log_table(precision=40):
    bins = 2**15
    table = []
    for j in range(bins+1):
        v = ((1+arb(j)/bins).log()*2**precision).floor().unique_fmpz()
        assert v is not None
        table.append(int(v))
    return table

def sieve(binary, edges, *, fine=False):
    numbers = [*log_table(58 if fine else 40), *edges]
    feed = "".join(f"{n}\n" for n in numbers)
    out = subprocess.run([binary, "--fine"] if fine else [binary], input=feed, text=True, capture_output=True, check=True)
    blocks = {}
    prefix = [0]*4
    for edge, line in zip(edges, out.stdout.splitlines(), strict=True):
        e, lo, hi, lr, hr = map(int,line.split())
        assert e == edge and lo <= hi and lr <= hr
        vals = (lo,hi,lr,hr)
        prefix = [v+d for v,d in zip(prefix,vals)]
        blocks[edge] = prefix
    return blocks

def box(lo,hi,bits):
    return arb(lo+hi,hi-lo)/2**(bits+1)

def score(y,blocks):
    def band(u,v):
        lo,hi,lr,hr = [b-a for a,b in zip(blocks[u],blocks[v])]
        return box(lo,hi,40),box(lr,hr,72)
    a,b = band(y,2*y)
    c,d = band(2*y,4*y)
    sqrt2=arb(2).sqrt()
    return 8*(2*a-2*y*b+4*y*d-c)/(2*sqrt2*(sqrt2-1)**2*arb(y).sqrt())

def cover(lower,upper,gains,losses):
    worst = (Fraction(0),0,0,0)
    for h in sorted({lower,upper,*(v for v in gains+losses if lower<=v<=upper)}):
        g=gains[bisect_right(gains,h)]
        l=losses[bisect_right(losses,h)]
        ratio=Fraction(2*max(g,l),h)
        if ratio>worst[0]:
            worst=ratio,h,g,l
    return worst

HEIGHT=1_000_000_000
TABLE=Path(__file__).resolve().parent/"evidence"/"windows.tsv"


def sample_grid():
    grid=[1]
    while grid[-1] < 4*HEIGHT:
        grid.append((33*grid[-1]+31)//32)
    return grid


def check_small(binary):
    # Independent primality tests and Euler's criterion; no congruence-class sieve here.
    edges=(3,4,16,31,547,727,5000,50000,250000)
    result=sieve(binary,edges)
    rows=[]
    for p in range(3,edges[-1]+1,2):
        if p in (547,727) or not fmpz(p).is_prime(): continue
        u,v=(1 if pow(p,(q-1)//2,q)==1 else -1 for q in (547,727))
        f=next(k for k in (1,2,4) if (u*v)**k==1 and pow(p,k,16) in (1,15))
        if p**f<=edges[-1]:
            rows.append((p**f,-(u**f)*arb(p).log()))
    for e in edges:
        a=sum((z for n,z in rows if n<=e),arb(0))
        b=sum((z/n for n,z in rows if n<=e),arb(0))
        lo,hi,lr,hr=result[e]
        assert arb(lo)<=a*2**40<=arb(hi) and arb(lr)<=b*2**72<=arb(hr)
    print("independent primality tests and direct Arb logarithms agree through norm 250000")


def generate(binary):
    ctx.prec=144
    check_small(binary)
    grid=sample_grid()
    sample=sorted(set(grid)|{2*y for y in grid})
    edges=sorted({k*y for y in sample for k in (1,2,4)})
    assert edges[-1]==32467011360
    blocks=sieve(binary,edges)
    scores={y:score(y,blocks) for y in sample}
    lines=["y\tslope_lower_times_10000\tslope_upper_times_10000"]
    for y in grid:
        slope=(scores[2*y]-scores[y])*10000/arb(2).log()
        lo=int(slope.lower().floor().unique_fmpz())
        hi=int(slope.upper().ceil().unique_fmpz())
        lines.append(f"{y}\t{lo}\t{hi}")
    return "\n".join(lines)+"\n"


def verify_coverage(data):
    rows=[tuple(map(int,line.split())) for line in data.splitlines()[1:]]
    assert [y for y,_,_ in rows]==sample_grid() and all(lo<=hi for _,lo,hi in rows)
    gains=[y for y,lo,_ in rows if lo>19735]
    losses=[y for y,_,hi in rows if hi< -19735]
    assert len(rows)==625
    for lower,cap in ((1,Fraction(18)),(10**6,Fraction(7))):
        ratio,h,g,l=cover(lower,HEIGHT,gains,losses)
        assert ratio<=cap
        print(f"H in [{lower},{HEIGHT}]: both endpoints <= ({ratio})H; worst grid start {h}")
    g,l=gains[bisect_right(gains,HEIGHT)],losses[bisect_right(losses,HEIGHT)]
    print(f"above one billion: gain {g} to {2*g}; loss {l} to {2*l}")
    print("625 integer slope enclosures; largest prime norm sieved 32467011360")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary",type=Path,nargs="?",help="compiled sieve executable")
    parser.add_argument("--write-table",type=Path,help="write regenerated certificate to a file")
    args=parser.parse_args()
    if not args.binary:
        if args.write_table: parser.error("--write-table requires the sieve executable")
        verify_coverage(TABLE.read_text())
        return
    data=generate(args.binary.resolve())
    if args.write_table:
        args.write_table.write_text(data)
    elif data!=TABLE.read_text():
        raise ArithmeticError("Recomputed table differs from evidence/windows.tsv")
    verify_coverage(data)


if __name__=="__main__":
    main()
