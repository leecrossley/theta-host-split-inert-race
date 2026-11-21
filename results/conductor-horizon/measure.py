"""Recover both omitted conductor masses from six independently enclosed finite scores."""
from argparse import ArgumentParser
from fractions import Fraction
from pathlib import Path
import subprocess

from flint import arb, ctx, fmpz

ROOT = Path(__file__).resolve().parent
PRIMES = (727, 547)
NORMS = (727**2, 547**4)
TICKS = tuple(k for n in NORMS for k in (n-1, n, n+1))
EDGES = sorted({t//d for t in TICKS for d in (4, 2, 1)})
BINS = 1 << 16


def boundaries(binary, edges):
    table = []
    for j in range(BINS+1):
        cell = ((1+arb(j)/BINS).log()*2**80).floor().unique_fmpz()
        assert cell is not None
        table.append(str(cell))
    data = "\n".join([*table, *map(str, edges)]) + "\n"
    result = subprocess.run([binary], input=data, text=True, capture_output=True, check=True)
    return result.stdout


def parse(text, edges):
    rows = {}
    for edge, line in zip(edges, text.splitlines(), strict=True):
        k, a, b, c, d = map(int, line.split())
        assert k == edge and a <= b and c <= d
        rows[k] = a,b,c,d
    return rows


def interval(a,b,precision):
    return arb(a+b,b-a)/2**(precision+1)


def band(rows, l, u):
    selected = [row for k,row in rows.items() if l < k <= u]
    a,b,c,d = (sum(row[i] for row in selected) for i in range(4))
    return interval(a,b,80),interval(c,d,112)


def raw(rows,k):
    a,b=band(rows,k//4,k//2)
    c,d=band(rows,k//2,k)
    return 8*(2*a-arb(k)*b/2+k*d-c)


def check_small(binary):
    edges=(3,4,16,31,547,727,5000,50000,250000)
    rows=parse(boundaries(binary,edges),edges)
    sums=[[arb(0),arb(0)] for _ in edges]
    for p in range(3,edges[-1]+1,2):
        if p in (547,727) or not fmpz(p).is_prime():
            continue
        s,t=(int(fmpz(p).jacobi(q)) for q in (547,727))
        f=next(i for i in (1,2,4) if pow(p,i,16) in (1,15) and (s*t)**i==1)
        n=p**f
        if n>edges[-1]:
            continue
        j=next(j for j,k in enumerate(edges) if n<=k)
        z=-(s**f)*arb(p).log()
        sums[j][0]+=z
        sums[j][1]+=z/n
    for k,(a,b) in zip(edges,sums):
        lo,hi,lr,hr=rows[k]
        assert interval(lo,hi,80).contains(a) and interval(lr,hr,112).contains(b)


def samples(text):
    rows=parse(text,EDGES)
    actual=[raw(rows,k) for k in TICKS]
    output=[]
    for k,v in zip(TICKS,actual):
        lo=int((v*10**12).lower().floor().unique_fmpz())
        hi=int((v*10**12).upper().ceil().unique_fmpz())
        assert lo < hi and hi-lo <= 5
        output.append(f"{k}\t{lo}\t{hi}\n")
    return "".join(output)


def verify(table):
    lines=[tuple(map(int,row.split())) for row in table.splitlines()]
    assert [k for k,_,_ in lines]==list(TICKS)
    for i,(p,n) in enumerate(zip(PRIMES,NORMS)):
        values=[arb(lo+hi,hi-lo)/(2*10**12) for _,lo,hi in lines[3*i:3*i+3]]
        mass=n*(values[0]-2*values[1]+values[2])/arb(p).log()
        # The unramified alternative is -8 in either of these residue degrees.
        assert mass.contains(0) and abs(mass)<1
        assert all((hi-lo)*n < 2*10**12*arb(p).log() for _,lo,hi in lines[3*i:3*i+3])
    print("Six actual raw scores certify absence at norms 528529 and 89526025681.")
    print("Each measured mass divided by log p is enclosed strictly between -1 and 1. An unramified prime would give -8.")


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("binary",type=Path,nargs="?")
    ap.add_argument("--generate",action="store_true")
    args=ap.parse_args()
    ctx.prec=224
    saved=ROOT/"evidence"/"bands.tsv"
    obs=ROOT/"evidence"/"raw-scores.tsv"
    if args.binary:
        exe=str(args.binary.resolve())
        check_small(exe)
        data=boundaries(exe,EDGES)
        if args.generate:
            saved.write_text(data)
            obs.write_text(samples(data))
        assert data==saved.read_text()
        print("Segmented sieve reproduced every integer enclosure through norm 89526025682.")
        print("Independent primality tests and direct Arb logarithms agree by band through norm 250000.")
    elif args.generate:
        ap.error("--generate requires the compiled sieve")
    text=samples(saved.read_text())
    assert text==obs.read_text()
    verify(text)


if __name__=="__main__":
    main()
