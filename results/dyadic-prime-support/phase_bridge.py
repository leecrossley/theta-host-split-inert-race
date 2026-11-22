"""Certify visibility at arbitrary phases, including a possible knot inside a prime cluster."""

import argparse
from collections import Counter, defaultdict
import csv
from pathlib import Path
import subprocess
import tempfile
from dyadic_witnesses import CONDUCTORS, OMIT, label, prime_divisors, trace

HERE = Path(__file__).resolve().parent
MODULUS = 65537
BITS = range(785,1461)


def direction(a,b):
    return (-a,-b) if a < 0 or (a == 0 and b < 0) else (a,b)


def initial_octaves():
    # The integer coordinate t = 8y admits an ordinary prime precisely for p < t < 2p.
    end = 16000
    events=defaultdict(list)
    for p in prime_divisors(end):
        if p in OMIT or p % 16 not in (1,15):
            continue
        a,b = trace(label(p))
        if a or b:
            v=direction(a,b)
            events[p].append((1,v))
            events[2*p].append((-1,v))
    active=Counter()
    failures=[]
    for t in sorted(set(events) | {8,end}):
        if t>end:
            break
        for sign,v in events[t]:
            if sign==-1:
                active[v]-=1
        if sum(active.values())-max(active.values(),default=0)<2:
            failures.append((t,"point"))
        for sign,v in events[t]:
            if sign==1:
                active[v]+=1
        if t<end and sum(active.values())-max(active.values(),default=0)<2:
            failures.append((t,"right"))
    assert failures[-1] == (271,"point")
    print("The last deficient ordinary-prime octave through y=2000 occurs exactly at y=271/8. Every larger real start through 2000 has two primes in every direction.")


def certified_groups(backend):
    groups=defaultdict(dict)
    with (HERE / "evidence/phase-proth-witnesses.tsv").open() as src:
        rows=csv.DictReader(src,delimiter="\t")
        assert rows.fieldnames==["bin","pattern","coefficient_offset","proth_base"]
        for row in rows:
            j,c,o,a=(int(row[k]) for k in rows.fieldnames)
            assert c not in groups[j]
            groups[j][c]=(o,a)
    codes={c for c in range(16) if any(trace(c))}
    assert len(codes)==12
    bins=[]
    primes=[]
    used=set()
    j,lower=0,1<<24
    while lower.bit_length()<=BITS.stop-1:
        if lower.bit_length() in BITS:
            assert set(groups[j])==codes
            used.add(j)
            e=lower.bit_length()-1
            step=1 << (e//2+1)
            cluster=[]
            for c,(offset,base) in groups[j].items():
                k=(lower-1)//step+offset
                p=1+k*step
                assert 0<k<step and k%2 and 2<=base<p
                assert lower<=p and 64*p<=65*lower and label(p)==c
                primes.append((p,base))
                if backend is None:
                    assert pow(base,(p-1)//2,p)==p-1
                a,b=trace(c)
                cluster.append((p,a,b))
            bins.append(sorted(cluster))
        lower=13*lower//10
        j+=1
    assert set(groups)==used
    assert len(primes)==len({p for p,_ in primes})==21432
    assert 169**2 < 2*128**2 and 13 < 2*10
    if backend is not None:
        with tempfile.TemporaryDirectory() as temporary:
            batch=Path(temporary)/"certificates.txt"
            batch.write_text("".join(f"{p} {a}\n" for p,a in primes))
            output=subprocess.check_output([str(backend),str(batch)],text=True)
            assert output==f"{len(primes)}\n"
    return bins


def augmented_rank(rows):
    basis={}
    for source in rows:
        row=list(source)
        for j in sorted(basis):
            v=row[j]
            if v:
                row[j:]=[(a-v*b)%MODULUS for a,b in zip(row[j:],basis[j][j:])]
        pivot=next((j for j,x in enumerate(row) if x),None)
        if pivot is not None:
            inverse=pow(row[pivot],-1,MODULUS)
            basis[pivot]=[x*inverse%MODULUS for x in row]
    return len(basis)


def certify_knots(bins):
    cases=0
    for primes in bins:
        records=[]
        for p,a,b in primes:
            x=pow(p,-1,MODULUS)
            square_sign=int(a==0)-int(b==0)
            ground=[t%MODULUS for t in (a,a*x,b,b*x)]
            tail=[square_sign%MODULUS,(-square_sign*x*x)%MODULUS]
            records.append((direction(a,b),ground+[0]*4+tail,[0]*4+ground+tail))
        for split in range(13):
            for delete in range(12):
                counts=(Counter(),Counter())
                equations=[]
                for k,(d,left,right) in enumerate(records):
                    if k!=delete:
                        side=int(k>=split)
                        counts[side][d]+=1
                        equations.append(right if side else left)
                upper=1+sum(min(4,sum(min(2,n) for n in c.values())) for c in counts)
                assert augmented_rank(equations)>upper
                cases+=1
    print(f"Verified {len(bins):,} full clusters and 21,432 prime certificates. All {cases:,} knot/deletion systems are inconsistent over the reals.")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proth",type=Path,help="exact GMP modular checker")
    options=parser.parse_args()
    initial_octaves()
    certify_knots(certified_groups(options.proth.resolve() if options.proth else None))
