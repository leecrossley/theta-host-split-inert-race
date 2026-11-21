"""Replay the two exact chains of prime witnesses up to the analytic rank bound."""
import csv
from fractions import Fraction
from itertools import product
from math import isqrt
from pathlib import Path

from verify_rank import FIELDS


def verify():
    omissions = {2,547,727,571,599,343057911569,343070636977,413718278401,603339156001}
    squares = {q:{k*k%q for k in range(1,q)} for field in FIELDS for q in field}
    evidence = Path(__file__).resolve().parent / "evidence/prime-bridges.tsv"
    last = [0,0]
    small = large = 0
    with evidence.open() as source:
        table = csv.DictReader(source,delimiter="\t")
        for e,j in product(range(18,532),range(5,10)):
            row = next(table)
            assert (e,j)==(int(row["octave"]),int(row["bin"]))
            n = 4 if e<34 else e//2+1
            step = 1<<n
            for h,host in enumerate(("A","B")):
                k = (j<<e)//(5*step)+int(row[host+"_offset"])
                p = 1+k*step
                assert (j<<e) < 5*p < ((j+1)<<e)
                assert p not in omissions and p%16==1
                a = int(row[host+"_base"])
                if e<34:
                    assert a==0 and all(p%d for d in range(2,isqrt(p)+1))
                    small += 1
                else:
                    assert 0<k<step and k%2==1 and 1<a<p
                    assert pow(a,(p-1)//2,p)==p-1
                    large += 1
                flags=[(p%u in squares[u])==(p%v in squares[v]) for u,v in FIELDS]
                assert flags[h] and not flags[1-h]
                if last[h]:
                    assert last[h]<p and p*p<2*last[h]*last[h]
                else:
                    assert p<1_000_000
                last[h]=p
        assert next(table,None) is None
    assert min(last)>2**531
    # log 2 > 2*(1/3 + (1/3)^3/3 + (1/3)^5/5) > 366/529.
    t=Fraction(1,3)
    assert 2*(t+t**3/3+t**5/5)*529>366
    assert (small,large)==(160,4980)
    print("Two overlapping prime chains cover the interval from one million to exp(366). The 5,140 witnesses passed 160 trial divisions and 4,980 Proth certificates, with every successive ratio below sqrt(2).")


if __name__ == "__main__":
    verify()
