"""Independent Arb/FLINT replay of the product separation certificates."""
import csv
from fractions import Fraction
from math import isqrt
from pathlib import Path

from flint import arb, ctx, fmpz, fmpz_mat
from crossing_recheck import TOP, field_data
from product_relations import PRIMES, ROWS

HERE = Path(__file__).resolve().parent


def read(name):
    with (HERE/"evidence"/name).open() as file:
        return list(csv.DictReader(file,delimiter="\t"))


def true_logs():
    ctx.prec = 1880
    ps = [p for p in range(3,TOP,2) if fmpz(p).is_prime()]
    data = [field_data(ps,h) for h in (0,1)]
    logs = []
    first = read("exact-crossings.tsv")
    deep = read("height-brackets-512.tsv")
    assert len(first)==len(deep)==62
    for row,record in zip(first,deep):
        assert all(row[k]==record[k] for k in ("host","sign","left_eighths"))
        _,move = data["AB".index(row["host"])]
        l,r = (int(row[k]) for k in ("left_eighths","right_eighths"))
        lo,hi = arb(l)/8,arb(r)/8
        gl,gr = lo.sqrt()*move(l,8),hi.sqrt()*move(r,8)
        u = (gr-gl)/(hi-lo)
        v = gl-u*lo
        target = int(row["sign"])*arb(3947)/2000
        radical = (target*target-4*u*v).sqrt()
        assert radical>0
        expected = int(record["height_floor_1e512"])
        early = int(row["lower_billionths"])
        found = []
        for numerator in (target+radical,target-radical):
            root = numerator/(2*u)
            if root>0 and early<root**2*10**9<early+1:
                found.append(root**2)
        assert len(found)==1 and expected<found[0]*10**512<expected+1
        logs.append(found[0].log())
    return logs+[arb(p).log() for p in PRIMES]


def independent_lattice(logs,row):
    exp,height,g,margin = row
    rounded = [int((v*10**exp+arb(1)/2).floor().unique_fmpz()) for v in logs]
    n = len(logs)
    basis = fmpz_mat([[int(i==j) for j in range(n)]+[rounded[i]] for i in range(n)])
    with (HERE/"evidence"/f"product-basis-{exp}.tsv").open() as file:
        transform = fmpz_mat([list(map(int,r)) for r in csv.reader(file,delimiter="\t")])
    assert abs(transform.det())==1
    reduced = transform*basis
    gram = (reduced*reduced.transpose()).tolist()
    dets = [1]
    for i in range(1,n+1):
        dets.append(int(fmpz_mat([r[:i] for r in gram[:i]]).det()))
    least = min(Fraction(b,a) for a,b in zip(dets,dets[1:]))
    assert isqrt(least.numerator//least.denominator)==g
    rounding = sum(abs(z-10**exp*v) for z,v in zip(rounded,logs))
    assert height*rounding + margin <= (arb(g)**2-n*height**2).sqrt()


def main():
    logs = true_logs()
    for row in ROWS:
        independent_lattice(logs,row)
    print("Arb at 1,880 bits reproduced all 62 positive quadratic roots in their 512-place brackets. FLINT determinants and direct log enclosures verified all four lattice bounds.")


if __name__=="__main__":
    main()
