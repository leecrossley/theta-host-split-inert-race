"""Exact cell counts for signed moves at starting-height ratio 200."""
import csv
from collections import Counter
from fractions import Fraction
from pathlib import Path

from crossing_cells import LIMIT, SQRT2, ZERO, events, primes_through, score
from height_ratios import HOSTS, prime_data

HERE = Path(__file__).resolve().parent
K = 200
SCALE = 8*K
ROOT = SQRT2.scaled(10)
UNIT = 10**12


def early_identities(primes):
    """Compare prime, norm, band and sign on every early open joint cell."""
    knots = {0,641}
    characters = {}
    for host in HOSTS:
        terms = {}
        for p in primes:
            if p > 641:
                break
            first = prime_data(host,p,True)
            assert first == prime_data(host,p,False)
            if first is None or first[0] > 641:
                continue
            norm,trace = first
            terms[p] = (norm,trace)
            knots.update(k for k in (norm,2*norm,4*norm,8*norm) if k < 641)
        characters[host] = terms
    empty = []
    ordered = sorted(knots)
    for left,right in zip(ordered,ordered[1:]):
        mid = Fraction(left+right,2)
        rows = []
        for host in HOSTS:
            rows.append({p:(n,sum(f*n<mid for f in (1,2,4)),trace)
                         for p,(n,trace) in characters[host].items() if n<mid<8*n})
        a,b = rows
        if not a and not b:
            empty.append((left,right))
            continue
        for eps in (1,-1):
            assert a != {p:(n,band,eps*trace) for p,(n,band,trace) in b.items()}
    assert len(ordered)-1 == 51 and empty == [(0,31)]
    print("51 joint cells through 641/8 checked for signed identities. The only common zero cell is (0,31/8).")


def contacts(early,late):
    fast = {n:(a*ROOT,b*ROOT.inv()) for n,(a,b) in late.items()}
    slow = {K*n:delta for n,delta in early.items() if K*n<=LIMIT}
    keys = sorted({SCALE,LIMIT,*(n for n in slow if n>=SCALE),*(n for n in fast if n>=SCALE)})
    u = v = p = q = ZERO
    for n,(x,z) in slow.items():
        if n < SCALE:
            u,v = u+x,v+z
    for n,(x,z) in fast.items():
        if n < SCALE:
            p,q = p+x,q+z
    rows = []
    for left,right in zip(keys,keys[1:]):
        if left in slow:
            x,z = slow[left]
            u,v = u+x,v+z
        if left in fast:
            x,z = fast[left]
            p,q = p+x,q+z
        for eps in (1,-1):
            a,b = (u-p,v-q) if eps==1 else (u+p,v+q)
            lhs = (a.scaled(left,SCALE)+b).sign()
            rhs = (a.scaled(right,SCALE)+b).sign()
            assert lhs and rhs
            if lhs != rhs:
                y = (-b)*a.inv()
                value = score(u,v,y)
                rows.append((eps,left,right,y.floor_at(UNIT),value.sign(),value.abs().floor_at(UNIT)))
    return len(keys)-1,rows


def main():
    primes = list(primes_through(LIMIT))
    early_identities(primes)
    by_host = [events(primes,i) for i in range(2)]
    rows = []
    cells = []
    for i,first in enumerate("AB"):
        for j,second in enumerate("AB"):
            count,found = contacts(by_host[i],by_host[j])
            cells.append(count)
            rows.extend((first,second,*row) for row in found)
    with (HERE/"evidence/factor-200-contacts.tsv").open() as file:
        actual = [(r[0],r[1],*map(int,r[2:])) for r in list(csv.reader(file,delimiter="\t"))[1:]]
    assert rows == actual
    positive = sorted((r for r in rows if r[6]),key=lambda r:r[7])
    assert len(rows)==70 and len(positive)==56 and len({r[7] for r in positive})==56
    counts = Counter(r[0]+r[1] for r in positive)
    assert [counts[i] for i in ("AA","AB","BA","BB")] == [14,16,15,11]
    below = max(r[7] for r in positive if r[7]<19735*10**8)
    above = min(r[7] for r in positive if r[7]>19735*10**8)
    assert (below,above)==(1745425318934,2057225896979)
    print(f"Factor 200, starts 1 through 5000. Joint cells in AA, AB, BA, BB order: {cells}.")
    print("All endpoints strict. Signed equations give 70 interior roots, 56 at different positive amplitudes.")
    print("Positive contacts in AA, AB, BA, BB order: 14, 16, 15, 11. Adjacent amplitude floors at 10^12: 1745425318934 and 2057225896979.")


if __name__=="__main__":
    main()
