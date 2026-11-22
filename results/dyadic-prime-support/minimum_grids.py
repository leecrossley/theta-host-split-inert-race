"""Certify the lists of prime weights obtainable from two target-specific algebraic grids."""

from bisect import bisect_left
from collections import defaultdict
from csv import DictReader
from fractions import Fraction as F
from itertools import combinations_with_replacement, pairwise
from pathlib import Path
from support_recheck import BAD, characters, hat, prime_data, small_primes

HERE=Path(__file__).resolve().parent
EXPECTED=({7,31,47,113,127,223,271}, {3,7,79,113,127,239,241})


def field_primes(h,limit):
    residue=characters(h)
    return sorted((m,p,sign) for p in small_primes(limit) if p not in BAD[h]
                  for m,sign in [prime_data(p,h,residue)] if m<=limit)


def last_octave(data):
    end=2*6554
    events=defaultdict(list)
    for m,_,_ in data:
        events[m].append(1)
        events[2*m].append(-1)
    count=0
    last=None
    for t in sorted(set(events)|{end}):
        if t>end:
            break
        count-=events[t].count(-1)
        if count<3:
            last=(t,"point")
        count+=events[t].count(1)
        if t<end and count<3:
            last=(t,"right")
    assert last==(289,"point")


def representative_phases(masses):
    # A grid knot can change its order relative to a norm only at the
    # corresponding phase m / 2^j. Include 289, the terminal support ceiling.
    bounds={F(1),F(2)}
    for m in masses | {289}:
        bounds.add(F(m,2**(m.bit_length()-1)))
    sample=bounds-{F(2)}
    for l,r in pairwise(sorted(bounds)):
        sample.update(((2*l+r)/3,(l+2*r)/3))
    return sorted(sample)


def possible_at(data,knots,target):
    """Force zero values at knots; single interior zeros propagate along a path."""
    forced={len(knots)-1}
    cells=defaultdict(int)
    for m,p,_ in data:
        if p==target:
            continue
        k=bisect_left(knots,m)
        if k==len(knots):
            continue
        if knots[k]==m:
            forced.add(k)
        else:
            assert k>0
            cells[k]+=1
    links=defaultdict(set)
    for k,count in cells.items():
        if count>=2:
            forced.update((k-1,k))
        else:
            links[k-1].add(k)
            links[k].add(k-1)
    queue=list(forced)
    while queue:
        for k in links[queue.pop()]:
            if k not in forced:
                forced.add(k)
                queue.append(k)
    norm,= [m for m,p,_ in data if p==target]
    j=bisect_left(knots,norm)
    return j<len(knots) and (j not in forced or (norm!=knots[j] and j-1 not in forced))


def classify(h,full):
    data=[entry for entry in full if entry[0]<289]
    phases=representative_phases({m for m,_,_ in data})
    possible=set()
    count=0
    for a,b in combinations_with_replacement(phases,2):
        knots=sorted({phase*2**j for phase in (a,b) for j in range(9) if phase*2**j<=289})
        for _,p,_ in data:
            if possible_at(data,knots,p):
                possible.add(p)
        count+=1
    assert possible==EXPECTED[h]
    assert count==(666,378)[h]
    print(f"Field {h}: {count} order representatives exhaust all two-grid placements; possible primes {sorted(possible)}.")


def formulas(fields):
    grouped=defaultdict(list)
    with (HERE/'evidence/minimum-grid-formulas.tsv').open() as stream:
        rows=DictReader(stream,delimiter="\t")
        assert rows.fieldnames==["host","prime","phase","level_exponent","coefficient"]
        for row in rows:
            h=int(row["host"]);p=int(row["prime"])
            phase=F(row["phase"]);j=int(row["level_exponent"]);c=F(row["coefficient"])
            assert h in (0,1) and 1<=phase<2 and 0<=j<=6 and c
            grouped[h,p].append((phase,j,c))
    for (h,p),terms in grouped.items():
        assert len({a for a,j,c in terms})<=(1 if p==(31,7)[h] else 2)
        assert len(terms)==len({(a,j) for a,j,c in terms})
        assert max(4*a*2**j for a,j,c in terms)<512
        result={q:sign*sum(c*hat(m,a*2**j) for a,j,c in terms) for m,q,sign in fields[h] if m<512}
        assert {q:v for q,v in result.items() if v}=={p:F(1)}
    assert set(grouped)=={(h,p) for h in range(2) for p in EXPECTED[h]}
    print(f"All {sum(map(len,grouped.values()))} nonzero coefficients in fourteen exact recovery formulas checked term by term.")


def three_point_check(fields):
    for data in fields:
        early=[v for v in data if v[0]<=512]
        for m,p,_ in early:
            z=F(m,4);d=F(1,4)
            for norm,q,_ in early:
                diff=hat(norm,z+d)-2*hat(norm,z)+hat(norm,z-d)
                assert diff==(F(1,m) if p==q else 0)
    print("The rational three-level second difference selects precisely the target norm on both fields through norm 512.")


def main():
    fields=[field_primes(h,2*6554) for h in (0,1)]
    for h,data in enumerate(fields):
        last_octave(data)
        classify(h,data)
    print("Both fields have at least three good norms in every (t,2t) from t>289/2 through 6554; t=289/2 has exactly two.")
    formulas(fields)
    three_point_check(fields)


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
