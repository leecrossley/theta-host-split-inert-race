"""Check all six crossings, two independent prime sums and the remaining grid gaps."""
from bisect import bisect_right
from fractions import Fraction
from pathlib import Path
import sys
from flint import arb, ctx
from band_values import values
from near_roots import NEIGHBORHOODS, all_crossings, exact_roots, higher_degrees, knots
from sharp_deadline import interval

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")

HEIGHT = 10**9
EXCEPTIONS = ((1,6145473,20405658),(-1,19787304,61780789),(-1,549154839,1714594397))


def check_global():
    rows=[tuple(map(int,line.split())) for line in
          (Path(__file__).parent/'evidence'/'windows.tsv').read_text().splitlines()[1:]]
    for sign in (1,-1):
        good=[y for y,l,u in rows if (l if sign>0 else -u)>19735]
        exceptional={(p,q) for s,p,q in EXCEPTIONS if s==sign}
        seen=set()
        for h in {10**6,HEIGHT,*(v for v in good if 10**6<=v<=HEIGHT)}:
            nxt=good[bisect_right(good,h)]
            if (h,nxt) in exceptional:
                seen.add((h,nxt))
            else:
                assert Fraction(2*nxt,h)<6
        assert seen==exceptional


def main(binary):
    ctx.prec=160
    keys,direct,roots=exact_roots()
    fine=values(binary,keys,8,fine=True)
    assert all(direct[k].overlaps(fine[k]) for k in keys)
    groups=knots(higher_degrees(8*sum(NEIGHBORHOODS[-1])))
    second=all_crossings(groups,fine)
    for (u,v,x),(i,j,z) in zip(roots,second,strict=True):
        assert (u,v)==(i,j) and x.overlaps(z) and z.rad()<arb('0.001')
    factors=[2*roots[i+1][2]/roots[i][2] for i in (0,2,4)]
    assert factors[0]>factors[1]>factors[2]>6
    corroborated=[2*second[i+1][2]/second[i][2] for i in (0,2,4)]
    assert all(x.overlaps(y) and y.rad()<arb("1e-11") for x,y in zip(factors,corroborated))
    for i,((sign,p,q),c) in enumerate(zip(EXCEPTIONS,factors,strict=True)):
        a,b=roots[2*i][2],roots[2*i+1][2]
        assert p<a<b<q and 2*a/p<3 and 2*q/b<3
        assert 10**6<a and 2*b/6<HEIGHT
        if i<2:
            assert 2*b/6<roots[2*i+2][2]
        lo,hi=roots[2*i][0:2];low,high=roots[2*i+1][0:2]
        print(f"{sign:+d} left in 8y cell [{lo},{hi}], y {interval(a,7)}")
        print(f"{sign:+d} right in 8y cell [{low},{high}], y {interval(b,7)}; factor {interval(c,12)}")
    integer_factors=[2*roots[i+1][2]/roots[i][2].ceil().unique_fmpz() for i in (0,2,4)]
    assert integer_factors[0]>integer_factors[1]>integer_factors[2]>6
    print("sharp integer-start factors for gain and loss:", interval(integer_factors[0],12), interval(integer_factors[1],12))
    counts=[]
    for i in (0,2,4):
        first=int(roots[i][2].ceil().unique_fmpz())
        last=int((roots[i+1][2]/3).floor().unique_fmpz())
        counts.append(last-first+1)
    assert counts==[480218,279415,4613140] and sum(counts)==5372773
    print("6H misses 480218 gains and 4892555 losses at 5372773 distinct integer starts")
    check_global()
    print(f"{len(keys)} boundary knots and endpoints agree in the independent 160-bit and integer sums")
    print("all remaining same-sign grid gaps on [1000000,1000000000] finish before 6H")


if __name__=='__main__':
    main(Path(sys.argv[1]).resolve())
