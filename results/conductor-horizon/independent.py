"""Check the lower triple by finite field factorization and both triples by rational interval arithmetic."""
from fractions import Fraction
from math import gcd
from pathlib import Path
from flint import arb,ctx,fmpz,nmod_poly

ROOT=Path(__file__).resolve().parent

def deg(coeff,p):
    factors=nmod_poly(coeff,p).factor()[1]
    degrees={poly.degree() for poly,_ in factors}
    assert len(degrees)==1 and all(m==1 for _,m in factors)
    return degrees.pop()


def direct_at_727():
    n=727**2
    raw=[arb(0),arb(0),arb(0)]
    for q in range(3,n+1,2):
        if q in (547,727) or not fmpz(q).is_prime():
            continue
        i=deg([2,0,-4,0,1],q)
        j=deg([-397669,0,1],q)
        f=i*j//gcd(i,j)
        norm=q**f
        if norm<n//4 or norm>n:
            continue
        last=deg([547,0,1],q)
        sign=1 if f%last==0 else -1
        for u,tick in enumerate((n-1,n,n+1)):
            w=max(Fraction(0),min(2-Fraction(tick,2*norm),Fraction(tick,norm)-1))
            raw[u]-=8*sign*arb(q).log()*w.numerator/w.denominator
    return raw


def check_archived_bands(rows):
    bands = [tuple(map(int, line.split())) for line in
             (ROOT / "evidence" / "bands.tsv").read_text().splitlines()]
    edges = [band[0] for band in bands]
    assert edges == sorted(set(edges))
    for tick, score_lo, score_hi in rows:
        assert {tick // 4, tick // 2, tick} <= set(edges)
        lower = upper = Fraction(0)
        for edge, log_lo, log_hi, ratio_lo, ratio_hi in bands:
            if tick // 4 < edge <= tick // 2:
                lower += Fraction(16 * log_lo, 2**80) - Fraction(4 * tick * ratio_hi, 2**112)
                upper += Fraction(16 * log_hi, 2**80) - Fraction(4 * tick * ratio_lo, 2**112)
            elif tick // 2 < edge <= tick:
                lower += Fraction(8 * tick * ratio_lo, 2**112) - Fraction(8 * log_hi, 2**80)
                upper += Fraction(8 * tick * ratio_hi, 2**112) - Fraction(8 * log_lo, 2**80)
        assert Fraction(score_lo, 10**12) <= lower <= upper <= Fraction(score_hi, 10**12), tick


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=192
    rows=[tuple(map(int,row.split())) for row in (ROOT/"evidence"/"raw-scores.tsv").read_text().splitlines()]
    check_archived_bands(rows)
    for x,(_,a,b) in zip(direct_at_727(),rows[:3],strict=True):
        assert arb(a)/10**12 <= x <= arb(b)/10**12
    for group,p in ((rows[:3],727),(rows[3:],547)):
        k=group[1][0]
        assert [t for t,_,_ in group]==[k-1,k,k+1]
        lower=Fraction(k*(group[0][1]-2*group[1][2]+group[2][1]),10**12)
        upper=Fraction(k*(group[0][2]-2*group[1][1]+group[2][2]),10**12)
        assert -1<lower<=0<=upper<1
        # Every contributing rational prime has a single possible residue degree.
        # The second-difference in this odd tick is a Kronecker delta at k.
        assert k==p**deg([2,0,-4,0,1],p)
    print("Direct factorization over residue fields and termwise sums independently enclose the three scores at 727 squared.")
    print("Exact rational sums of the archived sieve bands enclose all six raw scores.")
    print("Exact rational differences exclude the alternative mass at both conductor norms.")


if __name__=="__main__":
    main()
