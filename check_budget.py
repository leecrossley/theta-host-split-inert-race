"""Verify alternative coefficient enclosures and audit the integer allocation by exact fractions."""
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 192
    path = Path(__file__).parent/"evidence"
    inputs = [tuple(map(int,line.split())) for line in (path/"zero-brackets.tsv").read_text().splitlines()]
    bins = [tuple(map(int,line.split())) for line in (path/"phase-bins.tsv").read_text().splitlines()]
    assert len(inputs) == len(bins) == 6000
    a = arb(2).log()
    root = arb(2).sqrt()
    costs = []
    for (_,_,lo,hi),(w,b) in zip(inputs,bins):
        assert w > 0 and b >= 1
        t = (arb(lo)/10**14).union(arb(hi)/10**14)
        s = (a*t/2).sin()
        # Half-angle form of the positive compact test.
        num = (root-1)**2+4*root*s*s
        den = 2*(root-1)**2*(t*t+arb(1)/4)
        assert abs(4*s*num/(a*den)) < arb(w)/2**32
        costs.append(Fraction(w,2**(32+b)))
    total = sum(costs)
    last_included = min(cost for cost,(_,b) in zip(costs,bins) if b>1)
    next_excluded = max(cost/2 for cost in costs)
    assert next_excluded <= last_included
    target = Fraction(1,10)
    upper = Fraction(44,7)*total
    previous = Fraction(44,7)*(total+last_included)
    assert upper < target <= previous
    assert 2*arb.pi() < arb(44)/7
    exponent = sum(b for _,b in bins)
    print("6000 independent half-angle coefficient bounds passed at 192 bits.")
    print(f"Shift error upper bound {upper} < 1/10; best bound with one fewer bit {previous} >= 1/10.")
    print(f"Uniform lower log density >= 2^-{exponent}; finite prefix fraction >= 2^-{exponent+1}.")


if __name__ == "__main__":
    main()
