"""Enclose the first actual factor-400 miss subject to the explicit finite region."""
from fractions import Fraction
from pathlib import Path
import os
import subprocess
import tempfile
from flint import arb, ctx

from early_400 import EARLIEST, make_terms, zero_tail
from make_grid import build

GRID_START = Fraction(12500)
LOWER = Fraction(1620542262999, 10**6)
UPPER = Fraction(1620542276, 1000)
SAFE_WIDTH = Fraction(5298, 1000)
BRIDGES = (208856045, 245587081, 1099489440, 1176832610, 1620542263)


def arithmetic_bounds():
    a = arb(2).log()
    x = arb(12494)
    q = 16*547*727
    b, c = arb(10).exp(), arb(25).exp()
    low = 2*10**6*(-x/(20*(q*b).log())).exp()
    mid = 2000*zero_tail(b,q)*(-x/(20*(q*c).log())).exp()
    high = 2000*zero_tail(c,q)
    exceptional = 20/(arb(1).exp()*x)
    trivial = 128*(-(x-2)/2).exp()/(1-(-(x-2)).exp())
    other = 10**6*((-x/10).exp()+(-x/3).exp())
    assert low+mid+high+exceptional+trivial+other < arb(1)/1000
    # Bounds on the finite odd region at all covered log centres, and
    # a smaller bound at the one supplementary point with a small margin.
    eps=arb(4)/10**8
    p=arb(49)/50
    B=arb(10)**12
    A=(2*11632*B).log()+1
    weighted=160/B*(A/(1-p)+p/(1-p)**2)
    for height, target in ((1620550,22),(1177000,12)):
        X=arb(height)
        band=(2000*eps*(eps*X).exp()+50*((eps*X).cosh()-1))*zero_tail(arb(8590),11632)
        # Two extra log units occur at the end of the grid. Charge the
        # slope at the actual upper coordinate in that last band.
        extended = (1000*(eps*X).exp()+50)*weighted
        if height>EARLIEST.ceiling:
            slope=p*X/EARLIEST.ceiling
            extended=(1000*(eps*X).exp()+50)*160/B*(A/(1-slope)+slope/(1-slope)**2)
        tail=(1000*(X/2).exp()+50)*zero_tail(arb(EARLIEST.log_top).exp(),11632)
        assert band+extended+tail < arb(target)/1000
    assert arb(200).log() > arb(5298)/1000
    assert arb(5298)/1000 > arb(42)/8
    assert arb(276)/100-arb(2)/1000-arb(751)/1000-arb(22)/1000-arb(1)/1000 > arb(19735)/10000
    ratio=(arb(UPPER.numerator)/UPPER.denominator-arb(LOWER.numerator)/LOWER.denominator).exp()
    assert ratio < arb(10131)/10000
    print("Uniform base and arithmetic cost from log centre 12,494 is below 0.001. The finite odd region costs less than 0.022 everywhere and 0.012 at the narrow bridge.")


def verify_points():
    a,r=arb(2).log(),arb(2).sqrt()
    terms=make_terms(a,r)
    def head(x):
        return sum((b*(g*x).sin() for g,_,_,b in terms),arb(0))
    for milli in BRIDGES:
        y=head(arb(milli)/1000)
        cost=12 if milli==1176832610 else 22
        assert y-arb(751+cost+1)/1000 > arb(19735)/10000
    # Join the earlier 5.315-wide closed interval across its left edge.
    def extend(lo,hi):
        middle=(lo+hi)/2
        x=arb(middle.numerator)/middle.denominator
        rad=(hi-lo)/2
        radius=arb(rad.numerator)/rad.denominator
        if head(x+arb(0,1)*radius) < arb(117)/100:
            return 1
        if rad < Fraction(1,10**9):
            raise ArithmeticError("Unresolved left endpoint")
        return extend(lo,middle)+extend(middle,hi)
    pieces=extend(UPPER, Fraction(1620542278,1000))
    assert arb(200).log()<arb(5315)/1000
    print(f"Five extra gain witnesses use direct true-frequency sines. {pieces} closed cells join the upper endpoint to the certified no-gain interval.")


def covers(intervals, low, high):
    right=None
    for lo,hi in sorted(intervals):
        if right is None:
            if lo<low<hi:right=hi
        elif lo<right and hi>right:
            right=hi
        if right is not None and right>high:return
    raise ArithmeticError(f"Open interval cover stopped at {right}")


def check_runs(output):
    runs={1:[], -1:[]}
    counts=None
    for line in output.splitlines():
        sign,first,last=map(int,line.split())
        if sign==0:
            if counts is not None:raise ValueError("Duplicate totals")
            counts=(first,last)
        else:
            if sign not in runs or first>last:raise ValueError("Invalid signed run")
            runs[sign].append((Fraction(first,8)-SAFE_WIDTH, Fraction(last,8)))
    assert counts is not None
    assert all(value>5_900_000 for value in counts)
    runs[1].extend((Fraction(m,1000)-SAFE_WIDTH,Fraction(m,1000)) for m in BRIDGES)
    covers(runs[1],GRID_START,LOWER)
    covers(runs[-1],GRID_START,UPPER)
    return counts


def full_grid(directory):
    root=Path(__file__).resolve().parent
    CXX=os.environ.get('CXX','c++')
    results=[]
    for size in (32768,16384):
        data=directory/f"grid-{size}.txt"
        binary=directory/f"grid-{size}"
        build(root/'evidence/thousand-brackets.tsv',data,size)
        subprocess.run([CXX,'-std=c++17','-O3',f'-DGRID_SIZE={size}',str(root/'first_grid.cpp'),'-o',str(binary)],check=True)
        with data.open() as src:
            output=subprocess.check_output([str(binary)],stdin=src,text=True)
        counts=check_runs(output)
        print(f"The {size:,}-point integer grid certifies {counts[0]:,} gain samples and {counts[1]:,} loss samples, with all gaps and strict endpoints checked.")
        results.append(output)
    assert results[0]==results[1]


def main():
    if not __debug__:raise RuntimeError("Run with assertions enabled.")
    ctx.prec=144
    arithmetic_bounds()
    verify_points()
    with tempfile.TemporaryDirectory() as d:
        full_grid(Path(d))
    print("The first missed move at factor 400 is a gain. The two bounding starting heights have ratio below 1.0131 if the finite odd exceptional region is empty.")


if __name__=="__main__":
    main()
