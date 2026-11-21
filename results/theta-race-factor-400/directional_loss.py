"""Enclose the first missed loss and give the smaller finite test for every earlier height."""
from fractions import Fraction
from pathlib import Path
import os
import subprocess
import tempfile
from flint import arb, ctx

from early_400 import make_terms, zero_tail
from cell_check import even_moments
from first_cover import covers
from dyadic_alias import verify as verify_alias
from loss_moments import filtered_tail
from make_grid import build
from paired_windows import full_mass
from thousand_window import grid_enclosures

END = 43380920
LOWER = Fraction(4338091247,100)
LEFT = Fraction(7715795044,100)
WIDTH = Fraction(555,100)
NARROW_TOP = 3_500_000_000
HIGH_CEILING = 77157957
POINTS = (3736873352,5287782048,21001988699,21152680551,25023214200,
          25813056312,27767858727,27767863907,36683185825,37238881656,
          39604668514,40127139341,40127144590,40277113564)
MEANS = (21893231280,32155058334)


def rational(x):
    return arb(x.numerator)/x.denominator


def budgets():
    a = arb(2).log()
    q = 11632
    mass = zero_tail(arb(8590),q)
    eps = arb(12)/10**10
    cap = arb(NARROW_TOP)
    extent = 3*a+arb(3)/500
    even = (2000*eps*(eps*extent).exp()+40*((eps*extent).cosh()-1))*mass
    even += (40*(extent/2).exp()+40)*zero_tail(cap,q)
    assert even < arb(15)/10**5
    max_n = (END+1)//2+2
    radius = arb(max_n).sqrt()/250000
    assert radius < arb(2)/100
    assert 750000*arb(max_n).sqrt() < cap
    assert (arb(1)/2/250000).cosh()/3 < arb(34)/100
    x0 = arb(12494)
    # The displayed exponent decreases in the centre, including at the
    # smallest grid centre with an averaging displacement above 0.006.
    high = 10**6*((x0+arb(6)/1000)/2+x0/2*(arb(34)/100).log()).exp()
    assert high < arb(1)/10**10
    x = arb(END+1)
    odd = (2000*eps*(eps*x).exp()+50*((eps*x).cosh()-1))*mass
    assert odd + high + 50*zero_tail(cap,q) < arb(14)/1000
    assert arb(3805)/(250000*arb(3).sqrt()) < arb(9)/1000
    assert arb(276)/100-arb(2+752+9+14+1)/1000 > arb(19735)/10000
    assert arb(200).log() > arb(5298)/1000
    assert arb(42)/8 + arb(4)/100 < arb(5298)/1000
    assert rational(LOWER) > arb(26)*arb(1620542276)/1000
    assert rational(LEFT) < arb(48)*arb(1620542262999)/10**6
    # The window with no loss uses the separate pointwise fourth-root formula.
    x = arb(HIGH_CEILING+1)
    eps = arb(5)/10**10
    cut = arb(10)**12
    top_log = 19289505
    assert arb(HIGH_CEILING)/4+arb(5000000).log() < top_log
    p = arb(49)/50*x/HIGH_CEILING
    log_scale = (2*q*cut).log()+1
    near = (2000*eps*(eps*x).exp()+50*((eps*x).cosh()-1))*mass
    far = (1000*(eps*x).exp()+50)*160/cut*(log_scale/(1-p)+p/(1-p)**2)
    assert near+far+50*zero_tail(arb(top_log).exp(),q)+arb(1)/10**6 < arb(1)/100
    assert rational(LEFT+WIDTH) < HIGH_CEILING
    assert arb(111)/100 + arb(752+10+1)/1000 < arb(19735)/10000
    assert 2*rational(WIDTH).exp() > arb(514475)/1000
    print("A missing loss through log centre 43,380,912.47 would force an odd zero above real part 0.5000000012 at absolute ordinate between 8,590 and 3,500,000,000. Box smoothing costs below 0.009 on the retained spectrum.")


def signed_points(terms,a,r):
    assert sum(abs(b)*g for g,_,_,b in terms) < 3805
    filtered = filtered_tail(terms,a,r)
    threshold = arb(19735)/10000
    for point in POINTS:
        x = arb(point)/1000
        val = -sum(b*(g*x).sin() for g,_,_,b in terms)
        assert val-arb(752+14+9+1)/1000 > threshold
    for point in MEANS:
        x = arb(point)/1000
        val = sum(b*(g*x).sin() for g,b in filtered)
        assert val-arb(20+14+9+1)/1000 > threshold
    print("Four filtered cosine moments on 2,048 interval pairs bound the triangular omitted mass below 0.020. Fourteen supplementary points and two positive triangular averages pass the strict loss inequality.")
    v0,v1 = even_moments(a,r)
    m1 = (v0-v1)/2-sum(2*even*sine**2 for _,sine,even,_ in terms)
    assert 2*(m1-arb(15)/10**5)/a > arb(59)/100
    x = arb(43380915948)/1000
    value = -sum(b*(g*x).sin() for g,_,_,b in terms)
    assert arb(19)/10 < value < arb(196)/100
    print("At log centre 43,380,915.948 the retained signed loss lies between 1.90 and 1.96. The unknown positive spectrum has dyadic-moment ambiguity exceeding 0.59 in each direction in the relaxed measure problem.")


def check_grid(output):
    width = Fraction(5298,1000)
    box = Fraction(2,100)
    tri = Fraction(6,1000)
    # A positive mean centred at c with radius s supplies a start inside
    # (h,h+width] whenever c+s-width < h < c-s.
    windows = []
    count = None
    for line in output.splitlines():
        sign,first,last = map(int,line.split())
        if sign == 0:
            assert first==0 and count is None
            count = last
        else:
            assert sign==-1 and first<=last
            windows.append((Fraction(first,8)+box-width,Fraction(last,8)-box))
    for points,rad in ((POINTS,box),(MEANS,box+tri)):
        windows.extend((Fraction(j,1000)+rad-width,Fraction(j,1000)-rad) for j in points)
    covers(windows,Fraction(12500),LOWER)
    assert count==159453395
    return count


def grid(root,tmp):
    results = []
    for size in (32768,16384):
        data = tmp/f"coeff-{size}"
        exe = tmp/f"grid-{size}"
        build(root/"evidence/thousand-brackets.tsv",data,size,END)
        subprocess.run([os.environ.get("CXX","c++"),"-std=c++17","-O3",
                        f"-DGRID_SIZE={size}","-DEXTENDED_LOSS_GRID",
                        str(root/"first_grid.cpp"),"-o",str(exe)],check=True)
        with data.open() as inp:
            result = subprocess.check_output([str(exe)],stdin=inp,text=True)
        results.append(result)
        print(f"The {size:,}-point integer grid certifies {check_grid(result):,} negative samples and the open real-height cover through 43,380,912.47.")
    assert results[0]==results[1]


def negative_window(terms,a,r):
    left,width = rational(LEFT),rational(WIDTH)
    nodes,trunc,interp,shift,removed = grid_enclosures(left,width)
    assert all(-v+trunc+interp+shift < arb(111)/100 for v in nodes)
    tail = full_mass(a,r,2*(r-1)**2)-removed
    assert tail+120*arb(15)/10**5/(54*a) < arb(752)/1000
    def cell(j,depth):
        x = left+arb(2*j+1,1)*width/(2*(1<<depth))
        if -sum(b*(g*x).sin() for g,_,_,b in terms) < arb(111)/100:
            return 1
        if depth==27:raise ArithmeticError("Unresolved no-loss cell")
        return cell(2*j,depth+1)+cell(2*j+1,depth+1)
    cells = sum(cell(j,9) for j in range(512))
    print(f"An interval DFT with 32,769 nodes and {cells:,} direct closed cells enclose the 5.55-wide no-loss interval beginning at 77,157,950.44. Its completion factor exceeds 514.475.")


def main():
    if not __debug__:raise RuntimeError("Run with assertions enabled.")
    ctx.prec=144
    budgets()
    a,r=arb(2).log(),arb(2).sqrt()
    terms=make_terms(a,r)
    signed_points(terms,a,r)
    verify_alias()
    negative_window(terms,a,r)
    root=Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory() as tmp:
        grid(root,Path(tmp))
    print("Under the stated finite exclusions the first missed loss lies after 26 times and before 48 times the first missed gain in centred logarithmic height.")


if __name__=="__main__":
    main()
