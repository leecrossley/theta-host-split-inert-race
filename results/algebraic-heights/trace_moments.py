"""Verify the rational inequalities for joint local trace projections in every late prime interval."""
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial, prod


def exp_upper(x):
    # Positive Taylor polynomial, followed by a geometric majorant for the tail.
    terms = [F(1)]
    for k in range(1,65):
        terms.append(terms[-1]*x/k)
    first = terms[-1]*x/65
    return sum(terms) + first/(1-x/66)


def verify():
    f = (547,727,571,599)
    assert max(16*p*q for p,q in combinations(f,2)) == 6_967_568 < 366**3
    assert 16*max(f) < 366**2 and exp_upper(F(59,10)) < 366
    assert exp_upper(F(16,1)/F(6397,1000)) < 16
    # In the prime-power bound, even with a second logarithm of saving
    # it suffices that exp(t/2)>400*t^4 for t>=366.
    assert F(183)**10/factorial(10) > 400*366**4
    assert 366 > 8
    for x,y,z,w in product((-1,1),repeat=4):
        a,b = F(x+y,2),F(z+w,2)
        assert a*a == F(1+x*y,2) and b*b == F(1+z*w,2)
        assert a*b == F(x*z+x*w+y*z+y*w,4)
    coarse = F(19,10)*7/366
    fine = F(73,1000)*7/366**2
    assert F(1,8)-2*coarse > F(1,20)
    assert F(3,2)*fine < F(1,1000) and F(1,4)+fine < F(1,3)
    assert (F(1,20)/F(3,2)-F(1,1000))/2 > F(1,64)
    assert (F(1,64)-F(1,3)/32)/F(3,2) == F(1,288)
    # The core of adjacent knots of ratio sqrt(2) has (U+L)/(U-L)<7.
    assert F(57,41)**2 < 2  # 8*(r+1)<49*(r-1)
    assert F(29,21)**2 < 2  # 7*(r-1)/8>1/3
    assert F(366)**3/6 > 480*(366+3)
    print("Checked the 16 trace patterns, the modulus bound 6,967,568 and the rational constants 1/20, 1/64 and 1/288 in the prime-mass estimate.")
    print("Checked the two comparisons with sqrt(2) and the numerical bound at 366 used in the proof of rank seven.")


if __name__ == "__main__":
    verify()
