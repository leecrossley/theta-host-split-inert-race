"""Verify the rational inequalities for joint local trace projections in every late prime interval."""
from fractions import Fraction as F
from itertools import combinations, product


def verify():
    conductors = (547, 727, 571, 599)
    largest_pair = max(16*p*q for p,q in combinations(conductors, 2))
    bound = 2_176_229_286_416
    assert largest_pair == 6_967_568 < bound
    assert bound * 3**29 < 8**29 and F(11, 4)**2 < 8  # M < e^29 and log 2 > 2/3
    assert F(11, 24) * 366**2 > F(29, 6) * 366 + 60  # Grenie-Molteni plus prime powers < 2 sqrt(x) t^2
    assert 2 * 10**7 * bound * 366**12 < 2**174 and 2**174 * 3**183 < 8**183
    assert F(1, 10**7 * 366**10) < min(F(19, 10*366), F(73, 1000*366**2))
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
