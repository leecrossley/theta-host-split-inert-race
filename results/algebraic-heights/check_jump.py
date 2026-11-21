"""Check the four prime coefficients that keep consecutive cell constants distinct."""
from fractions import Fraction as F

# Pairs denote x + y/sqrt(2). DB is the jump in the constant part of
# sqrt(y)D(y). DA is the jump in its linear part times the prime norm.
DA = ((0,8),(-4,-12),(6,4),(-2,0))
DB = ((0,-1),(1,3),(-3,-2),(2,0))


def main():
    coefficients=[]
    for j in range(4):
        ratio=F(2**j,8)
        protected=tuple(ratio*sum(row[h]for row in DA[:j+1])-sum(row[h]for row in DB[:j]) for h in (0,1))
        assert protected!=(0,0) and DA[j]!=(0,0)
        assert all(ratio*DA[j][h]+DB[j][h]==0 for h in (0,1))
        coefficients.append(protected)
    assert coefficients==[(0,1),(-1,0),(0,-2),(2,0)]
    print("At the four knot types the protected coefficients in the derivative sum are 1/sqrt(2), -1, -sqrt(2) and 2. Each derivative jump is nonzero.")


if __name__=="__main__":
    main()
