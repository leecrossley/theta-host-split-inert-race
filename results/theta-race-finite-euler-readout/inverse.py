"""Exhaust the discrete kernel and validate the explicit dyadic inverse in rational arithmetic."""
from fractions import Fraction


def w(k,n):
    return max(Fraction(0),min(2-Fraction(k,2*n),Fraction(k,n)-1))


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    cap=256
    for n in range(1,cap+1):
        for k in range(1,4*n+2):
            got=k*(w(k+1,n)-2*w(k,n)+w(k-1,n))
            exp=(1 if k==n else 0)-(3 if k==2*n else 0)+(2 if k==4*n else 0)
            assert got==exp
    weights=[Fraction(0)]+[Fraction((-1)**i*(3*i%29-12),i%7+1) for i in range(1,cap+1)]
    samples=[sum(weights[n]*w(k,n) for n in range(1,cap+1)) for k in range(cap+2)]
    for k in range(1,cap+1):
        total,j,m=Fraction(0),0,k
        while True:
            total += (2**(j+1)-1)*m*(samples[m+1]-2*samples[m]+samples[m-1])
            if m%2: break
            m,j=m//2,j+1
        assert total==weights[k]
    print("All single-mass quarter-grid stencils through norm 256 match the three prescribed jumps exactly.")
    print("Exact rational inversion recovers 256 signed masses, including every dyadic overlap.")


if __name__ == "__main__":
    main()
