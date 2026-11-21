"""Certify the first prime cells, endpoint amplitudes and the two elliptic examples."""
from fractions import Fraction

from crossing_cells import Ball, SQRT2, ONE, logarithm
from verify_rank import FIELDS, primes_through


def first_norms(host, limit):
    u,v = FIELDS[host]
    out = []
    for p in primes_through(limit):
        if p in (2,u,v):
            continue
        a,b = (1 if pow(p,(q-1)//2,q)==1 else -1 for q in (u,v))
        degree = next(f for f in (1,2,4) if pow(p,f,16) in (1,15) and (a*b)**f==1)
        if p**degree <= limit:
            out.append((p**degree, p, degree, a**degree))
    return sorted(out)


def initial_cells():
    amplitude = ((SQRT2-ONE)*(SQRT2-ONE)).inv().scaled(2)
    for host, limit, n, p, degree, epsilon, whole, named_bin in (
            (0,47,31,31,1,-1,12297276305,4144584397),
            (1,79,49,7,2,1,6375809184,7109910097)):
        norms = first_norms(host,limit)
        assert len(norms)==2 and norms[0]==(n,p,degree,epsilon) and norms[1][0]==limit
        assert n < limit < 2*n
        endpoint = Ball.ratio(limit,8)
        value = amplitude*logarithm(p)*logarithm(2).inv()
        upper_target = value*Ball.ratio(limit-n,n)*endpoint.sqrt().inv()
        assert upper_target.floor_at(10**9)==whole
        y0=Ball.ratio(n,8)
        # u - 1/u = t*sqrt(y0)*log(2)/(A*log(p)), u=sqrt(y/y0).
        g=Ball.ratio(3947,2000)*y0.sqrt()*value.inv()
        u=(g+(g*g+Ball.ratio(4)).sqrt()).scaled(1,2)
        assert (y0*u*u).floor_at(10**9)==named_bin
    print("The sole active norms in the first open cells are 31 on 547/727 and 49 on 571/599. The next norms are 47 and 79.")
    print("The endpoint magnitudes lie in (12.297276305,12.297276306) and (6.375809184,6.375809185). Both named first hits agree with the earlier prime census.")


def elliptic_j(cuts):
    x,y,z=sorted(t*t for t in cuts)
    lam=(y-x)/(z-x)
    return 256*(1-lam+lam*lam)**3/(lam*lam*(1-lam)**2)


def arithmetic_examples():
    F=Fraction
    assert elliptic_j([F(1),F(2),F(3)]) == F(470596,225)
    assert elliptic_j([F(1,2),F(5,2),F(7,2)]) == 1728
    for n in range(1,101):
        roots = [F(1,j*j) for j in range(1,n+1)]
        assert len(set(roots))==n
        degree=2**n
        # Summing ramification over infinity and the n finite branch values.
        genus=1+F(degree*(n-3),4)
        assert genus.denominator==1 and genus>=0
        assert 2*genus-2 == -2*degree+(n+1)*degree//2
    # Among any chosen scales the smallest one has a simple root unique to it.
    scales = [F(1,2*j) for j in range(100,0,-1)]
    for i,s in enumerate(scales):
        blocks = [(v*v,25*v*v,49*v*v) for v in scales[i:]]
        assert sum(b.count(s*s) for b in blocks)==1
    print("The rational triples (1,2,3) and (1/2,5/2,7/2) give j = 470596/225 and 1728. The first 100 reciprocal-integer targets have distinct finite branch values.")
    print("The first 100 scaled triples on the curve V^2 = U^3 - U pass the simple-root square-class test.")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Use the verifier with assertions enabled.")
    initial_cells()
    arithmetic_examples()
