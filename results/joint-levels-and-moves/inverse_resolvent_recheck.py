"""Recover the two elliptic invariants from independent prime bands and branch cross-ratios."""

from itertools import combinations
from flint import acb_poly, arb, ctx
from inverse_branch_recheck import directions, rational_map


def pair_test(s,p,a,b):
    # Reduce the two quartics modulo z^2-s*z+p by Horner's rule.
    def remainder(poly):
        u,v=0,0
        for i in range(4,-1,-1):
            u,v = u*s+v, poly[i]-u*p
        return u,v
    u,v=remainder(a)
    w,t=remainder(b)
    return u*t-v*w


def partition_test(s,v,a,b):
    # Complete the fourth-degree polynomial to z^4+P*z^2+Q*z+R.
    # If a pair of depressed roots sums to d, its complementary pair sums to -d.
    c=[s*a[i]-b[i] for i in range(5)]
    e,d,k,m=[c[i]/c[4] for i in (0,1,2,3)]
    p=k-3*m**2/8
    q=d-m*k/2+m**3/8
    r=e-m*d/4+m**2*k/16-3*m**4/256
    u=m**2/4-k+v/c[4]
    return c[4]**3*(u*(u+p)**2-4*r*u-q**2)


def cross_ratio_invariant(test,origin,a,b):
    # On each line through the rational point interpolate the residual quadratic
    # from three original cubic evaluations, then interpolate its branch quartic.
    variable=acb_poly([0,1])
    quartic=acb_poly([0])
    slopes=(3,5,7,9,11)
    for slope in slopes:
        vals=[test(arb(s),origin+slope*s,a,b)/s for s in (1,2,3)]
        quadratic=(vals[2]-2*vals[1]+vals[0])/2
        linear=vals[1]-vals[0]-3*quadratic
        constant=vals[0]-linear-quadratic
        disc=linear**2-4*quadratic*constant
        basis=acb_poly([1])
        for k in slopes:
            if k!=slope:
                basis=basis*(variable-k)*(arb(1)/(slope-k))
        quartic+=basis*disc
    assert quartic.degree()==4 and not quartic[4].contains(0)
    roots=quartic.roots(tol=arb(1)/10**38,maxprec=640)
    assert len(roots)==4 and all(not u.overlaps(v) for u,v in combinations(roots,2))
    r,s,t,u=roots
    lam=(r-t)*(s-u)/((r-u)*(s-t))
    return 256*(1-lam+lam**2)**3/(lam**2*(1-lam)**2)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=448
    g=directions()
    a,b=rational_map(g)
    p=-g[3]/g[2]
    q=-(g[3]-g[1])/(g[2]-g[1])
    # Use the other rational pole pair as the origin for projection on the pair cubic.
    x=cross_ratio_invariant(pair_test,q,a,b)
    y=cross_ratio_invariant(partition_test,-b[4]*(p+q),a,b)
    for value,n in zip((x,y),(2159399670638032703450,1372295119267165244835)):
        assert value.imag.contains(0)
        assert arb(n)/10**18 < value.real < arb(n+1)/10**18
    print("Independent polynomial remainders, depressed quartics and four certified branch roots recover both 18-place j intervals.")


if __name__=="__main__":
    main()
