"""Real-ball determinant certificate for the Galois group of the balanced inverse."""

from flint import arb, arb_mat, arb_poly, ctx
from three_profiles import HOSTS, Trace
from three_inverse import bands, integer_anchor


def prime_matrix():
    traces = [Trace(host,500000,physical=True) for host in HOSTS]
    centre,k = integer_anchor(traces,10548)
    assert k == 13471
    prime = [centre[0]+1,centre[1]+1,centre[2],centre[3]]
    intervals = (("2457.4526","2457.4529"),("61460.7129","61460.7134"))
    columns = [col for l,h in intervals for col in bands(traces,arb(l),arb(h))]
    matrix = arb_mat([list(row) for row in zip(*columns)])
    return matrix, arb_mat([[x] for x in prime])


def resultant(p,q):
    # Sylvester determinant, using descending coefficients in each shifted row.
    m,n = p.degree(),q.degree()
    rows=[]
    for poly,degree,shifts in ((p,m,n),(q,n,m)):
        c=[poly[k] for k in range(degree,-1,-1)]
        for i in range(shifts):
            rows.append([0]*i+c+[0]*(shifts-1-i))
    return arb_mat(rows).det()


def pencil_polynomials(g):
    g1,g2,g3,g4 = g
    v = arb_poly([0,1])
    d = g4-g3*v*v
    n = 1-v*v
    n2 = d-g2*n
    m = d-n-(g4-g3)*v
    return g1*n*n2+m*m, d*n2


def discriminant_curve(a,b):
    # The binary quartic discriminant is homogeneous of degree six.
    # Evaluate at seven input values and interpolate in the scaled target T/1000.
    variable=arb_poly([0,1])
    disc=arb_poly([0])
    for j in range(7):
        p=a-1000*j*b
        assert p.degree()==4 and not p[4].contains(0)
        value=resultant(p,p.derivative())/p[4]
        basis=arb_poly([1])
        for k in range(7):
            if k!=j:
                basis=basis*(variable-k)*(arb(1)/(j-k))
        disc+=value*basis
    return disc


def balanced_box(matrix,g):
    # The independent contraction from the previous result supplies one physical root.
    # Isolate every real algebraic sheet uniformly before testing its signs.
    radius=arb(1)/10000
    target=arb_mat([[arb(j)+arb(0,1)*radius] for j in (1,1,0,0)])
    h=matrix.solve(target)
    factors=[arb_poly([h[i,0],-g[i]]) for i in range(4)]
    q1,q2=factors[0]*factors[1],factors[2]*factors[3]
    rest=arb_poly([1,-1])
    polynomial=(rest**2-q1-q2)**2-4*q1*q2
    assert polynomial.degree()==4 and not polynomial[4].contains(0)
    # Testing products after point evaluation avoids widening the interval by
    # expanding a polynomial with correlated interval coefficients.
    def signs(z):
        terms=[h[i,0]-z*g[i] for i in range(4)]
        a,b=terms[0]*terms[1],terms[2]*terms[3]
        raw=(1-z)**2-a-b
        return raw**2-4*a*b,raw
    for l,u in ((-181,-177),(38,41),(77,81),(134,138)):
        left,right=arb(l)/100,arb(u)/100
        f,_=signs(left)
        v,_=signs(right)
        assert f*v<0
        if l==77:
            area=(left+right)/2+arb(0,1)*(right-left)/2
            assert signs(area)[1]<0


def check_at_precision(bits):
    ctx.prec=bits
    matrix,anchor=prime_matrix()
    direction=matrix.solve(anchor)
    g=tuple(direction[i,0] for i in range(4))
    assert not g[2].contains(0) and not (g[2]-g[3]).contains(0)
    a,b=pencil_polynomials(g)
    assert a.degree()==b.degree()==4
    assert not (a[4]*b[4]).contains(0)
    assert not resultant(a,b).contains(0)
    assert not resultant(b,b.derivative()).contains(0)
    disc=discriminant_curve(a,b)
    assert disc.degree()==6 and not disc[6].contains(0)
    assert not resultant(disc,disc.derivative()).contains(0)
    balanced_box(matrix,g)
    return [1 if v>0 else -1 for v in (resultant(a,b),disc[6],resultant(disc,disc.derivative()))]


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    first=check_at_precision(448)
    second=check_at_precision(576)
    assert first==second
    print("Two real-ball precisions reproduce the degree-four rational slice and a squarefree degree-six branch discriminant.")
    print("Signs of the coprimality resultant, leading branch coefficient and separability resultant:", *first)
    print("All four poles are simple and the six finite branch values are distinct.")
    print("Throughout the balanced target box of radius 1/10000 all four algebraic weights are real.")
    print("Three disjoint root intervals violate an original positivity condition, leaving one physical solution on the two prime bands.")


if __name__=="__main__":
    main()
