"""Elliptic invariants of the pair and partition curves of the certified quartic slice."""

from flint import arb, arb_poly, ctx

from inverse_monodromy import pencil_polynomials, prime_matrix


def pair_cubic(s, p, a, b):
    # (A(z)B(w)-A(w)B(z))/(z-w), in s=z+w and p=zw.
    divided_powers = (1, s, s**2-p, s**3-2*s*p)
    return sum((a[i]*b[j]-a[j]*b[i])*p**j*divided_powers[i-j-1]
               for i in range(1,5) for j in range(i))


def partition_cubic(s, v, a, b):
    # In the coordinate s=1/T the coefficients of the binary quartic are s*A-B.
    c0,c1,c2,c3,c4 = [s*a[i]-b[i] for i in range(5)]
    return (v**3-c2*v**2+(c3*c1-4*c4*c0)*v
            +4*c4*c2*c0-c4*c1**2-c3**2*c0)


def projected_invariants(cubic, base):
    variable = arb_poly([0, 1])
    quartic = arb_poly([0])
    for slope in range(-2, 3):
        # The rational point (0,base) lies on the cubic exactly.
        line = cubic(variable, base+slope*variable)
        assert line.degree() <= 3 and line[0].contains(0)
        disc = line[2]**2-4*line[1]*line[3]
        basis = arb_poly([1])
        for j in range(-2, 3):
            if slope != j:
                basis *= (variable-j)*(arb(1)/(slope-j))
        quartic += disc*basis
    e,d,c,b,a = [quartic[i] for i in range(5)]
    invariant_i = 12*a*e-3*b*d+c**2
    invariant_j = 72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
    disc = 4*invariant_i**3-invariant_j**2
    assert not disc.contains(0)
    return 6912*invariant_i**3/disc


def elliptic_invariants(bits):
    ctx.prec = bits
    matrix,anchor = prime_matrix()
    direction = matrix.solve(anchor)
    g = [direction[i,0] for i in range(4)]
    a,b = pencil_polynomials(g)
    product1 = -g[3]/g[2]
    product2 = -(g[3]-g[1])/(g[2]-g[1])
    j_pair = projected_invariants(lambda s,p: pair_cubic(s,p,a,b), product1)
    v0 = -b[4]*(product1+product2)
    j_partition = projected_invariants(lambda s,v: partition_cubic(s,v,a,b), v0)
    return j_pair,j_partition


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    first = elliptic_invariants(448)
    second = elliptic_invariants(576)
    # Rational endpoints make the published enclosures independent of formatting.
    floors = (2159399670638032703450, 1372295119267165244835)
    scale = 10**18
    for x,y,n in zip(first,second,floors):
        assert x.overlaps(y) and arb(n)/scale < x and x < arb(n+1)/scale
        assert arb(n)/scale < y and y < arb(n+1)/scale
    x,y = second
    phi2 = (x**3+y**3-x**2*y**2+1488*x*y*(x+y)-162000*(x**2+y**2)
            +40773375*x*y+8748000000*(x+y)-157464000000000)
    assert phi2.contains(0) and abs(phi2)<arb(1)/10**90
    print("The unordered-pair j-invariant lies between 2159.399670638032703450 and 2159.399670638032703451.")
    print("The partition j-invariant lies between 1372.295119267165244835 and 1372.295119267165244836.")
    print("Both plane cubics are nonsingular. Their invariants enclose the level-two modular equation to within 10^-90.")


if __name__ == "__main__":
    main()
