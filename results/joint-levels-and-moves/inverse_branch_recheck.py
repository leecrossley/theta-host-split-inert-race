"""Independent complex discs for all six critical values of the rational inverse slice."""

from math import factorial
from itertools import combinations
from flint import acb, acb_poly, arb, arb_mat, ctx
from three_inverse_recheck import FIELDS, Reference, conic, original

CENTRES = (
    ("-0.0137941315920537","0"),("0.85591355713268","0"),
    ("-0.0032038874547672","-0.0019222008482537"),("-0.0032038874547672","0.0019222008482537"),
    ("0.0021722391060978","-0.0037721135455298"),("0.0021722391060978","0.0037721135455298"),
)


def directions():
    refs=[Reference(host,True,500000) for host in FIELDS]
    anchor=original(refs,arb(13471))
    # Two original tent values per physical band determine its two columns.
    columns=[]
    for t in ("2457.45276163","61460.71312352"):
        coeff=conic(refs,arb(t),arb(1)/100000)
        columns.extend(zip(*coeff))
    matrix=arb_mat([list(row) for row in zip(*columns)])
    vec=matrix.solve(arb_mat([[x] for x in anchor]))
    return tuple(vec[i,0] for i in range(4))


def rational_map(g):
    g1,g2,g3,g4=g
    v=acb_poly([0,1])
    den=(g4-g3*v*v)*(g4-g2+(g2-g3)*v*v)
    num=acb_poly([0])
    for i in range(5):
        s=arb(i)
        h=(1-s*s)/(g4-g3*s*s)
        second_weight=s*(1-g3*h)
        first_weight=1-h-second_weight
        target=g1*h+first_weight**2/(1-g2*h)
        basis=acb_poly([1])
        for j in range(5):
            if i!=j:
                basis=basis*(v-j)*(arb(1)/(i-j))
        num+=target*den(s)*basis
    return num,den


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=288
    g=directions()
    num,den=rational_map(g)
    assert num.degree()==den.degree()==4
    raw=num.derivative()*den-num*den.derivative()
    # The seventh-degree terms of the Wronskian cancel identically.
    critical=acb_poly([raw[k] for k in range(7)])
    assert not critical[6].contains(0)
    derivatives=[critical]
    for _ in range(6):
        derivatives.append(derivatives[-1].derivative())
    rad=arb(1)/10**10
    discs=[]
    values=[]
    for re,im in CENTRES:
        centre=acb(arb(re),arb(im))
        rect=acb(arb(re)+arb(0,1)*rad,arb(im)+arb(0,1)*rad)
        main_term=abs(derivatives[1](centre))*rad
        rest=abs(critical(centre))
        for j in range(2,7):
            rest+=abs(derivatives[j](centre))*rad**j/factorial(j)
        assert rest<main_term and not den(rect).contains(0)
        discs.append(rect)
        values.append(num(rect)/den(rect))
    assert all(not p.overlaps(q) for p,q in combinations(discs,2))
    assert all(not p.overlaps(q) for p,q in combinations(values,2))
    assert 1375<values[0].real<1376 and -arb(25)/1000<values[1].real<-arb(24)/1000
    assert all(abs(val.imag)>10000 for val in values[2:])
    # Four real simple sheets above zero are witnessed by opposite signs,
    # including a tiny bracket around the rational parameter 1.
    for lo,hi in (("-0.009","-0.008"),("-0.0021","-0.0018"),("0.74","0.76"),("0.99","1.01")):
        assert num(arb(lo)).real*num(arb(hi)).real<0
    # At infinity the two denominator quadratics contribute two real poles.
    assert g[3]/g[2] > 0 and (g[3]-g[1])/(g[2]-g[1]) < 0
    print("Independent original-tent interpolation and six disjoint Rouche discs of radius 10^-10.")
    print("All six critical images are separated. The two real branch values lie in (-0.025,-0.024) and (1375,1376).")
    print("The other four branch values have imaginary part of absolute value above 10000.")
    print("The slice above zero has four distinct real parameters; two of the four simple poles are real.")


if __name__=="__main__":
    main()
