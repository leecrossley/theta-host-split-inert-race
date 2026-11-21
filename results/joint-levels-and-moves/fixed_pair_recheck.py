"""Recheck the sharp two-anchor error by direct equations and nodal secants."""

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from flint import arb, arb_mat, ctx

from three_inverse_recheck import FIELDS, Reference, conic, curve, derivative, original

DATA = Path(__file__).resolve().parent/"evidence"


def centred(values):
    return (values[0]-1,values[1]-1,values[2],values[3])


def local_solution(refs):
    first,second = [centred(original(refs,arb(y))) for y in (2457,13471)]
    floors=(61455609822534600169779,246827861901412719,394992237138323180,2274350938622389)
    sign=(-1,1,1,-1)
    half=arb(1)/(2*10**18)
    middle=[(arb(n)+arb(1)/2)/10**18 for n in floors]
    boxes=[x+arb(0,1)*half for x in middle]
    coeff=conic(refs,middle[0],arb("0.011"))
    def data(vec):
        y,a,b,e=vec
        V=centred(curve(coeff,y)); deriv=derivative(coeff,y)
        equations=arb_mat([[a*A+b*B+(1-a-b)*v-e*s] for A,B,v,s in zip(first,second,V,sign)])
        jac=arb_mat([[(1-a-b)*t,A-v,B-v,-s]
                     for t,A,B,v,s in zip(deriv,first,second,V,sign)])
        return equations,jac
    f,J=data(middle)
    numerical=J.inv()
    inverse=arb_mat([[numerical[i,j].mid() for j in range(4)] for i in range(4)])
    assert not inverse.det().contains(0)
    _,wide=data(boxes)
    identity=arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    defect=identity-inverse*wide; shift=inverse*f
    refined=[]
    for i in range(4):
        norm=sum(abs(defect[i,j]) for j in range(4))
        assert norm<arb(1)/1000000 and abs(shift[i,0])+norm*half<half
        refined.append(middle[i]-shift[i,0]+arb(0,1)*norm*half)
    free=1-refined[1]-refined[2]
    assert arb(358179900960264099)/10**18<free<arb(358179900960264100)/10**18

    # Find the two supporting covectors as cofactors, independent of the quadratic.
    V=centred(curve(coeff,refined[0]))
    tangent=derivative(coeff,arb("61455.61")+arb(0,arb("0.01")))
    for omitted, direction in ((0,1),(1,-1)):
        ids=[i for i in range(4) if i!=omitted]
        v=[first[i]-second[i] for i in ids]
        u=[first[i]-V[i] for i in ids]
        n=[v[(j+1)%3]*u[(j+2)%3]-v[(j+2)%3]*u[(j+1)%3] for j in range(3)]
        mass=sum(sign[i]*x for i,x in zip(ids,n))
        n=[x/mass for x in n]
        assert all(sign[i]*x>0 for i,x in zip(ids,n))
        assert direction*sum(tangent[i]*x for i,x in zip(ids,n))>0
    print("Independent original equations contract all four 18-place boxes; the fifth weight is recovered.")
    print("Both supporting-plane derivatives have the required strict signs throughout the local band.")
    broad=[arb(1)/10000,arb(1)/1000000,arb(1)/1000000,arb(1)/1000000]
    domain=[x+arb(0,h) for x,h in zip(middle,broad)]
    _,jac=data(domain)
    wide_defect=identity-inverse*jac
    for i in range(4):
        delta=sum(abs(wide_defect[i,j])*broad[j] for j in range(4))
        displacement=sum(abs(inverse[i,j]) for j in range(4))*4/10**9
        assert delta<broad[i]/2 and abs(shift[i,0])+delta+displacement<broad[i]
    diff=[b-a for a,b in zip(first,second)]
    coeffs=[]
    for chosen in (0,1):
        unit=[int(i==chosen) for i in range(4)]
        def value(column):
            return arb_mat([list(row) for row in zip(diff,sign,unit,column)]).det()
        coeffs.append((value([p for p,q in coeff]),value([d-a for d,a in zip((-1,-1,0,0),first)]),
                       value([q for p,q in coeff])))
    (a,b,c),(d,e,f)=coeffs
    assert not (a*d).contains(0) and not ((a*f-c*d)**2-(a*e-b*d)*(b*f-c*e)).contains(0)
    assert not b.contains(0) or not e.contains(0)
    print("A uniform interval contraction covers target perturbations 4/10^9 in all four coordinates.")
    print("Independent projected minors and the quadratic resultant give a degree-two rational map.")

    return first,second


def secant_error(lo,hi):
    """Bound the prime part's smooth second derivative and its atomic second derivative."""
    def real(x):
        return arb(x.numerator)/x.denominator
    z0,z1=real(lo).exp(),real(hi).exp()
    h=real(hi-lo);a=arb(2).log();r=arb(2).sqrt();B=2*r*(r-1)**2
    def bound(mult):
        L,U=z0*mult,z1*mult
        logweight=8*(2*U).log()
        smooth=(2*U+1)*logweight*4*(a/2).sinh()/(B*(L/2).sqrt())
        jumps=arb(0)
        for k,coefficient in ((arb(1)/2,2/B),(1,4*(a/2).cosh()/B),(2,2/B)):
            jumps+=coefficient*logweight*(1+k*U*h)/(k*L).sqrt()
        return smooth*h*h/32 + jumps*h/4
    level=bound(1)
    move=(bound(r)+bound(1/r))/a
    return (level,level,move,move)


def global_check(refs,first,second):
    substrate=((1,"1228.5"),(4914,"6735.5"),(26942,"61455.601"),("61455.619",300000))
    trees=defaultdict(list)
    with (DATA/"fixed-pair-separators.tsv").open() as source:
        next(source)
        for line in source:
            k,sequence,*row=line.split(); path="" if sequence=="-" else sequence
            d=2**len(path); numerator=int(path,2) if path else 0
            trees[int(k)].append((Fraction(numerator,d),Fraction(numerator+1,d),tuple(map(int,row))))
    assert set(trees)==set(range(4))
    @lru_cache(maxsize=None)
    def evaluate(x):
        xx=arb(x.numerator)/x.denominator
        (a,b),(c,d)=[ref.observation(xx) for ref in refs]
        return (a,c,b,d)

    count=0
    for region,(start,end) in enumerate(substrate):
        grid=10**10
        s,t=[(arb(n)*2*arb(2).sqrt()).log()*grid for n in (start,end)]
        u=s.floor().unique_fmpz(); v=t.ceil().unique_fmpz()
        assert u is not None and v is not None
        lower=Fraction(int(u),grid);width=Fraction(int(v-u),grid)
        last=Fraction(0)
        for x,z,row in sorted(trees[region]):
            assert x==last and x<z
            last=z
            lin=lambda V: sum(n*q for n,q in zip(row,V))
            rhs=sum(map(abs,row))*arb("0.00227436")
            assert min(lin(first),lin(second))>rhs
            stack=[(lower+x*width,lower+z*width,0)]
            while stack:
                l,h,depth=stack.pop()
                err=secant_error(l,h)
                loss=sum(abs(n)*t for n,t in zip(row,err))
                if min(lin(evaluate(l)),lin(evaluate(h)))-loss>rhs:
                    count+=1
                    continue
                assert depth<36
                mid=(l+h)/2
                stack.extend(((l,mid,depth+1),(mid,h,depth+1)))
        assert last==1
    return count


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=264
    refs=[Reference(host,True) for host in FIELDS]
    a,b=local_solution(refs)
    print("Independent polynomial splitting and",global_check(refs,a,b),
          "closed secants validate every global separator.")


if __name__=="__main__":
    main()
