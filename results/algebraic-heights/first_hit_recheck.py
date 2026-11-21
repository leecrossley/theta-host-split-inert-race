"""Independently recover the first knots and elliptic invariants using Arb and FLINT."""
from math import isqrt
from flint import arb, ctx, fmpq, fmpq_poly


def prime(p):
    return p>1 and all(p%d for d in range(2,isqrt(p)+1))


def original_tent(n,y):
    x=arb(n)/y
    if 1<x<=2:
        return 2*(1-1/x)
    if 2<x<4:
        return 4/x-1
    return arb(0)


def prime_score(field,y):
    squares=[{j*j%q for j in range(1,q)} for q in field]
    B=2*arb(2).sqrt()*(arb(2).sqrt()-1)**2
    levels=[]
    for start in (y,2*y):
        total=arb(0)
        for p in range(3,81,2):
            if not prime(p):
                continue
            su,sv=(1 if p%q in residues else -1 for q,residues in zip(field,squares))
            if p%16 in (1,15) and su==sv:
                f=1
            elif p*p%16 in (1,15):
                f=2
            else:
                f=4
            total -= 8*su**f*arb(p).log()*original_tent(p**f,start)
        levels.append(total/(B*start.sqrt()))
    return (levels[1]-levels[0])/arb(2).log()


def verify():
    ctx.prec=256
    for field, left, right, sign, bound in (((547,727),31,47,1,12),((571,599),49,79,-1,6)):
        assert prime_score(field,arb(left)/8).contains(0)
        for j in range(1,17):
            y=arb(16*left+(right-left)*j)/128
            p=31 if left==31 else 7
            A=2/(arb(2).sqrt()-1)**2
            v=A*arb(p).log()/arb(2).log()*(8*y/left-1)/y.sqrt()
            actual=sign*prime_score(field,y)
            assert actual.overlaps(v) and max(actual.rad(),v.rad())<arb(10)**-60
        assert actual > bound
    x=fmpq_poly([0,1])
    for targets,jvalue in (([fmpq(1),fmpq(2),fmpq(3)],fmpq(470596,225)),
                           ([fmpq(1,2),fmpq(5,2),fmpq(7,2)],fmpq(1728))):
        cubic=fmpq_poly([1])
        for t in targets:
            cubic*=x-t*t
        # Short Weierstrass invariants after removing the quadratic coefficient.
        a,b,c=cubic[2],cubic[1],cubic[0]
        P=b-a*a/3
        Q=2*a**3/27-a*b/3+c
        assert -16*(4*P**3+27*Q**2)!=0
        assert 1728*4*P**3/(4*P**3+27*Q**2)==jvalue
    for j in (1,2,3,5,7,100):
        s=fmpq(1,2*j)
        mid=25*s*s
        step=24*s*s
        cubic=(s*s-x)*(mid-x)*(49*s*s-x)
        assert cubic==(mid-x)**3-step**2*(mid-x)
    print("Direct prime-ideal tents at 256 bits reproduce both first branches. Two translated Weierstrass cubics give the same exact j-invariants.")
    print("Polynomial identities recover the scaled points on V^2 = U^3 - U, including overlapping target triples.")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Use the verifier with assertions enabled.")
    verify()
