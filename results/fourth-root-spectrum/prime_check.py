"""Compare complete prime sums with the rational correction at the tent breakpoints."""
from fractions import Fraction
from math import isqrt
from flint import arb, ctx

FIELDS = ((547,727,343057911569,343070636977),
          (571,599,413718278401,603339156001))


def sieve(stop):
    flags = bytearray([1])*(stop+1)
    flags[:2] = b"\x00\x00"
    for p in range(2,isqrt(stop)+1):
        if flags[p]:
            count = (stop-p*p)//p+1
            flags[p*p::p] = bytes(count)
    return [p for p in range(2,stop+1) if flags[p]]


def symbol(p, q):
    val = pow(p,(q-1)//2,q)
    if val == 1:
        return 1
    if val == q-1:
        return -1
    raise ValueError("Ramified prime in a good-prime calculation")


def norms(field, primes):
    q1,q2,u,v = field
    rows = []
    for p in primes:
        if p in (2,q1,q2,u,v):
            continue
        c,d = symbol(p,q1),symbol(p,q2)
        f = next(j for j in (1,2,4) if pow(p,j,16) in (1,15) and (c*d)**j==1)
        rows.append((p**f, -8*c**f*arb(p).log()))
    return rows


def tent(t):
    if 1 < t <= 2:
        return 2*(1-1/t)
    if 2 < t < 4:
        return 4/t-1
    return Fraction(0)


def real(q):
    return arb(q.numerator)/q.denominator


def full_level(rows,y):
    value=sum((coef*real(tent(Fraction(n)/y)) for n,coef in rows if y<n<4*y),arb(0))
    return value/real(y).sqrt()


def sparse_defect(rows,y,v):
    total=arb(0)
    active=0
    for n,coef in rows:
        t=Fraction(n)/y
        change=tent(t/v)+v*tent(t*v)-(v+1)*tent(t)
        crosses=any(k/v <= t <= k*v for k in (1,2,4))
        if not crosses:
            assert change==0
        if change:
            total+=coef*real(change)
            active+=1
    return total/(2*real(v*y).sqrt()),active


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=160
    a,r=arb(2).log(),arb(2).sqrt()
    B=2*r*(r-1)**2
    v=Fraction(101,100)
    c=real(v).sqrt()
    m=(c+1/c)/2
    starts=list(map(Fraction,(401,547,599,997,1009,1723)))
    starts.extend((Fraction(1237,2),Fraction(2019,2)))
    primes=sieve(int(8*max(starts)*v)+2)
    for field in FIELDS:
        rows=norms(field,primes)
        active=0
        signs=set()
        for y in starts:
            def move(t):
                return (full_level(rows,2*t)-full_level(rows,t))/(a*B)
            full=(move(y*v)+move(y/v))/2-m*move(y)
            front,n0=sparse_defect(rows,y,v)
            back,n1=sparse_defect(rows,2*y,v)
            sparse=(back-front)/(a*B)
            assert full.overlaps(sparse) and abs(full-sparse)<arb(1)/10**38
            active+=n0+n1
            if sparse>0: signs.add(1)
            if sparse<0: signs.add(-1)
        assert signs=={-1,1} and active>0
        print(f"{field[0]}/{field[1]}: eight full prime-move defects agree with the exact rational breakpoint corrections within 10^-38; {active} nonzero norm contributions.")


if __name__=="__main__":
    main()
