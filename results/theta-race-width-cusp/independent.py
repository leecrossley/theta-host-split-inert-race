"""A second local sum, using primitive Conrey labels and a separate gamma evaluation."""
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt
from flint import arb, ctx, dirichlet_char
from width_kernel import EVIDENCE, PrimeTable, Window, ball, check_inputs

LABELS = ((547,546),(4376,1093),(8752,1093),(8752,5469),
          (727,726),(5816,1453),(11632,7269),(11632,1453))


@lru_cache(maxsize=None)
def power_trace(p,j):
    t = 0
    imag = 0
    for q,label in LABELS:
        char = dirichlet_char(q,label)
        assert char.conductor() == q and char.is_primitive()
        if gcd(q,p) != 1:
            continue
        v, d = 4*j*int(char.chi_exponent(p)), int(char.group().exponent())
        assert v % d == 0
        k = (v//d)%4
        t += (1,0,-1,0)[k]
        imag += (0,1,0,-1)[k]
    assert imag == 0
    return t


@lru_cache(maxsize=None)
def tvalue(scale):
    v = min(scale,1/scale)
    if v == 1:
        return arb(2).log()
    x = ball(v)
    if v <= Fraction(1,2):
        term = x*x.sqrt()
        total = arb(0)
        for k in range(180):
            total += term/((2*k+1)*(2*k+2))
            term *= x*x
        tail = term/(361*362*(1-x*x))
        assert tail < arb(10)**-100
        return total+arb(0,1)*tail
    t = -x.log()/2
    return t.cosh()*(2*t.cosh()).log()+t.sinh()*(2*t.sinh()).log()-t*t.exp()


def arch(b,z):
    k, c = ball(b).sqrt(), 2*(ball(b).sqrt()-1)**2
    positive = ((ball(b).log()-abs(ball(z).log()))/2).max(0)
    g = 2*k/c*positive.sinh()
    constant = sum((q/arb.pi()).log() for q,_ in LABELS)-8*(arb.const_euler()+arb(2).log())
    return constant*g-8/c*(k*(tvalue(b*z)+tvalue(b/z))-(ball(b)+1)*tvalue(z))


def arithmetic(b,z):
    k, c = ball(b).sqrt(), 2*(ball(b).sqrt()-1)**2
    total = arb(0)
    count = 0
    for p in range(2,int(b*z)+1):
        if any(p%d == 0 for d in range(2,isqrt(p)+1)):
            continue
        n,j = p,1
        while n <= b*z:
            distance = max(b/n,n/b)
            if distance < z:
                beta = arb(p).log()/arb(n).sqrt()*power_trace(p,j)
                total += beta*((ball(z/distance).log()/2).sinh())
                count += 1
            n,j=n*p,j+1
    return 4*k/c*total,count


def main():
    check_inputs()
    ctx.prec = 352
    rows = {}
    for line in (EVIDENCE/"weighted-extrema.tsv").read_text().splitlines():
        p,j,l,u = map(int,line.split())
        assert (p,j) not in rows
        rows[p,j] = (l,u)
    assert set(rows) == {(p,j) for p in (2,3,5,7,11,17,31,47,547,727) for j in (1,2)}
    for (p,j),(l,u) in rows.items():
        n = p**j
        b,z = Fraction(n),Fraction(2*n+1,2*n)
        h = ball(z).log()/2
        d = h.cosh()
        background = 2*(d-1)+2*(d*arch(b,Fraction(1))-arch(b,z))
        corrected = arb(2*l+1,1)/(2*10**32)-background
        c = 2*(arb(n).sqrt()-1)**2
        result = corrected*c/(4*arb(p).log()*h.sinh())
        assert result.unique_fmpz() == power_trace(p,j) == u
    table = PrimeTable(2000)
    counts = []
    for b,z in ((Fraction(13,8),Fraction(11,10)), (Fraction(4),Fraction(7,4)),
                (Fraction(15),Fraction(3,2)),(Fraction(31),Fraction(9,8)),
                (Fraction(95,2),Fraction(25,16)),(Fraction(257),Fraction(7,4)),
                (Fraction(727),Fraction(8,5))):
        _,value,_ = Window(b,table).measurements(z)
        d = (ball(z).log()/2).cosh()
        back = 2*(d-1)+2*(d*arch(b,Fraction(1))-arch(b,z))
        direct,count = arithmetic(b,z)
        assert (value-back-direct).contains(0) and abs(value-back-direct)<arb(10)**-83
        counts.append(count)
    assert min(counts) == 0 and max(counts)>80
    print("Independent gamma series and hyperbolic logarithms recover all 20 traces using the eight primitive Conrey characters.")
    print("Seven additional bands match direct local sums to 83 decimals, from an empty band to one containing "+str(max(counts))+" prime powers.")


if __name__ == "__main__":
    main()
