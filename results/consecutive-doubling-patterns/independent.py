"""Repeat the sign tests by direct summation and local polynomial factorizations."""
from fractions import Fraction
from math import lcm
from pathlib import Path
from flint import acb, arb, ctx, dirichlet_char, fmpz, nmod_poly

ROOT = Path(__file__).resolve().parent
PAIRS = ((547,727),(571,599))


def value(x):
    return arb(x.numerator)/x.denominator


def gamma_series(z):
    z = min(z, 1/z)
    if z == 1:
        return arb(2).log()
    r = value(z)
    limit = 128
    total = sum((r**(2*k)*r*r.sqrt()/((2*k+1)*(2*k+2)) for k in range(limit)), arb(0))
    error = r**(2*limit+1)*r.sqrt()/((2*limit+1)*(2*limit+2)*(1-r*r))
    return total + arb(0,error)


def tent(z):
    r = min(z, 1/z)
    if r <= Fraction(1,2):
        return arb(0)
    a = (2*value(r)).sqrt()
    return 2*(a-1/a)/(2*arb(2).sqrt()*(arb(2).sqrt()-1)**2)


def factor_degrees(coefficients,p):
    _, factors = nmod_poly(list(coefficients),p).factor()
    assert all(m == 1 for _,m in factors)
    return [int(f.degree()) for f,_ in factors]


def verify(pair,data):
    u,v = pair
    characters = []
    for p in pair:
        for multiplier,residue in ((1,1),(8,5),(16,5),(16,13)):
            q = multiplier*p
            n = next(x for x in range(p-1,q,p) if x%multiplier == residue%multiplier)
            ch = dirichlet_char(q,n)
            assert ch.is_primitive() and ch.parity()==1
            characters.append(ch)
    centers = {max(t,1/t) for _,z,_ in data for t in (z,2*z,4*z)} | {Fraction(1)}
    starts = {s for _,_,y in data for s in (y,2*y,4*y)}
    sparse = {z:arb(0) for z in centers}
    raw = {y:arb(0) for y in starts}
    spec_bound = int(max(centers)*2)
    max_norm = int(4*max(starts))
    for p in range(2,max(spec_bound,max_norm)+1):
        if not fmpz(p).is_prime():
            continue
        logp = arb(p).log()
        if p <= spec_bound:
            roots = [ch(p) for ch in characters]
            n, j = p, 1
            while n <= spec_bound:
                total = sum((z**j for z in roots),acb(0))
                tr = total.real.unique_fmpz()
                assert tr is not None and abs(total-int(tr)) < arb(1)/10
                if tr:
                    mass = int(tr)*logp/arb(n).sqrt()
                    for z in centers:
                        if z/2 < n < 2*z:
                            sparse[z] += mass*tent(z/n)
                j, n = j+1,n*p
        if p in (2,u,v) or p>max_norm:
            continue
        quartic = factor_degrees((2,0,-4,0,1),p)
        real = factor_degrees((-u*v,0,1),p)
        f, = {lcm(g,h) for g in quartic for h in real}
        n = p**f
        if n > max_norm:
            continue
        imag = factor_degrees((u,0,1),p)
        g, = {lcm(f,h) for h in imag}
        mass = (8 if g != f else -8)*logp
        for y in starts:
            if y < n < 4*y:
                weight = min(2-2*y/n, 4*y/n-1)
                raw[y] += value(weight)*mass
    r = arb(2).sqrt()
    c = 2*(r-1)**2
    a = arb(2).log()
    constant = sum((int(ch.modulus())/arb.pi()).log() for ch in characters)-8*(arb.const_euler()+a)
    for z in centers:
        arch = constant*tent(z) - 8/c*(r*(gamma_series(z/2)+gamma_series(2*z))-3*gamma_series(z))
        sparse[z] = arch - sparse[z]
    for y in starts:
        raw[y] /= r*c*value(y).sqrt()
    assert 0 < sparse[Fraction(1)] < 180
    i = PAIRS.index(pair)
    for code,z,y in data:
        for k in (0,1):
            weight = 1 if (code >> (3-i-2*k))%2 else -1
            s = y*2**k
            p0,p1 = z*2**k,z*2**(k+1)
            assert weight*(raw[2*s]-raw[s]) > a*16/10
            assert weight*(sparse[max(p1,1/p1)]-sparse[max(p0,1/p0)]) > 91*a/10


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 144
    data=[]
    for row in (ROOT/"evidence/witnesses.tsv").read_text().splitlines():
        code,p,q,a,b=map(int,row.split())
        data.append((code,Fraction(p,q),Fraction(a,b)))
    assert sorted(x[0] for x in data)==list(range(16))
    for pair in PAIRS:
        verify(pair,data)
    print("Polynomial splitting, primitive Conrey traces, direct prime sums and the convergent gamma series independently confirm all 128 strict inequalities.")


if __name__ == "__main__":
    main()
