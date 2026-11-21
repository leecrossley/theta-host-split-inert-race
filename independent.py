"""Check all denominators against primitive Conrey labels and selected finite scores term by term."""
from fractions import Fraction
from hashlib import sha256
from math import gcd,isqrt
from flint import arb,ctx,dirichlet_char,nmod_poly
from decode import EVIDENCE,LIMIT,TYPE_PRIMES,EXAMPLES

LABELS=((547,546),(4376,1093),(8752,1093),(8752,5469),
        (727,726),(5816,1453),(11632,7269),(11632,1453))


def sieve(limit):
    composite=set()
    for k in range(2,isqrt(limit)+1):
        if k not in composite:
            composite.update(range(k*k,limit+1,k))
    return [j for j in range(2,limit+1) if j not in composite]


def direct_factor(p,chars):
    z=[0,0,0,0]
    for q,char,exponent in chars:
        if gcd(q,p)>1:
            continue
        a=4*int(char.chi_exponent(p))
        assert a%exponent==0
        z[a//exponent%4]+=1
    assert z[1]==z[3]
    nd=sum(z)
    if z[1]:
        assert len(set(z))==1
        return 4,1,nd//4
    if z[0] and z[2]:
        assert z[0]==z[2]
        return 2,1,nd//2
    return 1,(1 if z[0] else -1),nd


def independent_scores(p,primes):
    acc=[arb(0),arb(0),arb(0)]
    for q in primes:
        if q>p: break
        if q in (2,547,727): continue
        # Factor over the residue field, independently of character exponents.
        def deg(coeff):
            fac=nmod_poly(coeff,q).factor()[1]
            d={poly.degree() for poly,m in fac}
            assert len(d)==1 and all(m==1 for _,m in fac)
            return d.pop()
        i,j,k=deg([2,0,-4,0,1]),deg([-397669,0,1]),deg([547,0,1])
        f=i*j//gcd(i,j)
        norm=q**f
        if norm>p: continue
        sign=1 if f%k==0 else -1
        for l,tick in enumerate((p-1,p,p+1)):
            weight=max(Fraction(0),min(2-Fraction(tick,2*norm),Fraction(tick,norm)-1))
            acc[l]-=8*sign*arb(q).log()*weight.numerator/weight.denominator
    return acc


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=224
    primes=sieve(LIMIT)
    chars=[]
    for q,n in LABELS:
        char=dirichlet_char(q,n)
        assert char.is_primitive() and char.conductor()==q
        chars.append((q,char,int(char.group().exponent())))
    result=sha256()
    for p in primes:
        f,s,d=direct_factor(p,chars)
        result.update(f"{p}\t{f}\t{s}\t{d}\n".encode())
    assert result.hexdigest()==(EVIDENCE/"factor-digest.txt").read_text().strip()
    for p in TYPE_PRIMES:
        assert direct_factor(p,chars)==(1,1,8)
    table={int(row.split()[0]):list(map(int,row.split()[1:])) for row in (EVIDENCE/"finite-samples.tsv").read_text().splitlines()}
    assert set(table)==EXAMPLES
    for p,values in table.items():
        samples=[arb(2*v,1)/(2*10**24) for v in values[:3]]
        f,s,d=direct_factor(p,chars)
        tr=d*s if f==1 else 0
        obs=-p*(samples[0]-2*samples[1]+samples[2])/arb(p).log()
        assert obs.unique_fmpz()==tr==values[3]
        if p<=1009:
            for got,box in zip(independent_scores(p,primes),samples):
                assert box.contains(got)
    print("Eight primitive Conrey characters independently reproduce every ordered local denominator and both omitted type factors.")
    print("Independent polynomial factorizations and direct field sums reproduce each archived triple through p=1009.")
    print("All archived three-sample intervals recover their integer trace.")


if __name__ == "__main__":
    main()
