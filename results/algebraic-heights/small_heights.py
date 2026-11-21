"""Prove the sharp six-score rank transition on the initial interval using exact minors."""
import csv
from pathlib import Path
from functools import cmp_to_key
from fractions import Fraction
from math import comb
from verify_rank import FIELDS, primes_through

# A knot is (a + b*sqrt(2))/8. Entries and determinants below lie in Z[sqrt(2)].
START = (0,79)
STOP = (8*1500,0)
ZERO = (0,0)


def add(x,y):
    return x[0]+y[0], x[1]+y[1]


def sub(x,y):
    return x[0]-y[0], x[1]-y[1]


def scale(x,k):
    return x[0]*k,x[1]*k


def mul(x,y):
    a,b=x; c,d=y
    return a*c+2*b*d, a*d+b*c


def divide(x,y):
    a,b=x; c,d=y
    norm=c*c-2*d*d
    u,v=divmod(a*c-2*b*d,norm)
    j,k=divmod(b*c-a*d,norm)
    assert v==k==0
    return u,j


def sign(x):
    a,b=x
    if a==0 or b==0 or (a>0)==(b>0):
        return (a+b>0)-(a+b<0)
    norm=a*a-2*b*b
    return ((norm>0)-(norm<0))*(1 if a>0 else -1)


def compare(x,y):
    return sign(sub(x,y))


def data():
    rows=[]; knots={START,STOP}
    for p in primes_through(STOP[0]):
        if p==2:
            continue
        features=[]
        for u,v in FIELDS:
            if p in (u,v):
                features.append(None)
                continue
            e=1 if pow(p,(u-1)//2,u)==1 else -1
            f=1 if pow(p,(v-1)//2,v)==1 else -1
            d=next(j for j in (1,2,4) if pow(p,j,16) in (1,15) and (e*f)**j==1)
            features.append((p**d,e**d))
        if not any(x and x[0]<=STOP[0] for x in features):
            continue
        rows.append((p,features))
        for spec in features:
            if spec is None:
                continue
            n=spec[0]
            for k in (1,2,4):
                for t in ((8*n//k,0),(0,4*n//k),(4*n//k,0)):
                    if compare(START,t)<0 and compare(t,STOP)<0:
                        knots.add(t)
    return rows,sorted(knots,key=cmp_to_key(compare))


def endpoint(row,t):
    # Multiply each row by 8 times its largest norm. All its norms divide that norm.
    _,specs=row
    common=max(s[0] for s in specs if s)
    entries=[]
    starts=(t,(2*t[1],t[0]),scale(t,2))
    for spec in specs:
        if spec is None:
            entries.extend([ZERO]*3)
            continue
        n,e=spec
        for s in starts:
            left=sub((16*n,0),scale(s,2))
            right=sub(scale(s,4),(8*n,0))
            weight=left if compare(left,right)<=0 else right
            entries.append(scale(weight,e*(common//n)) if sign(weight)>0 else ZERO)
    return entries


def determinant(matrix):
    rows=[list(row) for row in matrix]
    previous=(1,0); polarity=1
    for k in range(5):
        pivot=next((j for j in range(k,6) if rows[j][k]!=ZERO),None)
        if pivot is None:
            return ZERO
        if pivot!=k:
            rows[k],rows[pivot]=rows[pivot],rows[k]
            polarity=-polarity
        diagonal=rows[k][k]
        for i in range(k+1,6):
            first=rows[i][k]
            for j in range(k+1,6):
                value=sub(mul(diagonal,rows[i][j]),mul(first,rows[k][j]))
                rows[i][j]=divide(value,previous)
        previous=diagonal
    return scale(rows[-1][-1],polarity)


def interpolation():
    # Invert the fixed evaluation matrix for degree-six Bernstein coefficients.
    a=[[Fraction(comb(6,k)*j**k*(6-j)**(6-k)) for k in range(7)]
       +[Fraction(j==k) for k in range(7)] for j in range(7)]
    for k in range(7):
        pivot=next(i for i in range(k,7) if a[i][k])
        a[k],a[pivot]=a[pivot],a[k]
        a[k]=[x/a[k][k] for x in a[k]]
        for i in range(7):
            if i!=k:
                a[i]=[u-a[i][k]*v for u,v in zip(a[i],a[k])]
    return [row[7:] for row in a]


def certificate(chosen,left,right,inverse):
    if len(chosen)!=6:
        return False
    a=[endpoint(row,left) for row in chosen]
    b=[endpoint(row,right) for row in chosen]
    values=[]
    for j in range(7):
        matrix=[[add(scale(x,6-j),scale(y,j)) for x,y in zip(l,h)] for l,h in zip(a,b)]
        values.append(determinant(matrix))
    signs=[]
    for row in inverse:
        coeff=tuple(sum(c*t[i] for c,t in zip(row,values)) for i in (0,1))
        signs.append(sign(coeff))
    return signs[-1]!=0 and all(s*signs[-1]>=0 for s in signs)


def verify():
    rows,knots=data()
    early=[(p,s[1][0]) for p,s in rows if s[1] and s[1][0]<=79]
    assert early==[(7,49),(79,79)],early
    inverse=interpolation()
    by_prime=dict(rows)
    path=Path(__file__).resolve().parent/"evidence/initial-minors.tsv"
    with path.open() as source:
        table=csv.DictReader(source,delimiter="\t")
        for left,right in zip(knots,knots[1:]):
            line=next(table)
            assert right==(int(line["right_a"]),int(line["right_b"]))
            primes=[int(line["p"+str(i)]) for i in range(1,7)]
            assert len(set(primes))==6
            chosen=[(p,by_prime[p]) for p in primes]
            assert certificate(chosen,left,right,inverse),(left,right)
        assert next(table,None) is None
    print(f"Exact degree-six Bernstein minors cover {len(knots)-1} knot intervals from 79*sqrt(2)/8 (open) through 1,500. At and below that transition the first two second-host readings use at most the logarithm of 7.")


if __name__ == "__main__":
    verify()
