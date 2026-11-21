"""Check nine sign brackets and every mixed frequency for both theta hosts."""
from itertools import product
from pathlib import Path
from flint import acb,arb,ctx,dirichlet_char

ROOT=Path(__file__).resolve().parent
SCALE=10**14
PARAMETERS=(547,727,571,599)
TERMS=12
GRID=128


def characters():
    chars=[]
    for p in PARAMETERS:
        for multiplier,frame in ((1,1),(8,5),(16,5),(16,13)):
            q=p*multiplier
            label=next(n for n in range(p-1,q,p) if n%multiplier==frame%multiplier)
            char=dirichlet_char(q,label)
            assert char.is_primitive() and char.conductor()==q and char.parity()==1
            chars.append((p,q,label,char))
    return chars


def targets():
    rows=[tuple(map(int,line.split())) for line in (ROOT/"evidence"/"brackets.tsv").read_text().splitlines()]
    assert len(rows)==9 and sum(p in (547,727) for p,_,_,_,_ in rows)==5
    return rows


def tests(rows):
    # Return (twice-scaled midpoint, twice-scaled radius, permitted conductor family).
    out=[]
    for p,_,_,u,v in rows:
        out.append((u+v,v-u,0 if p in (547,727) else 1))
    for coeff in product((-1,0,1),repeat=len(rows)):
        nz=[e for e in coeff if e]
        if len(nz)<2 or nz[0]<0:
            continue
        mid=abs(sum(s*(u+v) for s,(_,_,_,u,v) in zip(coeff,rows)))
        rad=sum(abs(s)*(v-u) for s,(_,_,_,u,v) in zip(coeff,rows))
        assert mid>rad
        out.append((mid,rad,-1))
    assert len(out)==9841
    return out


def check_character(task):
    i,rows=task
    ctx.prec=144
    p,q,label,char=characters()[i]
    from math import factorial
    checked=0
    samples={}
    polys={}
    # The error holds componentwise on the nearest twelve dyadic nodes.
    error=arb(48*q*factorial(6)**2)/32**TERMS
    for mid,rad,owner in tests(rows):
        if owner==i//8:
            continue
        t=arb(mid,rad)/(2*SCALE)
        assert abs(t)<4
        cell=(GRID*mid)//(2*SCALE)
        z=GRID*t-cell
        if z>=0 and z<=1:
            if cell not in polys:
                coeff=[]
                for n in range(cell-5,cell+7):
                    if n not in samples:
                        samples[n]=char.l_function(acb(arb(1)/2,arb(n)/GRID))
                differences=[samples[n] for n in range(cell-5,cell+7)]
                for k in range(TERMS):
                    coeff.append(differences[0]/factorial(k))
                    differences=[b-a for a,b in zip(differences,differences[1:])]
                polys[cell]=coeff
            val=acb(0)
            for k in reversed(range(TERMS)):
                val=val*(z+5-k)+polys[cell][k]
            if abs(val)>error:
                checked+=1
                continue
        # Direct ball evaluation also covers tubes meeting a cell boundary.
        assert abs(char.l_function(acb(arb(1)/2,t)))>0,(q,label,mid,rad)
        checked+=1
    return checked


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    from multiprocessing import Pool
    from os import environ
    ctx.prec=144
    rows=targets()
    labels={(q,l):ch for _,q,l,ch in characters()}
    for _,q,label,u,v in rows:
        ch=labels[q,label]
        ends=[ch.hardy_z(arb(k)/SCALE) for k in (u,v)]
        assert all(z.imag.contains(0) for z in ends) and ends[0].real*ends[1].real<0
    with Pool(min(8,int(environ.get("THETA_WORKERS","4")))) as pool:
        found=pool.map(check_character,[(i,rows) for i in range(16)])
    count=sum(found)
    assert count==9832*16+72
    print("Nine sign changes certified in the two theta hosts. All 9,832 mixed positive frequencies miss all sixteen spectra.")
    print(f"{count:,} rigorous frequency exclusions include 72 singleton tests against the other host.")


if __name__=="__main__":
    main()
