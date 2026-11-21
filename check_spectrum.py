"""Certify that every odd-family nontrivial zero of absolute ordinate at most 48 is listed."""
from hashlib import sha256
from pathlib import Path
from flint import acb, acb_series, arb, ctx, dirichlet_char

ROOT = Path(__file__).resolve().parent
TABLE_HASH = "e5483b464160fbce34313daa263f18dae8cc22035044a2053784a274e13a10f3"
FAMILY = ((547,546),(4376,1093),(8752,1093),(8752,5469),
          (727,726),(5816,1453),(11632,7269),(11632,1453))


def read_brackets():
    path = ROOT / "evidence" / "zero-brackets.tsv"
    if sha256(path.read_bytes()).hexdigest() != TABLE_HASH:
        raise ValueError("Incorrect signed bracket table")
    rows = [tuple(map(int,line.split())) for line in path.read_text().splitlines()]
    return [(q, label, arb(lo+hi, hi-lo)/(2*10**14)) for q,label,lo,hi in rows]


def check(displacement, centres, cap):
    brackets = read_brackets()
    for centre in centres:
        point = acb(arb(1)/2+displacement, centre)
        full = arb(0)
        for q, label in FAMILY:
            vals = acb_series([point,1],2).dirichlet_l(dirichlet_char(q,label))
            completed = ((q/arb.pi()).log() + ((point+1)/2).digamma())/2 + vals[1]/vals[0]
            full += completed.real/displacement
        used = sum(1/(displacement**2+(v-centre)**2) + 1/(displacement**2+(v+centre)**2)
                   for _,_,v in brackets)
        omitted = full-used
        if not (0 < omitted < cap):
            raise ArithmeticError(f"Poisson remainder inconclusive at {centre}: {omitted}")


def main():
    ctx.prec = 160
    sigma = arb(1)/2
    cap = arb(4)/25
    check(sigma, range(0,49,4), cap)
    if not cap < 1/(sigma**2+4):
        raise ArithmeticError("Insufficient covering radius")
    print("Thirteen positive Poisson remainders below 4/25 certify completeness through ordinate 48.")


if __name__ == "__main__":
    main()
