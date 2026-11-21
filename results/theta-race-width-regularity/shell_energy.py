"""Certify a finite mean-square lower bound for a balanced two-width mixture."""
from hashlib import sha256
from pathlib import Path
from flint import arb, ctx, dirichlet_char

TABLE = Path(__file__).resolve().parent/"evidence"/"selected-brackets.tsv"
SCALE = 10**14
SHELLS = (1,2,4,8,16,32)
COUNTS = (8,20,42,91,197,423)
FLOORS = (79,51,73,35,31,28)
FAMILY = {(547,546),(4376,1093),(8752,1093),(8752,5469),
          (727,726),(5816,1453),(11632,7269),(11632,1453)}


def brackets():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    data = TABLE.read_bytes()
    assert sha256(data).hexdigest() == "9b0e17cf840efa27a7e7f3444233d1212b4225edf8112daa38e301dcabf9bc1d"
    return [tuple(map(int,line.split())) for line in data.splitlines()]


def main():
    ctx.prec = 144
    a, r = arb(2).log(), arb(2).sqrt()
    totals, counts = [arb(0) for _ in SHELLS], [0 for _ in SHELLS]
    chars = {key:dirichlet_char(*key) for key in FAMILY}
    prev = {}
    low = 0
    for q,k,lo,hi in brackets():
        assert prev.get((q,k),0) < lo < hi
        prev[q,k] = hi
        char = chars[q,k]
        assert char.is_primitive() and char.parity() == 1 and char.conductor() == q
        left,right = char.hardy_z(arb(lo)/SCALE),char.hardy_z(arb(hi)/SCALE)
        assert left.imag.contains(0) and right.imag.contains(0) and left.real*right.real < 0
        t = arb(lo+hi,hi-lo)/(2*SCALE)
        if 0 < t < arb(1)/20:
            assert (2*t).cos() > arb(99)/100
            low += 1
            continue
        for i,n in enumerate(SHELLS):
            if n < t < 2*n:
                f = -2*r*(a*t).cos()+arb(12)/5*(2*a*t).cos()
                e = 2*f*(1-(t/n).cos())/(t*t+arb(1)/4)
                totals[i] += 2*e*e
                counts[i] += 1
                break
        else:
            raise ValueError("Straddled shell endpoint")
    assert tuple(counts)==COUNTS and low==1 and set(prev)==FAMILY
    for n,s,b in zip(SHELLS,totals,FLOORS):
        assert b < s*n**3/arb(2*n).log() < b+1
    print("782 true-ordinate Hardy sign changes rechecked for eight primitive characters (1564 endpoints).")
    print("At h=1/T for T=1,2,4,8,16,32 the corresponding T^3/log(2T)-scaled mean squares exceed 79,51,73,35,31,28.")
    print("The positive cosine coefficient of the smooth example exceeds 1/10 at a certified zero below 1/20.")


if __name__ == "__main__":
    main()
