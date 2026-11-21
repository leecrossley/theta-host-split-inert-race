"""Check the finite witnesses, their perturbations, the buffered deadline and the geometric grid."""
from hashlib import sha256
from pathlib import Path
from flint import arb, ctx

ROOT = Path(__file__).resolve().parent


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 160
    evidence = ROOT/"evidence"
    files = {"finite-windows.tsv":"046a5c2fc45054e9159b307be8fbd6d5a86ff135a211c75cf68beb4e72f33e24",
             "deadline-check.txt":"00847fa1cf5a53a20c34042d8d54115e81d11b502eb419d2a0e345df88344194"}
    for name,digest in files.items():
        assert sha256((evidence/name).read_bytes()).hexdigest() == digest
    rows = {int(y):(int(lo),int(hi)) for y,lo,hi in
            (line.split() for line in (evidence/"finite-windows.tsv").read_text().splitlines()[1:])}
    delta, a, root = arb(1)/2000000, arb(2).log(), arb(2).sqrt()
    b = 2*root*(root-1)**2
    for y,direction,cap in ((1016187975,1,14),(2893002095,-1,32)):
        low,high = rows[y]
        signed = arb(low if direction == 1 else -high)/10000
        slope = arb(0)
        for z in (y,2*y):
            k = 0
            while not 4*z*delta.exp() < 2**k:
                k += 1
            slope += 20*2**(k+1)/(b*(z*(-delta).exp()).sqrt())
        assert slope*delta < cap and signed-cap > arb(19735)/10000
        assert y*(-delta).exp() > 10**9 and 2*y*delta.exp() < 10**18
    e58 = 2000000*(-arb(58)/4).exp() + 128*(-arb(56)/2).exp()/(1-(-arb(56)).exp()) + 22000*(-3*arb(58)/20).exp()
    assert e58 < 6 and 20-2*e58-arb(141)/10000 > arb(19735)/10000
    rows_mid = ((9*arb(10).log(),55,(6050,5840)),(55,80,(8366,8800)),(80,100,(10910,11144)))
    for bot,top,pair in rows_mid:
        for y in pair:
            x = arb(y)/100
            assert x-delta > 58 and x-delta-3*a/2 > top and x+delta-a/2 < bot+18*arb(10).log()
    gamma_lower = arb(3287)/1000
    width = 18*arb(10).log()-a
    assert 2*delta + arb(155)/4 + 2*arb.pi()/gamma_lower < width
    assert int((width*10**6).floor().unique_fmpz())+1 == 40753385
    print("Finite prime counts and piecewise slope bounds preserve both signs for every H through 10^9.")
    print("All six intermediate intervals and the expanded large-height filter fit strictly within the deadline.")
    print("A log grid of step 1/1000000 has at most 40753385 eligible starts at any height.")


if __name__ == "__main__":
    main()
