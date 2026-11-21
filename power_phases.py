"""Certify the positive kernel mass and the low ordinates used in the finite power test."""
from hashlib import sha256
from pathlib import Path
import sys
from flint import arb, ctx

EVIDENCE = Path(__file__).resolve().parent / "evidence"
TABLE_HASH = "e5483b464160fbce34313daa263f18dae8cc22035044a2053784a274e13a10f3"
CONDUCTORS = (547, 4376, 8752, 8752, 727, 5816, 11632, 11632)
LIMIT = 10**9


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 160
    table = EVIDENCE / "zero-brackets.tsv"
    assert sha256(table.read_bytes()).hexdigest() == TABLE_HASH
    a = arb(2).log()
    root = arb(2).sqrt()
    c = 2 * (root - 1)**2
    g0 = 4 * (a / 2).sinh() / (root * c)
    conductor = sum((q / arb.pi()).log() for q in CONDUCTORS)
    constant = conductor - 8 * (arb.const_euler() + a)
    half = arb(1) / 2
    moment = ((1+half)*(1+half).log() + (1-half)*(1-half).log()) / (2*half.sqrt())
    full = constant*g0 - 8*(2*root*moment - 3*a)/c
    selected = arb(0)
    phases = []
    low_brackets = []
    for row, line in enumerate(table.read_text().splitlines(), start=1):
        _, _, lo, hi = map(int, line.split())
        v = arb(lo+hi, hi-lo) / (2*10**14)
        kernel = (3 - 2*root*(a*v).cos()) / (c*(v*v + arb(1)/4))
        selected += 2*kernel
        if kernel > arb(7)/5:
            theta = a*v / (2*arb.pi())
            p = int((theta.mid()*2**64).floor().unique_fmpz())
            assert 0 < p < 2**64
            assert LIMIT*abs(theta - arb(p)/2**64) < arb(1)/1000
            phases.append((row, p))
            low_brackets.append((lo, hi))
    low_brackets.sort()
    assert all(left[1] < right[0] for left, right in zip(low_brackets, low_brackets[1:]))
    assert len(phases) == 25 and 0 < full-selected < arb(647)/1000 and selected < 175
    assert 11*arb(1)/64 - arb(1)/1000 > arb(1)/6
    assert arb(7)/5 > arb(19735)/10000*a > arb(4)/3
    # Sum all dyadic prime-only errors from the first log centre onward.
    x = arb(1000)
    summed = a*(2*10**6*(-x/4).exp()/(1-(-a/4).exp())
                + 128*(-(x-2)/2).exp()/((1-(-(x-2)).exp())*(1-(-a/2).exp()))
                + 22000*(-3*x/20).exp()/(1-(-3*a/20).exp()))
    assert summed < arb(10)**-19
    assert 2*arb.pi()*175/32768 + 2*arb(647)/1000 + summed < arb(4)/3
    text = "".join(f"{row}\t{p}\n" for row,p in phases)
    if len(sys.argv) == 2 and sys.argv[1] == "--print-table":
        print(text, end="")
        return
    assert text == (EVIDENCE / "low-phases.tsv").read_text()
    print("6000 signed pairs have even kernel mass below 175; omitted positive mass is below 0.647.")
    print("25 disjoint low ordinates have kernel > 7/5 and certified 64-bit phases accurate to 1/1000 through 10^9.")
    print("The sum of all dyadic prime-only errors is below 10^-19 from log centre 1000.")
    print("The 32768-arc partition certifies < 4/3, strictly below 1.9735 log 2.")


if __name__ == "__main__":
    main()
