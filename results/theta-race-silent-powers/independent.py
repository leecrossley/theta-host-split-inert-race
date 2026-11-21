"""Check the omitted even mass by gamma quadrature and rescan in independent decimal arithmetic."""
import os
from pathlib import Path
import subprocess
import tempfile
from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
CONDUCTORS = (547, 4376, 8752, 8752, 727, 5816, 11632, 11632)


def integral_moment(a, root):
    b = 2*root*(root-1)**2
    origin = 4*(a/2).sinh()/b
    constant = sum((q/arb.pi()).log() for q in CONDUCTORS) + 8*(arb(3)/4).digamma()
    integral = arb(0)
    pieces = 65536
    for j in range(-20, 5):
        low = a*arb(2)**j
        width = low/pieces
        for i in range(pieces):
            t = low + arb(2*i+1, 1)*width/2
            g = 4*((a-t).max(0)/2).sinh()/b
            integral += width*(-t/2).exp()*(g-origin)/t.sinh()
    lip = 2*(a/2).cosh()/b
    end = 32*a
    missed = lip*a/2**20 + 4*origin*(-3*end/2).exp()/((arb(3)/2)*(1-(-2*end).exp()))
    return constant*origin - 8*(integral + arb(0,1)*missed)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 176
    a, root = arb(2).log(), arb(2).sqrt()
    v0 = integral_moment(a, root)
    selected = arb(0)
    phases = []
    selected_intervals = []
    for line in (ROOT/"evidence"/"zero-brackets.tsv").read_text().splitlines():
        _, _, lo, hi = map(int, line.split())
        v = (arb(lo)/10**14).union(arb(hi)/10**14)
        s = (a*v/2).sin()
        weight = ((root-1)**2 + 4*root*s*s)/(2*(root-1)**2*(v*v+arb(1)/4))
        selected += 2*weight
        if weight > arb(7)/5:
            theta = a*v/(2*arb.pi())
            p = int((theta.mid()*10**16 + arb(1)/2).floor().unique_fmpz())
            assert 10**9*abs(theta-arb(p)/10**16) < arb(1)/1000
            phases.append(p)
            assert all(hi < old_lo or old_hi < lo for old_lo,old_hi in selected_intervals)
            selected_intervals.append((lo,hi))
    assert 0 < v0-selected < arb(649)/1000 and selected < 175
    assert 2*arb.pi()*175/32768 + 2*arb(649)/1000 + arb(10)**-19 < arb(4)/3
    assert len(phases) == 25 and arb(17)/100 - arb(1)/1000 > arb(1)/6
    print("Independent real gamma integral gives omitted mass < 0.649, sufficient for 4/3.", flush=True)
    with tempfile.TemporaryDirectory() as directory:
        binary = str(Path(directory)/"decimal_check")
        subprocess.run([os.environ.get("CC","cc"), "-std=c11", "-O3", "-Wall", "-Wextra", "-Werror",
                        str(ROOT/"decimal_check.c"), "-o", binary], check=True)
        subprocess.run([binary], input="\n".join(map(str,phases))+"\n", text=True, check=True)


if __name__ == "__main__":
    main()
