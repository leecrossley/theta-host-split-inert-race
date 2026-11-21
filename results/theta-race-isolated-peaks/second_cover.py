"""Independently certify a staggered cover using a rational lower Taylor bound."""
from flint import arb, ctx
from inputs import brackets, run_scan


def main():
    ctx.prec = 192
    a, r, tau = arb(2).log(), arb(2).sqrt(), 2*arb.pi()
    # All 128 lowest ordinates, irrespective of amplitude.
    low = sorted(brackets(), key=lambda row: row[2])[:128]
    modes = []
    for _, _, l, u in low:
        t = (arb(l)/10**14).union(arb(u)/10**14)
        s = (a*t/2).sin()
        g = ((r-1)**2+4*r*s**2)/(2*(r-1)**2*(arb(1)/4+t**2))
        modes.append((t, g))
    # Interleave different parts of the spectrum for the early exit in the integer scan.
    modes = [m for block in zip(modes[:32], modes[32:64], modes[64:96], modes[96:]) for m in block]
    modulus = 2**61-1
    lines = []
    for t, g in modes:
        f = t/(256*tau)
        n = int((f.mid()*modulus+arb(1)/2).floor().unique_fmpz())
        error = arb(1)/250
        assert 0 < n < modulus and 256*(10**9+2)*abs(f-arb(n)/modulus) < error
        lines.append(str(n))
        for depth in range(6):
            radius = arb(3)/(8*2**depth)
            weights = []
            for k in range(129):
                # The Mersenne denominator is smaller than 2^61, so k/256 is a lower bound.
                d = (arb(k)/256-t*radius/tau-error).lower().max(0)
                assert 0 <= d and d < arb(1)/2
                z = (tau*d)**2
                # Lower Taylor bound for 1-cos at angles between 0 and pi.
                p = z*(arb(1)/2-z*(arb(1)/24-z*(arb(1)/720-z/40320)))
                weights.append(max(0, int((200*g*p).lower().floor().unique_fmpz())))
            lines.append(" ".join(map(str, weights)))
    print("128 low ordinates and 99072 eighth-degree lower bounds pass at 192 bits.", flush=True)
    run_scan("scan_staggered.c", lines)


if __name__ == "__main__":
    main()
