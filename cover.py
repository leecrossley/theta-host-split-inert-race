"""Construct certified whole-cell cosine bounds on a dyadic grid with binary phases."""
from flint import arb, ctx
from inputs import brackets, run_scan


def main():
    ctx.prec = 160
    a, root, tau = arb(2).log(), arb(2).sqrt(), 2*arb.pi()
    modes = []
    for _, _, lo, hi in brackets():
        t = arb(lo+hi, hi-lo)/(2*10**14)
        g = (3-2*root*(a*t).cos())/(2*(root-1)**2*(t*t+arb(1)/4))
        # The ranking only chooses a subset. All subsequent inequalities use full balls.
        modes.append((float(g.mid()), t, g))
    modes.sort(key=lambda row: row[0], reverse=True)
    lines = []
    for _, t, g in modes[:100]:
        phase = t/(128*tau)
        p = int((phase.mid()*2**64).floor().unique_fmpz())
        error = arb(1)/10000
        assert 0 < p < 2**64 and 128*(10**9+1)*abs(phase-arb(p)/2**64) < error
        lines.append(str(p))
        for level in range(7):
            radius = arb(1)/2**(level+1)
            costs = []
            for k in range(256):
                lower = (arb(min(k, 255-k))/256-radius*t/tau-error).lower().max(0)
                assert 0 <= lower and lower < arb(1)/2
                charge = int((200*g*(1-(tau*lower).cos())).lower().floor().unique_fmpz())
                costs.append(max(0, charge))
            lines.append(" ".join(map(str, costs)))
    print("100 true-ordinate binary phases and all 179200 directed cosine bounds pass.", flush=True)
    run_scan("scan_cells.c", lines)


if __name__ == "__main__":
    main()
