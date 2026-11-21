"""Test the two-parameter asymptotic directly against the whole compact formula."""
from fractions import Fraction
from flint import arb, ctx
from width_kernel import PrimeTable, Window, ball, check_inputs, trace

CENTRES = ((Fraction(3,2),None),(Fraction(2), (2,1)),(Fraction(4),(2,2)),
           (Fraction(6),None),(Fraction(7),(7,1)),(Fraction(31),(31,1)),
           (Fraction(47),(47,1)),(Fraction(49),(7,2)),(Fraction(961),(31,2)))


def main():
    check_inputs()
    ctx.prec = 384
    table = PrimeTable(1000)
    cases = 0
    for k in (6,12,18):
        unit = Fraction(10**k+1,10**k)
        z = unit**2
        h = ball(z).log()
        for centre, prime in CENTRES:
            for j in range(-4,5):
                b = centre*unit**j
                window = Window(b, table)
                c,bp = window.c,ball(b)+1
                curv, _, _ = window.measurements(z)
                remainder = c*curv/(2*h)-4*bp*(1/h).log()-bp*(8+window.constant)/2
                target = arb(0) if prime is None else (arb(prime[0]).log()*trace(*prime)*
                                                        ball(max(0,1-abs(Fraction(j,2)))))
                assert abs(remainder-target) < 100000*h
                cases += 1
    print(f"{cases} complete-profile cusp evaluations agree with the triangular limit within 100000 times the displacement.")
    print("Nine centres include empty, invisible, positive and negative prime-power weights. Nine detunings include both edges; three step sizes reach 2*log(1+10^-18).")


if __name__ == "__main__":
    main()
