"""Verify the isolated moving frequency and the constants in the waiting bound."""
from fractions import Fraction

from flint import acb, acb_series, arb, ctx, dirichlet_char

if not __debug__:
    raise RuntimeError("Run these checks with assertions enabled.")

LOW, HIGH, SCALE = 3287178298086, 3287178298087, 10**12


def odd_factors():
    for p in (547, 727):
        for m, idx in ((1, 1), (8, 5), (16, 5), (16, 13)):
            q = p*m
            label = next(a for a in range(p-1, q, p) if a % m == idx % m)
            chi = dirichlet_char(q, label)
            assert chi.is_primitive() and int(chi.conductor()) == q and chi.parity() == 1
            yield q, label, chi


def main():
    ctx.prec = 112
    g = arb(LOW + HIGH, HIGH - LOW)/(2*SCALE)
    mid = arb(LOW + HIGH)/(2*SCALE)
    selected = dirichlet_char(547, 546)
    for n, direction in ((LOW, 1), (HIGH, -1)):
        v = selected.hardy_z(arb(n)/SCALE)
        assert v.imag.contains(0) and direction*v.real > 0
    print(f"547.546: opposite Hardy signs at {LOW}/{SCALE} and {HIGH}/{SCALE}")
    radius = arb(1)/1000 + arb(1)/10**11
    # Subdivision counts give strict interval margins on every closed subinterval.
    meshes = ((4,1), (4,4), (4,1), (4,4), (1,4), (1,4), (64,16), (64,4))
    for (q, lab, c), mesh in zip(odd_factors(), meshes, strict=True):
        for at_gamma, (centre, n) in enumerate(zip((arb(0), mid), mesh)):
            for j in range(n):
                t = centre + (2*j+1-n)*radius/n + arb(0, radius/n)
                s = acb(arb(1)/2, t)
                if (q,lab) == (547,546) and at_gamma:
                    derivative = acb_series([s,1],2).dirichlet_l(c)[1]
                    assert derivative.real > 1
                else:
                    assert abs(c.l_function(s)) > 0
        print(f"{q}.{lab}: both critical-line neighborhoods certified")
    b, a = arb(2).sqrt(), arb(2).log()
    v = (3-2*b*(a*g).cos())/(2*(b-1)**2*(arb(1)/4+g**2)) * 2*(a*g/2).sin()/a
    assert arb(3342843)/10**6 < v < arb(3342844)/10**6
    print("selected residue magnitude lies in (3.342843, 3.342844)")

    assert 8500*11632 < 10**8
    assert arb(1).exp() < arb(8500-4)/3000
    assert 4_000_000*(1/arb(2**40) + (-arb(999)/2).exp()) < arb(1)/1000
    assert 128*(-arb(997)/2).exp()/(1-arb(-997).exp()) < arb(1)/1000
    assert 10**6*(-arb(999)/10).exp() < arb(1)/1000
    assert (1-arb(2)**(-arb(2)/5)) > arb(1)/5
    assert 16*4*13 < 1000
    alpha = Fraction(3001,2999)
    assert alpha*81000 + 80000 + Fraction(1,3000) < 162000
    assert (Fraction(1001,1000)-alpha)*500_000_000 > 162000
    print("height 8500, spectrum, trivial zeros, prime powers and endpoint constants verified")


if __name__ == "__main__":
    main()
