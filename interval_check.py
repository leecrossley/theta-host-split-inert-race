"""Cross-check the spectral mass by real quadrature and the profile by direct interval sines."""
from pathlib import Path
from flint import arb, ctx


def mass_from_integral(a, r):
    b = 2*r*(r-1)**2

    def g(s):
        return 4*((a-abs(s)).max(0)/2).sinh()/b

    constant = sum((q/arb.pi()).log() for q in (547,4376,8752,8752,727,5816,11632,11632))
    constant += 8*(arb(3)/4).digamma()
    integral0, integral1 = arb(0), arb(0)
    pieces = 8192
    origin = g(0)
    for j in range(-16,5):
        low = a*arb(2)**j
        width = low/pieces
        for i in range(pieces):
            t = low + arb(2*i+1,1)*width/2
            weight = width*(-t/2).exp()/(2*t.sinh())
            integral0 += weight*2*(g(t)-origin)
            integral1 += weight*(g(a-t)+g(a+t))
    lip = 2*(a/2).cosh()/b
    end = 32*a
    error = lip*a/2**16 + 4*origin*(-3*end/2).exp()/((arb(3)/2)*(1-(-2*end).exp()))
    enclosure = arb(0,1)*error
    v0 = constant*origin - 8*(integral0+enclosure)
    va = -8*(integral1+enclosure)
    return (3*v0-va)/(2*a)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 144
    a, r = arb(2).log(), arb(2).sqrt()
    mass = mass_from_integral(a,r)
    assert mass < arb(38676)/100 and mass > arb(38671)/100
    head = []
    removed = arb(0)
    table = Path(__file__).parent / "evidence" / "zero-brackets.tsv"
    for line in table.read_text().splitlines():
        _, _, lo, hi = map(int, line.split())
        t = (arb(lo)/10**11).union(arb(hi)/10**11)
        s = (a*t/2).sin()
        # This identity uses only the half-angle sine.
        g = ((r-1)**2+4*r*s*s)/(2*(r-1)**2*(t*t+arb(1)/4))
        removed += 2*g*(1+s*s)/a
        head.append((t,4*g*s/a))
    missing = mass-removed
    assert missing < arb(174)/100
    ceiling = arb(19)/10

    def covers(index, depth):
        count = 1 << depth
        x = (arb(389568623) + arb(2*index+1,1)*528/(2*count))/100
        neg_head = sum((weight*(t*x).sin() for t,weight in head), arb(0))
        if missing-neg_head < ceiling:
            return 1
        if depth >= 24:
            raise ArithmeticError(f"Unresolved cell {index}/2^{depth}")
        return covers(2*index,depth+1) + covers(2*index+1,depth+1)

    cells = sum(covers(i,8) for i in range(256))
    print("Real gamma integral encloses the full mass in (386.71, 386.76); omitted mass < 1.74.")
    print(f"Direct interval sines on {cells} closed cells independently give U(x) < 1.9 throughout.")


if __name__ == "__main__":
    main()
