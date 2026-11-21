"""Verify both signs by real gamma quadrature, Cauchy-Schwarz and full-cell sine intervals."""
from pathlib import Path
from flint import arb, ctx


def even_moments(a,r):
    b = 2*r*(r-1)**2

    def g(s):
        return 4*((a-abs(s)).max(0)/2).sinh()/b

    constant = sum((q/arb.pi()).log() for q in (547,4376,8752,8752,727,5816,11632,11632))
    constant += 8*(arb(3)/4).digamma()
    integral0, integral1 = arb(0), arb(0)
    pieces = 16384
    origin = g(0)
    for j in range(-16,5):
        low = a*arb(2)**j
        width = low/pieces
        for i in range(pieces):
            t = low+arb(2*i+1,1)*width/2
            weight = width*(-t/2).exp()/(2*t.sinh())
            integral0 += 2*weight*(g(t)-origin)
            integral1 += weight*(g(a-t)+g(a+t))
    lip = 2*(a/2).cosh()/b
    end = 32*a
    error = lip*a/2**16 + 4*origin*(-3*end/2).exp()/((arb(3)/2)*(1-(-2*end).exp()))
    uncertainty = arb(0,1)*error
    return constant*origin-8*(integral0+uncertainty), -8*(integral1+uncertainty)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 144
    a,r = arb(2).log(),arb(2).sqrt()
    v0,va = even_moments(a,r)
    brack = Path(__file__).parent / "evidence" / "zero-brackets.tsv"
    mass0,mass1 = v0,(v0-va)/2
    terms = []
    for line in brack.read_text().splitlines():
        _,_,lo,hi = map(int,line.split())
        t = (arb(lo)/10**14).union(arb(hi)/10**14)
        s = (a*t/2).sin()
        g = ((r-1)**2+4*r*s*s)/(2*(r-1)**2*(t*t+arb(1)/4))
        mass0 -= 2*g
        mass1 -= 2*g*s*s
        terms.append((t,-4*g*s/a))
    assert mass0 > 0 and mass1 > 0
    tail = 2*(mass0*mass1).sqrt()/a
    assert tail < arb(152)/100
    print("Two real gamma integrals and Cauchy-Schwarz independently bound omitted absolute mass by 1.52.")

    def verify(left,right,sgn):
        def covers(i,depth):
            n = 1 << depth
            x = (arb(left) + arb(2*i+1,1)*(right-left)/(2*n))/1000
            value = sum((weight*(t*x).sin() for t,weight in terms),arb(0))
            if sgn*value-tail > arb(1)/20:
                return 1
            if depth >= 26:
                raise ArithmeticError(f"Inconclusive cell {i}/2^{depth} at {left}")
            return covers(2*i,depth+1)+covers(2*i+1,depth+1)
        return sum(covers(i,8) for i in range(256))

    up = verify(469417229340,469417235340,1)
    down = verify(484947092225,484947097850,-1)
    print(f"Direct true-frequency sines give U > 1/20 on {up} closed gain cells and U < -1/20 on {down} closed loss cells.")


if __name__ == "__main__":
    main()
