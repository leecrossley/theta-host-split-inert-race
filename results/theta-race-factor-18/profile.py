"""Enclose the odd zero sum on a full log interval of length 11/5."""
from flint import arb, ctx

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")


def segment(left, right, denominator):
    return arb(left + right, right - left)/(2*denominator)


def prime_powers():
    # Each row is n, its underlying rational prime, and the coefficient divided by log p.
    rows = []
    for p in range(2,72):
        if any(p%d == 0 for d in range(2,p) if d*d <= p):
            continue
        signs = [1 if pow(p,(q-1)//2,q) == 1 else -1 for q in (547,727)]
        n, k = p, 1
        while n <= 71:
            multiplicity = sum(s**k for s in signs)
            if p != 2:
                multiplicity *= 4 if pow(p,k,16) in (1,15) else 0
            if multiplicity:
                rows.append((n,p,-multiplicity))
            n, k = n*p, k+1
    return sorted(rows)


class Profile:
    def __init__(self):
        self.a = arb(2).log()
        self.b = arb(2).sqrt()
        self.drift = 2*self.b*(self.b-1)**2
        self.rows = [(n,c*arb(p).log()) for n,p,c in prime_powers()]
        assert prime_powers() == [(4,2,-2),(16,2,-2),(31,31,8),(47,47,-8),(49,7,-8),(64,2,-2)]
        assert (arb(3214)/1000+3*self.a/2).exp() < 71
        assert (-arb(1014)/1000+3*self.a/2).exp() < 2
        conductors = [p*m for p in (547,727) for m in (1,8,16,16)]
        self.constant = sum((q/arb.pi()).log() for q in conductors) + 8*(arb(3)/4).digamma()

    def j(self,v):
        def g(w):
            return 4/self.drift*((self.a-abs(w)).max(0)/2).sinh()
        return (g(v+self.a/2)-g(v-self.a/2))/self.a

    def raw(self,x):
        def score(y):
            val = arb(0)
            for n,c in self.rows:
                z = y/n
                w = (2*(1-z)).min(4*z-1).max(0)
                val += c*w
            return val/(self.drift*y.sqrt())
        y = x.exp()/(2*self.b)
        return (score(2*y)-score(y))/self.a


def kink_neighborhood(prof):
    # Kernel K(t)=e^(-t/2)/(2 sinh t); directed rectangle sums for the integral.
    grid = []
    for start, stop, step in ((1,128,1),(128,1000,2),(1000,3000,5),(3000,8000,20)):
        for lo in range(start,stop,step):
            hi = min(stop,lo+step)
            t = segment(lo,hi,1000)
            weight = (hi-lo)*(-t/2).exp()/(2000*t.sinh())
            grid.append((t,weight))
    a,b = prof.a,prof.drift
    assert 4*(a/2).cosh()/(b*a) < 14  # Lipschitz constant of j
    assert 8*(a/2).sinh()/(b*a) < 9
    tail = 36*arb(-12).exp()/((arb(3)/2)*(1-arb(-16).exp()))
    allowance = 8*(arb(14)/1000+tail)
    for left in range(2028,2080):
        x = segment(left,left+1,2000)
        jx = prof.j(x)
        integral = sum(w*(prof.j(x-t)+prof.j(x+t)-2*jx) for t,w in grid)
        lower = prof.raw(x)+prof.constant*jx-8*integral-allowance
        assert lower > -arb(9)/5, (left,lower)
    print("52 cells through the kernel breakpoint: odd zero sum > -9/5")


def positive_trivial_tail(prof):
    u = arb(3)/2
    a,b = prof.a,prof.b
    G = (2*b*(a*u).cosh()-3)/(2*(b-1)**2*(u*u-arb(1)/4))
    first = 8*G*2*(a*u/2).sinh()/a
    assert first > 0 and 3*a/2 < arb(1040)/1000
    for left in range(2080,6428):
        x = segment(left,left+1,2000)
        assert prof.raw(x) + first*(-u*x).exp() > -arb(9)/5, left
    print("4348 further cells: primes and the first trivial-zero term alone exceed -9/5")


def main():
    ctx.prec=112
    prof=Profile()
    kink_neighborhood(prof)
    positive_trivial_tail(prof)
    assert (arb(11)/5).exp()*2 > arb(1805)/100
    print("certified profile [1.014, 3.214]; twice exp(11/5) > 18.05")


if __name__ == "__main__":
    main()
