"""Independent secant-error check using Conrey characters and polynomial splitting."""

from bisect import bisect_right
from functools import lru_cache
from itertools import product
from math import isqrt, lcm
from pathlib import Path

from flint import arb, arb_mat, ctx, dirichlet_char, nmod_poly

ROOT = Path(__file__).resolve().parent
FIELDS = ((547,727),(571,599))
SCALE = 16000


def primelist(n):
    # Odd-index sieve. Index i represents the integer 2*i+1.
    flags = bytearray([1])*((n+1)//2)
    flags[0] = 0
    for i in range(1,(isqrt(n)-1)//2+1):
        if flags[i]:
            p = 2*i+1
            start = p*p//2
            flags[start::p] = bytes(len(range(start,len(flags),p)))
    yield 2
    yield from (2*i+1 for i in range(1,len(flags)) if flags[i])


def primitive_rows(conditions, limit):
    characters = []
    for p in conditions:
        for m, frame in ((1,1),(8,5),(16,5),(16,13)):
            q = p*m
            label, = [n for n in range(p-1,q,p) if n%m == frame%m]
            chi = dirichlet_char(q,label)
            assert chi.order() in (2,4) and chi.is_primitive() and chi.parity() == 1
            values = []
            for n in range(q):
                z = chi(n)
                re, im = z.real.unique_fmpz(), z.imag.unique_fmpz()
                assert re is not None and im is not None
                values.append((int(re),int(im)))
            characters.append((q,values))
    for p in primelist(limit):
        n = p
        while n <= limit:
            components = [v[n%q] for q,v in characters]
            assert sum(im for _,im in components) == 0
            coefficient = -sum(re for re,_ in components)
            if coefficient:
                yield n, p, coefficient
            n *= p


def degree(polynomial, p):
    _, factors = nmod_poly(polynomial,p).factor()
    assert all(mult == 1 for _,mult in factors)
    d, = {int(f.degree()) for f,_ in factors}
    return d


def physical_rows(conditions, limit):
    u, v = conditions
    for p in primelist(limit):
        if p in (2,u,v):
            continue
        # In the real cyclotomic quartic, degree one requires p = +1 or -1 modulo 16.
        if p*p > limit and p%16 not in (1,15):
            continue
        f = lcm(degree([2,0,-4,0,1],p),degree([-u*v,0,1],p))
        if p**f <= limit:
            g = lcm(f,degree([u,0,1],p))
            assert g in (f,2*f)
            yield p**f, p, -8 if g == f else 8


class Reference:
    def __init__(self, host, physical, limit=3_000_000):
        self.physical = physical
        self.limit = limit
        rows = (physical_rows if physical else primitive_rows)(host,self.limit)
        self.n, self.p, self.q = [], [arb(0)], [arb(0)]
        for n,p,c in sorted(rows):
            w = c*arb(p).log()
            self.n.append(n)
            self.p.append(self.p[-1]+w)
            self.q.append(self.q[-1]+w/n)
        self.a = arb(2).log()
        self.r = arb(2).sqrt()
        self.B = 2*self.r*(self.r-1)**2
        self.correction = sum((arb(p*m)/arb.pi()).log() for p in host for m in (1,8,16,16))
        self.correction -= 8*(self.a+arb.const_euler())

    def gamma(self,z):
        def kernel(t):
            q = t if t <= 1 else 1/t
            if q == 1:
                return self.a
            return ((1+q)*(1+q).log()+(1-q)*(1-q).log())/(2*q.sqrt())
        x = z.log()
        tent = max(arb(0),self.a-abs(x))/2
        arch = 3*kernel(z)-self.r*(kernel(2*z)+kernel(z/2))
        return 4/self.B*tent.sinh()*self.correction + 8*self.r/self.B*arch

    def level(self,z):
        assert z >= 1 and 2*z <= self.limit
        idx = []
        for k in (z/2,z,2*z):
            n = k.floor().unique_fmpz()
            assert n is not None
            idx.append(bisect_right(self.n,int(n)))
        lo,mid,hi = idx
        p,q = self.p,self.q
        left = 2*(p[mid]-p[lo])-z*(q[mid]-q[lo])
        right = 2*z*(q[hi]-q[mid])-(p[hi]-p[mid])
        bias = -1 if self.physical else self.gamma(z)
        return bias+(left+right)/(self.B*(z/2).sqrt())

    def physical_at(self,k):
        root = arb(2).sqrt()
        raw = self.level(2*k*root)
        change = (self.level(arb(4*k))-self.level(arb(2*k)))/self.a
        return raw,change

    def observation(self, x):
        if x == 0:
            return self.level(arb(1)),arb(0)
        z = abs(x).exp()
        d = (self.level(z*self.r)-self.level(z/self.r))/self.a
        return self.level(z), (d if x>0 else -d)


def secant_errors(left,right):
    """Componentwise interpolation errors over a grid of step 1/SCALE."""
    h = arb(1)/SCALE
    r = arb(2).sqrt()
    B = 2*r*(r-1)**2
    a = arb(2).log()
    small = arb(min(abs(left),abs(right)))/1000
    large = arb(max(abs(left),abs(right)))/1000
    assert left*right > 0 and small-a/2-a > 4.5
    def bound(mult):
        z0,z1 = small.exp()*mult,large.exp()*mult
        logweight = 8*(2*z1).log()
        max_prime = (2*z1+1)*logweight*4*(a/2).sinh()/(B*(z0/2).sqrt())
        jumps = arb(0)
        for k,c in ((arb(1)/2,2/B),(arb(1),4*(a/2).cosh()/B),(arb(2),2/B)):
            jumps += c*logweight*(1+k*z1*h)/(k*z0).sqrt()
        # Smooth part: prime'' = prime/4. On these arcs |gamma''| < 100.
        return (max_prime/4+100)*h*h/8 + jumps*h/4
    level = bound(arb(1))
    move = (bound(1/r)+bound(r))/a
    return (level,level,move,move)


def evaluate(physical):
    with (ROOT/"evidence"/"three-start-charts.tsv").open() as fh:
        next(fh)
        charts = [tuple(map(int,s.split())) for s in fh]
    rows_by_chart = []
    refs = [Reference(host,physical) for host in FIELDS]
    @lru_cache(maxsize=None)
    def point(k):
        x = arb(k)/SCALE
        s0,d0 = refs[0].observation(x)
        s1,d1 = refs[1].observation(x)
        return (s0,s1,d0,d1,arb(1))
    assert SCALE%1000 == 0
    factor = SCALE//1000
    identity = arb_mat([[int(i==j) for j in range(5)] for i in range(5)])
    cells = [(i,j) for i in range(-10,10) for j in range(-8,8)]
    for chart in charts:
        vertices = [list(point(t*factor)) for t in chart]
        if physical:
            target = (arb(chart[0])/1000).exp()/(2*arb(2).sqrt())
            k = int((target+arb(1)/2).floor().unique_fmpz())
            assert abs((arb(k)/target).log()) < arb(1)/2000
            v,d = refs[0].physical_at(k)
            u,e = refs[1].physical_at(k)
            vertices[0] = [v,u,d,e,arb(1)]
        vertices = arb_mat(vertices)
        basis = vertices.solve(identity)
        columns = [tuple(basis[i,j] for i in range(5)) for j in range(5)]
        ceiling = [arb(0) for _ in range(5)]
        for u,v,which in ((chart[1],chart[2],(0,3,4)),(chart[3],chart[4],(0,1,2))):
            err = secant_errors(u,v)
            nodes = [point(k) for k in range(min(u,v)*factor,max(u,v)*factor+1)]
            for j,c in enumerate(columns):
                vals = (sum(a*b for a,b in zip(z,c)).upper() for z in nodes)
                term = max(vals) + sum(abs(b)*e for b,e in zip(c,err))
                assert sum(abs(x) for x in c[:4]) < arb(3)/20 and term < arb(3)/2
                if j in which:
                    ceiling[j] = ceiling[j].max(term)
        row = []
        for i,j in cells:
            low = arb(1000)
            for s,t,m,n in product((i,i+1),(j,j+1),(-1,1),(-1,1)):
                target = (s,t,m*arb(7)/2,n*arb(11)/5,1)
                for c,top in zip(columns,ceiling):
                    value = sum(q*v for q,v in zip(c,target))-top
                    low = low.min(value)
            row.append(low.lower())
        rows_by_chart.append(row)
    print("Polynomial splitting" if physical else "Full Conrey table", "and independent secant errors:")
    print(len(charts),"configurations, 320 closed rectangular cells, 16 corners per cell")
    return rows_by_chart


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 148
    finite = evaluate(True)
    limiting = evaluate(False)
    for k in range(320):
        assert any(f[k] > arb(1)/250 and l[k] > arb(1)/250 for f,l in zip(finite,limiting))
    print("Same chart covers both traces on each cell with clearance > 1/250.")
    print("Every arc coordinate < 3/2 and every row norm < 3/20 by the independent secant bound.")
