"""Enclose two consecutive dyadic moves from finite prime and prime-power sums."""
from bisect import bisect_right
from fractions import Fraction
from math import isqrt
from pathlib import Path
from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
HOSTS = ((547, 727), (571, 599))


def ball(r):
    return arb(r.numerator) / r.denominator


def witnesses():
    result = {}
    for line in (ROOT / "evidence/witnesses.tsv").read_text().splitlines():
        code, p, q, a, b = map(int, line.split())
        assert code not in result and min(p, q, a, b) > 0
        result[code] = Fraction(p, q), Fraction(a, b)
    assert set(result) == set(range(16))
    return result


def primes(n):
    sieve = bytearray([1]) * (n+1)
    sieve[:2] = bytes(2)
    for p in range(2, isqrt(n)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes((n-p*p)//p+1)
    return [p for p in range(n+1) if sieve[p]]


def sign_table(p):
    table = [-1] * p
    table[0] = 0
    for n in range(1, p):
        table[n*n % p] = 1
    return table


class Prefix:
    """Two directed prefix sums evaluate every translate of the same tent."""
    def __init__(self, rows):
        self.ns = []
        self.a = [arb(0)]
        self.b = [arb(0)]
        for n, w in sorted(rows, key=lambda row: row[0]):
            self.ns.append(n)
            self.a.append(self.a[-1]+w)
            self.b.append(self.b[-1]+w/n)

    def value(self, z, root, norm):
        l, m, h = (bisect_right(self.ns, u) for u in (z/2, z, 2*z))
        x = ball(z).sqrt()
        a, b = self.a, self.b
        return 2/norm * (root/x*(a[m]-a[l])-x/root*(b[m]-b[l])
                         +root*x*(b[h]-b[m])-(a[h]-a[m])/(root*x))


class Host:
    def __init__(self, pair, spectral_limit, raw_limit):
        u, v = pair
        self.a = arb(2).log()
        self.root = arb(2).sqrt()
        self.c = 2*(self.root-1)**2
        self.norm = self.root*self.c
        self.constant = sum((f*q/arb.pi()).log() for q in pair for f in (1,8,16,16))
        self.constant -= 8*(arb.const_euler()+self.a)
        table = (sign_table(u), sign_table(v))
        spec, raw = [], []
        for p in primes(max(spectral_limit, raw_limit)):
            ls, rs = table[0][p%u], table[1][p%v]
            logp = arb(p).log()
            if p <= spectral_limit:
                m, n = 1, p
                while n <= spectral_limit:
                    k = 1 if p == 2 else 4 if n % 16 in (1,15) else 0
                    tr = k*(ls**m+rs**m)
                    if tr:
                        spec.append((n, tr*logp))
                    m, n = m+1, n*p
            if p not in (2,u,v):
                f = next(f for f in (1,2,4) if pow(p,f,16) in (1,15) and (ls*rs)**f == 1)
                if p**f <= raw_limit:
                    raw.append((p**f, -8*ls**f*logp))
        self.spectral, self.raw = Prefix(spec), Prefix(raw)

    def arch_series(self, x):
        ratio = min(x, 1/x)
        if ratio == 1:
            return self.a
        z = ball(ratio)
        return ((1+z)*(1+z).log()+(1-z)*(1-z).log())/(2*z.sqrt())

    def profile(self, z):
        z = max(z, 1/z)
        s = (2/ball(z)).sqrt()
        tent = 2*(s-1/s)/self.norm if z < 2 else arb(0)
        arch = self.root*(self.arch_series(2*z)+self.arch_series(z/2))-3*self.arch_series(z)
        return self.constant*tent-8*arch/self.c-self.spectral.value(z,self.root,self.norm)

    def score(self, y):
        return self.raw.value(2*y,self.root,self.norm)


def signs(code):
    return tuple(1 if code & (1 << bit) else -1 for bit in (3,2,1,0))


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 144
    rows = witnesses()
    zmax = max(2*max(4*z,1/z) for z,_ in rows.values())
    ymax = max(16*y for _,y in rows.values())
    hosts = [Host(pair,int(zmax)+1,int(ymax)+1) for pair in HOSTS]
    a = hosts[0].a
    for h in hosts:
        assert 0 < h.profile(Fraction(1)) < 180
    for word,(z,y) in rows.items():
        for h,host in enumerate(hosts):
            for j in (0,1):
                sign = signs(word)[2*j+h]
                phase = z*2**j
                start = y*2**j
                measured = (host.score(2*start)-host.score(start))/a
                limiting = (host.profile(2*phase)-host.profile(phase))/a
                assert sign*measured > arb(16)/10, (word,h,j,measured)
                assert sign*limiting > arb(91)/10, (word,h,j,limiting)
    print("All 16 pairs of consecutive joint moves certified: 64 finite moves exceed 1.6 and 64 profile moves exceed 9.1 in their specified directions.")
    print(f"Finite sums use norms at most {int(ymax)} and primitive prime powers at most {int(zmax)}. Both spectral masses are below 180.")


if __name__ == "__main__":
    main()
