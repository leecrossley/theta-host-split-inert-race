"""Interval evaluation of physical and limiting dyadic prime traces on compact log intervals."""

from bisect import bisect_left, bisect_right
from math import isqrt

from flint import arb

HOSTS = ((547, 727), (571, 599))


def primes_through(limit):
    flags = bytearray([1])*(limit+1)
    flags[:2] = bytes(2)
    for p in range(2, isqrt(limit)+1):
        if flags[p]:
            flags[p*p::p] = bytes(1+(limit-p*p)//p)
    return (p for p, flag in enumerate(flags) if flag)


def legendre_table(q):
    squares = {n*n % q for n in range(1, q)}
    return [0] + [1 if n in squares else -1 for n in range(1, q)]


def integer_floor(x):
    value = x.floor().unique_fmpz()
    if value is None:
        raise ArithmeticError(f"Cannot place the positive scale {x} between integers")
    return int(value)


class Trace:
    def __init__(self, host, limit, physical=False):
        self.physical = physical
        self.a = arb(2).log()
        self.root = arb(2).sqrt()
        self.c = 2*(self.root-1)**2
        self.B = self.root*self.c
        self.limit = limit
        p0, p1 = host
        factors = [m*p for p in host for m in (1, 8, 16, 16)]
        self.conductor_term = sum((arb(q)/arb.pi()).log() for q in factors)
        self.conductor_term -= 8*(arb.const_euler()+self.a)
        residues = [legendre_table(q) for q in host]
        rows = []
        for p in primes_through(limit):
            u, v = residues[0][p % p0], residues[1][p % p1]
            if physical:
                if p in (2, p0, p1):
                    continue
                degree = next(f for f in (1, 2, 4) if pow(p, f, 16) in (1, 15)
                              and (u*v)**f == 1)
                norm = p**degree
                if norm <= limit:
                    rows.append((norm, p, 8*u**degree))
                continue
            n, k = p, 1
            while n <= limit:
                twist = (1 if p == 2 else 4 if n % 16 in (1, 15) else 0)
                coefficient = twist*(u**k+v**k)
                if coefficient:
                    rows.append((n, p, coefficient))
                n *= p
                k += 1
        self.norms = []
        self.plain, self.divided, self.absolute = [arb(0)], [arb(0)], [arb(0)]
        for n, p, coefficient in sorted(rows):
            w = coefficient*arb(p).log()
            self.norms.append(n)
            self.plain.append(self.plain[-1]+w)
            self.divided.append(self.divided[-1]+w/n)
            self.absolute.append(self.absolute[-1]+abs(w))

    def gamma_series(self, q):
        if q == 1:
            return self.a
        if q < arb(1)/4:
            value = sum((q**(2*k)/((2*k+1)*(2*k+2)) for k in range(7)), arb(0))
            rest = q**14/(15*16*(1-q*q))
            return q*q.sqrt()*(value + arb(0,1)*rest)
        if not q < 1:
            raise ArithmeticError(f"Unresolved gamma kernel at {q}")
        return ((1+q)*(1+q).log()+(1-q)*(1-q).log())/(2*q.sqrt())

    def gamma(self, z):
        def ratio(x):
            if x == 1:
                return arb(1)
            if x < 1:
                return x
            if x > 1:
                return 1/x
            raise ArithmeticError(f"Unresolved reflection at scale {x}")
        q = ratio(z)
        tent = arb(0)
        if q > arb(1)/2:
            u = (2*q).sqrt()
            tent = 2*(u-1/u)/self.B
        elif not q <= arb(1)/2:
            raise ArithmeticError(f"Unresolved tent boundary at scale {z}")
        terms = self.root*(self.gamma_series(ratio(2*z))+self.gamma_series(ratio(z/2)))
        return self.conductor_term*tent-8/self.c*(terms-3*self.gamma_series(q))

    def shell(self, z, d, k):
        scale = z*k
        lo = scale*(-d).exp()
        hi = scale*d.exp()
        # Possible knots at an integer endpoint are included on both sides.
        low = bisect_left(self.norms, integer_floor(lo))
        high = bisect_right(self.norms, integer_floor(hi)+1)
        return (self.absolute[high]-self.absolute[low])/lo.sqrt()

    def level(self, z, d):
        """Enclose V(log(z)+u) for |u| <= d, with z >= 1."""
        if not (z >= 1 and 2*z*d.exp() < self.limit):
            raise ValueError("Finite prime-power table does not cover this support.")
        low, mid, high = [bisect_right(self.norms, integer_floor(v)) for v in (z/2,z,2*z)]
        p, q, r = self.plain, self.divided, self.root
        first = 2/self.B*(r*(p[mid]-p[low])-(p[high]-p[mid])/r)/z.sqrt()
        second = 2/self.B*(r*(q[high]-q[mid])-(q[mid]-q[low])/r)*z.sqrt()
        prime = first+second
        if d == 0:
            return (arb(-1) if self.physical else self.gamma(z))-prime
        # Away from knots the extrapolation is first*exp(-u/2)+second*exp(u/2).
        # A slope jump contributes 2*jump*sinh((u-u_k)/2) after its knot.
        jumps = 2/self.B*(self.shell(z,d,arb(1)/2)+self.shell(z,d,2))
        jumps += 4*(self.a/2).cosh()/self.B*self.shell(z,d,1)
        error = abs(prime)*((d/2).cosh()-1) + (abs(second-first)+2*jumps)*(d/2).sinh()
        zr = z*arb(0,d).exp()
        return (arb(-1) if self.physical else self.gamma(zr))-prime+arb(0,1)*error

    def pair(self, x, d=0):
        """Positive x gives the even level and the odd normalized move."""
        d = arb(d)
        if x == 0:
            assert d == 0
            return self.level(arb(1), d), arb(0)
        z = x.exp()
        assert z > self.root
        level = self.level(z, d)
        move = (self.level(z*self.root, d)-self.level(z/self.root, d))/self.a
        return level, move


def vector(traces, x, d=0):
    pairs = [trace.pair(x, d) for trace in traces]
    return [pairs[0][0], pairs[1][0], pairs[0][1], pairs[1][1]]
