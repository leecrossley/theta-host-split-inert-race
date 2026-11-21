"""Generate the example curvature enclosures from complete finite prime-power sums."""
from bisect import bisect_left, bisect_right
from fractions import Fraction
from math import isqrt
import sys
from flint import arb, ctx
from curvature import CASES, DEN, EVIDENCE, GammaTerms, ball, check_inputs, locations


def sieve(bound):
    primes = bytearray(b"\1")*(bound+1)
    primes[:2] = b"\0\0"
    for p in range(2, isqrt(bound)+1):
        if primes[p]:
            k = p*p
            primes[k::p] = bytes((bound-k)//p+1)
    return primes


def quadratic(p, q):
    v = pow(p, (q-1)//2, q)
    return v if v <= 1 else -1


def prime_records(bound):
    primality = sieve(bound)
    rows = []
    for p, isprime in enumerate(primality):
        if not isprime:
            continue
        first = quadratic(p, 547)
        second = quadratic(p, 727)
        n, m = p, 1
        while n <= bound:
            if p == 2:
                mult = first**m + second**m
            elif pow(p, m, 16) in (1, 15):
                mult = 4*(first**m + second**m)
            else:
                mult = 0
            if mult:
                rows.append((n, p, mult))
            n *= p
            m += 1
    return sorted(rows)


class PrimeSums:
    def __init__(self, bound):
        self.norms = []
        self.plain = [arb(0)]
        self.divided = [arb(0)]
        for n, p, mult in prime_records(bound):
            w = mult*arb(p).log()
            self.norms.append(n)
            self.plain.append(self.plain[-1]+w)
            self.divided.append(self.divided[-1]+w/n)

    def at(self, z, gamma):
        l, m, h = bisect_left(self.norms, z/2), bisect_right(self.norms, z), bisect_right(self.norms, 2*z)
        x = ball(z).sqrt()
        a, b = self.plain, self.divided
        return 2/gamma.b*(gamma.root/x*(a[m]-a[l])-x/gamma.root*(b[m]-b[l])
                         +gamma.root*x*(b[h]-b[m])-(a[h]-a[m])/(gamma.root*x))


def main():
    check_inputs()
    ctx.prec = 256
    gamma = GammaTerms()
    keys = [(p, m) for p in CASES for m in ((1,2) if p in (2,547,727) else (1,2,4))]
    largest = max(locations(p**m)[-1] for p, m in keys)
    primes = PrimeSums(int(2*largest)+1)
    records = []
    for p, m in keys:
        enc = []
        for z in locations(p**m):
            value = 2*(gamma.origin-gamma.at(z)+primes.at(z,gamma))
            rounded = (value*DEN).floor().unique_fmpz()
            assert rounded is not None and arb(rounded)/DEN < value < arb(rounded+1)/DEN
            enc.append(int(rounded))
        records.append("\t".join(map(str, (p, m, *enc))))
    body = "\n".join(records)+"\n"
    if sys.argv[1:] == ["--print-table"]:
        print(body, end="")
        return
    assert not sys.argv[1:] and body == (EVIDENCE/"curvatures.tsv").read_text()
    print(f"{3*len(records)} curvature samples in intervals of width 10^-36 match complete finite prime-power sums.")
    print("All three values use the full arithmetic and Archimedean terms, including the sample at log 2.")


if __name__ == "__main__":
    main()
