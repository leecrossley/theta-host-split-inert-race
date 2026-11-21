"""Verify finite arithmetic witnesses and the refinement calculation using exact fractions."""
from fractions import Fraction
from math import gcd, isqrt, prod

FIELDS = ((547, 727), (571, 599))
BANDS = ((11, 14), (15, 19), (21, 26))  # endpoints in tenths
LOW, HIGH = 1_500, 1_000_000
MODULI = (16, 547, 727, 571, 599)


def primes_through(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[:2] = bytes(2)
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit-p*p)//p + 1)
    return (p for p in range(2, limit + 1) if sieve[p])


def degree_one(p, host):
    u, v = FIELDS[host]
    if p % 16 not in (1,15) or gcd(p,u*v) != 1:
        return False
    square_u = pow(p, (u-1)//2, u) == 1
    square_v = pow(p, (v-1)//2, v) == 1
    return square_u == square_v


def crt(residues):
    modulus = prod(MODULI)
    return sum(a * (modulus//q) * pow(modulus//q, -1, q)
               for q, a in zip(MODULI, residues)) % modulus


def separating_classes():
    for host, nonsquare_index in ((0, 3), (1, 1)):
        residues = [1]*len(MODULI)
        q = MODULI[nonsquare_index]
        residues[nonsquare_index] = next(a for a in range(2,q) if pow(a,(q-1)//2,q) == q-1)
        representative = crt(residues)
        assert gcd(representative, prod(MODULI)) == 1
        assert all(representative % q == a for a,q in zip(residues,MODULI))
        assert degree_one(representative,host)
        assert not degree_one(representative,1-host)
    print(f"Two reduced classes separate the fields modulo {prod(MODULI):,}.")


def cover(primes, left, right):
    # The start interval of this prime is (10p/right, 10p/left).
    # Use the greatest reachable right endpoint at each step.
    end = Fraction(LOW)
    j = 0
    chain = []
    while end <= HIGH:
        first = j
        while j < len(primes) and 10*primes[j] < right*end:
            j += 1
        assert j > first
        p = primes[j-1]
        assert left*end < 10*p
        chain.append(p)
        end = Fraction(10*p,left)
    return chain


def independent_check(chain, host, left, right, squares):
    # Each chain member is trial-divided anew and its character tested by
    # explicit square-residue sets, independently of the sieve and Euler tests.
    previous = None
    for p in chain:
        assert all(p % n for n in range(2,isqrt(p)+1))
        assert p % 16 in (1,15)
        good = []
        for u,v in FIELDS:
            assert p % u and p % v
            good.append((p%u in squares[u]) == (p%v in squares[v]))
        assert good[host] and not good[1-host]
        if previous is None:
            assert left*LOW < 10*p < right*LOW
        else:
            assert previous < p and left*p < right*previous
        previous = p
    assert 10*chain[-1] > left*HIGH


def finite_cover():
    limit = HIGH * BANDS[-1][1]//10
    omissions = {2,547,727,571,599,343057911569,343070636977,413718278401,603339156001}
    small = set(MODULI[1:])
    exclusive = ([],[])
    for p in primes_through(limit):
        if p in omissions or p%16 not in (1,15):
            continue
        flags = [degree_one(p,h) for h in (0,1)]
        if flags[0] != flags[1]:
            exclusive[0 if flags[0] else 1].append(p)
    squares = {q:{n*n%q for n in range(1,q)} for q in small}
    counts = []
    for host in (0,1):
        for left,right in BANDS:
            chain = cover(exclusive[host],left,right)
            independent_check(chain,host,left,right,squares)
            counts.append(len(chain))
    # First two bands end before the next scale and all three lie in
    # the positive support of their respective tent. Also p^2 > 8y.
    assert 10 < 11 < 14 and 14**2 < 200 < 15**2 and 15 < 19 < 20 < 21 < 26 < 40
    assert LOW > 8
    print("Six strict interval covers give six independent readings at every algebraic start from 1,500 through 1,000,000.")
    print("Trial division and quadratic-residue tables verified the six prime chains of lengths " + ", ".join(map(str,counts)) + ".")


def convolution(a,b):
    result = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            result[i+j] += x*y
    return result


def hat(z,b,u):
    return max(0,min(b*(z-u), b*b*u-z))


def refinements():
    examples = 0
    for b in (Fraction(6,5),Fraction(3,2),Fraction(2),Fraction(3)):
        for m in range(1,10):
            powers = [b**(m-1-k) for k in range(m)]
            coefficients = convolution([1]*m,powers)
            wide = [0]*(2*m+1)
            wide[0],wide[m],wide[-1] = b**m,-(b**m+1),1
            assert convolution([b,-(b+1),1],coefficients) == wide
            u = Fraction(7,5)
            knots = [u*b**k for k in range(2*m+1)]
            points = [knots[0]/2,*knots,*(sum(pair)/2 for pair in zip(knots,knots[1:])),2*knots[-1]]
            for z in points:
                assert hat(z,b**m,u) == sum(c*hat(z,b,u*b**j) for j,c in enumerate(coefficients))
            examples += 1
    print(f"Exact slope jumps and direct piecewise-linear evaluations agree for {examples} geometric refinements.")


if __name__ == "__main__":
    separating_classes()
    finite_cover()
    refinements()
