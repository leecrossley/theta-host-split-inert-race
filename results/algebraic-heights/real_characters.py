"""Certify real-axis nonvanishing for the quadratic characters in the pair moduli."""
from itertools import combinations
from math import isqrt, prod

PRIMES = (547, 727, 571, 599)
GRID = 1 << 40
CUTOFF = 1 << 16


def small_character(d, n):
    if d == 1:
        return 1
    if n % 2 == 0:
        return 0
    if d == -4:
        return 1 if n%4 == 1 else -1
    sign = 1 if n%8 in (1,7) else -1
    return sign if d == 8 else sign * small_character(-4,n)


def outward_roots(n):
    # Enclose n^(1/256) on the fixed grid by eight directed square roots.
    lo = hi = n * GRID
    for _ in range(8):
        lo = isqrt(lo * GRID)
        radicand = hi * GRID
        hi = isqrt(radicand)
        if hi * hi < radicand:
            hi += 1
    lower = [GRID]
    upper = [GRID]
    for _ in range(16):
        lower.append(lower[-1]*lo//GRID)
        upper.append((upper[-1]*hi+GRID-1)//GRID)
    return lower,upper


def partitions():
    cut = [1]
    while cut[-1] < CUTOFF:
        left = cut[-1]
        cut.append(min(CUTOFF,left+max(1,left//128)))
    return cut


def lower_bounds(tables, cond, cuts, powers):
    # The complete character period gives an absolute bound for all partial sums.
    moduli = [len(rows) for rows in tables]
    total = 0
    maximum = 0
    minima = []
    block_min = 0
    block = 0
    for n in range(1, max(cond,CUTOFF-1)+1):
        v = 1
        for table,mod in zip(tables,moduli):
            v *= table[n%mod]
        total += v
        if n == cond:
            assert total == 0
        maximum = max(maximum, abs(total))
        if n >= CUTOFF:
            continue
        if n == cuts[block]:
            block_min = total
        else:
            block_min = min(block_min,total)
        if n+1 == cuts[block+1]:
            minima.append(block_min)
            block += 1
    assert len(minima)==len(cuts)-1
    lower = []
    for k in range(16):
        acc = 0
        for left,right,value in zip(cuts,cuts[1:],minima):
            if value >= 0:
                scale = (255-k)*powers[left][0][k]
            else:
                scale = (256-k)*powers[right][1][k+1]
            acc += value*scale*(right-left)//(256*left*right)
        acc -= (maximum*powers[CUTOFF][1][k+1]+CUTOFF-1)//CUTOFF
        lower.append(acc)
    return min(lower)


def verify():
    cuts = partitions()
    powers = {x:outward_roots(x) for x in cuts}
    tables = {}
    for p in PRIMES:
        squares = {n*n%p for n in range(1,p)}
        tables[p] = [0] + [1 if a in squares else -1 for a in range(1,p)]
    bounds = []
    for count in range(3):
        for subset in combinations(PRIMES,count):
            for d in (1,-4,8,-8):
                cond = abs(d)*prod(subset)
                if cond==1:
                    continue
                chars = [tables[p] for p in subset]
                if d != 1:
                    chars.append([small_character(d,n) for n in range(abs(d))])
                bounds.append((lower_bounds(chars,cond,cuts,powers),d*prod(-p for p in subset)))
    assert len(bounds)==43
    low,discriminant = min(bounds)
    assert 1000*low > 363*GRID
    assert discriminant == -547 and len(cuts)-1 == 1002
    print("All 43 nonprincipal real primitive characters exceed 0.363 on [15/16,1] by sixteen directed Abel bounds with 1,002 integration blocks per cell.")


if __name__ == "__main__":
    verify()
