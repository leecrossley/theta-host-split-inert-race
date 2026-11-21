"""Rebuild the ordinate-dependent dyadic partition and its exact shift-error budget."""
from hashlib import sha256
from heapq import heapify, heappop, heappush
from pathlib import Path
from flint import arb, ctx

ROOT = Path(__file__).parent
TABLE = ROOT/"evidence"/"zero-brackets.tsv"
BINS = ROOT/"evidence"/"phase-bins.tsv"
TABLE_HASH = "e5483b464160fbce34313daa263f18dae8cc22035044a2053784a274e13a10f3"
WEIGHT_BITS = 32
ALIGN = 64


def rounded_weights():
    if sha256(TABLE.read_bytes()).hexdigest() != TABLE_HASH:
        raise ValueError("Wrong predecessor bracket table")
    ctx.prec = 160
    a,r = arb(2).log(),arb(2).sqrt()
    c = 2*(r-1)**2
    result = []
    for line in TABLE.read_text().splitlines():
        _,_,lo,hi = map(int,line.split())
        t = arb(lo+hi,hi-lo)/(2*10**14)
        sinc = (a*t/2).sin()
        kernel = (3-2*r*(a*t).cos())/(c*(t*t+arb(1)/4))
        weight = abs(4*kernel*sinc/a)
        numerator = int((weight.abs_upper()*2**WEIGHT_BITS).ceil().unique_fmpz())+2
        result.append(numerator)
    if len(result) != 6000:
        raise ValueError("Wrong number of phase coordinates")
    return result


def allocate(weights):
    # All bounds share denominator 2^(WEIGHT_BITS+ALIGN) before the factor 44/7.
    bits = [1]*len(weights)
    queue = [(-(w << (ALIGN-1)),j) for j,w in enumerate(weights)]
    heapify(queue)
    error = -sum(cost for cost,_ in queue)
    unit = 1 << (WEIGHT_BITS+ALIGN)
    while 440*error >= 7*unit:
        negative, j = heappop(queue)
        if bits[j] >= ALIGN:
            raise ArithmeticError("Exhausted the exact integer scale")
        half = (-negative)//2
        error -= half
        bits[j] += 1
        heappush(queue,(-half,j))
    return bits


def main():
    weights = rounded_weights()
    bits = allocate(weights)
    expected = [tuple(map(int,line.split())) for line in BINS.read_text().splitlines()]
    if list(zip(weights,bits)) != expected:
        raise ValueError("The published weight and bin table differs from the recomputation")
    print(f"6000 certified coefficient ceilings reproduced at 160 bits; total phase bits {sum(bits)}.")
    print(f"Dyadic exponents range from {min(bits)} through {max(bits)}; the greedy table matches exactly.")


if __name__ == "__main__":
    main()
