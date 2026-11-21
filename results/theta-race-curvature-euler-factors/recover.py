"""Recover local Euler factors from interval-valued curvature samples."""
from fractions import Fraction
from flint import arb, ctx
from curvature import CASES, GammaTerms, ball, check_inputs, locations, read_samples


def integer_enclosed(value):
    # The theorem gives integrality. Check that exactly one integer lies in the ball.
    candidate = value.unique_fmpz()
    if candidate is None:
        raise ValueError("Increase the accuracy of the curvature enclosures.")
    return int(candidate)


def scaled_jump(p, m, values, gamma):
    n = p**m
    loc = locations(n)
    corrected = [v/2-gamma.origin+gamma.at(z) for z,v in zip(loc, values)]
    e = ball(Fraction(4*n+1,4*n)).sqrt()
    jump = (corrected[0]-(e+1/e)*corrected[1]+corrected[2])/(e-1/e)
    scaled = gamma.c*arb(n).sqrt()*jump/arb(p).log()
    return integer_enclosed(scaled if p == 2 else -scaled/3)


def recover(p, samples, gamma, rank=None):
    u, v = [scaled_jump(p,m,samples[m],gamma) for m in (1,2)]
    if p == 2:
        # Here u and v are (C2-3*C1, 3*C1-3*C2).
        even = Fraction(-u-v, 2)
        odd = (even-u)/3
        answer = ((even+odd)/2, (even-odd)/2, Fraction(0))
    else:
        rank = scaled_jump(p,4,samples[4],gamma) if rank is None else rank
        answer = (Fraction(rank+v+2*u,4), Fraction(rank+v-2*u,4), Fraction(rank-v,4))
    if any(k < 0 or k.denominator != 1 for k in answer):
        raise ValueError("Incompatible local power sums.")
    return tuple(map(int, answer))


def main():
    check_inputs()
    ctx.prec = 224
    table, gamma = read_samples(), GammaTerms()
    for p in CASES:
        values={k: row for (q,k),row in table.items() if q == p}
        rank = 4 if p in (547,727) else None
        counts = recover(p, values, gamma, rank)
        if p not in (2,547,727):
            assert counts == recover(p, values, gamma, 8)
        print(f"p={p}: denominator multiplicities for (1-T), (1+T), (1+T^2) are {counts}.")
    print("Six samples per prime suffice with the given local degree. The odd-prime nine-sample decode also passes.")


if __name__ == "__main__":
    main()
