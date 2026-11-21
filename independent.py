"""Rebuild signed frequency sums and count their exact midpoint distances with binary searches."""
import numpy as np
from flint import arb, ctx
from certificate import (PAIRS, UNIT, charge, diagonal_checks, read_inputs, sorted_nodes,
                         verify_proportions)


def make_spectrum():
    rows, head = read_inputs()
    ctx.prec = 224
    a, root = arb(2).log(), arb(2).sqrt()
    table = [[], []]
    half_squared = fourth_powers = pair_squares = remainder = arb(0)
    for j, (_, _, left, right) in enumerate(rows):
        # Enclose the full zero interval using midpoint and radius in integer ticks.
        t = arb(left+right, right-left)/(2*10**14)
        s = (t*a/2).sin()
        positive = (root-1)**2+4*root*s**2
        magnitude = abs(positive*s)/(a*(root-1)**2*(t**2+arb(1)/4))
        if j in head:
            b = magnitude/27
            pair_squares += half_squared*b**2
            half_squared += b**2
            fourth_powers += b**4
            group = 0
        else:
            b = magnitude
            remainder += 2*b**2
            group = 1
        table[group].append((left+right, right-left, charge(b)))
    # Diagonal convolution energy, classified by equal and different index pairs.
    diagonal = diagonal_checks(2*half_squared, 6*fourth_powers+24*pair_squares, remainder)
    return [np.array(x, dtype=np.int64).T for x in table], diagonal


def signed_sums(table):
    mid, error, cap = table
    m, b = np.r_[mid, -mid], np.r_[cap, cap]
    assert int(max(b)) < UNIT
    assert sum(int(v)**2 for v in b) < 2**63-UNIT
    first, second = np.triu_indices(len(m))
    points = m[first]+m[second]
    positive = points > 0
    i, j = first[positive], second[positive]
    new = points[positive]
    weights = ((1+(i != j).astype(np.int64))*b[i]*b[j]+UNIT-1)//UNIT
    freq = np.r_[m, new, -new, 0]
    baseline, linear, quadratic = [np.zeros_like(freq) for _ in range(3)]
    linear[:len(m)] = b
    quadratic[len(m):-1] = np.r_[weights, weights]
    quadratic[-1] = (sum(int(v)**2 for v in b)+UNIT-1)//UNIT
    baseline[-1] = UNIT
    assert int((points == 0).sum()) == 800 and len(freq) == 1281601
    return freq, [baseline, linear, quadratic], 4*int(max(error))


def exact_cumulative(position, columns, uncertainty, pairs):
    p, cols, powers = sorted_nodes(position, columns, uncertainty)
    pref = [np.r_[np.int64(0), np.cumsum(c)] for c in cols]
    levels = []
    for power in powers:
        h = 2**power
        # Search for the first index outside each open midpoint neighborhood.
        upper = np.searchsorted(p, p+2*h, side="left")
        nearby = [prefix[upper]-prefix[1:] for prefix in pref]
        counts = [int(cols[i]@nearby[j])+int(cols[j]@nearby[i]) for i, j in pairs]
        levels.append((h, counts))
    return levels, uncertainty


def main():
    (selected, other), diagonal = make_spectrum()
    head = exact_cumulative(*signed_sums(selected), PAIRS)
    m, e, b = other
    rest = exact_cumulative(np.r_[m, -m], [np.r_[b, b]], 2*int(max(e)), [(0, 0)])
    print("Half-angle coefficients and grouped diagonal energies pass at 224 bits.")
    print("Independent signed-pair enumeration and exact cumulative-distance counts pass.")
    verify_proportions(head, rest, diagonal)


if __name__ == "__main__":
    main()
