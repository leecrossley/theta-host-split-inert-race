"""Bound the four finite-window moments by counting neighboring frequency bins."""
import numpy as np
from flint import arb, ctx
from certificate import (PAIRS, UNIT, charge, diagonal_checks, read_inputs,
                         sorted_nodes, verify_proportions)


def spectrum():
    rows, head = read_inputs()
    ctx.prec = 192
    root, a = arb(2).sqrt(), arb(2).log()
    positions, weights, radii = [[], []], [[], []], [0, 0]
    second = fourth_sum = residual = arb(0)
    for index, (_, _, lo, hi) in enumerate(rows):
        t = (arb(lo)/10**14).union(arb(hi)/10**14)
        g = (3-2*root*(a*t).cos())/(2*(root-1)**2*(t*t+arb(1)/4))
        amplitude = abs(4*g*(a*t/2).sin()/a)
        group = 0 if index in head else 1
        value = amplitude/27 if group == 0 else amplitude
        positions[group].append(lo+hi)
        radii[group] = max(radii[group], hi-lo)
        weights[group].append(charge(value/2))
        if group == 0:
            second += value**2/2
            fourth_sum += value**4
        else:
            residual += value**2/2
    diagonal = diagonal_checks(second, 3*second**2-3*fourth_sum/8, residual)
    return ([np.array(p, dtype=np.int64) for p in positions],
            [np.array(w, dtype=np.int64) for w in weights], radii, diagonal)


def head_nodes(m, b):
    n = len(m)
    assert int(max(b)) < UNIT and n == 800
    assert 2*sum(int(v)**2 for v in b) < 2**63-UNIT
    i, j = np.triu_indices(n)
    sums = m[i]+m[j]
    # A pair with different indices occurs twice in the convolution.
    plus = ((1+(i != j).astype(np.int64))*b[i]*b[j]+UNIT-1)//UNIT
    i, j = np.triu_indices(n, 1)
    both = np.r_[sums, m[i]-m[j]]
    pair_weights = np.r_[plus, (2*b[i]*b[j]+UNIT-1)//UNIT]
    frequencies = np.r_[0, m, -m, both, -both]
    zero, one, two = [np.zeros_like(frequencies) for _ in range(3)]
    zero[0] = UNIT
    one[1:1+2*n] = np.r_[b, b]
    two[0] = (2*int(b@b)+UNIT-1)//UNIT
    two[1+2*n:] = np.r_[pair_weights, pair_weights]
    assert len(frequencies) == 1281601
    return frequencies, [zero, one, two]


def count_bins(position, columns, indices, uncertainty):
    p, cols, powers = sorted_nodes(position, columns, uncertainty)
    diagonal = [int(cols[i]@cols[j]) for i, j in indices]
    levels = []
    for power in powers:
        h = 1 << power
        # If centre distance is < 2h, these half-h bin indices differ by at most 4.
        bins = (p-p[0])//(h//2)
        first = np.r_[0, np.flatnonzero(np.diff(bins))+1]
        occupied = bins[first]
        masses = [np.add.reduceat(c, first) for c in cols]
        after = np.searchsorted(occupied, occupied+5)
        prefixes = [np.r_[np.int64(0), np.cumsum(c)] for c in masses]
        neighbors = [prefix[after]-prefix[1:] for prefix in prefixes]
        counts = [int(masses[i]@masses[j])-same + int(masses[i]@neighbors[j])
                  + int(masses[j]@neighbors[i]) for (i, j), same in zip(indices, diagonal)]
        levels.append((h, counts))
    return levels, uncertainty


def main():
    positions, weights, radii, diagonal = spectrum()
    m, p = positions
    b, c = weights
    head = count_bins(*head_nodes(m, b), PAIRS, 4*radii[0])
    rest = count_bins(np.r_[p, -p], [np.r_[c, c]], [(0, 0)], 2*radii[1])
    print("800 ordinates give 1281601 separated Fourier nodes through degree two.")
    print("Grouped integer weights certify all four moment errors and the residual variance.")
    verify_proportions(head, rest, diagonal)


if __name__ == "__main__":
    main()
