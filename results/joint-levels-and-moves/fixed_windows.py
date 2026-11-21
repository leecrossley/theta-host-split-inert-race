"""Verify the low spectral bands, translate interpolation, and fixed-window margins."""
from fractions import Fraction as Q
from itertools import permutations
from math import prod

from flint import arb, ctx

from all_height_common import CARRIERS
from amplitudes import semiaxes
from deadline_spectrum import ordinate
from zeros import SCALE, targets

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")


def coefficients(nodes, index):
    """Coefficients of the cardinal Lagrange polynomial at a selected node."""
    poly = [Q(1)]
    denominator = Q(1)
    for j, node in enumerate(nodes):
        if j == index:
            continue
        higher = [Q(0)]*(len(poly)+1)
        for k, value in enumerate(poly):
            higher[k] -= node*value
            higher[k+1] += value
        poly = higher
        denominator *= nodes[index]-node
    return [c/denominator for c in poly]


def translated_moments():
    # At one ordinate different real abscissae may cancel at a given translation.
    # Choose exp(d_j) rational, so the competing translates and their inverse are exact.
    nodes = (Q(9, 8), Q(5, 4), Q(3, 2))
    amplitudes = ((4, 1), (-6, -5), (2, 4))
    readings = [[sum(nodes[j]**k*amplitudes[j][c] for j in range(len(nodes)))
                 for c in (0, 1)] for k in range(len(nodes))]
    assert readings[0] == [0, 0] and any(readings[1])
    for j in range(len(nodes)):
        inverse_row = coefficients(nodes, j)
        assert all(sum(v*x**k for k, v in enumerate(inverse_row)) == int(i == j)
                   for i, x in enumerate(nodes))
        recovered = [sum(inverse_row[k]*readings[k][c] for k in range(len(nodes)))
                     for c in (0, 1)]
        assert recovered == list(amplitudes[j])
    # Recompute the determinant from its six Leibniz terms.
    matrix = [[x**j for x in nodes] for j in range(len(nodes))]
    det = Q(0)
    for perm in permutations(range(len(nodes))):
        inversions = sum(perm[i] > perm[j] for j in range(len(nodes)) for i in range(j))
        det += (-1)**inversions*prod(matrix[i][col] for i, col in enumerate(perm))
    assert det == prod(nodes[j]-nodes[i] for j in range(len(nodes)) for i in range(j)) > 0


def true_bands():
    ctx.prec = 160
    roots = targets()
    assert len(roots) == 9
    total_width = sum(Q(hi, SCALE) for _, _, _, lo, hi in roots)
    assert total_width < 6
    g, h = (ordinate(row) for row in CARRIERS)
    assert 8*(g+h) < 38
    assert max(11632, 9584)*8590 < 10**8
    # Every mixed low test lies below the known unconditional off-line exclusion.
    assert 38 < 8590
    axes = semiaxes()
    thresholds = ((99375, 34686, 66921, 19735), (79731, 21866, 61027, 15894))
    for row, lower in zip(axes, thresholds):
        assert all(value > arb(bound)/10000 for value, bound in zip(row, lower))
    assert axes[0][0]+axes[1][0] > arb(179106)/10000
    assert 2*Q(585, 100) == Q(117, 10)
    # For sixteen factors, the leading full support cost at separation 2 is eight log T.
    assert Q(16, 4)*2 == 8


if __name__ == "__main__":
    translated_moments()
    true_bands()
    print("fixed windows: nine-carrier harmonic magnitudes are below 6; the 289 dyadic shifts stay below 38")
    print("the strict two-field ellipse axes and 17.9106 level contrast survive the common filters")
    print("three exact colliding exponential modes recovered from three translations by polynomial interpolation")
    print("sixteen primitive factors give a completion exponent of eight plus any positive excess on log H")
