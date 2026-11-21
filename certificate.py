"""Shared inputs and exact rational inequalities for the two spectral calculations."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import numpy as np
from flint import arb

BASE = Path(__file__).resolve().parent / "evidence"
UNIT = 10**7
SCALE = 2 * 10**14
PAIRS = ((0, 1), (1, 1), (1, 2), (2, 2))
COEFF = tuple(F(n, 325057991399) for n in
              (-22196224000, 99943833600, 108661670400, -2970592000, -16515360000))
CASES = ((20000, F(1, 10), F(1, 100)), (50000, F(1, 20), F(1, 16)),
         (100000, F(1, 50), F(1, 12)), (250000, F(1, 50), F(1, 10)))


def read_inputs():
    if not __debug__:
        raise RuntimeError("Run with Python assertions enabled.")
    for name, digest in (
        ("zero-brackets.tsv", "e5483b464160fbce34313daa263f18dae8cc22035044a2053784a274e13a10f3"),
        ("bracket-check.txt", "c454eaa9ebff709db4df3668def47cc4bcba2f4879b7118a2dacc1ffbc093212"),
        ("deadline-check.txt", "00847fa1cf5a53a20c34042d8d54115e81d11b502eb419d2a0e345df88344194")):
        assert sha256((BASE/name).read_bytes()).hexdigest() == digest
    rows = [tuple(map(int, line.split())) for line in (BASE/"zero-brackets.tsv").read_text().splitlines()]
    head = [int(s)-1 for s in (BASE/"head-rows.txt").read_text().split()]
    assert len(rows) == 6000 and len(head) == 800 and head == sorted(set(head))
    assert min(head) >= 0 and max(head) < len(rows)
    assert all(0 < lo < hi < 10**17 for _, _, lo, hi in rows)
    return rows, set(head)


def charge(x):
    """Upper-bound a nonnegative Fourier magnitude in units of 10^-7."""
    cap = int((x*UNIT).upper().ceil().unique_fmpz())
    assert cap >= 0 and x <= arb(cap)/UNIT and cap < 10*UNIT
    return cap


def diagonal_checks(second, fourth, residual):
    assert second > arb(10065)/10000 and second < arb(10066)/10000
    assert fourth < arb(30034)/10000 and residual < arb(923)/100000
    return COEFF[0] + COEFF[2]*F(10065, 10000) + COEFF[4]*F(30034, 10000)


def sorted_nodes(position, columns, uncertainty):
    order = np.argsort(position)
    position = position[order]
    columns = [column[order] for column in columns]
    assert all(np.all(column >= 0) for column in columns)
    totals = [sum(map(int, c)) for c in columns]
    assert max(totals)**2 < 2**63
    # This bound protects signed positions, their differences and the search keys.
    assert max(abs(int(position[0])), abs(int(position[-1]))) < 2**58
    gap = int(np.diff(position).min())
    assert gap > 0
    first = gap.bit_length()-1
    assert 2**first > uncertainty and first >= 1
    stop = int(position[-1]-position[0]).bit_length()
    return position, columns, range(first, stop)


def error_bounds(levels, uncertainty, length, edge):
    """Integrate a decreasing dyadic majorant against cumulative close-pair masses."""
    def kernel(gap):
        d = F(gap-uncertainty, SCALE)
        return min(1, 2/(d*length*(1-edge))) * min(1, 2/(d*length*edge))
    result = [F(0)] * len(levels[0][1])
    for h, counts in levels:
        chi = kernel(h)-kernel(2*h)
        assert chi >= 0
        for j, count in enumerate(counts):
            assert count >= 0
            result[j] += chi*count
    last_h, last_counts = levels[-1]
    return [(value+kernel(2*last_h)*last)/UNIT**2 for value, last in zip(result, last_counts)]


def verify_proportions(head, rest, diagonal):
    for length, edge, target in CASES:
        moment_errors = error_bounds(*head, length, edge)
        tail_error, = error_bounds(*rest, length, edge)
        average = diagonal-sum(abs(c)*e for c, e in zip(COEFF[1:], moment_errors))
        bad = (F(923, 100000)+tail_error) / F(149, 100)**2
        proportion = (1-edge)*(average-bad) - F(100, length)
        assert proportion > target
        # Display only certified lower decimals rounded towards zero.
        lower = (1000000*proportion).__floor__()
        assert F(lower, 1000000) < proportion
        print(f"Log length >= {length}: signed fraction > {lower}/1000000 > {target}.")
    asymptotic = (F(10065, 10000)-F(25, 729))**2/(2*F(30034, 10000))
    asymptotic -= F(923, 100000)/F(149, 100)**2
    assert asymptotic > F(3, 20)
    assert 5-F(149, 100)-F(152, 100)-F(7, 1000) > F(19735, 10000)
    print("Each signed excursion set has lower uniform log density greater than 3/20.")
