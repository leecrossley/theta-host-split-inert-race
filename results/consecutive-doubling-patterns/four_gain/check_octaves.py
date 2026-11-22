"""Reconstruct directed prime scores from octave cell moments and the tagged dyadic tables."""

import gzip
from pathlib import Path

from flint import arb, ctx, fmpz


PAPER = Path(__file__).resolve().parents[1]
EVIDENCE = PAPER / "evidence"
OCTAVES = (24, 37, 38, 39, 40, 41, 42, 43)
# Endpoints in hundredths of a normalized move.
WITNESSES = (
    (0, 36, -5206, -5203),
    (0, 38, -995, -993),
    (1, 40, -3166, -3162),
    (0, 41, -2267, -2263),
)
PAIRS = ((547, 727), (571, 599))
EXCEPTIONAL = ((343057911569, 343070636977), (413718278401, 603339156001))


def tagged(path, mode):
    shifts = (58, 96) if mode == "primary" else (62, 102)
    result = {}
    for field, (q, _) in enumerate(PAIRS):
        rows = (path / f"dyadic-{mode}-{q}.txt").read_text().splitlines()
        assert len(rows) == 39
        for k, line in enumerate(rows):
            end, p0, p1, q0, q1 = map(int, line.split())
            assert end == 2**k and p0 <= p1 and q0 <= q1
            result[k - 1, field] = (
                arb(p0 + p1, p1 - p0) / 2 ** (shifts[0] + 1),
                arb(q0 + q1, q1 - q0) / 2 ** (shifts[1] + 1),
            )
    return result


def read_moments(path):
    with gzip.open(path, "rt") as stream:
        prefix, marker, k, name, count = stream.readline().split()
        assert (prefix, marker, name) == ("#", "octave", "cells")
        k, count = int(k), int(count)
        assert 15 <= k <= 44 and count == 2**15
        width = 2**k // count
        values = [[arb(0), arb(0)] for _ in PAIRS]
        primes = 0
        for j in range(count):
            idx, num, *moments = map(int, stream.readline().split())
            assert idx == j and len(moments) == 6 and 0 <= num <= width
            primes += num
            start = arb(2**k + width * j)
            centre = start + arb(width) / 2
            ln = centre.log()
            inv = 1 / centre
            linear_q = (1 - ln) * inv**2
            # Taylor remainder per unit absolute weight. The raw integer moments
            # have weight 4*sign/f, so division by four follows summation.
            mag_p = (arb(width) / start) ** 2 / 8
            mag_q = mag_p * (2 * (start + width).log() + 3) / start
            for field in range(2):
                sign, offset, absolute = moments[field * 3 : field * 3 + 3]
                assert abs(sign) <= absolute <= 4 * width
                assert abs(offset) <= absolute * (width - 1)
                displacement = arb(offset) - arb(sign * width) / 2
                values[field][0] += (
                    sign * ln + displacement * inv + arb(0, (absolute * mag_p).upper())
                ) / 4
                values[field][1] += (
                    sign * ln * inv
                    + displacement * linear_q
                    + arb(0, (absolute * mag_q).upper())
                ) / 4
        assert not stream.read()
    return k, primes, values


def exceptional(octaves):
    """Restore the optional exceptional rational primes for the conductor-only convention."""
    result = dict(octaves)
    for field, ((u, v), primes) in enumerate(zip(PAIRS, EXCEPTIONAL, strict=True)):
        for p in primes:
            assert fmpz(p).is_prime() and p % u and p % v
            su = pow(p, (u - 1) // 2, u)
            sv = pow(p, (v - 1) // 2, v)
            s = 1 if su == 1 else -1
            f = (
                1
                if p % 16 in (1, 15) and (su == 1) == (sv == 1)
                else 2
                if p * p % 16 in (1, 15)
                else 4
            )
            sign = -s if f == 1 else -1
            norm = p**f
            key = norm.bit_length() - 1, field
            if key in result:
                initial_p, initial_q = result[key]
                mass = sign * arb(p).log()
                result[key] = initial_p + mass, initial_q + mass / norm
    return result


def moves(octaves, top):
    rt = arb(2).sqrt()
    fac = 8 / (2 * rt * (rt - 1) ** 2 * arb(2).log())
    vals = {}
    for k in range(top + 1):
        y = 1 << k
        for field in range(2):
            (p, q), (s, t), (v, w) = [octaves[k + j, field] for j in range(3)]
            first = -2 * p + 2 * y * q
            middle = (1 + rt) * s - (4 + 2 * rt) * y * t
            last = -v / rt + 4 * rt * y * w
            vals[field, k] = fac / arb(y).sqrt() * (first + middle + last)
    return vals


def main():
    ctx.prec = 160
    first, second = (tagged(EVIDENCE, m) for m in ("primary", "independent"))
    assert first.keys() == second.keys()
    for key in first:
        assert all(x.overlaps(y) for x, y in zip(first[key], second[key], strict=True))
    updated = dict(second)
    for expected in OCTAVES:
        path = EVIDENCE / "four-gain" / f"octave-{expected}.txt.gz"
        k, total, vals = read_moments(path)
        assert k == expected
        for field, row in enumerate(vals):
            key = k, field
            if key in second:
                for observed in (first[key], second[key]):
                    assert all(
                        x.overlaps(y) for x, y in zip(row, observed, strict=True)
                    ), key
            else:
                updated[key] = tuple(row)
        print(
            f"Octave [{2**k}, {2 ** (k + 1)}) contains {total} rational primes.",
            flush=True,
        )
    for convention, data in (
        ("trace omissions", updated),
        ("conductors only", exceptional(updated)),
    ):
        values = moves(data, 41)
        for n in range(42):
            witnesses = (
                value for (field, k), value in values.items() if k in (n, n + 1)
            )
            assert any(value < arb(3947) / 2000 for value in witnesses), (convention, n)
        for field, k, low, high in WITNESSES:
            value = 100 * values[field, k]
            assert low < value < high, (convention, field, k, value)
        print(
            f"{convention}: all 42 starts excluded; four signed move intervals verified."
        )


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
