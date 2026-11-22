"""Certify the target's contributions from all prime ideals of residue degree two or four."""

import argparse
from fractions import Fraction
from pathlib import Path
import sys
from flint import arb, ctx


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("primary", type=Path)
    ap.add_argument("independent", type=Path)
    args = ap.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from word_cover import Race, prefix_tables

    ctx.prec = 200
    edges = [1 << i for i in range(59)]
    results = {}
    targets = (
        (0, 54, -1080, -1079),
        (1, 54, 101, 102),
        (0, 55, 1078, 1079),
        (1, 55, -25, -24),
    )
    micro_bounds = {(f, k): (lo, hi) for f, k, lo, hi in targets}
    for mode, binary in (("primary", args.primary), ("independent", args.independent)):
        fields, bits, recip = prefix_tables(binary.resolve(), mode, edges)
        group = [Race(f, bits, recip) for f in fields]
        for i, race in enumerate(group):
            for k in range(30, 56):
                if mode == "primary":
                    y = Fraction(1 << k)
                    val = race.move(y, 1)
                    results[i, k] = val
                else:
                    y = 1 << k
                    rows = []
                    for j in range(k, k + 3):
                        lo, hi = 1 << j, 2 << j
                        log_low, log_high, reciprocal_low, reciprocal_high = [
                            v - u
                            for u, v in zip(fields[i][lo], fields[i][hi], strict=True)
                        ]
                        rows.append(
                            (
                                arb(log_low + log_high, log_high - log_low)
                                / 2 ** (bits + 1),
                                arb(
                                    reciprocal_low + reciprocal_high,
                                    reciprocal_high - reciprocal_low,
                                )
                                / 2 ** (recip + 1),
                            )
                        )
                    (p, q), (t, u), (v, w) = rows
                    rt = arb(2).sqrt()
                    val = (
                        8
                        / ((2 * rt * (rt - 1) ** 2) * arb(2).log() * arb(y).sqrt())
                        * (
                            -2 * p
                            + 2 * y * q
                            + (1 + rt) * t
                            - (4 + 2 * rt) * y * u
                            - v / rt
                            + 4 * rt * y * w
                        )
                    )
                    assert results[i, k].overlaps(val), (i, k)
                    results[i, k] = results[i, k].intersection(val)
                if (i, k) in micro_bounds:
                    lo, hi = micro_bounds[i, k]
                    assert lo < 10**6 * val < hi, (mode, i, k, val)
    for field, exponent, low, high in targets:
        name = ("547/727", "571/599")[field]
        print(
            f"{name} at exponent {exponent}: higher-degree move lies strictly between {low} and {high} millionths."
        )
    print(
        "All 52 higher-degree moves agree across both full sieves and both score formulas."
    )


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
