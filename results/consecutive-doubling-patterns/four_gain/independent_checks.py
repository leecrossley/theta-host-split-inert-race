"""Check ordinary prime totals, independently sieved cells, and octave secant enclosures."""

import argparse
import gzip
import subprocess
from pathlib import Path

from flint import arb, ctx

from check_octaves import (
    EVIDENCE,
    EXCEPTIONAL,
    OCTAVES,
    PAIRS,
    WITNESSES,
    exceptional,
    tagged,
)

# Numbers of primes at the exact powers of two. OEIS A007053, table by David Baugh.
PI2 = {
    24: 1077871,
    25: 2063689,
    37: 5586502348,
    38: 10866266172,
    39: 21151907950,
    40: 41203088796,
    41: 80316571436,
    42: 156661034233,
    43: 305761713237,
    44: 597116381732,
    45: 1166746786182,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("independent_sieve", type=Path)
    args = ap.parse_args()
    ctx.prec = 180
    base = tagged(EVIDENCE, "primary")
    samples = {}
    for k in OCTAVES:
        file = EVIDENCE / "four-gain" / f"octave-{k}.txt.gz"
        width = 2 ** (k - 15)
        values = [[arb(0), arb(0)] for _ in PAIRS]
        selected = {0, 997, 16383, 25309, 32767}
        for primes in EXCEPTIONAL:
            for p in primes:
                if 2**k <= p < 2 ** (k + 1):
                    selected.add((p - 2**k) // width)
        prime_count = 0
        with gzip.open(file, "rt") as stream:
            assert stream.readline().strip() == f"# octave {k} cells 32768"
            for j in range(32768):
                idx, count, *moment = map(int, stream.readline().split())
                assert idx == j
                prime_count += count
                if j in selected:
                    samples[k, j] = (count, *moment)
                lo = arb(2**k + j * width)
                hi = lo + width
                l0, l1 = lo.log(), hi.log()
                for field in range(2):
                    total, first, mag = moment[3 * field : 3 * field + 3]
                    m = arb(first) / width
                    # Linear interpolant between endpoint weights, with curvature remainder.
                    limit = arb(mag) * (arb(width) / lo) ** 2 / 32
                    values[field][0] += ((total - m) * l0 + m * l1) / 4 + arb(
                        0, limit.upper()
                    )
                    qlimit = limit * (2 * l1 + 3) / lo
                    values[field][1] += ((total - m) * l0 / lo + m * l1 / hi) / 4 + arb(
                        0, qlimit.upper()
                    )
            assert not stream.read()
        assert prime_count == PI2[k + 1] - PI2[k], k
        for f, row in enumerate(values):
            if (k, f) in base:
                assert all(
                    x.overlaps(y) for x, y in zip(row, base[k, f], strict=True)
                ), (k, f)
            else:
                base[k, f] = row
        print(
            f"Octave {k}: ordinary census agrees with the independent pi(2^k) table.",
            flush=True,
        )
    keys = sorted(samples)
    run = subprocess.run(
        [str(args.independent_sieve.resolve())],
        input="".join(f"{k} {j}\n" for k, j in keys),
        capture_output=True,
        text=True,
        check=True,
    )
    observed = {}
    for line in run.stdout.splitlines():
        k, j, *row = map(int, line.split())
        assert (k, j) not in observed
        observed[k, j] = tuple(row)
    assert len(samples) == 44 and observed == samples
    print(
        f"All {len(samples)} independently sieved cells agree on every field moment and ordinary count."
    )
    root = arb(2).sqrt()
    B = 2 * root * (root - 1) ** 2
    a = arb(2).log()

    def verify(data):
        def level(f, k):
            y = 2**k
            (P, Q), (R, T) = data[k, f], data[k + 1, f]
            return 8 * ((2 * P - 2 * y * Q) + (4 * y * T - R)) / (B * arb(y).sqrt())

        def move(f, k):
            return (level(f, k + 1) - level(f, k)) / a

        for n in range(42):
            assert any(
                move(f, k) < arb(3947) / 2000
                for f in range(2)
                for k in (n, n + 1)
                if (k + 2, f) in data
            ), n
        for f, k, lo, hi in WITNESSES:
            assert lo < 100 * move(f, k) < hi, (f, k)

    verify(base)
    verify(exceptional(base))
    print(
        "Independent endpoint interpolation confirms the 42 exclusions and all four signed intervals under both conventions."
    )


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
