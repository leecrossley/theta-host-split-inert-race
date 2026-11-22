"""Verify rational-prime witnesses for every finite octave used in the two-field proof."""

import argparse
import csv
from bisect import bisect_left
from math import isqrt
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent
CONDUCTORS = (547, 727, 571, 599)
OMIT = {2, *CONDUCTORS, 343057911569, 343070636977, 413718278401, 603339156001}
RESIDUES = [{x*x % q for x in range(1, q)} for q in CONDUCTORS]
DIRECT = (7, 13, 15)
MIXED = (1, 2, 13, 14, 4, 8, 7, 11)


def prime_divisors(limit):
    sieve = bytearray([1]) * (limit+1)
    sieve[:2] = bytes(2)
    for d in range(2, isqrt(limit)+1):
        if sieve[d]:
            sieve[d*d::d] = bytes((limit-d*d)//d+1)
    return [d for d, flag in enumerate(sieve) if flag]


def label(p):
    assert p % 16 in (1, 15) and p not in OMIT
    return sum(1 << j for j, q in enumerate(CONDUCTORS) if p % q in RESIDUES[j])


def trace(code):
    return (1-((code & 1) + ((code >> 1) & 1)), 1-(((code >> 2) & 1) + ((code >> 3) & 1)))


def check_initial():
    # Only the first primes in three different directions are needed.
    prime_list = prime_divisors(1 << 15)
    sizes = []
    for n in range(8, 22):
        start = 1 << (n+2)
        found = set()
        for p in range(start+1, 2*start, 2):
            if p % 16 not in (1, 15) or p in OMIT:
                continue
            if any(p % d == 0 for d in prime_list if d*d <= p):
                continue
            u, v = trace(label(p))
            if (u, v) == (0, 0):
                continue
            if u < 0 or (u == 0 and v < 0):
                u, v = -u, -v
            found.add((u, v))
            if len(found) == 3:
                break
        assert len(found) == 3
        sizes.append(len(found))
    print("Three distinct directions at ordinary primes in all 14 octaves for exponents 8 through 21.")


def check_archive(backend):
    with (ROOT / "evidence/dyadic-proth-witnesses.tsv").open() as source:
        reader = csv.DictReader(source, delimiter="\t")
        assert reader.fieldnames == ["bin", "pattern", "coefficient_offset", "proth_base"]
        records = {}
        for row in reader:
            j, code, offset, base = (int(row[name]) for name in reader.fieldnames)
            assert (j, code) not in records and j >= 0 and 0 <= code < 16
            records[j, code] = (offset, base)
    edges = [1 << 24]
    for _ in range(max(j for j, _ in records)):
        edges.append(13*edges[-1]//10)
    # Beneath 2^40 all trial candidates lie below 2^41.
    divisors = prime_divisors(isqrt(1 << 41))
    checked = {}
    proth = []
    for (j, code), (offset, base) in records.items():
        lower = edges[j]
        octave = lower.bit_length()-1
        step = 1 << (4 if octave < 40 else octave//2+1)
        coefficient = (lower-1)//step + offset
        p = 1 + coefficient*step
        assert lower <= p and 64*p <= 65*lower and label(p) == code
        if octave < 40:
            assert base == 0 and p < 1 << 41
            for d in divisors:
                if d*d > p:
                    break
                assert p % d
        else:
            assert 0 < coefficient < step and coefficient % 2 and 2 <= base < p
            if backend is None:
                assert pow(base, (p-1)//2, p) == p-1
            proth.append((p, base))
        checked[j, code] = p
    if backend is not None:
        with tempfile.TemporaryDirectory() as temporary:
            batch = Path(temporary) / "certificates.txt"
            batch.write_text("".join(f"{p} {b}\n" for p, b in proth))
            output = subprocess.check_output([str(backend), str(batch)], text=True)
            assert output == f"{len(proth)}\n"
    assert len(set(checked.values())) == len(checked)
    used = set()
    for n in range(22, 2900):
        y = 1 << n
        direct = n <= 1580
        lower = 4*y if direct else isqrt(4*y-1)+1
        j = bisect_left(edges, lower)
        for code in (DIRECT if direct else MIXED):
            p = checked[j, code]
            used.add((j, code))
            assert 4*y < (p if direct else p*p) < 8*y
            if not direct:
                assert p.bit_length()-1 == (n+2)//2
    assert used == set(checked)
    assert {trace(c) for c in DIRECT} == {(-1, 0), (0, -1), (-1, -1)}
    assert [trace(c) for c in MIXED] == [(0,1)]*2 + [(0,-1)]*2 + [(1,0)]*2 + [(-1,0)]*2
    print(f"All {len(checked):,} archived primes certified ({len(proth):,} Proth tests).")
    print("Three terminal directions verified for exponents 22 through 1580; eight signed square-norm witnesses in one rational-prime octave for every exponent 1581 through 2899.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proth", type=Path, help="optional exact GMP modular checker")
    options = parser.parse_args()
    check_initial()
    check_archive(options.proth.resolve() if options.proth else None)
