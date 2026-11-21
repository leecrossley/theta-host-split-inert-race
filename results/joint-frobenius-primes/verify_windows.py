"""Replay the small prime gaps, certified chains and all sharp endpoints."""
from fractions import Fraction
from itertools import groupby, product
from math import isqrt
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
QS = (547, 727, 571, 599)
OMIT = set(QS) | {2, 343057911569, 343070636977, 413718278401, 603339156001}
SPLIT = tuple({a*a % q for a in range(1, q)} for q in QS)
SIEVE_LIMIT = 1 << 25


def pattern(p):
    assert p not in OMIT and p % 16 in (1, 15)
    return sum(1 << j for j, q in enumerate(QS) if p % q in SPLIT[j])


def pair(code):
    return ((code & 1) + ((code >> 1) & 1) - 1,
            ((code >> 2) & 1) + ((code >> 3) & 1) - 1)


def initial_primes():
    flags = bytearray(b'\1') * (SIEVE_LIMIT + 1)
    flags[:2] = b'\0\0'
    for d in range(2, isqrt(SIEVE_LIMIT) + 1):
        if flags[d]:
            flags[d*d::d] = bytes((SIEVE_LIMIT-d*d)//d + 1)
    patterns = [[] for _ in range(16)]
    for residue in (1, 15):
        for p in range(residue, SIEVE_LIMIT+1, 16):
            if flags[p] and p not in OMIT:
                patterns[pattern(p)].append(p)
    return [sorted(ps) for ps in patterns]


def summaries(sequences):
    rows = []
    for label, ps in sequences:
        # Lower bound 1 represents the initial interval before this state's first prime.
        rows.append((label, 1, ps[0]))
        for p, q in zip(ps, ps[1:]):
            if 3*q > 4*p:
                rows.append((label, p, q))
    return rows


def check_map(patterns):
    pairs = {t: [] for t in product((-1, 0, 1), repeat=2)}
    for code, ps in enumerate(patterns):
        pairs[pair(code)].extend(ps)
    patterns_map = summaries([(f"signs-{code:02d}", ps) for code, ps in enumerate(patterns)])
    pairs_map = summaries([(f"trace-{a},{b}", sorted(ps)) for (a, b), ps in pairs.items()])
    corners = [r for r in pairs_map if '0' not in r[0]]
    all_rows = patterns_map + pairs_map
    table = "state\tleft\tnext_prime\n" + "".join(f"{name}\t{p}\t{q}\n" for name,p,q in all_rows)
    assert (HERE / "evidence/wide-gaps.tsv").read_text() == table
    for group, right, worst, last in ((corners, 4831, 337, ("trace-1,1", 2879, 4831)),
                                       (pairs_map, 4831, 431, ("trace-1,1", 2879, 4831)),
                                       (patterns_map, 8737, 769, ("signs-14", 6047, 8737))):
        assert max(group, key=lambda row: row[2]) == last and last[2] == right
        assert max(Fraction(q,p) for _,p,q in group) == worst
    for label, l, u in all_rows:
        # Independent trial divisions and Euler's criterion at EVERY published deficient gap.
        code_label = label.startswith("signs-")
        def matches(p):
            if p in OMIT or p % 16 not in (1, 15) or p < 2:
                return False
            if any(p % d == 0 for d in range(2, isqrt(p)+1)):
                return False
            symbols = [pow(p, (q-1)//2, q) for q in QS]
            code = sum(1 << j for j, value in enumerate(symbols) if value == 1)
            a, b = pair(code)
            return label == (f"signs-{code:02d}" if code_label else f"trace-{a},{b}")
        assert (l == 1 or matches(l)) and matches(u) and not any(matches(n) for n in range(l+1, u))
    print(f"Exact sieve to {SIEVE_LIMIT:,}. The complete exceptional-gap table has "
          f"{len(patterns_map)} sign rows and {len(pairs_map)} trace rows including initial gaps.")
    print("Below 2^1585 the final thresholds are 14,493/4 for corners and trace pairs and 26,211/4 for all signs. The three optimal factors on this finite range are 337, 431 and 769.")


def trial_prime(p):
    return p > 2 and p % 2 and all(p % d for d in range(3, isqrt(p)+1, 2))


def check_chains(patterns, backend):
    path = HERE / "evidence/prime-chains.tsv"
    lower = 1 << 24
    previous = [max(p for p in ps if p < lower) for ps in patterns]
    counts = [0, 0]
    largest = 0
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / "proth-input"
        with work.open("w") as certs, path.open() as lines:
            assert next(lines).rstrip() == "bin\tpattern\tcoefficient_offset\tbase"
            blocks = groupby((tuple(map(int, line.split())) for line in lines), key=lambda t: t[0])
            bins = 0
            for j, rows in blocks:
                assert j == bins and (not j or prior_lower < 1 << 1585)
                e = lower.bit_length() - 1
                n = e//2 + 1 if e >= 40 else 4
                step = 1 << n
                origin = (lower - 1)//step
                rows = list(rows)
                assert [row[1] for row in rows] == list(range(16))
                for _, code, offset, a in rows:
                    k = origin + offset
                    p = 1 + k*step
                    assert lower <= p and 64*p <= 65*lower
                    assert previous[code] < p and 3*p <= 4*previous[code]
                    assert pattern(p) == code
                    if e < 40:
                        assert a == 0 and trial_prime(p)
                        counts[0] += 1
                    else:
                        assert k % 2 == 1 and 0 < k < step and 2 <= a < p
                        if backend is None:
                            assert pow(a, (p-1)//2, p) == p-1
                        else:
                            certs.write(f"{p} {a}\n")
                        counts[1] += 1
                    previous[code] = p
                    largest = max(largest, p)
                prior_lower = lower
                lower = 13*lower//10
                bins += 1
        assert min(previous) >= 1 << 1585
        if backend is not None:
            answer = subprocess.check_output([str(backend), str(work)], text=True)
            assert answer == f"{counts[1]}\n"
    print(f"Verified {bins:,} complete bins, {sum(counts):,} labelled primes "
          f"({counts[0]:,} trial and {counts[1]:,} Proth), and every connecting gap. "
          f"Largest witness has {largest.bit_length():,} bits.")


if __name__ == "__main__":
    assert len(sys.argv) == 2
    ps = initial_primes()
    check_map(ps)
    check_chains(ps, None if sys.argv[1] == "--python" else Path(sys.argv[1]).resolve())
