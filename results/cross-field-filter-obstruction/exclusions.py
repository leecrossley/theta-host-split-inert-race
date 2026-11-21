"""Exclude every opposite-field zero on each of the indicated closed brackets."""
from collections import defaultdict
from multiprocessing import Pool
from os import environ
from pathlib import Path
from sys import argv
from flint import acb, arb, ctx, dirichlet_char

ROOT = Path(__file__).resolve().parent
SCALE = 10**14


def factors():
    result = []
    for field, primes in enumerate(((547, 727), (571, 599))):
        for prime in primes:
            for modulus, residue in ((1, 1), (8, 5), (16, 5), (16, 13)):
                conductor = prime * modulus
                label = next(k for k in range(prime-1, conductor, prime)
                             if k % modulus == residue % modulus)
                result.append((field, conductor, label))
    return result


def check(job):
    q, label, rows, method = job
    ctx.prec = 104 if method == 'interval' else 128
    char = dirichlet_char(q, label)
    assert char.parity() == 1 and char.is_primitive() and char.conductor() == q
    half = arb(1) / 2
    for left, right in rows:
        if method == 'interval':
            t = arb(left + right, right - left) / (2*SCALE)
            value = char.l_function(acb(half, t))
            assert not value.contains(0), (q, label, left, right, value)
        else:
            midpoint = arb(left + right) / (2*SCALE)
            radius = arb(right - left) / (2*SCALE)
            tmax = arb(right) / SCALE
            strength = abs(char.hardy_z(midpoint))
            # |Z(t)| = |L(1/2+it)|. Abel summation bounds |L'(1/2+it)|.
            derivative = q * (2 + 4 * (half + tmax))
            assert strength > radius * derivative, (q, label, left, right, strength)
    return len(rows)


def main(method):
    if not __debug__ or method not in ('interval', 'midpoint'):
        raise RuntimeError('Run with assertions and choose interval or midpoint.')
    data = [tuple(map(int, row.split()))
            for row in (ROOT/'evidence/brackets.tsv').read_text().splitlines()]
    chars = factors()
    assert len(data) == 710 and set(row[:3] for row in data) == set(chars)
    intervals = defaultdict(list)
    for field, _, _, left, right in data:
        assert 0 < left < right
        intervals[field].append((left, right))
    for rows in intervals.values():
        rows.sort()
        assert all(r < s for (_, r), (s, _) in zip(rows, rows[1:]))
    jobs = []
    for field, q, label in chars:
        rows = intervals[1-field]
        jobs.extend((q, label, rows[j:j+24], method) for j in range(0, len(rows), 24))
    workers = max(1, min(8, int(environ.get('THETA_WORKERS', '4'))))
    with Pool(workers) as pool:
        total = sum(pool.map(check, jobs))
    assert total == 5680
    name = 'Closed-interval Dirichlet L' if method == 'interval' else 'Midpoint Hardy modulus with derivative enclosure'
    print(f'{name} excludes all {total:,} opposite-field factors on 710 certified intervals (381 and 329).')


if __name__ == '__main__':
    main(argv[1])
