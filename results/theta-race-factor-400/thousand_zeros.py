"""Replay the extended disjoint family, reusing exact endpoint pairs checked earlier."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import os
from check_zeros import FAMILY, check_batch, read_table

EXTENDED = Path(__file__).parent/'evidence/thousand-brackets.tsv'


def extended_table():
    rows = [tuple(map(int,line.split())) for line in EXTENDED.read_text().splitlines()]
    allowed = {(q,n) for q,n,_,_ in FAMILY}
    if rows != sorted(rows) or {row[:2] for row in rows} != allowed:
        raise ValueError('Unexpected factor family or bracket order')
    last = {}
    for q,n,lo,hi in rows:
        if not last.get((q,n),0) < lo < hi:
            raise ValueError(f'Invalid positive bracket {q,n,lo,hi}')
        last[q,n] = hi
    return rows


def main():
    previous = set(read_table())
    rows = extended_table()
    additional = [row for row in rows if row not in previous]
    jobs = [additional[i:i+16] for i in range(0,len(additional),16)]
    workers = int(os.environ.get('THETA_WORKERS',min(12,os.cpu_count() or 1)))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        if sum(pool.map(check_batch,jobs)) != len(additional):
            raise ArithmeticError('Incomplete extended zero checks')
    print(f'The extended table contains {len(rows):,} disjoint-within-factor positive brackets. The earlier replay supplies {len(rows)-len(additional):,} identical pairs and this replay checks {len(additional):,} further Hardy sign changes.')
    print('Extended counts by modulus and character: '+', '.join(f'{q}/{n}={k}' for (q,n),k in sorted(Counter(v[:2] for v in rows).items()))+'.')


if __name__ == '__main__':
    main()
