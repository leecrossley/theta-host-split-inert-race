"""Certified lower bounds on logarithmic root-mean-square error, in thousandths."""
from pathlib import Path
from flint import arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parent
# Each row gives level and normalized move for 547/727, then for 571/599.
BOUNDS = {
    1: (14419, 26294, 14586, 26690),
    2: (14124, 24533, 14304, 25005),
    4: (14087, 22590, 14245, 22959),
    8: (13966, 21535, 14220, 21613),
    16: (13338, 20198, 13313, 19995),
    32: (11311, 16265, 12686, 18542),
    64: (4967, 7754, 5732, 8325),
    128: (1500, 2180, 1750, 2500),
}


def brackets(field):
    data = [tuple(map(int,line.split())) for line in (ROOT/'evidence/brackets.tsv').read_text().splitlines()]
    return [(left, right) for f, _, _, left, right in data if f == field]


def main():
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    ctx.prec = 160
    size = max(BOUNDS)
    a, r = arb(2).log(), arb(2).sqrt()
    for field in (0, 1):
        moments = [[arb(0) for _ in range(size)] for _ in range(2)]
        for left, right in brackets(field):
            t = arb(left + right, right - left)/(2*10**14)
            theta = t*a
            g = (3-2*r*theta.cos())/(2*(r-1)**2*(arb(1)/4+t*t))
            weights = [2*g*g, 2*(g*2*(theta/2).sin()/a)**2]
            for j in range(size):
                elt = (j*theta).cos()
                for kind in (0, 1):
                    moments[kind][j] += weights[kind]*elt
        for kind, row in enumerate(moments):
            for m, thresholds in BOUNDS.items():
                matrix = arb_mat([[row[abs(i-j)] for j in range(m)] for i in range(m)])
                rhs = arb_mat(m, 1, [0]*(m-1) + [1])
                inverse_entry = matrix.solve(rhs, algorithm='precond')[m-1, 0]
                assert inverse_entry > 0
                bound = arb(thresholds[2*field+kind])/1000
                assert 1 > bound**2*inverse_entry, (field,kind,m,inverse_entry)
    print('Rigorous preconditioned linear systems certify 32 lower bounds from the two disjoint sets of field-specific frequencies.')
    for m, bounds in BOUNDS.items():
        print(f'{m} coordinates, RMS bounds in thousandths (547/727 level and move, 571/599 level and move): '+', '.join(map(str,bounds))+'.')


if __name__ == '__main__':
    main()
