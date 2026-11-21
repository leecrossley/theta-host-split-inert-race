"""Certify uniform coverage for every prescribed start from 13,450 through 13,490."""

from collections import defaultdict
from functools import lru_cache
from itertools import product
import argparse
from pathlib import Path

from flint import arb, ctx
from three_profiles import HOSTS, Trace, vector
from three_starts import affine, integer_anchor, inverse

DATA = Path(__file__).resolve().parent / "evidence"
SIGNS = tuple(product((-1, 1), repeat=4))


def load_cells(count, prefix, max_depth):
    """Read one complete depth-first subdivision of every closed root box."""
    table = (DATA / (prefix+"-tree.tsv")).read_text().splitlines()
    assert table[0] == "cell\tpreorder" and len(table) == 1281
    cells = defaultdict(list)
    radius = (arb(1)/2, arb(1)/2, arb(7)/4, arb(11)/10)
    roots = product(range(-10, 10), range(-8, 8), (-1, 1), (-1, 1))
    for number, ((i, j, u, v), record) in enumerate(zip(roots, table[1:])):
        label, codes = record.split("\t")
        assert int(label) == number
        tokens = iter(map(int, codes.split()))

        def visit(centre, widths, depth):
            chart = next(tokens)
            if chart:
                assert 1 <= chart <= count
                cells[chart].append((centre, widths))
                return
            assert depth < max_depth
            half = tuple(x/2 for x in widths)
            for signs in SIGNS:
                moved = tuple(x+s*h for x,s,h in zip(centre,signs,half))
                visit(moved, half, depth+1)

        visit((i+radius[0], j+radius[1], u*radius[2], v*radius[3]), radius, 0)
        assert next(tokens, None) is None
    assert set(cells) == set(range(1, count+1))
    return cells


def geometry(charts, max_log, finish):
    for fixed, b, c, d, e in charts:
        assert fixed == 10548 and b < c and d < e
        intervals = sorted(((10546, 10550), (b,c), (d,e)))
        assert intervals[0][0] >= 8092 and intervals[-1][1] <= max_log
        assert all(hi[0]-lo[1] >= 702 for lo,hi in zip(intervals,intervals[1:]))
        assert intervals[-1][1]-intervals[0][0] <= 4304
    assert (arb(702)/1000).exp() > arb(2017)/1000
    assert 2*(arb(4304)/1000).exp() < 148
    assert (arb(8092)/1000).exp()/(2*arb(2).sqrt()) > 1155
    assert (arb(max_log)/1000).exp()/arb(2).sqrt() < finish


def target_caps(basis, cells):
    """Reserve a strict 1/200 clearance even at every closed-cell corner."""
    result = []
    for row in basis:
        minima = [affine(row, centre) - sum(abs(c)*h for c,h in zip(row,widths))
                  for centre,widths in cells]
        result.append(min(x.lower() for x in minima)-arb(1)/200)
    assert all(x > 0 for x in result)
    return result


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 112
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace",choices=("physical","limiting"))
    physical = parser.parse_args().trace == "physical"
    prefix = "prescribed-anchor" if physical else "limiting-anchor"
    count, depth, max_log, finish = (41,2,13624,583866) if physical else (35,4,13628,586206)
    lines = (DATA / (prefix+"-charts.tsv")).read_text().splitlines()
    assert lines[0].split() == ["singleton","left1","right1","left2","right2"]
    charts = [tuple(map(int,line.split())) for line in lines[1:]]
    assert len(charts) == len(set(charts)) == count
    cells = load_cells(count,prefix,depth)
    geometry(charts,max_log,finish)
    traces = [Trace(host, 3_000_000, physical=physical) for host in HOSTS]
    fixed, reference_height = integer_anchor(traces, 10548)
    assert reference_height == 13471
    exact = [(2*arb(2).sqrt()*k).log() for k in (13450,13490)]
    # Outward log endpoints also avoid placing a shell endpoint on an integer knot.
    scale = 10**9
    indices = [(exact[0]*scale).floor().unique_fmpz(), (exact[1]*scale).ceil().unique_fmpz()]
    assert all(n is not None for n in indices)
    anchor = [arb(int(n))/scale for n in indices]
    assert anchor[0] < exact[0] and exact[1] < anchor[1]
    assert arb(10546)/1000 < anchor[0] < anchor[1] < arb(10550)/1000

    @lru_cache(maxsize=None)
    def endpoint(n):
        return vector(traces, arb(n)/1000)

    calls = 0
    for number, (_,b,c,d,e) in enumerate(charts,1):
        basis = inverse([fixed]+[endpoint(n) for n in (b,c,d,e)])
        assert all(sum(abs(v) for v in row[:4]) < arb(3)/20 for row in basis)
        caps = target_caps(basis,cells[number])
        pieces = [(anchor,(1,2,3,4)), ([arb(b)/1000, arb(c)/1000],(0,3,4)),
                  ([arb(d)/1000, arb(e)/1000],(0,1,2))]
        for (left,right), affected in pieces:
            bound = [min(arb(2),cap) if i in affected else arb(2)
                     for i,cap in enumerate(caps)]
            stack = [(left,right,0)]
            while stack:
                lo, hi, depth = stack.pop()
                mid = (lo+hi)/2
                measurement = vector(traces, mid, ((hi-lo)/2).upper())
                calls += 1
                if all(affine(row,measurement) < cap for row,cap in zip(basis,bound)):
                    continue
                if depth >= 27:
                    raise ArithmeticError((number,affected,lo,hi))
                stack.extend(((lo,mid,depth+1), (mid,hi,depth+1)))

    print("Direct", "physical" if physical else "limiting", "trace:",count,"charts, 1280 closed root boxes and",
          sum(map(len,cells.values())),"certified leaves.")
    print("All real anchors in [13450,13490] pass all five face tests with clearance > 1/200.")
    print("Every affine coordinate < 2 on its complete chart, gradient norm < 3/20.")
    print(f"Template starts > 1155, every finish < {finish}, adjacent start ratio > 2.017, completion factor < 148.")
    print("Adaptive prime-norm interval evaluations:",calls)


if __name__ == "__main__":
    main()
