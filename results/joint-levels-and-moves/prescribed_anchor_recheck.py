"""Check the prescribed-anchor certificate with polynomial splitting and secants."""

import argparse
from functools import lru_cache
from itertools import product
from pathlib import Path

from flint import arb, arb_mat, ctx
from three_recheck import FIELDS, Reference, SCALE, secant_errors

DATA = Path(__file__).resolve().parent/"evidence"


def read_rectangle_corners(chart_count,prefix,den_limit):
    """Decode each tree by tracking dyadic integer addresses on its four axes."""
    trees = (DATA/(prefix+"-tree.tsv")).read_text().splitlines()[1:]
    assert len(trees) == 1280
    boxes = [[] for _ in range(chart_count)]
    for n,record in enumerate(trees):
        first, encoded = record.split("\t")
        assert int(first) == n
        a, rem = divmod(n,64)
        b, rem = divmod(rem,4)
        u,v = divmod(rem,2)
        origins = (arb(a-10), arb(b-8), (u-1)*arb(7)/2, (v-1)*arb(11)/5)
        size = (arb(1), arb(1), arb(7)/2, arb(11)/5)
        stack = [(1,(0,0,0,0))]
        units = 0
        for token in map(int,encoded.split()):
            denominator,address = stack.pop()
            if token == 0:
                assert denominator < den_limit
                for bits in reversed(list(product((0,1),repeat=4))):
                    stack.append((2*denominator,tuple(2*j+t for j,t in zip(address,bits))))
                continue
            assert 1 <= token <= chart_count and den_limit%denominator == 0
            units += (den_limit//denominator)**4
            intervals = [tuple(o+h*(k+t)/denominator for t in (0,1))
                         for o,h,k in zip(origins,size,address)]
            boxes[token-1].append(tuple(product(*intervals)))
        assert not stack and units == den_limit**4
    return boxes


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 148
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace",choices=("physical","limiting"))
    physical=parser.parse_args().trace == "physical"
    prefix="prescribed-anchor" if physical else "limiting-anchor"
    den_limit,leaves,count=(4,1415,41) if physical else (16,1430,35)
    refs = [Reference(f,physical) for f in FIELDS]
    charts = [tuple(map(int,s.split())) for s in
              (DATA/(prefix+"-charts.tsv")).read_text().splitlines()[1:]]
    boxes = read_rectangle_corners(len(charts),prefix,den_limit)
    grid_scale = SCALE if physical else 2*SCALE
    factor = grid_scale//1000
    assert factor*1000 == grid_scale
    identity = arb_mat([[int(i==j) for j in range(5)] for i in range(5)])

    def exact_anchor(k):
        (a,b),(c,d) = [r.physical_at(k) for r in refs]
        return (a,c,b,d,arb(1))

    @lru_cache(maxsize=None)
    def node(n):
        (a,b),(c,d) = [r.observation(arb(n)/grid_scale) for r in refs]
        return (a,c,b,d,arb(1))

    loglo,loghi = [(2*arb(2).sqrt()*k).log() for k in (13450,13490)]
    low_index = (loglo*grid_scale).ceil().unique_fmpz()
    high_index = (loghi*grid_scale).floor().unique_fmpz()
    assert low_index is not None and high_index is not None
    anchor_nodes = [exact_anchor(13450)] + [node(n) for n in range(int(low_index),int(high_index)+1)] + [exact_anchor(13490)]
    for chart,rectangles in zip(charts,boxes):
        _,b,c,d,e = chart
        vertices = [exact_anchor(13471)] + [node(n*factor) for n in (b,c,d,e)]
        columns = arb_mat(vertices).solve(identity)
        dual = [tuple(columns[i,j] for i in range(5)) for j in range(5)]
        ceiling = [arb(0) for _ in range(5)]
        pieces = [(10546,10550,(1,2,3,4),anchor_nodes)]
        for lo,hi,faces in ((b,c,(0,3,4)),(d,e,(0,1,2))):
            pieces.append((lo,hi,faces,[node(i) for i in range(lo*factor,hi*factor+1)]))
        for lo,hi,faces,points in pieces:
            errors = secant_errors(lo,hi,grid_scale)
            for j,row in enumerate(dual):
                assert sum(abs(t) for t in row[:4]) < arb(3)/20
                maximum = max(sum(x*y for x,y in zip(value,row)).upper() for value in points)
                upper = maximum + sum(abs(x)*y for x,y in zip(row,errors))
                assert upper < 2
                if j in faces:
                    ceiling[j] = ceiling[j].max(upper)
        for corners in rectangles:
            for target in corners:
                for row,bound in zip(dual,ceiling):
                    coordinate = row[-1] + sum(x*y for x,y in zip(row,target))
                    assert coordinate > bound + arb(1)/200, (chart,target,coordinate-bound)

    method="Polynomial factorization" if physical else "Primitive Conrey tables"
    print(f"{method} and independent secants: all {leaves} closed leaves pass all 16 corners.")
    print("Both moving arcs and the full prescribed-anchor interval have uniform face clearance > 1/200.")
    assert sum(map(len,boxes)) == leaves and len(charts) == count


if __name__ == "__main__":
    main()
