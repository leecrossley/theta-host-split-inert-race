"""Certify the three-start covering of a four-dimensional target box."""

from flint import arb, arb_mat, ctx
from three_profiles import HOSTS, Trace, vector

from pathlib import Path

ROOT = Path(__file__).resolve().parent
CELLS = tuple((i, j) for i in range(-10, 10) for j in range(-8, 8))


def charts():
    raw = (ROOT/"evidence"/"three-start-charts.tsv").read_text().splitlines()
    assert raw[0].split() == ["singleton", "left1", "right1", "left2", "right2"]
    return [tuple(map(int, row.split())) for row in raw[1:]]


def inverse(vertices):
    matrix = arb_mat([[v[k] for v in vertices] for k in range(4)] + [[1]*5])
    dual = matrix.inv()
    return [[dual[i, k] for k in range(5)] for i in range(5)]


def affine(row, v):
    return row[4] + sum((u*t for u,t in zip(row[:4],v)), arb(0))


def margins(rows, faces):
    lower = []
    for i, j in CELLS:
        centre = (i+arb(1)/2, j+arb(1)/2, 0, 0)
        vals = []
        for row, face in zip(rows, faces):
            cost = (abs(row[0])+abs(row[1]))/2 + 7*abs(row[2])/2 + 11*abs(row[3])/5
            vals.append(affine(row, centre)-cost-face)
        lower.append(min(v.lower() for v in vals))
    return lower


def check_geometry(configurations):
    scales = [t for chart in configurations for t in chart]
    assert min(scales) == 8092 and max(scales) == 13808
    assert (arb(8092)/1000).exp()/(2*arb(2).sqrt()) > 1155
    assert (arb(13808)/1000).exp()/arb(2).sqrt() < 701816
    assert 2*(arb(13808-8092)/1000).exp() < 608
    for a,b,c,d,e in configurations:
        blocks = sorted(((a,a),(b,c),(d,e)))
        assert b < c and d < e
        gaps = [blocks[k+1][0]-blocks[k][1] for k in (0,1)]
        assert min(gaps) >= 720
        assert blocks[-1][1]-blocks[0][0] <= 4004
    assert (arb(720)/1000-arb(1)/2000).exp() > arb(41)/20
    assert 2*(arb(4004)/1000+arb(1)/2000).exp() < 110


def integer_anchor(traces, t):
    root = arb(2).sqrt()
    y = (arb(t)/1000).exp()/(2*root)
    k = (y+arb(1)/2).floor().unique_fmpz()
    assert k is not None
    k = int(k)
    assert 1500 < k < 330000 and abs((arb(k)/y).log()) < arb(1)/2000
    z = 2*k*root
    pairs = []
    for tr in traces:
        level = tr.level(z,arb(0))
        move = (tr.level(arb(4*k),arb(0))-tr.level(arb(2*k),arb(0)))/tr.a
        pairs.append((level,move))
    return [pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1]], k


def table_margins(configurations, physical):
    traces = [Trace(host, 3_000_000, physical=physical) for host in HOSTS]
    points = {t: vector(traces,arb(t)/1000) for chart in configurations for t in chart}
    cache = {}
    per_unit = 12
    half = arb(1)/(2000*per_unit)
    margins_by_chart = []
    anchor_integers = []
    for chart in configurations:
        vertices = [points[t] for t in chart]
        if physical:
            vertices[0], k = integer_anchor(traces,chart[0])
            anchor_integers.append(k)
        rows = inverse(vertices)
        assert all(sum(abs(u) for u in row[:4]) < arb(3)/20 for row in rows)
        faces = [arb(0) for row in rows]
        for left, right, active in ((chart[1],chart[2],(0,3,4)),(chart[3],chart[4],(0,1,2))):
            for k in range(left*per_unit,right*per_unit):
                if k not in cache:
                    cache[k] = vector(traces, (2*k+1)*half, half)
                for j, row in enumerate(rows):
                    value = affine(row,cache[k]).upper()
                    assert value < arb(3)/2
                    if j in active:
                        faces[j] = faces[j].max(value)
        margins_by_chart.append(margins(rows,faces))
    print("physical" if physical else "limiting", "trace:", len(cache),
          "whole-arc enclosures, row norms < 3/20, all arc coordinates < 3/2")
    if physical:
        print("integer singleton for each timing row:", " ".join(map(str,anchor_integers)))
    return margins_by_chart


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 112
    atlas = charts()
    check_geometry(atlas)
    finite = table_margins(atlas, True)
    limiting = table_margins(atlas, False)
    selections = []
    for k in range(len(CELLS)):
        clearances = [min(a[k],b[k]) for a,b in zip(finite,limiting)]
        viable = [i+1 for i,v in enumerate(clearances) if v > arb(1)/250]
        assert viable, CELLS[k]
        selections.append(viable[0])
    print(len(atlas),"configurations cover all",len(CELLS),"closed boxes on both traces")
    print("one shared configuration per level cell; all five face clearances > 1/250")
    print("starts > 1155, finishes < 701816, adjacent ratios > 41/20, group completion factor < 110")
    print("positive weights > 1/400 persist under sup-norm perturbation at most 1/1000")
    print("chart indices, A-level cells -10 through 9; B-level cells -8 through 7:")
    for j in range(0, len(selections), 16):
        print(" ".join(map(str,selections[j:j+16])))


if __name__ == "__main__":
    main()
