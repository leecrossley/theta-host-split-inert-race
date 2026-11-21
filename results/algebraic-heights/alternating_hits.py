"""Count disjoint alternating excursions in three tagged signed-support tables."""
from bisect import bisect_left
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import tarfile

HERE = Path(__file__).resolve().parent
DIGESTS = {
    "finite-1000.tsv": "71fd24544f0aba5a313d0670e2a52638bb10b0ea224fb28648eceecbd7f75476",
    "finite-2000.tsv": "d49ec3e9de025e45f584ae8b8d44f050c69287a00ec6db8293d7295749583394",
    "sphere-cover.tsv": "8b9af3d0b7d257d4f8c535f69979e09f2b0b1dc1e31caeafcf8b28ee2ec9d74c",
}


def supports():
    selected = []
    with tarfile.open(HERE/"evidence/signed-supports.tar.xz","r:xz") as tables:
        assert set(tables.getnames()) == DIGESTS.keys()
        for name,digest in DIGESTS.items():
            rows = tables.extractfile(name).read()
            assert sha256(rows).hexdigest() == digest
            for line in rows.splitlines():
                values = tuple(map(int,line.split()))
                if name == "sphere-cover.tsv":
                    top,sign,m,k = values
                    centre = Fraction(m,64)
                    radius = Fraction(top+82,16984)+Fraction(2,k)
                else:
                    m,sign,_ = values
                    centre = Fraction(m,100)
                    n,k = (560,16) if name == "finite-1000.tsv" else (1070,20)
                    radius = Fraction(n,3000)+Fraction(2,k)
                left,right = centre-radius,centre+radius
                if 15 < left and right <= 12500:
                    selected.append((right,left,sign))
    return sorted(selected)


def alternating(intervals,first):
    last = 0
    wanted = first
    size = 0
    for right,left,sign in intervals:
        if sign == wanted and last < left:
            size += 1
            last = right
            wanted = -wanted
    return size


def dynamic_length(intervals):
    # Each prefix stores the largest chain ending in each sign. Binary
    # search gives the prefix ending strictly before the new left endpoint.
    ends = [r for r,l,s in intervals]
    best = [[0],[0]]
    for i,(right,left,sign) in enumerate(intervals):
        before = bisect_left(ends,left,0,i)
        index = int(sign==1)
        new = 1+best[1-index][before]
        best[index].append(max(best[index][-1],new))
        best[1-index].append(best[1-index][-1])
    return max(seq[-1] for seq in best)


def verify():
    intervals = supports()
    sizes = tuple(alternating(intervals,sign) for sign in (-1,1))
    assert sizes == (4893,4894) and dynamic_length(intervals)==max(sizes)
    # log(2 sqrt(2) * 10^6) < 15, using lower Taylor sums for exp(15)
    # and the rational upper bound sqrt(2) < 3/2.
    term = total = Fraction(1)
    for n in range(1,31):
        term *= Fraction(15,n)
        total += term
    assert total > 3_000_000
    traversals = max(sizes)-1
    positives,negatives = traversals+17,traversals+18
    assert (positives,negatives)==(4910,4911)
    assert (positives+negatives+2)//3 == 3274
    assert (positives+2)//3 == (negatives+2)//3 == 1637
    assert (traversals+1)//2 == 2447
    early = [row for row in intervals if row[0] <= 1000]
    assert [alternating(early, sign) for sign in (-1, 1)] == [505, 506]
    assert dynamic_length(early) == 506
    assert 2 * 505 + 35 == 1045 and (1045 + 2) // 3 == 349
    print("Two exact interval algorithms find 4894 alternating disjoint signed supports. "
          "The 4893 intervening traversals give that many crossings of each signed target and of zero.")
    print("Including the early crossings gives 9821 signed hits before completion height exp(12500). Under Dirichlet GRH these give at least 3274 classes, 1637 per sign, and the zero hits give 2447 classes.")
    print("The finite prime chains alone give 349 classes among 1045 signed hits before completion height exp(1000).")


if __name__ == "__main__":
    verify()
