"""Compute which prime-logarithm coefficients disappear from the three grid thresholds."""

from math import isqrt

HEIGHTS = (46656, 46691, 46695)
LAST = 46993
FIELDS = ((547, 727, 343057911569, 343070636977), (571, 599, 413718278401, 603339156001))


def primes_to(limit):
    marked = bytearray([1]) * (limit + 1)
    marked[:2] = bytes(2)
    for d in range(2, isqrt(limit) + 1):
        if marked[d]:
            marked[d * d :: d] = bytes((limit - d * d) // d + 1)
    return (p for p in range(3, limit + 1, 2) if marked[p])


def norm(p, field):
    if p in field:
        return 0
    u, v = field[:2]
    common = (pow(p, (u - 1) // 2, u) == 1) == (pow(p, (v - 1) // 2, v) == 1)
    degree = next(d for d in (1, 2, 4) if pow(p, d, 16) in (1, 15) and (common or d % 2 == 0))
    return p**degree


def determinant(a, b, c):
    return sum(
        a[i] * (b[(i + 1) % 3] * c[(i + 2) % 3] - b[(i + 2) % 3] * c[(i + 1) % 3])
        for i in range(3)
    )


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    assert LAST < 2 * min(HEIGHTS)
    lost = []
    for p in primes_to(LAST):
        n = norm(p, FIELDS[1])
        last = max(n, norm(p, FIELDS[0]))
        if min(HEIGHTS) < n and last <= LAST:
            # One common nonzero scalar per row; the three thresholds are scaled by sqrt(y).
            lost.append((last, p, n, tuple(max(n - y, 0) for y in HEIGHTS)))
    lost.sort()
    assert lost == [
        (46687, 46687, 46687, (31, 0, 0)),
        (46703, 46703, 46703, (47, 12, 8)),
        (46993, 46993, 46993, (337, 302, 298)),
    ]
    coeffs = [item[3] for item in lost]
    assert determinant(*coeffs) == 35960
    survivor = (0, -2, 3)
    dots = [sum(x * y for x, y in zip(row, survivor, strict=True)) for row in coeffs]
    assert dots == [0, 0, 290]
    for last, p, _, row in lost:
        print(f"Cutoff {last}: rational prime {p}, scaled coefficient row {row}.")
    print(
        "The last three rows have determinant 35960; the surviving combination after 46703 has coefficients (0, -2, 3)."
    )


if __name__ == "__main__":
    main()
