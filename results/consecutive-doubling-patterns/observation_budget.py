"""Verify the small norm intervals and large rational-prime witnesses for each observation budget."""

from bisect import bisect_right
from fractions import Fraction
from math import isqrt

from recover_thresholds import FIELDS, norm, primes_to

TARGET = 46769
LOW = 46447
HIGH = 46559
NEAR = (17258590369, 17258590801)
FAR = (17498715103, 17498715233)
LAST = 17341923439


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def delayed_on_b(p):
    return (
        p % 16 in (1, 15)
        and is_prime(p)
        and norm(p, FIELDS[0]) == p
        and norm(p, FIELDS[1]) == p * p
    )


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    n1, n2, n = LOW**2, HIGH**2, TARGET**2
    delayed = [
        p for p in primes_to(46992) if norm(p, FIELDS[0]) == p * p and norm(p, FIELDS[1]) == p
    ]
    assert delayed[-4:] == [46337, LOW, HIGH, TARGET]
    small = [p for p in delayed if p < 46656]
    assert len(small) == 306 and small[-2:] == [LOW, HIGH]
    norms = [p * p for p in small]
    z = Fraction(n, 8)
    start = bisect_right(norms, z) - 1
    stop = norms.index(n1)
    # Between norms[k] and norms[k+1] the two preceding expiries force a second
    # observation to be at least norms[k-1]. The next two norms are on their last branch.
    ratios = [Fraction(norms[k + 2], norms[k - 1]) for k in range(start, stop)]
    assert len(ratios) == 188 and max(ratios) == Fraction(377952481, 325766401) < 2
    assert n1 < 2 * norms[stop - 1] and 2 * n1 > n2
    for p in (571, 599, *FIELDS[1][2:]):
        a = norm(p, FIELDS[0])
        assert a < n1 or a > 8 * n
    assert all(delayed_on_b(p) for p in (*NEAR, *FAR, LAST))
    p, q = NEAR
    for last, edge in ((p, 8 * n1), (LAST, 8 * n2)):
        assert last < edge and not any(delayed_on_b(k) for k in range(last + 1, edge + 1))
    assert 4 * n < p < 8 * n1 < q < 8 * n2 < 8 * n < min(FAR) < max(FAR) < 2 * p
    assert 2 * n2 < p and q < 2 * p and n2 < n < p // 4
    assert Fraction(p * p, 8) == 37232367690617444520 + Fraction(1, 8)
    assert Fraction(LAST**2, 8) == 37592788570517198340 + Fraction(1, 8)
    assert 136174087494 * n < p * p < 136174087495 * n
    print(
        "306 delayed nuisance norms; 188 adjacent-interval tests have maximum ratio 377952481/325766401 < 2."
    )
    print(
        "The A-height infima for one, two and three A readings are 2167740481 (attained), 2157323809 (open), and 2187339361/8 (open)."
    )
    print(
        f"The largest delayed primes below eight times the two relevant norms are {LAST} and {p}."
    )
    print(
        f"Additional obstruction primes {q}, {FAR[0]} and {FAR[1]} have degree one on A and two on B."
    )
    print(
        "With one or two A readings the common upper endpoint must strictly exceed 37592788570517198340 + 1/8 or 37232367690617444520 + 1/8, respectively."
    )
    print(
        "A third A reading reduces the infimum by a factor strictly between 136174087494 and 136174087495."
    )


if __name__ == "__main__":
    main()
