"""Compare each single-field observation space and certify the sharp final observation height."""

from recover_thresholds import FIELDS, HEIGHTS, norm, primes_to


def scaled_weight(n, y):
    return max(0, min(2 * n - 2 * y, 4 * y - n))


def combined_row_parts(n):
    # Coefficients of 1 and sqrt(2) after multiplying by a nonzero prime-dependent factor.
    return (
        tuple(-2 * scaled_weight(n, y) for y in HEIGHTS),
        tuple(scaled_weight(n, 2 * y) for y in HEIGHTS),
    )


def pairing(row, linear_form):
    return sum(c * w for c, w in zip(row, linear_form, strict=True))


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    support = []
    for p in primes_to(8 * max(HEIGHTS)):
        n = norm(p, FIELDS[1])
        if not min(HEIGHTS) < n < 8 * max(HEIGHTS):
            continue
        parts = combined_row_parts(n)
        if any(any(row) for row in parts):
            support.append((p, norm(p, FIELDS[0]), n, parts))
    assert len(support) == 3327 and all(a for _, a, _, _ in support)
    a_rows = sorted((a, p, n, parts) for p, a, n, parts in support if a <= 46687)
    assert [(a, p, n) for a, p, n, _ in a_rows] == [
        (223, 223, 49729),
        (271, 271, 73441),
        (353, 353, 124609),
        (367, 367, 134689),
        (463, 463, 214369),
        (46687, 46687, 46687),
    ]
    first, second = (a_rows[k][3][0] for k in (0, 1))
    assert not any(a_rows[0][3][1]) and not any(a_rows[1][3][1])
    assert first[0] * second[1] != first[1] * second[0]
    a_survivor = (4, -39, 35)
    for *_, parts in a_rows[:-1]:
        assert all(pairing(row, a_survivor) == 0 for row in parts)
    assert pairing(a_rows[-1][3][0], a_survivor) != 0
    b_small = sorted((n, p, a, parts) for p, a, n, parts in support if n <= 46993)
    assert [(n, p, a) for n, p, a, _ in b_small] == [
        (46687, 46687, 46687),
        (46703, 46703, 46703),
        (46769, 46769, 2187339361),
        (46993, 46993, 46993),
    ]
    mixed_survivor = (0, -2, 3)
    assert [
        (pairing(parts[0], mixed_survivor), pairing(parts[1], mixed_survivor))
        for *_, parts in b_small
    ] == [(0, 0), (0, 0), (-264, 0), (-1160, 0)]
    a_norm = b_small[2][2]
    assert a_norm == 46769**2 == 8 * 273417420 + 1 and max(HEIGHTS) * 8 < a_norm // 8
    for p in (46447, 46559):
        assert norm(p, FIELDS[0]) == p * p and norm(p, FIELDS[1]) == p and p < min(HEIGHTS)
    assert 46447**2 == 4 * 539330952 + 1 and 46559**2 < a_norm < 2 * 46447**2
    print(
        "Degree-two primes 46447 and 46559 force three A readings through the upper endpoint 539330952 + 1/4."
    )
    print(
        "The 3327 contributing rational primes give single-field rank loss at 223, 271, 46687 on A and 46687, 46703, 46769 on B."
    )
    print("Before 46687 the last A-only combination has scaled coefficients (4, -39, 35).")
    print(
        "For 46769 <= H < 46993 the joint survivor needs prime 46769 on A at norm 2187339361; the sharp upper observation endpoint is strictly above 273417420 + 1/8."
    )


if __name__ == "__main__":
    main()
