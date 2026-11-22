"""Certify optimum additive and octave sampling counts at the early four-gain window."""

from fractions import Fraction as F

from flint import arb, ctx

from patterns import Host, ball


def bisect_value(move, threshold, left, right, rising):
    orientation = 1 if rising else -1
    assert orientation * (move(left) - threshold) < 0
    assert orientation * (move(right) - threshold) > 0
    for _ in range(102):
        middle = (left + right) / 2
        sign = orientation * (move(middle) - threshold)
        if sign < 0:
            left = middle
        elif sign > 0:
            right = middle
        else:
            break
    return ball(left).union(ball(right))


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 218
    host = Host((571, 599), 0, 16 * 46696)

    def move(k, y):
        return (host.score(2 * k * y) - host.score(k * y)) / host.a

    left, right = F("46695.21"), F("46695.22")

    def line(k):
        f = move(k, left) * ball(left).sqrt()
        g = move(k, right) * ball(right).sqrt()
        slope = (g - f) / ball(right - left)
        return slope, f - slope * ball(left)

    u, v = line(1)
    s, t = line(2)
    peak = (t - v) / (u - s)
    margin = (u * peak + v) / peak.sqrt()
    integer = F(46695)
    integer_margin = move(1, integer)
    inverse_slope = 2 * peak.sqrt() / (u - v / peak) - 2 * peak.sqrt() / (s - t / peak)
    assert arb("1472.01872576754168") < inverse_slope < arb("1472.01872576754169")
    for label, q, pinned, expected in (
        ("Integer threshold", integer_margin, True, (4, 97121)),
        (
            "One millionth below the real threshold",
            margin - arb(1) / 10**6,
            False,
            (680, 21987939),
        ),
    ):
        x = (
            ball(integer)
            if pinned
            else bisect_value(lambda y: move(1, y), q, integer, F("46695.216167402620"), True)
        )
        z = bisect_value(lambda y: move(2, y), q, F("46695.216167402622"), F(46696), False)
        width = z - x
        octave = (z / x).log() / host.a
        for length, count in zip((width, octave), expected, strict=True):
            assert count - 1 < 1 / length < count
        print(f"{label}: {expected[0]} phases per unit and {expected[1]} phases per octave.")
        for title, val in (("left", x), ("right", z), ("width", width), ("octave", octave)):
            n = int((val * 10**18).floor().unique_fmpz())
            print(f"  {title} floor at eighteen decimal places: {n}")
        if pinned:
            missed = 1 - 3 * width
            assert arb("0.00021166213190490") < missed < arb("0.00021166213190491")
    print(
        "Limiting width divided by the margin deficit is in (1472.01872576754168, 1472.01872576754169)."
    )


if __name__ == "__main__":
    main()
