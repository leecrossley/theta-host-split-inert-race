"""Recover the first tangency at higher precision from four finite score evaluations."""

from fractions import Fraction
from flint import ctx

from patterns import Host, ball

LEFT = Fraction(4669521, 100)
RIGHT = Fraction(4669522, 100)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 230
    host = Host((571, 599), 0, int(16 * RIGHT) + 1)

    def line(k):
        def sample(y):
            delta = host.score(2 * k * y) - host.score(k * y)
            return ball(y).sqrt() * delta / host.a

        left, right = sample(LEFT), sample(RIGHT)
        slope = (right - left) / ball(RIGHT - LEFT)
        return slope, left - slope * ball(LEFT)

    slope, intercept = line(1)
    other_slope, other_intercept = line(2)
    height = (other_intercept - intercept) / (slope - other_slope)
    margin = (slope * height + intercept) / height.sqrt()
    integer_height = Fraction(46695)
    integer_margin = (host.score(2 * integer_height) - host.score(integer_height)) / host.a
    for label, value, expected in (
        ("Peak height", height, 46695216167402621654621),
        ("Peak margin", margin, 1682449287659191709),
        ("Integer margin", integer_margin, 1682222888957865575),
    ):
        assert expected < value * 10**18 < expected + 1
        print(f"{label} is in ({expected}, {expected + 1}) / 10^18 by direct prime counting.")


if __name__ == "__main__":
    main()
