"""Check the usable frequency band and the margins after log smoothing."""
from flint import arb, ctx

from quadrature import read, verify_frames


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 144
    brackets = read("brackets.tsv")
    unit = 10**14
    assert len(brackets) == 710
    assert all(0 < left < right < 133 * unit for _, _, _, left, right in brackets)
    assert max(q for _, q, _, _, _ in brackets) == 11632
    assert 11632 * 8590 < 10**8
    # An even probability supported on [-epsilon,epsilon] has real
    # Fourier transform at least 1 - (epsilon*t)^2/2 at ordinate t.
    epsilon = arb(1) / 10**6
    retention = 1 - (epsilon * 133)**2 / 2
    verify_frames(retention)
    print("All 710 bracket ordinates are below 133 and lie inside Platt's zero-verification range for their conductors.")
    print("All 28 strict frame bounds survive even log smoothing of radius one millionth by a second-moment estimate.")


if __name__ == "__main__":
    main()
