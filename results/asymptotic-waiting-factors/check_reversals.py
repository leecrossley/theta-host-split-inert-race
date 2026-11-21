"""Enclose the critical-line variation over log(1601/1600) by two evaluations."""
from hashlib import sha256
from pathlib import Path
import lzma

from flint import arb, ctx

TABLE = Path(__file__).parent/'evidence'/'ordinate-brackets.tsv.xz'
DIGEST = '9a26df7437ee2b66c5d89fc458c0327c98f54b171e38bc575004f2b9a7c4eb6c'


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 144
    data = lzma.decompress(TABLE.read_bytes())
    assert sha256(data).hexdigest() == DIGEST
    rows = [tuple(map(int, row.split())) for row in data.splitlines()]
    assert len(rows) == 13946 and rows == sorted(rows)
    a = arb(2).log()
    r = arb(2).sqrt()
    delta = (arb(1601)/1600).log()
    first = second = arb(0)
    for _, _, left, right in rows:
        g = arb(left+right, right-left)/(2*10**14)
        # Twice the magnitude of the bilateral exponential coefficient.
        s = (a*g/2).sin()
        amplitude = 2*(3-2*r*(a*g).cos())*abs(s)/(a*(r-1)**2*(g*g+arb(1)/4))
        first += amplitude*g
        # Bound the full sine difference directly. Use 3-2 sqrt(2) cos(v)
        # = (sqrt(2)-1)^2 + 4 sqrt(2) sin^2(v/2).
        v = a*g/2
        amplitude2 = 2*((r-1)**2+4*r*v.sin()**2)*abs(v.sin())/(a*((r-1)*g)**2+a*((r-1)/2)**2)
        upper = (g*delta/2).upper()
        sine_bound = upper.sin() if upper < arb.pi()/2 else arb(1)
        second += amplitude2*2*sine_bound
    assert first < 3805
    assert second < 3805*delta
    assert delta < arb(1)/1600
    assert 3805*delta+arb(3)/2 < arb(389)/100
    assert second+arb(3)/2 < arb(389)/100
    assert arb(389)/100 < 2*arb(19735)/10000
    print("13,946 archived frequency brackets match the original SHA-256.")
    print("The retained derivative sum is below 3805; an independent sine-difference bound verifies the variation.")
    print("At C = 1601/800 the critical-line profile changes by less than 3.89, below twice 1.9735.")


if __name__ == "__main__":
    main()
