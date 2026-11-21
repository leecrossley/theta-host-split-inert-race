"""Certify the sharp finite deadline and the first jump in its amplitude curve."""
from bisect import bisect_right
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
from prime_knots import event_sweep, direct_prefix

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")

CENTER = 804028    # k = 8y
P, Q = 20068, 214710
HEIGHT = 10**9


def interval(value, digits):
    den = 10**digits
    lo = (value*den).floor().unique_fmpz()
    assert lo is not None, value
    def fmt(n):
        a, b = divmod(int(n), den)
        return f"{a}.{b:0{digits}d}"
    return f"({fmt(lo)}, {fmt(lo+1)})"


def check_grid():
    path = Path(__file__).parent/"evidence"/"windows.tsv"
    rows = [tuple(map(int, line.split())) for line in path.read_text().splitlines()[1:]]
    gains = [y for y, lo, hi in rows if lo > 50000]
    losses = [y for y, lo, hi in rows if hi < -50000]
    assert P in losses and losses[bisect_right(losses, P)] == Q
    def deadlines(witnesses):
        for h in {2, HEIGHT, *(y for y in witnesses if 2 <= y <= HEIGHT)}:
            nxt = witnesses[bisect_right(witnesses, h)]
            yield h, nxt, 2*Fraction(nxt, h)
    assert max(c for _, _, c in deadlines(gains)) == Fraction(41, 4)
    assert max(c for h, _, c in deadlines(losses) if h != P) <= 10
    bridged = sorted([*losses, Fraction(CENTER, 8)])
    assert max(c for _, _, c in deadlines(bridged)) <= Fraction(41, 4)
    print("625 signed grid enclosures: gains by 41H/4 from H=2; losses by 10H except one witness gap")
    print("strict loss witnesses 20068 and 214710 flank the gap; the central bridge gives 41H/4")


def locate(pieces, keys, t):
    """Enclose the two crossings at the ends of the central gap."""
    left = [k for k in keys if 8*P <= k < CENTER]
    right = [k for k in keys if CENTER < k <= 8*Q]
    for k in left+right:
        assert pieces[k][2] < -t or pieces[k][2] > -t
    i = max(k for k in left if pieces[k][2] < -t)
    j = min(k for k in right if pieces[k][2] < -t)
    def cross(u, v, increasing):
        a, b, _ = pieces[u]
        lo, hi = arb(u)/8, arb(v)/8
        derivative = b*lo.union(hi)-a
        assert derivative > 0 if increasing else derivative < 0
        for _ in range(64):
            mid = (lo+hi)/2
            val = (a+b*mid)/mid.sqrt()
            assert val < -t or val > -t
            if (val < -t) == increasing:
                lo = mid
            else:
                hi = mid
        return lo.union(hi)
    a = cross(i, keys[bisect_right(keys, i)], True)
    b = cross(keys[bisect_right(keys, j)-2], j, False)
    assert 20068 < a < 21100 and 207000 < b < 214710
    return a, b


def main():
    ctx.prec = 176
    pieces = event_sweep()
    direct = direct_prefix()
    keys = sorted(pieces)
    assert keys == sorted(direct)
    for k in keys:
        assert pieces[k][2].overlaps(direct[k]) and max(pieces[k][2].rad(), direct[k].rad()) < arb("1e-30")
    gridpath = Path(__file__).parent/"evidence"/"windows.tsv"
    checked = 0
    for row in gridpath.read_text().splitlines()[1:]:
        y, lo, hi = map(int,row.split())
        if y > Q:
            break
        k = keys[bisect_right(keys, 8*y)-1]
        a, b, _ = pieces[k]
        v = (a+b*y)/arb(y).sqrt()*10000
        assert arb(lo) <= v <= arb(hi)
        checked += 1
    print(f"{len(keys)} knot and endpoint values agree in two independent enumerations through y={Q}; {checked} grid checks")

    # Before the decreasing cell every normalized slope is nonnegative.
    assert [k for k in keys if k <= 79] == [8, 31, 47, 49, 62, 79]
    assert all(pieces[k][2] >= 0 for k in (8, 31, 47, 49, 62) if k != 31)
    assert pieces[31][2].contains(0)  # all weights at their first support boundary equal zero
    A, negative_B, _ = pieces[62]
    B = -negative_B
    root = arb(2).sqrt()
    scale = 8/(2*root*(root-1)**2*arb(2).log())
    l31, l47, l7 = (arb(p).log() for p in (31, 47, 7))
    intercept = scale*((1+2/root)*l31+(l47+l7)/root)
    slope = scale*(4*(1+1/root)*l31/31+8*(l47/47+l7/49)/root)
    assert A.overlaps(intercept) and B.overlaps(slope)
    assert A > 0 and B > 0 and pieces[79][2] < -5
    assert (pieces[31][0]+pieces[31][1]*5)/arb(5).sqrt() > 5
    def early(t):
        return 2*((t+(t*t+4*A*B).sqrt())/(2*B))**2

    tau = -pieces[CENTER][2]
    assert arb("3.99703") < tau < arb("3.99704")
    core = [k for k in keys if 168286 <= k <= 1659407]
    assert core[0] == 168286 and core[-1] == 1659407 and CENTER in core
    assert all(pieces[k][2] > -tau+arb("0.000006") for k in core if k != CENTER)
    assert early(arb(0)) > 17 and early(arb(5)) < arb("18.1")
    print(f"{len(core)} central knots: unique minimum at y=201007/2; next knot values exceed it by 0.000006")
    print("threshold tau in", interval(tau, 12))
    print("sharp factor at 1.9735 in", interval(early(arb("1.9735")), 12))
    print("left limit at tau in", interval(early(tau), 9))
    for label, t in (("tau", tau), ("4", arb(4)), ("4.5", arb("4.5")), ("5", arb(5))):
        a, b = locate(pieces, keys, t)
        print(f"{label}: left {interval(a, 7)}, right {interval(b, 7)}, sharp factor {interval(2*b/a, 9)}")
    check_grid()


if __name__ == "__main__":
    main()
