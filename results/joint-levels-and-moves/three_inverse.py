"""Isolate a real quartic root that gives three exact balanced common readings."""

from bisect import bisect_right
from flint import arb, arb_mat, arb_poly, ctx
from three_profiles import HOSTS, Trace, integer_floor
from three_starts import integer_anchor


def branch(trace, k, lo, hi):
    cuts = []
    for multiplier in (k/2,k,2*k):
        left = bisect_right(trace.norms, integer_floor(multiplier*lo))
        right = bisect_right(trace.norms, integer_floor(multiplier*hi))
        assert left == right
        cuts.append(left)
    i,j,m = cuts
    r,B = trace.root,trace.B
    p,q = trace.plain,trace.divided
    first = -2*k.sqrt()/B*(r*(q[m]-q[j])-(q[j]-q[i])/r)
    second = -2/(B*k.sqrt())*(r*(p[j]-p[i])-(p[m]-p[j])/r)
    return first,second


def bands(traces, lo, hi):
    levels, moves = [], []
    for trace in traces:
        levels.append(branch(trace,2*trace.root,lo,hi))
        low,high = (branch(trace,arb(k),lo,hi) for k in (2,4))
        moves.append(tuple((q-p)/trace.a for p,q in zip(low,high)))
    return list(zip(*(levels+moves)))


def hull(lo,hi):
    return (lo+hi)/2 + arb(0,1)*(hi-lo)/2


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 224
    traces = [Trace(host,500_000,physical=True) for host in HOSTS]
    anchor,k = integer_anchor(traces,10548)
    assert k == 13471
    raw = [anchor[0]+1,anchor[1]+1,anchor[2],anchor[3]]
    windows = ((arb("2457.4526"),arb("2457.4529")),(arb("61460.7129"),arb("61460.7134")))
    columns = [v for lo,hi in windows for v in bands(traces,lo,hi)]
    matrix = arb_mat([list(row) for row in zip(*columns)])
    base = matrix.solve(arb_mat([[1],[1],[0],[0]]))
    slope = matrix.solve(arb_mat([[z] for z in raw]))
    affine = [arb_poly([base[j,0],-slope[j,0]]) for j in range(4)]
    u,v,s,t = affine
    q1,q2 = u*v,s*t
    remainder = arb_poly([1,-1])
    polynomial = (remainder**2-q1-q2)**2-4*q1*q2
    assert polynomial.degree() == 4 and not polynomial[4].contains(0)
    lo,hi = arb("0.394907"),arb("0.394909")
    dy = polynomial.derivative()(hull(lo,hi))
    orientation = 1 if dy > 0 else -1
    assert dy*orientation > 0 and polynomial(lo)*orientation < 0 < polynomial(hi)*orientation
    for _ in range(110):
        middle = (lo+hi)/2
        value = orientation*polynomial(middle)
        if value < 0:
            lo = middle
        else:
            assert value > 0
            hi = middle
    weight = hull(lo,hi)
    uv,st = q1(weight),q2(weight)
    assert all(f(weight)>0 for f in affine)
    assert (remainder**2-q1-q2)(weight)>0 and remainder(weight)>0
    y1,y2 = u(weight)/v(weight),s(weight)/t(weight)
    mu = (remainder(weight)**2+uv-st)/(2*remainder(weight))
    nu = remainder(weight)-mu
    assert all(l < y < h for y,(l,h) in zip((y1,y2),windows))
    assert weight>arb(1)/5 and mu>arb(1)/5 and nu>arb(1)/5
    assert (y1*2 < k and 2*k < y2)
    for label,value in zip(("first movable start","second movable start","integer-start weight",
                            "first movable weight","second movable weight"),(y1,y2,weight,mu,nu)):
        scale = 10**18
        low = (value*scale).floor().unique_fmpz()
        assert low is not None and arb(low)/scale < value < arb(low+1)/scale
        print(f"{label} in ({low}, {low+1}) / 10^18")
    print("Integer start 13471; one simple quartic root gives both levels exactly one and both moves exactly zero.")


if __name__ == "__main__":
    main()
