"""Isolate the physical quartic at every vertex of the four-dimensional target box."""

from csv import DictReader
from itertools import product
from pathlib import Path

from flint import arb, arb_mat, arb_poly, ctx
from three_inverse import bands, hull
from three_profiles import HOSTS, Trace
from three_starts import charts as read_charts, integer_anchor

ROOT = Path(__file__).resolve().parent
SCALE = 10**18


def isolate(traces, row, charts):
    target = [arb(row[key]) for key in ("level_A", "level_B", "move_A", "move_B")]
    data = charts[int(row["row"])-1]
    centre,k = integer_anchor(traces,data[0])
    assert k == int(row["anchor"])
    raw = [centre[0]+1, centre[1]+1, centre[2], centre[3]]
    starts = [arb(row[key]) for key in ("start_1", "start_2")]
    windows = [(x-arb(1)/10000, x+arb(1)/10000) for x in starts]
    assert all(8*hi < traces[0].limit for _,hi in windows)
    columns = [v for lo,hi in windows for v in bands(traces,lo,hi)]
    matrix = arb_mat([list(row) for row in zip(*columns)])
    direction = matrix.solve(arb_mat([[x] for x in raw]))
    base = matrix.solve(arb_mat([[x] for x in target]))
    factors = [arb_poly([base[i,0],-direction[i,0]]) for i in range(4)]
    u,v,s,t = factors
    q1,q2 = u*v,s*t
    rest = arb_poly([1,-1])
    f = (rest**2-q1-q2)**2-4*q1*q2
    assert f.degree()==4 and not f[4].contains(0)
    seed = arb(row["anchor_weight"])
    lo,hi = seed-arb(1)/10**6, seed+arb(1)/10**6
    derivative = f.derivative()(hull(lo,hi))
    sign = 1 if derivative>0 else -1
    assert sign*derivative>0 and sign*f(lo)<0<sign*f(hi)
    for _ in range(115):
        mid = (lo+hi)/2
        val = sign*f(mid)
        if val>0:
            hi=mid
        else:
            assert val<0
            lo=mid
    h = hull(lo,hi)
    assert all(term(h)>0 for term in factors)
    assert rest(h)>0 and (rest**2-q1-q2)(h)>0
    y,z = u(h)/v(h),s(h)/t(h)
    b = (rest(h)**2+q1(h)-q2(h))/(2*rest(h))
    c = rest(h)-b
    assert all(l<x<r for x,(l,r) in zip((y,z),windows))
    for x,(l,r) in zip((y,z),((data[1],data[2]),(data[3],data[4]))):
        centred_log = (2*arb(2).sqrt()*x).log()*1000
        assert l < centred_log < r
    assert min(b,c,h)>arb(1)/400
    ordered=sorted((y,arb(k),z))
    assert ordered[0]>1155 and 2*ordered[-1]<701816 and 2*ordered[-1]<110*ordered[0]
    assert all(20*y>41*x for x,y in zip(ordered,ordered[1:]))
    floors=[]
    for x in (y,z,h,b,c):
        lower=(x*SCALE).floor().unique_fmpz()
        assert lower is not None and arb(lower)/SCALE<x<arb(lower+1)/SCALE
        floors.append(str(lower))
    return floors


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=320
    traces=[Trace(host,2_200_000,physical=True) for host in HOSTS]
    with (ROOT/"evidence/corner-seeds.tsv").open() as source:
        seeds=list(DictReader(source,delimiter="\t"))
    columns=("level_A","level_B","move_A","move_B")
    targets=[tuple(row[key] for key in columns) for row in seeds]
    assert len(targets)==16 and set(targets)==set(product(("-9","11"),("-7","9"),("-3.5","3.5"),("-2.2","2.2")))
    charts=read_charts()
    print("level_A\tlevel_B\tmove_A\tmove_B\tanchor\tstart_1_floor\tstart_2_floor\tanchor_weight_floor\tweight_1_floor\tweight_2_floor")
    for row,t in zip(seeds,targets):
        print("\t".join([*t,row["anchor"],*isolate(traces,row,charts)]))


if __name__=="__main__":
    main()
