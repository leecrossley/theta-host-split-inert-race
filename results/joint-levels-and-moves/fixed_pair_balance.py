"""Certify the best simultaneous balance with the starts 2457 and 13471 fixed."""

from pathlib import Path
from flint import arb, arb_mat, arb_poly, ctx
from three_profiles import HOSTS, Trace, vector
from three_inverse import bands, hull

DATA = Path(__file__).resolve().parent/"evidence"
SIGNS = (-1,1,1,-1)
OUTSIDE = ((1,"1228.5"), (4914,"6735.5"), (26942,"61455.601"), ("61455.619",300000))
UPPER = "0.00227436"


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def observation(traces,y):
    z = 2*arb(2).sqrt()*y
    pairs = [(t.level(z,arb(0)), (t.level(4*y,arb(0))-t.level(2*y,arb(0)))/t.a)
             for t in traces]
    return (pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1])


def optimal_pencil(A,B,u,v):
    """The two projected band quadratics and a whole target neighbourhood."""
    from inverse_monodromy import resultant
    difference=[b-a for a,b in zip(A,B)]
    shift=[c-a for c,a in zip((-1,-1,0,0),A)]
    coordinates=[]
    base=arb_mat([[difference[i] for i in (2,3)],[SIGNS[i] for i in (2,3)]])
    for i in (0,1):
        free=base.solve(arb_mat([[-difference[i]],[-SIGNS[i]]]))
        form=[arb(0)]*4
        form[i]=1
        form[2],form[3]=free[0,0],free[1,0]
        coordinates.append(arb_poly([dot(form,v),dot(form,shift),dot(form,u)]))
    assert all(p.degree()==2 and not p[2].contains(0) for p in coordinates)
    assert not resultant(*coordinates).contains(0)
    assert any(not p[1].contains(0) for p in coordinates)

    radius=arb(4)/10**9
    targets=[arb(0,radius) for _ in range(4)]
    near=hull(arb("61455.60975"),arb("61455.60990"))
    band=hull(arb("61455.60"),arb("61455.62"))
    centre=arb("61455.60982")
    const=(-1,-1,0,0)
    def tangent(y):
        return [(p-q/y)/(2*y.sqrt()) for p,q in zip(u,v)]
    def profile(y):
        # Centre the slope formula to avoid dependency of the two powers of y.
        mid=[c+p*centre.sqrt()+q/centre.sqrt() for c,p,q in zip(const,u,v)]
        return [m+t*(y-centre) for m,t in zip(mid,tangent(near))]
    def det(X):
        return arb_mat([list(row) for row in zip([a-p for a,p in zip(A,targets)],difference,X,SIGNS)]).det()
    sense=1 if det(tangent(near))>0 else -1
    assert sense*det(tangent(near))>0
    for y,sgn in (("61455.60975",-1),("61455.60990",1)):
        assert sgn*sense*det([c-a for c,a in zip(profile(arb(y)),A)])>0
    C=profile(near)
    matrix=arb_mat([[A[i]-C[i],B[i]-C[i],-SIGNS[i]] for i in range(3)])
    sol=matrix.solve(arb_mat([[targets[i]-C[i]] for i in range(3)]))
    assert min(sol[0,0],sol[1,0],1-sol[0,0]-sol[1,0])>arb(1)/5 and sol[2,0]>0
    for omit,direction in ((1,-1),(0,1)):
        idx=[i for i in range(4) if i!=omit]
        matrix=arb_mat([[A[i]-B[i] for i in idx],[A[i]-C[i] for i in idx],[SIGNS[i] for i in idx]])
        n=matrix.solve(arb_mat([[0],[0],[1]]))
        assert all(SIGNS[i]*n[j,0]>0 for j,i in enumerate(idx))
        assert direction*sum(n[j,0]*tangent(band)[i] for j,i in enumerate(idx))>0
    assert arb(2274350938622390)/10**18+2*radius<arb(UPPER)
    print("The unique optimum persists on the closed target box of radius 4/10^9.")
    print("The projected quadratics have nonzero resultant and a nonzero linear coefficient.")



def crossing(traces,A,B):
    lo,hi = arb("61455.60"),arb("61455.62")
    u,v = bands(traces,lo,hi)
    constant = (-1,-1,0,0)
    def determinant(column):
        return arb_mat([list(row) for row in zip(A,B,column,SIGNS)]).det()
    polynomial = arb_poly([determinant(v),determinant(constant),determinant(u)])
    assert polynomial.degree() == 2
    l,h=lo.sqrt(),hi.sqrt()
    sign = 1 if polynomial(l)<0 else -1
    assert sign*polynomial(l)<0<sign*polynomial(h)
    assert sign*polynomial.derivative()(hull(l,h))>0
    for _ in range(135):
        mid = (l+h)/2
        if sign*polynomial(mid) < 0:
            l=mid
        else:
            assert sign*polynomial(mid)>0
            h=mid
    q = hull(l,h)
    y = q*q
    C = [c+a*q+b/q for c,a,b in zip(constant,u,v)]
    matrix = arb_mat([[A[i]-C[i],B[i]-C[i],-SIGNS[i]] for i in range(3)])
    weights = matrix.solve(arb_mat([[-C[i]] for i in range(3)]))
    a,b,error = (weights[i,0] for i in range(3))
    assert min(a,b,1-a-b)>arb(1)/5 and 0<error<arb(UPPER)
    assert max(A[2],B[2]) < -10
    for omitted,orientation in ((1,-1),(0,1)):
        ids=[i for i in range(4) if i!=omitted]
        matrix=arb_mat([[A[i]-B[i] for i in ids], [A[i]-C[i] for i in ids],
                        [SIGNS[i] for i in ids]])
        sol=matrix.solve(arb_mat([[0],[0],[1]]))
        normal=[arb(0) for _ in SIGNS]
        for j,i in enumerate(ids):
            normal[i]=sol[j,0]
            assert SIGNS[i]*normal[i]>0
        band=hull(lo,hi)
        tangent=[(s-t/band)/(2*band.sqrt()) for s,t in zip(u,v)]
        assert orientation*dot(normal,tangent)>0
    optimal_pencil(A,B,u,v)
    return y,error,a,b,1-a-b


def outside_logs():
    scale = 10**10
    for lo,hi in OUTSIDE:
        x,z = [(2*arb(2).sqrt()*arb(t)).log()*scale for t in (lo,hi)]
        l,h = x.floor().unique_fmpz(),z.ceil().unique_fmpz()
        assert l is not None and h is not None
        yield arb(l)/scale,arb(h)/scale


def verify_separators(traces,A,B):
    roots = list(outside_logs())
    lines = (DATA/"fixed-pair-separators.tsv").read_text().splitlines()
    assert lines[0].split() == ["root","path","n0","n1","n2","n3"]
    records = {}
    for line in lines[1:]:
        region,path,*values = line.split()
        key = (int(region),"" if path=="-" else path)
        assert key not in records and 0<=key[0]<4 and set(key[1])<={"0","1"}
        records[key] = tuple(map(int,values))
        assert len(records[key])==4
    checked=0
    for j,(start,end) in enumerate(roots):
        stack=[("",start,end)]
        while stack:
            path,lo,hi=stack.pop()
            mid=(lo+hi)/2
            normal=records.pop((j,path),None)
            if normal is None:
                assert len(path)<48
                stack.extend(((path+"0",lo,mid),(path+"1",mid,hi)))
                continue
            threshold=arb(UPPER)*sum(map(abs,normal))
            assert threshold>0 and min(dot(normal,A),dot(normal,B))>threshold
            whole=vector(traces,mid,((hi-lo)/2).upper())
            assert dot(normal,whole)>threshold
            checked+=1
    assert not records
    return checked


def decimal_cell(label, value):
    scale=10**18
    n=(value*scale).floor().unique_fmpz()
    assert n is not None and arb(n)/scale<value<arb(n+1)/scale
    print(f"{label} in ({n}, {n+1}) / 10^18")


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=256
    traces=[Trace(host,3_000_000,True) for host in HOSTS]
    A,B=[observation(traces,arb(y)) for y in (2457,13471)]
    values=crossing(traces,A,B)
    for label,v in zip(("optimal third start","minimum worst error","weight 2457","weight 13471",
                         "weight at the third start"),values):
        decimal_cell(label,v)
    print("Direct prime count:",verify_separators(traces,A,B),
          "closed separating intervals outside the local band.")
    print("Unique minimizer over disjoint third doublings with starts from 1 through 300000.")


if __name__=="__main__":
    main()
