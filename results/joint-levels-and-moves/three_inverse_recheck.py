"""A four-variable interval contraction for the balanced physical witness."""

from bisect import bisect_right
from flint import arb, arb_mat, ctx
from three_recheck import FIELDS, Reference


def original(refs,y):
    pairs = [ref.physical_at(y) for ref in refs]
    return [pairs[0][0]+1,pairs[1][0]+1,pairs[0][1],pairs[1][1]]


def conic(refs,y,extent):
    # Check a prime-free interval of moving cuts, then interpolate sqrt(y) times the scores.
    l,h = y-arb(1)/1000000,y+arb(1)/1000000
    for ref in refs:
        for scale in (arb(2),2*arb(2).sqrt(),arb(4)):
            for cut in (scale/2,scale,2*scale):
                lo,hi = (t.floor().unique_fmpz() for t in ((y-extent)*cut,(y+extent)*cut))
                assert lo is not None and hi is not None
                assert bisect_right(ref.n,int(lo)) == bisect_right(ref.n,int(hi))
    a,b = original(refs,l),original(refs,h)
    slope = [(v*h.sqrt()-u*l.sqrt())/(h-l) for u,v in zip(a,b)]
    intercept = [u*l.sqrt()-m*l for u,m in zip(a,slope)]
    return list(zip(slope,intercept))


def curve(c,y):
    t = y.sqrt()
    return [a*t+b/t for a,b in c]


def derivative(c,y):
    t = y.sqrt()
    return [(a-b/y)/(2*t) for a,b in c]


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 224
    refs = [Reference(host,True,500_000) for host in FIELDS]
    anchor = original(refs,arb(13471))
    bottom = (2457452761630441562944,61460713123523079975916,246867572221531610,358224372896804637)
    scale = 10**18
    radius = arb(1)/(2*scale)
    centres = [(arb(k)+arb(1)/2)/scale for k in bottom]
    boxes = [u+arb(0,1)*radius for u in centres]
    broad = [arb(1)/20,arb(2)/5,arb(1)/100000,arb(1)/100000]
    c,d = (conic(refs,t,extent) for t,extent in zip(centres[:2],broad[:2]))
    def data(v):
        y,z,s,t=v
        a,b=curve(c,y),curve(d,z)
        da,db=derivative(c,y),derivative(d,z)
        A=[q-u for q,u in zip(a,anchor)]
        B=[q-u for q,u in zip(b,anchor)]
        f=arb_mat([[u+s*m+t*n-k] for u,m,n,k in zip(anchor,A,B,(1,1,0,0))])
        jac=arb_mat([[s*x,t*y,m,n] for x,y,m,n in zip(da,db,A,B)])
        return f,jac
    f,jac=data(centres)
    # Use an exact dyadic midpoint preconditioner and verify it is nonsingular.
    inverse=jac.inv()
    inverse=arb_mat([[inverse[i,j].mid() for j in range(4)] for i in range(4)])
    assert not inverse.det().contains(0)
    _,over_box=data(boxes)
    identity=arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    defect=identity-inverse*over_box
    correction=inverse*f
    for i in range(4):
        lipschitz=sum(abs(defect[i,j]) for j in range(4))
        assert lipschitz < arb(1)/1000000
        assert abs(correction[i,0])+lipschitz*radius < radius
    fourth=1-(centres[2]-correction[2,0])-(centres[3]-correction[3,0])
    error=2*radius/1000000
    assert arb(394908054881663751)/scale < fourth-error
    assert fourth+error < arb(394908054881663752)/scale
    _,wide_jac = data([x+arb(0,1)*rad for x,rad in zip(centres,broad)])
    broad_defect = identity-inverse*wide_jac
    for i in range(4):
        variation = sum(abs(broad_defect[i,j])*broad[j] for j in range(4))
        shift = sum(abs(inverse[i,j]) for j in range(4))/10000
        assert variation < broad[i]/2
        assert abs(correction[i,0])+variation+shift < broad[i]
    print("Every target within 1/10000 of the balanced vector has a unique inverse in the larger specified box.")
    print("Finite-field splitting and two original tent evaluations per prime band:")
    print("the four physical equations give a strict contraction of the 18-place box, with norm < 1/1000000.")
    print("The integer-start weight agrees in its 18-place bracket; the zero-move and level-one witness is unique in this box.")


if __name__ == "__main__":
    main()
