"""Check the 16 corner witnesses by contractions of the original four equations."""

from csv import DictReader
from pathlib import Path
from flint import arb, arb_mat, ctx
from three_inverse_recheck import FIELDS, Reference, conic, curve, derivative, original

ROOT=Path(__file__).resolve().parent


def check(refs,row):
    scale=10**18
    one=arb(1)
    rad=one/(2*scale)
    names=("start_1_floor","start_2_floor","weight_1_floor","weight_2_floor")
    centres=[(arb(row[name])+one/2)/scale for name in names]
    target=[arb(row[name]) for name in ("level_A","level_B","move_A","move_B")]
    anchor=original(refs,arb(row["anchor"]))
    bands=[conic(refs,y,one/1000000) for y in centres[:2]]
    def eval_map(params):
        y,z,b,c=params
        p,q=(curve(band,x) for band,x in zip(bands,(y,z)))
        dp,dq=(derivative(band,x) for band,x in zip(bands,(y,z)))
        u,v=([x-a for x,a in zip(t,anchor)] for t in (p,q))
        f=arb_mat([[a+b*i+c*j-t] for a,i,j,t in zip(anchor,u,v,target)])
        jac=arb_mat([[b*i,c*j,k,l] for i,j,k,l in zip(dp,dq,u,v)])
        return f,jac
    f,jac=eval_map(centres)
    approx=jac.inv()
    pre=arb_mat([[approx[i,j].mid() for j in range(4)] for i in range(4)])
    assert not pre.det().contains(0)
    box=[x+arb(0,1)*rad for x in centres]
    _,dbox=eval_map(box)
    identity=arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    error=identity-pre*dbox
    correction=pre*f
    outs=[]
    for i in range(4):
        norm=sum(abs(error[i,j]) for j in range(4))
        assert norm<one/1000000 and abs(correction[i,0])+rad*norm<rad
        outs.append(centres[i]-correction[i,0]+arb(0,1)*rad*norm)
    h=1-outs[2]-outs[3]
    expected=arb(row["anchor_weight_floor"])
    assert expected/scale<h<(expected+1)/scale
    assert min(h,outs[2],outs[3])>one/100


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec=320
    refs=[Reference(host,True,2_200_000) for host in FIELDS]
    with (ROOT/"evidence/corner-enclosures.tsv").open() as source:
        rows=list(DictReader(source,delimiter="\t"))
    assert len(rows)==16
    for row in rows:
        check(refs,row)
    print("All 16 rational vertices have their own three-height solution, with every weight greater than 1/100.")
    print("At each vertex independent polynomial splitting and original-tent contraction certify all five 18-place brackets.")
    print("The contraction norm is below 1/1000000 at every corner.")


if __name__=="__main__":
    main()
