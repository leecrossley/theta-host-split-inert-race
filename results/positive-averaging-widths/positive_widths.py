"""Cosine separation and verified positive quadrature at the actual zero ordinates."""
import json
from flint import arb, arb_mat, ctx
from band_zeros import ROOT, TICKS, ball, rows

TARGETS = (('547/727',88436,88440,9),('571/599',47698,47702,7),
           ('both',103945,103949,16))


def designs():
    raw,_ = rows()
    data = json.loads((ROOT/'evidence/design.json').read_text())
    assert len(data)==len(TARGETS)
    for field,((name,lo,hi,count),entry) in enumerate(zip(TARGETS,data)):
        assert entry['field']==name and entry['width_milli']==[lo,hi]
        intervals = sorted((a,b) for f,_,_,a,b in raw if b<TICKS and (field==2 or f==field))
        assert len(intervals)==count and len(entry['dual_pico'])==count
        ts = [ball(it) for it in intervals]
        dual = [arb(n)/10**12 for n in entry['dual_pico']]
        ticks = entry['atom_nano']
        assert ticks==sorted(set(ticks)) and len(ticks)==count+1
        assert 0<=ticks[0] and 2*ticks[-1]<=hi*10**6
        atoms = [arb(n)/10**9 for n in ticks]
        yield name,lo,hi,ts,dual,atoms


def matrix(ts,atoms):
    return arb_mat([[arb(1) for _ in atoms]]+
                   [[(t*x).cos() for x in atoms] for t in ts])


def lower(lo,ts,dual):
    r = arb(lo)/2000
    grid = 16384
    step = r/grid
    curvature = sum(abs(c)*t*t for c,t in zip(dual,ts))
    # The direct sine/cosine values below bound every mesh vertex separately.
    allowance = curvature*step**2/8
    floor = arb(1)/10**6
    for i in range(grid+1):
        x = i*step
        y = sum(c*(t*x).cos() for c,t in zip(dual,ts))
        assert y-allowance > floor
    assert sum(abs(c) for c in dual) < arb(1000001)/1000000


def upper(ts,atoms):
    n = len(atoms)
    weights = matrix(ts,atoms).solve(arb_mat([[int(i==0)] for i in range(n)]), algorithm='precond')
    assert all(weights[i,0] > arb(1)/1000 for i in range(n))


def main():
    ctx.prec = 176
    for name,lo,hi,ts,dual,atoms in designs():
        lower(lo,ts,dual)
        upper(ts,atoms)
        print(f'{name}: a signed cosine sum exceeds 1/1,000,000 throughout half the lower width {lo}/1000. At width {hi}/1000 each of {len(atoms)} exact symmetric-pair weights exceeds 1/1000.')
        if name != 'both':
            level,move = {'547/727':(25000,1214),'571/599':(9460,1162)}[name]
            t = ts[0]
            a,r = arb(2).log(),arb(2).sqrt()
            g = (3-2*r*(a*t).cos())/(2*(r-1)**2*(t*t+arb(1)/4))
            assert 0 < 10*t < arb.pi()/2
            radius = r*g*(10*t).cos()
            assert radius > arb(level)/10000 and radius*2*(a*t/2).sin()/a > arb(move)/10000
            print(f'{name}: at log width 20 the private first mode alone forces level and move RMS above {level}/10000 and {move}/10000.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    main()
