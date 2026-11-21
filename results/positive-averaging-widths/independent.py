"""Check separation on adaptive closed cells and solve by rational preconditioning."""
from flint import arb, arb_mat, fmpq, fmpq_mat, ctx
from positive_widths import designs


def separating_cover(radius,ts,c):
    bound = sum(abs(v)*t*t for v,t in zip(c,ts))
    stack = [(arb(0),radius,0)]
    leaves = 0
    while stack:
        a,b,depth = stack.pop()
        m,h = (a+b)/2,(b-a)/2
        trig = [(t*m).sin_cos() for t in ts]
        val = sum(coeff*sc[1] for coeff,sc in zip(c,trig))
        slope = sum(-coeff*t*sc[0] for coeff,t,sc in zip(c,ts,trig))
        if val-abs(slope)*h-bound*h*h/2 > arb(1)/10**6:
            leaves += 1
            continue
        assert depth < 24
        stack.extend(((a,m,depth+1),(m,b,depth+1)))
    return leaves


def rational_midpoint(x):
    m,e = x.mid().man_exp()
    return fmpq(m << e) if e>=0 else fmpq(m,1 << -e)


def weights_by_defect(ts,atoms):
    # The inverse of a rational point matrix provides a fixed preconditioner.
    n = len(atoms)
    A = arb_mat([[arb(1)]*n]+[[(x*t).cos() for x in atoms] for t in ts])
    B = fmpq_mat([[rational_midpoint(A[i,j]) for j in range(n)] for i in range(n)])
    P = B.inv()
    inv = arb_mat([[arb(P[i,j]) for j in range(n)] for i in range(n)])
    pre = inv*A
    err = [sum(abs(pre[i,j]-int(i==j)) for j in range(n)) for i in range(n)]
    assert all(x < arb(1)/10**6 for x in err)
    # If PA = I+E, the solution differs from P e0 by at most
    # ||E||/(1-||E||) times ||P e0|| in the maximum norm.
    assert all(abs(inv[i,0]) < 1 for i in range(n))
    assert all(inv[i,0]-arb(1)/999999 > arb(1)/1000 for i in range(n))


def main():
    ctx.prec = 192
    for name,lo,hi,ts,dual,atoms in designs():
        count = separating_cover(arb(lo)/2000,ts,dual)
        weights_by_defect(ts,atoms)
        print(f'{name}: {count} Taylor cells independently separate every shorter filter. The rational inverse has row-sum defect below one millionth and gives strictly positive exact weights.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    main()
