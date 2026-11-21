"""A second continuous check using two real gamma integrals and true-frequency cells."""
from flint import arb, ctx
from cell_check import even_moments
from thousand_zeros import extended_table

LEFT = 3855250213808
RIGHT = 3855250220043
DEN = 10**14


def main():
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    ctx.prec = 112
    a,r = arb(2).log(),arb(2).sqrt()
    m0, first = even_moments(a,r)
    m1 = (m0-first)/2
    terms = []
    for _,_,u,v in extended_table():
        t = (arb(u)/DEN).union(arb(v)/DEN)
        sine = (a*t/2).sin()
        g = ((r-1)**2+4*r*sine**2)/(2*(r-1)**2*(arb(1)/4+t*t))
        m0 -= 2*g
        m1 -= 2*g*sine**2
        terms.append((t,-4*g*sine/a))
    assert m0>0 and m1>0
    tail = 2*(m0*m1).sqrt()/a
    assert tail < arb(83)/100

    def head(x):
        return sum((weight*(t*x).sin() for t,weight in terms),arb(0))

    def cell(i,depth):
        x = (arb(LEFT)+arb(2*i+1,1)*(RIGHT-LEFT)/(2*(1<<depth)))/1000
        if head(x)-tail > -arb(197)/100:
            return 1
        if depth==28:
            raise ArithmeticError(f'Inconclusive interval {i}/2^{depth}')
        return cell(2*i,depth+1)+cell(2*i+1,depth+1)

    count = sum(cell(i,9) for i in range(512))
    high = head(arb(3855250215346)/1000)
    low = head(arb(3855250215954)/1000)
    assert high-tail > arb(121)/10
    assert low+tail < -arb(5)/100
    print(f'Two gamma integrals give a separate Cauchy-Schwarz tail below 0.83. Direct interval sines certify the profile above -1.97 on {count:,} closed cells.')
    print('Inside that same interval the profile exceeds 12.1 at log centre 3855250215.346 and falls below -0.05 at 3855250215.954.')


if __name__ == '__main__':
    main()
