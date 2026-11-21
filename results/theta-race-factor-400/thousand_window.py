"""Enclose the full interval using a Taylor expansion in frequency and an interval DFT."""
from math import factorial
from flint import acb, arb, ctx
from paired_windows import full_mass, majorant
from thousand_zeros import extended_table

DEN = 10**14
LEFT = 3855250213808
WIDTH = 6235
GRID = 32768
DEGREE = 16


def grid_enclosures(left, width):
    """Return rigorous signed nodes and separate errors for a specified interval."""
    a,r = arb(2).log(),arb(2).sqrt()
    c = 2*(r-1)**2
    centre = left+width/2
    right = left+width
    removed, shift, curvature, mass = (arb(0) for _ in range(4))
    # All ordinate bins belong to the same equally spaced frequency lattice.
    bins = [{} for _ in range(DEGREE+1)]
    for _,_,lo,hi in extended_table():
        midpoint = arb(lo+hi)/(2*DEN)
        radius = arb(hi-lo)/(2*DEN)
        box = arb(lo+hi,hi-lo)/(2*DEN)
        def kernel(t):
            return (3-2*r*(a*t).cos())/(c*(arb(1)/4+t*t))
        sin_box = (a*box/2).sin()
        actual = -4*kernel(box)*sin_box/a
        beta = -4*kernel(midpoint)*(a*midpoint/2).sin()/a
        removed += 4*kernel(box)*majorant(sin_box**2)/a
        shift += abs(actual-beta)+abs(beta)*right*radius
        curvature += abs(beta)*midpoint**2
        mass += abs(beta)
        j = round(float(midpoint*width/(2*arb.pi())))
        delta = midpoint*width/2-j*arb.pi()
        assert abs(delta) < arb(8)/5
        i = j % GRID
        coefficient = beta*acb(0,midpoint*centre).exp()
        for order, terms in enumerate(bins):
            terms[i] = terms.get(i,acb(0))+coefficient
            coefficient *= acb(0,delta)/(order+1)
    # With negative DFT exponent, entry -k represents the required positive phase.
    values = [arb(0) for _ in range(GRID+1)]
    half = GRID//2
    for terms in reversed(bins):
        vector = [acb(0)]*GRID
        for index, coefficient in terms.items():
            vector[index] = coefficient
        transformed = acb.dft(vector)
        for i in range(GRID+1):
            k = i-half
            values[i] = values[i]*(arb(k)/half)+transformed[(-k)%GRID].imag
    taylor = mass*(arb(8)/5)**(DEGREE+1)/factorial(DEGREE+1)
    interpolation = curvature*(width/GRID)**2/8
    return values,taylor,interpolation,shift,removed


def main():
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    ctx.prec = 144
    a,r = arb(2).log(),arb(2).sqrt()
    c = 2*(r-1)**2
    width = arb(WIDTH)/1000
    values,taylor,interpolation,shift,removed = grid_enclosures(arb(LEFT)/1000,width)
    total = full_mass(a,r,c)
    tail = total-removed
    # Check every ball individually, including the two different boundary nodes.
    assert taylor < arb(1)/10**7
    assert all(v-taylor > -arb(106)/100 for v in values)
    assert interpolation < arb(1)/100
    assert shift < arb(14)/100
    assert tail < arb(75)/100
    assert all(v-taylor-interpolation-shift-tail > -arb(194)/100 for v in values)
    assert 2*width.exp()>1020 and width>arb(500).log()+arb(1)/50
    assert total < 346 and 346*2*arb.pi()/2**22 < arb(1)/1000
    assert 22*len(extended_table()) == 306812
    print('Four exact spectral moments bound the omitted absolute mass below 0.75. Ordinate uncertainty contributes less than 0.14 on the whole interval.')
    print('32,769 signed nodes exceed -1.06 after the frequency Taylor remainder. Interpolation of the midpoint sum costs less than 0.01 and the continuous profile exceeds -1.94.')
    print('A 22-bit phase mesh in 13,946 coordinates changes the retained sum by less than 0.001. The factor-1000 interval has more than 0.02 log units of slack.')
    print('The rational log width 6.235 exceeds log(500), and twice its exponential exceeds 1,020.')


if __name__ == '__main__':
    main()
