"""Cover the candidate loss-free window and bound every uncomputed contribution."""

from math import factorial
import sys
from flint import acb, arb, ctx
from onset_spectrum import DEN, FAMILY, table, terms
from splitting_square_bound import CAP, correction

GRID=8192
DEGREE=16


def majorant(u):
    return (((16*u-48)*u+75)*u+11)/54


def full_mass(a, r, c):
    g0 = 4*(a/2).sinh()/(r*c)
    constant = sum((q/arb.pi()).log() for q,_ in FAMILY)-8*(arb.const_euler()+a)

    def gamma_sum(k):
        if k == 0:
            return a
        z = arb(1)/2**abs(k)
        return ((1+z)*(1+z).log()+(1-z)*(1-z).log())/(2*z.sqrt())

    def cosine_moment(k):
        # The sole nonzero interior prime-power term occurs at n=4 and k=2.
        local = constant*g0 if k == 0 else -a*g0 if k == 2 else arb(0)
        arch = r*(gamma_sum(k+1)+gamma_sum(k-1))-3*gamma_sum(k)
        return local-8*arch/c

    v0, v1, v2, v3 = (cosine_moment(k) for k in range(4))
    m1 = (v0-v1)/2
    m2 = (3*v0-4*v1+v2)/8
    m3 = (10*v0-15*v1+6*v2-v3)/32
    return 2*(11*v0+75*m1-48*m2+16*m3)/(54*a)


def grid_enclosures(left, width):
    """Return rigorous signed nodes and separate errors for a specified interval."""
    a,r = arb(2).log(),arb(2).sqrt()
    c = 2*(r-1)**2
    centre = left+width/2
    right = left+width
    removed, shift, curvature, mass = (arb(0) for _ in range(4))
    # All ordinate bins belong to the same equally spaced frequency lattice.
    bins = [{} for _ in range(DEGREE+1)]
    for _,_,lo,hi in table():
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


def endpoints():
    return arb(400936)/10000, arb(21973)/10000


def whole_cells():
    """Use the true ordinate boxes directly throughout each closed dyadic cell."""
    ctx.prec = 112
    left,width = endpoints()
    modes = list(terms('split'))
    def visit(i, depth):
        x = left + arb(2*i+1,1)*width/(2 << depth)
        lower = sum((b*(t*x).sin() for t,b in modes),arb(0))
        if lower > -arb(1015)/1000:
            return 1
        if depth >= 25:
            raise ArithmeticError(f"Unresolved cell {i} at depth {depth}")
        return visit(2*i,depth+1) + visit(2*i+1,depth+1)
    assert sum(visit(i,8) for i in range(256)) >= 256
    print("Direct closed-cell arithmetic puts the retained true-ordinate sum above -1.015 throughout the proposed loss-free interval.")


def compare_with_actual():
    ctx.prec = 144
    a,r = arb(2).log(),arb(2).sqrt()
    B = 2*r*(r-1)**2
    left,width = endpoints()
    right = left+width
    assert width > arb(9).log()
    nodes,taylor,interp,shift,removed = grid_enclosures(left,width)
    assert all(s-taylor-interp-shift > -arb(1015)/1000 for s in nodes)
    total = full_mass(a,r,2*(r-1)**2)
    tail = total-removed
    assert tail < arb(75)/100

    top = arb(5000000)
    carrier = 197*top/200
    length = 25/carrier
    variance = length/carrier/(carrier*length).tanh()-1/carrier**2
    assert 0 < variance < arb(1)/10**12
    ymin,ymax = (left-3*a/2).exp(),(right-3*a/2).exp()
    lognorm = (8*ymax).log()+length
    # Count either sign of the degree-one coefficient. Taking one sign at a
    # time bounds either side of a signed smooth term or a signed slope jump.
    wheel = (3,5,7,11,13,17,19,23,29,31,37)
    classes = 2*273*363
    density = arb(classes)/(16*547*727)
    for p in wheel:
        density *= arb(p-1)/p
    assert ymin*(-length).exp() > 37**4
    # All remaining prime norms are integer squares or fourth powers.
    density += 1/(2*(ymin*(-length).exp()).sqrt())
    edge = classes*2**len(wheel)+1
    leading = 52*density*(1+r)*(3*length/2).exp()*ymax.sqrt()*variance
    boundary = edge*(1+1/r)*(length/2).exp()/ymin.sqrt()
    raw = lognorm/(B*a)*(leading+boundary*(variance+24*variance.sqrt()))
    assert raw < arb(142)/1000
    assert left-length-3*a/2 > 39 and right+length-3*a/2 < 42
    assert 2*r*((right+length-3*a/2)/2).exp()<CAP
    arithmetic = correction(39)
    assert arithmetic < arb(42)/1000
    # Integrate the positive exterior Poisson sum over a four-unit shoulder
    # on both sides of each counting endpoint (16 half-lines in the family).
    count = (1+4/top)**2*8*((2*11632*top).log()+10)/((top-4)*2*arb(2).atan())
    Gstrip = (2*r*(a/2).cosh()+3)/(2*(r-1)**2)*(1+1/top**2)
    moment_weight = Gstrip*(71+42*(a/2).cosh()+6*a.cosh()+(3*a/2).cosh())/(54*a)
    assert moment_weight < 65 and 65*count < arb(17)/10000
    barrier = (top**2-carrier**2-arb(1)/4).sqrt()
    # The strip is |Re z| <= 1/2.  For |Im z| >= top the bound on
    # |Im z|/(2 sqrt((Im z)^2-carrier^2-1/4)) decreases with |Im z|.
    assert top/(2*barrier)<3
    sphere = carrier/barrier*(3*length).cosh()/(carrier*length).sinh()
    high = 1000*(right/2).exp()*sphere*count
    trivial = 48*(-3*(left-length-2)/2).exp()/(1-(-2*(left-length-2)).exp())
    residual = sum((abs(b)*t*t*variance/2 for t,b in terms('direct')),arb(0))
    assert high+trivial+residual < arb(63)/10000
    budget = arb(1015)/1000+tail+raw+arithmetic+65*count+high+trivial+residual
    assert budget < arb(196)/100 < arb(19735)/10000
    for h,lo,hi in ((arb(400507)/10000,875,876),(left,913,915),(arb(4133)/100,3146,3147)):
        start = (h-3*a/2).exp()/10**14
        assert lo < start < hi
    assert (left-arb(400507)/10000).exp()<arb(1044)/1000
    assert left-arb(40041)/1000 < (arb(33)/4).log()
    assert arb(180013)/10000 < 2*width.exp() < arb(180014)/10000
    assert arb(42351)/1000-(arb(33)/4).log() < arb(4133)/100
    first_floor = (arb(400507)/10000-3*a/2).exp()
    assert (arb(399)/20)/first_floor < arb(228)/10**18
    assert 2*width.exp()-arb(228)/10**18 > arb(180013)/10000
    first = (left-3*a/2).exp().ceil().unique_fmpz()
    last = ((right-3*a/2).exp()/9).floor().unique_fmpz()
    assert first is not None and last is not None and last-first+1 > 6800000000000
    print("The interval transform also proves the 13,946 retained modes exceed -1.015 on the whole width 2.1973, which exceeds log 9.")
    print("If the eight factors have no off-line zero between ordinates 8,590 and 5,000,000 in absolute value, the actual normalized score exceeds -1.96 throughout the window. Zeros above 5,000,000 are unrestricted.")
    print(f"Under that finite condition all integer starts from {first} through {last} miss the loss at factor 18, a run of {last-first+1}.")
    print("The possible log onsets span less than the shortest tested gap. All factors above the sharp height-one factor through 2 exp(2.1973), between 18.0013 and 18.0014, have the same first missed height under the finite condition.")


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    if sys.argv[1:] == ['cells']:
        whole_cells()
    elif sys.argv[1:] == []:
        compare_with_actual()
    else:
        raise SystemExit("Usage: onset_window.py [cells]")
