"""Apply the four even spectral moments to positive triangular averages."""
from flint import arb

from early_400 import make_terms


DELTA = arb(3)/500
CONDUCTORS = (547,4376,8752,8752,727,5816,11632,11632)


def gamma_at(t, a):
    if t == 0:
        return a
    assert t > 0
    z = (-t).exp()
    m = -(-t).expm1()
    return ((1+z)*(1+z).log()+m*m.log())/(2*z.sqrt())


def filtered_tail(terms,a,r):
    c = 2*(r-1)**2
    kappa = sum((arb(q)/arb.pi()).log() for q in CONDUCTORS)-8*(arb.const_euler()+a)

    def gamma(u):
        # The positive exponential series decreases with |u|, even at zero.
        d = abs(u)
        low = d.lower()
        if not low > 0:
            low = arb(0)
        return gamma_at(d.upper(),a).union(gamma_at(low,a))

    def tent(u):
        return 4*((a-abs(u)).max(0)/2).sinh()/(r*c)

    def moment(u):
        result = kappa*tent(u)-8*(r*(gamma(u+a)+gamma(u-a))-3*gamma(u))/c
        for j in (2,4):
            result -= 2*a/(2**(j//2))*(tent(u-j*a)+tent(u+j*a))
        return result

    cells = 2048
    values = [arb(0) for _ in range(4)]
    for i in range(cells):
        u = DELTA*arb(2*i+1,1)/(2*cells)
        weight = arb(2*cells-2*i-1)/(2*cells*cells)
        for k in range(4):
            values[k] += weight*(moment(k*a+u)+moment(k*a-u))
    tail = (71*values[0]-42*values[1]-6*values[2]-values[3])/(54*a)
    filtered = []
    for g,sine,even,b in terms:
        weight = (g*DELTA/2).sinc()**2
        u = sine**2
        major = (16*u**3-48*u**2+75*u+11)/54
        tail -= 4*even*major*weight/a
        filtered.append((g,-b*weight))
    assert tail + 120*arb(15)/10**5/(54*a) < arb(2)/100
    return filtered
