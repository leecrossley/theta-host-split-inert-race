"""A finite zero-location alternative for the earlier 400H obstruction."""
from flint import arb, ctx

from cell_check import even_moments
from thousand_window import grid_enclosures
from thousand_zeros import extended_table

LEFT = 3895686221
WIDTH = 5299
DEN = 1000
BRACKET_DEN = 10**14


def zero_tail(height, modulus):
    return 80*((2*modulus*height).log()+1)/(height-2)


def finite_strip_errors():
    a = arb(2).log()
    x = arb(3895692)
    eps = arb(2)/10**8
    band = zero_tail(arb(8590),11632)
    top = arb(1947900).exp()
    high = zero_tail(top,11632)
    assert band < arb(19)/100
    assert (1000*(x/2).exp()+50)*high < arb(1)/10**8
    middle = (2000*eps*(eps*x).exp()+50*((eps*x).cosh()-1))*band
    # Beyond cutoff the allowed displacement increases with log ordinate.
    cutoff = arb(10)**12
    power = arb(49)/50
    log_scale = (2*11632*cutoff).log()+1
    weighted = 160/cutoff*(log_scale/(1-power)+power/(1-power)**2)
    upper_band = (1000*(eps*x).exp()+50)*weighted
    assert middle+upper_band+arb(1)/10**8 < arb(3)/100
    assert arb(1)/2+eps+power*(arb(10)**6-cutoff.log())/x > arb(3)/4
    # At log shifts zero and log(2) the even moments have much smaller errors.
    moment = (2000*eps*(eps*a).exp()+40*((eps*a).cosh()-1))*band
    assert moment + (40*(a/2).exp()+40)*zero_tail(cutoff,11632) < arb(1)/10**5
    base_q = 16*547*727
    base_t = arb(40).exp()
    delta = 1/(20*(base_q*base_t).log())
    # This includes either sign of the ordinates in each zero count.
    base_low = 2*10**6*(-delta*arb(LEFT)/DEN).exp()
    base_high = 2000*zero_tail(base_t,base_q)
    base_real = 8*20/(arb(1).exp()*(arb(LEFT)/DEN))
    higher = 10**6*(-arb(LEFT)/(DEN*10)).exp()
    deleted = 10**6*(-arb(LEFT)/(DEN*3)).exp()
    trivial = 128*(-(arb(LEFT)/DEN-2)/2).exp()/(1-(-(arb(LEFT)/DEN-2)).exp())
    assert base_low+base_high+base_real+higher+deleted+trivial < arb(1)/1000
    assert 2*(arb(WIDTH)/DEN).exp() > arb(400273)/1000
    print("The finite region's odd displacement costs less than 0.03, the individual moment errors are below 0.00001, and the base and remaining arithmetic terms together cost less than 0.001.")
    print("Above exp(1000000) a zero capable of overturning the gain obstruction must have real part greater than 3/4, up to the cutoff exp(1947900).")


def direct_enclosures():
    a,r = arb(2).log(),arb(2).sqrt()
    v0,va = even_moments(a,r)
    m0,m1 = v0,(v0-va)/2
    data = []
    for _,_,lo,hi in extended_table():
        g = arb(lo+hi,hi-lo)/(2*BRACKET_DEN)
        sine = (a*g/2).sin()
        even = ((r-1)**2+4*r*sine**2)/(2*(r-1)**2*(arb(1)/4+g*g))
        m0 -= 2*even
        m1 -= 2*even*sine**2
        data.append((g,-4*even*sine/a))
    error = arb(1)/10**5
    assert m0 > 0 and m1 > 0
    assert 2*((m0+error)*(m1+error)).sqrt()/a < arb(831)/1000

    def cover(i,depth):
        x = (arb(LEFT)+arb(2*i+1,1)*WIDTH/(2*(1<<depth)))/DEN
        value = sum((b*(t*x).sin() for t,b in data), arb(0))
        if value < arb(11)/10:
            return 1
        if depth == 30:
            raise ArithmeticError(f"Unresolved cell {i}/2^{depth}")
        return cover(2*i,depth+1)+cover(2*i+1,depth+1)

    count = sum(cover(i,9) for i in range(512))
    print(f"True-frequency interval sines bound the retained sum below 1.1 on {count:,} closed cells; the omitted projected mass is below 0.831.")


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 144
    finite_strip_errors()
    nodes,taylor,interpolation,shift,_ = grid_enclosures(arb(LEFT)/DEN,arb(WIDTH)/DEN)
    assert all(v+taylor+interpolation+shift < arb(11)/10 for v in nodes)
    print("The independent interval DFT covers all 32,769 nodes, their connecting intervals and the ordinate errors with the same 1.1 bound.")
    ctx.prec = 120
    direct_enclosures()


if __name__ == "__main__":
    main()
