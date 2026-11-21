"""Certify two prime-race windows and the finite zero regions that could overturn them."""
from dataclasses import dataclass

from flint import arb, ctx

from cell_check import even_moments
from paired_windows import full_mass
from thousand_window import grid_enclosures
from thousand_zeros import extended_table


@dataclass(frozen=True)
class Window:
    left_milli: int
    width_milli: int
    ceiling: int
    epsilon_units: int  # units of 10^-8
    log_top: int
    log_three_quarters: int
    head_percent: int
    odd_cost_milli: int


EARLIEST = Window(1620542278, 5315, 1620548, 4, 810330, 420000, 117, 22)
EXTENDED = Window(3895686221, 5299, 3895692, 2, 1947900, 1000000, 110, 30)
BRACKET_DEN = 10**14


def zero_tail(height, modulus):
    return 80*((2*modulus*height).log()+1)/(height-2)


def finite_region_errors(window):
    a = arb(2).log()
    x = arb(window.ceiling)
    left = arb(window.left_milli)/1000
    width = arb(window.width_milli)/1000
    assert left+width < x
    eps = arb(window.epsilon_units)/10**8
    cutoff = arb(10)**12
    power = arb(49)/50
    band = zero_tail(arb(8590), 11632)
    top = arb(window.log_top).exp()
    high = zero_tail(top, 11632)
    assert band < arb(19)/100
    unrestricted = (1000*(x/2).exp()+50)*high
    assert unrestricted < arb(1)/10**8
    middle = (2000*eps*(eps*x).exp()+50*((eps*x).cosh()-1))*band
    log_scale = (2*11632*cutoff).log()+1
    weighted = 160/cutoff*(log_scale/(1-power)+power/(1-power)**2)
    upper_band = (1000*(eps*x).exp()+50)*weighted
    assert middle+upper_band+unrestricted < arb(window.odd_cost_milli)/1000
    assert arb(1)/2+eps+power*(window.log_three_quarters-cutoff.log())/x > arb(3)/4

    extent = 3*a if window == EARLIEST else a
    moment = (2000*eps*(eps*extent).exp()+40*((eps*extent).cosh()-1))*band
    moment += (40*(extent/2).exp()+40)*zero_tail(cutoff, 11632)
    assert moment < arb(window.epsilon_units)/(2*10**5)

    base_q = 16*547*727
    base_t = arb(40).exp()
    delta = 1/(20*(base_q*base_t).log())
    base_low = 2*10**6*(-delta*left).exp()
    base_high = 2000*zero_tail(base_t, base_q)
    base_real = 8*20/(arb(1).exp()*left)
    higher = 10**6*(-left/10).exp()
    deleted = 10**6*(-left/3).exp()
    trivial = 128*(-(left-2)/2).exp()/(1-(-(left-2)).exp())
    assert base_low+base_high+base_real+higher+deleted+trivial < arb(1)/1000
    factor = 2*width.exp()
    lower = 406729 if window == EARLIEST else 400273
    assert arb(lower)/1000 < factor < arb(lower+1)/1000
    print(f"At log centre {window.left_milli}/1000 the odd displacement is below {window.odd_cost_milli}/1000, the remaining arithmetic costs less than 0.001, and the factor lies strictly between {lower}/1000 and {lower+1}/1000.")


def higher_even_moments(a, r):
    # At 2a and 3a the gamma integrand has support separated from zero.
    c = 2*(r-1)**2
    divisor = r*c
    g0 = 4*(a/2).sinh()/divisor
    steps = 32768
    dx = 2*a/steps
    values = []
    for k in (2, 3):
        integral = arb(0)
        for j in range(steps):
            t = (k-1)*a + arb(2*j+1, 1)*dx/2
            g = 4*((a-abs(k*a-t)).max(0)/2).sinh()/divisor
            integral += dx*(-t/2).exp()*g/(2*t.sinh())
        prime = -a*g0 if k == 2 else arb(0)
        values.append(prime-8*integral)
    return values


def make_terms(a, r):
    terms = []
    for _, _, lo, hi in extended_table():
        g = arb(lo+hi, hi-lo)/(2*BRACKET_DEN)
        sine = (a*g/2).sin()
        even = ((r-1)**2+4*r*sine**2)/(2*(r-1)**2*(arb(1)/4+g*g))
        terms.append((g, sine, even, -4*even*sine/a))
    return terms


def real_integral_budgets(terms, a, r):
    v0, v1 = even_moments(a, r)
    m0, m1 = v0, (v0-v1)/2
    v2, v3 = higher_even_moments(a, r)
    # Power reduction expanded into a cosine polynomial: 71,-42,-6,-1 over 2.
    cubic = (71*v0-42*v1-6*v2-v3)/(54*a)
    # Recompute using power moments before relying on the expanded coefficients.
    second = (3*v0-4*v1+v2)/8
    third = (10*v0-15*v1+6*v2-v3)/32
    power_form = 2*(11*v0+75*m1-48*second+16*third)/(54*a)
    assert cubic.overlaps(power_form)
    for _, sine, even, _ in terms:
        u = sine**2
        major = (16*u**3-48*u**2+75*u+11)/54
        cubic -= 4*even*major/a
        m0 -= 2*even
        m1 -= 2*even*u
    small = arb(1)/10**5
    assert m0 > 0 and m1 > 0
    assert 2*((m0+small)*(m1+small)).sqrt()/a < arb(831)/1000
    assert cubic + (71+42+6+1)*2*small/(54*a) < arb(751)/1000
    print("Real gamma integrals at four shifts give the earlier projected tail below 0.751. Two shifts independently give the other window's tail below 0.831.")


def direct_cover(window, terms):
    cap = arb(window.head_percent)/100
    def cover(i, depth):
        x = (arb(window.left_milli)+arb(2*i+1,1)*window.width_milli/(2*(1<<depth)))/1000
        value = sum((b*(g*x).sin() for g, _, _, b in terms), arb(0))
        if value < cap:
            return 1
        if depth == 30:
            raise ArithmeticError(f"Unresolved cell {i}/2^{depth}")
        return cover(2*i, depth+1)+cover(2*i+1, depth+1)
    count = sum(cover(i, 9) for i in range(512))
    print(f"The retained true frequencies are below {window.head_percent}/100 on all {count:,} closed cells starting at {window.left_milli}/1000.")


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    ctx.prec = 144
    a, r = arb(2).log(), arb(2).sqrt()
    for window in (EARLIEST, EXTENDED):
        finite_region_errors(window)
        left, width = arb(window.left_milli)/1000, arb(window.width_milli)/1000
        nodes, taylor, interp, shift, removed = grid_enclosures(left, width)
        assert all(v+taylor+interp+shift < arb(window.head_percent)/100 for v in nodes)
        if window == EARLIEST:
            closed_tail = full_mass(a, r, 2*(r-1)**2)-removed
            assert closed_tail+(71+42+6+1)*arb(2)/10**5/(54*a) < arb(751)/1000
            print("Closed gamma terms and the projected cubic majorant give the earlier tail below 0.751.")
        print(f"The interval Fourier transform covers 32,769 nodes, their connecting intervals and frequency uncertainty at {window.left_milli}/1000.")
    ctx.prec = 120
    terms = make_terms(a, r)
    real_integral_budgets(terms, a, r)
    for window in (EARLIEST, EXTENDED):
        direct_cover(window, terms)


if __name__ == "__main__":
    main()
