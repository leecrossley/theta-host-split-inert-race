"""Enclose both true signed profiles using a cubic majorant and grid interpolation."""
from flint import arb, ctx
from check_zeros import DEN, FAMILY, read_table

# Left endpoint and width in thousandths of a log unit. A positive sign means gains.
WINDOWS = ((469417229340, 6000, 1), (484947092225, 5625, -1))
CELLS = 8192


def majorant(u):
    return (((16*u-48)*u+75)*u+11)/54


def full_mass(a, r, c):
    g0 = 4*(a/2).sinh()/(r*c)
    constant = sum((q/arb.pi()).log() for q,_,_,_ in FAMILY)-8*(arb.const_euler()+a)

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


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 160
    a, r = arb(2).log(), arb(2).sqrt()
    c = 2*(r-1)**2
    total = full_mass(a,r,c)
    removed, move, curvature = arb(0), arb(0), arb(0)
    xmax = arb(max(left+width for left,width,_ in WINDOWS))/1000
    terms = []

    def kernel(t):
        return (3-2*r*(a*t).cos())/(c*(t*t+arb(1)/4))

    for _, _, lo, hi in read_table():
        t = arb(lo+hi)/(2*DEN)
        radius = arb(hi-lo)/(2*DEN)
        bracket = arb(lo+hi,hi-lo)/(2*DEN)
        s = (a*bracket/2).sin()
        exact_weight = -4*kernel(bracket)*s/a
        weight = -4*kernel(t)*(a*t/2).sin()/a
        removed += 4*kernel(bracket)*majorant(s*s)/a
        move += abs(exact_weight-weight)+abs(weight)*xmax*radius
        curvature += abs(weight)*t*t
        terms.append((t,weight))
    assert total-removed < arb(144)/100
    assert move < arb(43)/1000 and curvature < 343000
    print("Four cosine moments and the cubic majorant give omitted absolute mass < 1.44.")
    print("Frequency-bracket error < 0.043 on both windows; midpoint curvature < 343000.")
    for left, width, direction in WINDOWS:
        x0 = arb(left)/1000
        step = arb(width)/(1000*CELLS)
        values = [arb(0) for _ in range(CELLS+1)]
        for t,weight in terms:
            oriented = direction*weight
            sine_step, cosine_step = (t*step).sin_cos()
            for k in range(0,CELLS+1,32):
                sine,cosine = (t*(x0+k*step)).sin_cos()
                for i in range(k,min(k+32,CELLS+1)):
                    values[i] += oriented*sine
                    sine,cosine = sine*cosine_step+cosine*sine_step, cosine*cosine_step-sine*sine_step
        assert all(p > arb(164)/100 for p in values)
        assert curvature*step*step/8 < arb(24)/1000
    assert arb(164)/100-arb(144)/100-arb(43)/1000-arb(24)/1000 > arb(13)/100
    print("8193 signed midpoint nodes per window exceed 1.64; interpolation costs < 0.024.")
    print("U(x) > 0.13 throughout the six-unit gain window and U(x) < -0.13 throughout the loss window.")

    error = 2000000*arb(-250).exp()+128*arb(-499).exp()/(1-arb(-998).exp())+22000*arb(-150).exp()
    assert error < arb(1)/10**20
    assert arb(200).log() < arb(45)/8 < 6
    big, small = 2*arb(6).exp(), 2*(arb(45)/8).exp()
    assert arb(806857586985470)/10**12 < big < arb(806857586985471)/10**12
    assert arb(554544569046267)/10**12 < small < arb(554544569046268)/10**12
    print("Prime-only error < 10^-20 for x >= 1000; log(200) < 45/8 < 6.")
    print("806.857586985470 < 2 exp(6) < 806.857586985471; 554.544569046267 < 2 exp(45/8) < 554.544569046268.")


if __name__ == "__main__":
    main()
