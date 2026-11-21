"""Certify the spectral remainder and enclose the whole interval by linear interpolation."""
from flint import arb, ctx
from check_zeros import DEN, FAMILY, read_table

LEFT = 389568623
RIGHT = 389569151
CELLS = 4096


def full_mass(a, r, c):
    conductor = sum((q/arb.pi()).log() for q,_,_,_ in FAMILY)
    constant = conductor - 8*(arb.const_euler()+a)

    def gamma_sum(z):
        return ((1+z)*(1+z).log() + (1-z)*(1-z).log())/(2*z.sqrt())

    t1 = gamma_sum(arb(1)/2)
    t2 = gamma_sum(arb(1)/4)
    v0 = constant*4*(a/2).sinh()/(r*c) - 8*(2*r*t1-3*a)/c
    va = -8*(r*(t2+a)-3*t1)/c
    return (3*v0-va)/(2*a)


def main():
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    ctx.prec = 160
    a, r = arb(2).log(), arb(2).sqrt()
    c = 2*(r-1)**2
    left, right = arb(LEFT)/100, arb(RIGHT)/100
    step = (right-left)/CELLS
    total = full_mass(a,r,c)
    assert arb(3867375649676)/10**10 < total < arb(3867375649677)/10**10

    def kernel(t):
        return (3-2*r*(a*t).cos())/(c*(t*t+arb(1)/4))

    # U contains -2*alpha*sin(t*x) for each chosen positive zero and its partner.
    def coefficient(t):
        return 4*kernel(t)*(a*t/2).sin()/a

    nodes = [arb(0) for _ in range(CELLS+1)]
    removed, movement, curvature = arb(0), arb(0), arb(0)
    for _, _, lo, hi in read_table():
        center = arb(lo+hi)/(2*DEN)
        radius = arb(hi-lo)/(2*DEN)
        interval = arb(lo+hi, hi-lo)/(2*DEN)
        amplitude = coefficient(center)
        removed += kernel(interval)*(3-(a*interval).cos())/a
        movement += abs(coefficient(interval)-amplitude) + abs(amplitude)*right*radius
        curvature += abs(amplitude)*center**2
        angle = center*step
        sine_step, cosine_step = angle.sin_cos()
        # Re-anchor every 32 nodes to keep rounding radii small.
        for k in range(0, CELLS+1, 32):
            sine, cosine = (center*(left+k*step)).sin_cos()
            for i in range(k, min(k+32,CELLS+1)):
                nodes[i] += amplitude*sine
                sine, cosine = sine*cosine_step+cosine*sine_step, cosine*cosine_step-sine*sine_step
    tail = total-removed
    interp = curvature*step**2/8
    assert tail < arb(172)/100 and movement < arb(3)/100
    assert curvature < 320000 and interp < arb(67)/1000
    assert all(val > arb(2)/100 for val in nodes)
    assert arb(172)/100 + arb(3)/100 + arb(67)/1000 - arb(2)/100 < arb(18)/10
    print("Full majorant mass lies between 386.7375649676 and 386.7375649677.")
    print("Omitted mass < 1.72; frequency-bracket error < 0.03; midpoint curvature < 320000.")
    print("4097 grid lower bounds for the negative midpoint head exceed 0.02; interpolation costs < 0.067.")
    print("Hence U(x) < 1.8 throughout [3895686.23, 3895691.51].")

    base = 2000000 * arb(-250).exp()
    trivial = 128*arb(-499).exp()/(1-arb(-998).exp())
    powers = 22000*arb(-150).exp()
    assert base+trivial+powers < arb(1)/10**20
    factor = 2*(arb(132)/25).exp()
    assert arb(392739750703596)/10**12 < factor < arb(392739750703597)/10**12
    print("Prime-only error < 10^-20 for x >= 1000; 392.739750703596 < 2 exp(132/25) < 392.739750703597.")


if __name__ == "__main__":
    main()
