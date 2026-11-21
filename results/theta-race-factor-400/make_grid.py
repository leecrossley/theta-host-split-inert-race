"""Build certified fixed-point coefficients directly from true zero brackets."""
import sys
from math import factorial
from pathlib import Path
from flint import acb, arb, ctx

SCALE=2**48
ORDER=16


def integer(x):
    guess=int(float(x))
    # The double estimate only selects a candidate integer. Check it with balls.
    assert abs(x-guess) < 1
    return guess


def coords(z):
    return (integer(z.real*SCALE), integer(z.imag*SCALE))


def build(table, path, size, ceiling=1620550):
    assert size in (16384,32768)
    width=size//8
    ctx.prec=144
    a, r, pi = arb(2).log(), arb(2).sqrt(), arb.pi()
    rows=[tuple(map(int,line.split())) for line in Path(table).read_text().splitlines()]
    assert len(rows)==13946
    mass, shift, all_orders = arb(0), arb(0), arb(0)
    def coefficient(t):
        return -2*(3-2*r*(a*t).cos())*(a*t/2).sin()/(a*(r-1)**2*(arb(1)/4+t*t))
    with open(path,"w") as out:
        print(len(rows),file=out)
        for j in range(size//2):
            print(*coords(acb(0,2*pi*j/size).exp()),file=out)
        for _,_,low,high in rows:
            t=arb(low+high)/(2*10**14)
            d=arb(high-low)/(2*10**14)
            v=arb(low+high,high-low)/(2*10**14)
            b=coefficient(t); mass+=abs(b)
            shift+=abs(coefficient(v)-b)+abs(b)*ceiling*d
            q=round(float(t*width/(2*pi)))
            delta=t*width/2-q*pi
            assert abs(delta) < arb(8)/5
            z=b*acb(0,t*width/2).exp()
            values=[q%size,*coords(acb(0,t*width).exp())]
            for k in range(ORDER+1):
                all_orders+=abs(z)
                values.extend(coords(z))
                z*=acb(0,delta)/(k+1)
            print(*values,file=out)
    assert mass < 346 and all_orders < 2000
    assert shift < arb(2)/1000
    trunc=mass*(arb(8)/5)**17/factorial(17)
    assert trunc < arb(1)/10**7
    unit = arb(1)/SCALE
    blocks = ceiling*8//size+1
    phase_error = 4*blocks*unit*(1+2*unit)**blocks
    input_error = phase_error*2000+6*13946*(ORDER+1)*unit
    stage_error = arb(0)
    for _ in range(size.bit_length()-1):
        stage_error = (2+4*unit)*stage_error+5000*unit
    total_fixed = input_error+(ORDER+1)*(stage_error+unit)
    assert total_fixed < arb(1)/10000
    assert 4004*SCALE < 2**61
    assert shift+trunc+total_fixed < arb(2)/1000


if __name__=="__main__":
    build(sys.argv[1],sys.argv[2],int(sys.argv[3]))
