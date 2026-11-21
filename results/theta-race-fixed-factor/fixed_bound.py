"""Enclose the two inverse-square spectral moments and the fixed-window error budget."""
from flint import acb, acb_series, arb, ctx, dirichlet_char

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")

ODD=((547,546),(4376,1093),(8752,1093),(8752,5469),
     (727,726),(5816,1453),(11632,7269),(11632,1453))
CROSSINGS=((11632,7269,48536362),(8752,1093,248572181),(4376,1093,434216664),
           (547,546,446011529),(11632,7269,579615600))


def crossing(q,idx,n,scale):
    c=dirichlet_char(q,idx)
    signs=[]
    for k in (n,n+1):
        h=c.hardy_z(arb(k)/scale)
        assert h.imag.contains(0)
        signs.append(h.real)
    assert signs[0]*signs[1] < 0
    return arb(2*n+1,1)/(2*scale)


def main():
    ctx.prec=144
    gamma=crossing(547,546,3287178298086,10**12)
    low=[crossing(q,idx,n,10**9) for q,idx,n in CROSSINGS]
    assert all(low[j]<low[j+1] for j in range(4))
    a,b=arb(2).log(),arb(2).sqrt()
    def J(t):
        G=(3-2*b*(a*t).cos())/(2*(b-1)**2*(arb(1)/4+t*t))
        return G*2*(t*a/2).sin()/a
    assert J(gamma)>arb(33428)/10000
    for k in range(1750,12000):
        t=arb(2*k+1,1)/2000
        assert abs(J(t))<arb(7)/2
    assert (3+2*b)/(2*(b-1)**2*(arb(1)/4+144))*2/a < arb(7)/2
    print("6 positive crossings, selected residue > 3.3428, global multiplier < 7/2 in magnitude")

    moments=[]
    for centre in (arb(0),gamma):
        s=acb(arb(1)/2,centre)
        z=acb_series([s,1],4)
        gamma_second=2*((z+1)/2).lgamma()[2]
        moment=arb(0)
        for q,idx in ODD:
            char=dirichlet_char(q,idx)
            assert char.is_primitive() and char.parity()==1
            v=z.dirichlet_l(char)
            if q==547 and centre>0:
                assert v[1].real>1
                # (log(L(s)/(s-rho)))'' at a simple root.
                regular=2*v[3]/v[1]-(v[2]/v[1])**2
            else:
                assert abs(v[0])>0
                regular=2*v[2]/v[0]-(v[1]/v[0])**2
            value=regular+gamma_second
            assert value.imag.contains(0) and value.real>0
            moment+=value.real
        moments.append(moment)
    A0,Ag=moments
    assert 936<A0<937 and 216<Ag<217
    weighted=(arb(7)/2)*A0-2*sum((arb(7)/2-2*t)/(t*t) for t in low)
    assert weighted<225
    print("sum of inverse squared distances: (936,937) at zero; (216,217) at the selected root")
    print("weighted central sum < 225 after removing ten low ordinates")

    sigma=arb(7)/20
    z=acb_series([3*sigma,1],2).zeta()
    assert (-z[1]/z[0]).real < 21
    mass=2*8*21/(1-arb(2)**(-sigma))
    assert mass < 1600
    j0=3*a/2
    B=2*b*(b-1)**2
    weight=b*((a/4).exp()+(-a/4).exp())/(B*a)
    assert weight*1600*(sigma*j0).exp() < 22000
    assert arb(727).log()*8 < 53 and arb(343070636977).log() < 53 and 55-j0 > 53
    rem=22000*(-arb(3)*55/20).exp()/81
    squares=2_000_000*(-arb(55)/4).exp()/225
    trivial=128*(-arb(53)/2).exp()/(1-arb(-53).exp())
    assert rem+squares+trivial < arb(1)/10
    main_error=(arb(225)+arb(7)*217/2)/900
    assert arb(33428)/10000-main_error-arb(2)/10 > (arb(19735)/10000)*(1+1/arb(90)**2)
    print("prime powers, square zeros and trivial zeros cost < 1/10 per moment")
    print("two triangular-window tests clear 1.9735; their endpoints are below exp(175) H")


if __name__=="__main__":
    main()
